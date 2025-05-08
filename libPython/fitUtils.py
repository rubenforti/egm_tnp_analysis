#!/usr/bin/env python3

import ROOT
import os
import re
import sys
from .binUtils import binMinPt
from .rootUtils import safeOpenFile, safeGetObject
from copy import copy

from config.fit_settings import fitParsAndShapes


class fitterHandler():

    """
    Generic class that manages the base operations to run the fits using the
    compiled tnpFitter class. The key operations for the fits and the plots are
    performed by the tnpFitter class, this one is just a wrapper to help in the
    operations handling.
    """

    def __init__(self, sample, tnpBin, typeflag, strategy_name, massbins=60, massmin=60, massmax=120):
        """
        """
        self.sample = sample
        self.tnpBin = tnpBin
        self.typeflag = typeflag
        self.strategy = strategy_name
        self.strategy_dir = strategy_name[0].lower() + strategy_name[1:] + "Fit"
        self.massbins = massbins
        self.massmin = massmin
        self.massmax = massmax
        self.histFitter = None
        self.parameters = fitParsAndShapes(self.typeflag)[f"tnpPar{self.strategy}"]
        self.pdfs = fitParsAndShapes(self.typeflag)[f"tnpShapes{self.strategy}"]
        self.constraints = fitParsAndShapes(self.typeflag)["parConstraints"]


    def _removeDuplicatePar(self, lst):
        """
        Remove duplicates in the given list, keeping the latest element
        """
        seen = {}
        for el in lst:
            key = el.split('[')[0]
            seen[key] = el
        return list(seen.values())


    def _removeDuplicatePdf(self, lst):
        """
        Remove duplicates in the given list, keeping the latest element
        """
        seen = {}
        for el in lst:
            key = el.split('::')[1].split("(")[0]
            seen[key] = el
        return list(seen.values())
    
    
    def _checkExistingPars(self):
        """
        """
        valid_pdfs = []
        valid_constraints = []
        parameter_names = {par.split("[")[0] for par in self.parameters}
        allowed_params = parameter_names | {"x", "hBkgPass", "hBkgFail", "paramHistP", "paramHistF"}

        for item in [*self.pdfs, *self.constraints]:
            target_list = valid_constraints if "constrain" in item else valid_pdfs
            params = item.split("(")[1].split(")")[0].split(",")
            required_params = {par.strip().replace("{", "").replace("}", "") 
                               for par in params if not par.strip().replace('.', '', 1).isdigit() and par != "x"}

            missing_params = required_params - allowed_params
            if missing_params:
                print(f"Warning: Missing parameters {missing_params} in {item}, skipping.")
            else:
                target_list.append(item)

        self.pdfs = valid_pdfs
        self.constraints = valid_constraints

    def _applyPrefitOnPar(self, flag, parName, parObj, constrainMode=""):
        """
        Method that applies the postfit value of the fit on pure MC on the
        parameters for the fit on data, by modifying its central value. Other
        constraints can be applied with the "constrainMode" argument:
          - "gauss" : a Gaussian constraint on the parameter is applied
          - "legacy": the parameter is allowed to vary in a range of 3 sigma
          - "fix"   : the parameter is fixed to the prefit value
        """
        self.parameters.remove(parName)
        pName = parName.split("[")[0]
        pVal = parObj.getVal()
        if constrainMode == "":
            self.parameters.append(pName + f"[{pVal}, {parObj.getRange()[0]}, {parObj.getRange()[1]}]")
        elif constrainMode == "gauss":
            self.parameters.append(pName + f"[{pVal}, {parObj.getRange()[0]}, {parObj.getRange()[1]}]")
            self.constraints.append(f"RooGaussian::constrain{flag}_{pName}({pName}, {pVal}, {parObj.getError()})")
        elif constrainMode == "legacy":
            self.parameters.append(pName + f"[{pVal}, {pVal-3.0*parObj.getError()}, {pVal+3.0*parObj.getError()}]")
        elif constrainMode == "fix":
            self.parameters.append(pName + f"[{pVal}]")
        else:
            print(f"Error: Unknown constrainMode '{constrainMode}' for parameter '{pName}'")
            sys.exit(1)


    def setPars(self, tnpWorkspaceParam=""):
        """
        """
        self.parameters.extend(tnpWorkspaceParam)
        self.parameters = self._removeDuplicatePar(self.parameters)
            
    
    def setPdfs(self, tnpWorkspaceShapes=""):
        """
        """
        self.pdfs.extend(tnpWorkspaceShapes)
        self.pdfs = self._removeDuplicatePdf(self.pdfs)


    def setConstraints(self, tnWorkspaceConstraints=""):
        """
        """
        self.constraints.extend(tnWorkspaceConstraints)
        self.constraints = self._removeDuplicatePdf(self.constraints)


    def doPrefitSigOnMC(self, opts):
        """
        Method that uses the fit on a control MC sample to modify the 
        parameters of the fit on data. The instructions are passed with the 
        keyword argument "--doPrefit=sample-flag", where sample={sig,bkg} and
        flag={pass,fail,both}.
        """
        if not opts.get("doPrefit", None):
            print("Error: doPrefit option not set.")
            return
        samplePrefit, flagPrefit = opts.get("doPrefit").split("-")
        constrainMode = opts.get("constrainMode", "")

        fileref = getattr(self.sample.mcRef, "altSigFit" if samplePrefit=="sig" else "altBkgFit", None)
        if fileref is None:
            print(f"Error: No file reference found for {samplePrefit} prefit.")
            return
        
        filemc  = safeOpenFile(fileref, mode='READ')

        if flagPrefit not in ["pass", "fail", "both"]:
            print(f"Error: Invalid flag '{flagPrefit}' for prefit.")
            return
        
        for flag in [flagPrefit.capitalize()] if flagPrefit != "both" else ["Pass", "Fail"]:
            fitres = safeGetObject(filemc, f'{self.tnpBin['name']}_res{flag[0]}', detach=False)
            listOfParam = copy([tnpPar for tnpPar in self.parameters if flag[0] in tnpPar])
            for par in fitres.floatParsFinal():
                for parNameExtended in listOfParam:
                    x = re.compile('%s\\[.*?' % par.GetName())
                    if x.match(parNameExtended):
                        self._applyPrefitOnPar(flag, parNameExtended, par, constrainMode)
                        
        filemc.Close()


    def initFitter(self, opts):
        """
        """
        infile = ROOT.TFile(self.sample.getOutputPath(), "READ")
        hP = infile.Get(f"{self.tnpBin['name']}_Pass")
        hF = infile.Get(f"{self.tnpBin['name']}_Fail")
        self.histFitter = ROOT.tnpFitter(hP, hF, self.tnpBin["name"], self.massbins, self.massmin, self.massmax)
        infile.Close()

        self.histFitter.setPassStrategy(opts.get("passStrategy", 2))
        self.histFitter.setFailStrategy(opts.get("failStrategy", 2))
        self.histFitter.setPrintLevel(opts.get("printLevel", -1))
        self.histFitter.isMC(opts.get("isMC", False))

        self.outfileName = opts.get("outputFileName", 
                                    f"{self.sample.nominalFit.replace("nominalFit.root", self.strategy_dir)}_bin_{self.tnpBin["name"]}.root")
        print(f"Output file: {self.outfileName}")
        self.histFitter.setOutputFile(self.outfileName)

        self.plotPath = opts.get("plotPath", 
                                 os.path.abspath(os.path.dirname(self.outfileName)) + f"/plots/{self.sample.getName()}/{self.strategy_dir}/")
        self.histFitter.setPlotOutputPath(self.plotPath)
        

    def setZLineShapes(self, opts):
        """
        Sets the DY templates in the tnpFitter class by searching in the
        apposite file.
        The keyword argument "useAllProbesForFail" is a further control over
        the failing probes template:
          - useAllProbesForFail==0: default template (only failing probes) 
          - otherwise: add to the default template the passing probes template
                       (using standalone variables), if useAllProbesForFail is
                       negative or the integral of failing probes is smaller
                       than its value
        """
        fileDY = ROOT.TFile(self.sample.mcRef.getOutputPath(), "READ")
        histZLineShapeP = fileDY.Get(f"{self.tnpBin['name']}_Pass")
        histZLineShapeF = fileDY.Get(f"{self.tnpBin['name']}_Fail")
    
        useAllProbesForFail = opts.get("useAllProbesForFail", 0)
        if (useAllProbesForFail>0 and histZLineShapeF.Integral()<useAllProbesForFail) or useAllProbesForFail<0:
            altPass = f"{self.tnpBin['name']}_Pass_alt"
            histZLineShapeP_alt = fileDY.Get(altPass) if altPass in [k.GetName() for k in fileDY.GetListOfKeys()] else None
            histZLineShapeF.Add(histZLineShapeP_alt if histZLineShapeP_alt else histZLineShapeP)
        
        self.histFitter.setZLineShapes(histZLineShapeP, histZLineShapeF)
        fileDY.Close()


    def setBkgTemplatePdf(self, opts):
        """
        """
        fileBkg = ROOT.TFile(self.sample.bkgRef.getOutputPath(),'read')
        histBkgP = fileBkg.Get(f"{self.tnpBin['name']}_Pass")
        histBkgF = fileBkg.Get(f"{self.tnpBin['name']}_Fail")

        self.histFitter.setTotalBkgShapes(histBkgP, histBkgF)

        if opts.get("useBBFail", True):
            self.histFitter.setBarlowBeestonBkgPdf(False)
            self.constraints.extend(["RooHistConstraint::constrainF_histConstr(paramHistF)"])
        else:
            self.pdfs.extend(["RooHistPdf::bkgFail(x, hBkgFail, 0)"])

        self.pdfs.extend(["RooHistPdf::bkgFailBackup(x, hBkgFail, 0)"])
        
        fileBkg.Close()


    def importConstraints(self):
        """
        """
        constrP = ",".join([x.split("::")[1].split("(")[0] for x in self.constraints if "constrainP" in x])
        constrF = ",".join([x.split("::")[1].split("(")[0] for x in self.constraints if "constrainF" in x])
        if len(constrP):
            self.histFitter.updateConstraints("constrainP", constrP)
        if len(constrF):    
            self.histFitter.updateConstraints("constrainF", constrF)


    def importAndFit(self, analyticPhysicsShape=False, modelFSR=False):
        """
        """
        self._checkExistingPars()
        workspace = ROOT.vector("string")()
        for iw in [*self.parameters, *self.pdfs, *self.constraints]:
            workspace.push_back(iw)
       
        self.histFitter.setWorkspace(workspace, self.sample.isMonteCarlo(), analyticPhysicsShape, modelFSR)
        self.importConstraints()

        title = self.tnpBin['title'].replace(';',' - ')
        self.histFitter.fits(title)


###############################################################################


def histFitterAllTemplate(sample, tnpBin, typeflag, strategy_name, 
                          massbins=60, massmin=60, massmax=120, 
                          tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[], **opts):
    """
    Signal: MC gaussian-smeared template;
    Background: exponential for Pass, MC template for Fail enhanced with bin
                content fluctuations (Barlow-Beeston method) if useBBFail=True.
    """
    fitter = fitterHandler(sample, tnpBin, typeflag, strategy_name, massbins, massmin, massmax)

    fitter.setPars(tnpWorkspaceParam)
    fitter.setPdfs(tnpWorkspaceFunc)
    fitter.setConstraints(constrainPars)

    fitter.initFitter(opts)
    fitter.setZLineShapes(opts)
    fitter.setBkgTemplatePdf(opts)
    fitter.importAndFit()



def histFitterAnalyticSig(sample, tnpBin, typeflag, strategy_name, 
                          massbins=60, massmin=60, massmax=120,
                          tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[], **opts):
    """
    Signal: fixed Breit-Wigner, smeared with a double-tailed Crystal-Ball (a
            Gaussian for Fail if symmConvSigFail=True), plus a Gaussian to 
            model the FSR effect if modelFSR=True; 
    Background: exponential for Pass, MC template for Fail enhanced with bin
                content fluctuations (Barlow-Beeston method) if useBBFail=True.
    """

    fitter = fitterHandler(sample, tnpBin, typeflag, strategy_name, massbins, massmin, massmax)

    if binMinPt(tnpBin) >= 35:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if x!="tailLeft[1]"]
    else:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if x!="tailLeft[-1]"]

    if opts.get("symmConvSigFail", False):
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "RooCBExGaussShape::sigResFail" not in x]
    else:
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "Gaussian::sigResFail" not in x]

    modelFSR = opts.get("modelFSR", True if typeflag in ["iso", "trigger", "isonotrig"] else False)

    if not modelFSR:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if "fsr" not in x]
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "Gaussian::sigFsrFail" not in x]

    fitter.setPars(tnpWorkspaceParam)
    fitter.setPdfs(tnpWorkspaceFunc)
    fitter.setConstraints(constrainPars)

    fitter.initFitter(opts)
    fitter.setZLineShapes(opts)
    fitter.setBkgTemplatePdf(opts)
    fitter.importAndFit(analyticPhysicsShape=True, modelFSR=modelFSR)


def histFitterAnalyticBkg(sample, tnpBin, typeflag, strategy_name,
                          massbins=60, massmin=60, massmax=120, 
                          tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[], **opts):
    """
    Signal: MC gaussian-smeared template;
    Background: exponential for Pass, exponential or RooCMSShape for Fail.
    """

    fitter = fitterHandler(sample, tnpBin, typeflag, strategy_name, massbins, massmin, massmax)

    fitter.setPars(tnpWorkspaceParam)
    fitter.setPdfs(tnpWorkspaceFunc)
    fitter.setConstraints(constrainPars)

    fitter.initFitter(opts)
    fitter.setZLineShapes(opts)
    fitter.importAndFit()



def histFitterAllAnalytic(sample, tnpBin, typeflag, strategy_name, 
                          massbins=60, massmin=60, massmax=120, 
                          tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[], **opts):

    fitter = fitterHandler(sample, tnpBin, typeflag, strategy_name, massbins, massmin, massmax)

    if binMinPt(tnpBin) >= 35:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if x!="tailLeft[1]"]
    else:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if x!="tailLeft[-1]"]

    if opts.get("symmConvSigFail", False):
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "RooCBExGaussShape::sigResFail" not in x]
    else:
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "Gaussian::sigResFail" not in x]

    modelFSR = opts.get("modelFSR", True if typeflag in ["iso", "trigger", "isonotrig"] else False)

    if not modelFSR:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if "fsr" not in x]
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "Gaussian::sigFsrFail" not in x]

    fitter.setPars(tnpWorkspaceParam)
    fitter.setPdfs(tnpWorkspaceFunc)
    fitter.setConstraints(constrainPars)
    if opts.get("isMC", False) is False:
        fitter.doPrefitSigOnMC(opts)
    fitter.initFitter(opts)
    fitter.setZLineShapes(opts)
    fitter.importAndFit(analyticPhysicsShape=True, modelFSR=modelFSR)
