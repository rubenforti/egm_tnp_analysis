#!/usr/bin/env python3

import ROOT
import os

from ..config.fit_settings import fitParsAndShapes

def ptMin( tnpBin ):
    ptmin = 1
    if tnpBin['name'].find('pt_') >= 0:
        ptmin = float(tnpBin['name'].split('pt_')[1].split('p')[0])
    elif tnpBin['name'].find('et_') >= 0:
        ptmin = float(tnpBin['name'].split('et_')[1].split('p')[0])
    return ptmin


def histFitterNominal(sample, tnpBin, massbins=60, massmin=60, massmax=120, 
                      tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[], 
                      useAllTemplateForFail=False, maxFailIntegralToUseAllProbe=-1):
    """
    Nominal fit: signal made by MC gaussian-smeared template; background is an
    exponential for Pass, and the MC template for Fail enhanced with bin content
    fluctuations (Barlow-Beeston technique) if useBBfail==True.
    """

    analyticPhysicsShape = False
    useBBfail = True

    if useBBfail:
        defaultBkgShapes = ["Exponential::bkgPass(x, expalphaP)",
                            "RooHistPdf::bkgFailBackup(x, hBkgFail, 0)"]
    else:
        defaultBkgShapes = ["Exponential::bkgPass(x, expalphaP)", 
                            "RooHistPdf::bkgFail(x, hBkgFail, 0)",
                            "RooHistPdf::bkgFailBackup(x, hBkgFail, 0)"]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)

    # Init fitter
    infile = ROOT.TFile(sample.getOutputPath(), "READ")
    hP = infile.Get(f"{tnpBin['name']}_Pass")
    hF = infile.Get(f"{tnpBin['name']}_Fail")
    fitter = ROOT.tnpFitter(hP, hF, tnpBin["name"], massbins, massmin, massmax)
    infile.Close()

    # Fitter setup
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    outFileName = sample.nominalFit.rstrip(".root") + "_bin_" + tnpBin["name"] + ".root"
    fitter.setOutputFile(outFileName)
    plotPath = os.path.abspath(os.path.dirname(outFileName)) + f"/plots/{sample.getName()}/nominalFit/"
    fitter.setPlotOutputPath(plotPath)
    fitter.isMC(sample.isMonteCarlo())

    # Z lineshape
    fileDY = ROOT.TFile(sample.mcRef.getOutputPath(), "READ")
    histZLineShapeP = fileDY.Get(f"{tnpBin['name']}_Pass")
    histZLineShapeF = fileDY.Get(f"{tnpBin['name']}_Fail")
    altPass = f"{tnpBin['name']}_Pass_alt"
    if altPass in [k.GetName() for k in fileDY.GetListOfKeys()]:
        histZLineShapeP_alt = fileDY.Get(altPass)
    else:
        histZLineShapeP_alt = None

    if useAllTemplateForFail:
        if maxFailIntegralToUseAllProbe < 0 or histZLineShapeF.Integral() < maxFailIntegralToUseAllProbe:
            histZLineShapeF.Add(histZLineShapeP_alt if histZLineShapeP_alt else histZLineShapeP)

    fitter.setZLineShapes(histZLineShapeP, histZLineShapeF)
    fileDY.Close()

    # Background templates
    fileBkg = ROOT.TFile(sample.bkgRef.getOutputPath(),'read')
    histBkgP = fileBkg.Get('%s_Pass'%tnpBin['name'])
    histBkgF = fileBkg.Get('%s_Fail'%tnpBin['name'])

    fitter.setTotalBkgShapes(histBkgP, histBkgF)

    if useBBfail:
      fitter.setBarlowBeestonBkgPdf(False)
      constrainPars.extend(["RooHistConstraint::constrainF_histConstr(paramHistF)"])
    
    fileBkg.Close()
    
    # Constrain parameters
    if len(constrainPars):
        tnpWorkspace.extend(constrainPars)
        constrP = ",".join([x.split("::")[1].split("(")[0] for x in constrainPars if "constrainP" in x])
        constrF = ",".join([x.split("::")[1].split("(")[0] for x in constrainPars if "constrainF" in x])
        if len(constrP):
            fitter.updateConstraints("constrainP", constrP)
        if len(constrF):    
            fitter.updateConstraints("constrainF", constrF)
                

    # Set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace(workspace, sample.isMonteCarlo(), analyticPhysicsShape, False)

    title = tnpBin['title'].replace(';',' - ')

    # Perform the fit
    fitter.fits(title)


def histFitterAltSig(sample, tnpBin, massbins=60, massmin=60, massmax=120,
                     tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[],
                     symmConvSigFail=False, modelFSR=False):
    """
    Alternative signal fit: signal made by a fixed Breit-Wigner, smeared with a
    double-tailed Crystal-Ball (a Gaussian for Fail if symmConvSigFail==True), 
    plus a Gaussian to model the FSR effect if modelFSR==True; background is an
    exponential for Pass, and the MC template for Fail enhanced with bin 
    content fluctuations (Barlow-Beeston technique) if useBBfail==True.
    """

    analyticPhysicsShape = True
    useBBfail = True

    if useBBfail:
        defaultBkgShapes = ["Exponential::bkgPass(x, expalphaP)",
                            "RooHistPdf::bkgFailBackup(x, hBkgFail, 0)"]
    else:
        defaultBkgShapes = ["Exponential::bkgPass(x, expalphaP)", 
                            "RooHistPdf::bkgFail(x, hBkgFail, 0)",
                            "RooHistPdf::bkgFailBackup(x, hBkgFail, 0)"]

    if ptMin(tnpBin) >= 35:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if x!="tailLeft[1]"]
    else:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if x!="tailLeft[-1]"]

    if symmConvSigFail:
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "RooCBExGaussShape::sigResFail" not in x]
    else:
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "Gaussian::sigResFail" not in x]

    if not modelFSR:
        tnpWorkspaceParam = [x for x in tnpWorkspaceParam if "fsr" not in x]
        tnpWorkspaceFunc = [x for x in tnpWorkspaceFunc if "Gaussian::sigFsrFail" not in x]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)

    print(tnpWorkspace)

   # Init fitter
    infile = ROOT.TFile(sample.getOutputPath(), "READ")
    hP = infile.Get(f"{tnpBin['name']}_Pass")
    hF = infile.Get(f"{tnpBin['name']}_Fail")
    fitter = ROOT.tnpFitter(hP, hF, tnpBin["name"], massbins, massmin, massmax)
    infile.Close()

    # Fitter setup
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    outFileName = sample.altSigFit.rstrip(".root") + "_bin_" + tnpBin["name"] + ".root"
    fitter.setOutputFile(outFileName)
    plotPath = os.path.abspath(os.path.dirname(outFileName)) + f"/plots/{sample.getName()}/altSigFit/"
    fitter.setPlotOutputPath(plotPath)
    fitter.isMC(sample.isMonteCarlo())

    # Background templates
    fileBkg = ROOT.TFile(sample.bkgRef.getOutputPath(),"READ")
    histBkgP = fileBkg.Get('%s_Pass'%tnpBin['name'])
    histBkgF = fileBkg.Get('%s_Fail'%tnpBin['name'])
    fitter.setTotalBkgShapes(histBkgP, histBkgF)

    if useBBfail:
      fitter.setBarlowBeestonBkgPdf(False)
      constrainPars.extend(["RooHistConstraint::constrainF_histConstr(paramHistF)"])
    
    fileBkg.Close()
    
    # Constrain parameters
    if len(constrainPars):
        tnpWorkspace.extend(constrainPars)
        constrP = ",".join([x.split("::")[1].split("(")[0] for x in constrainPars if "constrainP" in x])
        constrF = ",".join([x.split("::")[1].split("(")[0] for x in constrainPars if "constrainF" in x])
        if len(constrP):
            fitter.updateConstraints("constrainP", constrP)
        if len(constrF):    
            fitter.updateConstraints("constrainF", constrF)
    
    # Set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace, sample.isMonteCarlo(), analyticPhysicsShape, modelFSR )

    title = tnpBin['title'].replace(';',' - ')
    fitter.fits(title)

    #rootfile.Close()


def histFitterAltBkg(sample, tnpBin, massbins=60, massmin=60, massmax=120, 
                     tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[],
                     useAllTemplateForFail=False, maxFailIntegralToUseAllProbe=-1):
    """
    Alternative background fit: signal made by MC gaussian-smeared template;
    background pdf is analytical, by default a RooCMSShape.
    """

    analyticPhysicsShape = False

    
    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)        

    
    # Init fitter
    infile = ROOT.TFile(sample.getOutputPath(), "READ")
    hP = infile.Get(f"{tnpBin['name']}_Pass")
    hF = infile.Get(f"{tnpBin['name']}_Fail")
    fitter = ROOT.tnpFitter(hP, hF, tnpBin["name"], massbins, massmin, massmax)
    infile.Close()

    # Fitter setup
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    outFileName = sample.altBkgFit.rstrip(".root") + "_bin_" + tnpBin["name"] + ".root"
    fitter.setOutputFile(outFileName)
    plotPath = os.path.abspath(os.path.dirname(outFileName)) + f"/plots/{sample.getName()}/altBkgFit/"
    fitter.setPlotOutputPath(plotPath)
    fitter.isMC(sample.isMonteCarlo())
    

    # Z lineshape
    fileDY = ROOT.TFile(sample.mcRef.getOutputPath(), "READ")
    histZLineShapeP = fileDY.Get(f"{tnpBin['name']}_Pass")
    histZLineShapeF = fileDY.Get(f"{tnpBin['name']}_Fail")
    altPass = f"{tnpBin['name']}_Pass_alt"
    if altPass in [k.GetName() for k in fileDY.GetListOfKeys()]:
        histZLineShapeP_alt = fileDY.Get(altPass)
    else:
        histZLineShapeP_alt = None

    if useAllTemplateForFail:
        if maxFailIntegralToUseAllProbe < 0 or histZLineShapeF.Integral() < maxFailIntegralToUseAllProbe:
            histZLineShapeF.Add(histZLineShapeP_alt if histZLineShapeP_alt else histZLineShapeP)

    fitter.setZLineShapes(histZLineShapeP, histZLineShapeF)
    fileDY.Close()

    # Constrain parameters
    if len(constrainPars):
        tnpWorkspace.extend(constrainPars)
        constrP = ",".join([x.split("::")[1].split("(")[0] for x in constrainPars if "constrainP" in x])
        constrF = ",".join([x.split("::")[1].split("(")[0] for x in constrainPars if "constrainF" in x])
        if len(constrP):
            fitter.updateConstraints("constrainP", constrP)
        if len(constrF):    
            fitter.updateConstraints("constrainF", constrF)
            
    ### set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace(workspace, sample.isMonteCarlo(), analyticPhysicsShape, False )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_eta','#eta')
    title = title.replace('probe_pt','p_{T}')
    fitter.fits(title)


def histFitterAllAnalytic(sample, tnpBin, massbins=60, massmin=60, massmax=120, 
                          tnpWorkspaceParam=[], tnpWorkspaceFunc=[], constrainPars=[],
                          altSignalFail=False, modelFSR=False, constrainSignalFailFromMC=False, zeroBackground=False):

    analyticPhysicsShape = True
    
    '''
    if sample.isMonteCarlo():
        tnpWorkspacePar = tnpWorkspaceParam
    else:
        tnpWorkspacePar = createWorkspaceForAltSig( sample,  tnpBin, tnpWorkspaceParam, constrainSignalFailFromMC=constrainSignalFailFromMC )
    '''

    ### tricky: use n < 0 for high pT bin (so need to remove param and add it back)
    ptmin = ptMin(tnpBin)        

    defaultBkgShapes = [
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
    ]
    
    if altSignalFail and not sample.isMonteCarlo():
        tnpWorkspaceFunc = [
            "tailLeft[%d]" % (-1 if ptmin >= 35 else 1),
            "RooCBExGaussShape::sigResPass(x,meanP,sigmaP,alphaP,nP,sigmaP_2,tailLeft)",
            "Gaussian::sigResFail(x,meanF,sigmaF)",
        ]
    else:
        tnpWorkspaceFunc = [
            "tailLeft[%d]" % (-1 if ptmin >= 35 else 1),
            "RooCBExGaussShape::sigResPass(x,meanP,sigmaP,alphaP,nP,sigmaP_2,tailLeft)",
            "RooCBExGaussShape::sigResFail(x,meanF,sigmaF,alphaF,nF,sigmaF_2,tailLeft)",
        ]

    tnpWorkspaceFunc.extend(bkgShapes if len(bkgShapes) else defaultBkgShapes)
        
    bwShapes = ["BreitWigner::sigPhysPass(x,91.1876,2.4952)",
                "BreitWigner::sigPhysFail(x,91.1876,2.4952)"]
    tnpWorkspaceFunc.extend(bwShapes)

    if modelFSR:
        tnpWorkspaceFunc.append("Gaussian::sigFsrFail(x,fsrMeanF,fsrSigmaF)")

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)

    ## init fitter
    infile = ROOT.TFile( sample.getOutputPath(), "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = ROOT.tnpFitter( hP, hF, tnpBin['name'], massbins, massmin, massmax )
    #    fitter.fixSigmaFtoSigmaP()
    infile.Close()
    
    ## setup
    ## make configurable from outside
    # fitter.useMinos()
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    outFileName = sample.altSigFit.rstrip(".root") + "_bin_" + tnpBin["name"] + ".root"
    fitter.setOutputFile(outFileName)
    plotPath = os.path.abspath(os.path.dirname(outFileName)) + f"/plots/{sample.getName()}/altSigFit/"
    fitter.setPlotOutputPath(plotPath)
    fitter.isMC(sample.isMonteCarlo())

    if zeroBackground:
       fitter.setZeroBackground()
    
    if len(constrainPars):
        tnpWorkspace.extend(constrainPars)
        # ugly, just to define constraints for passing or failing
        constraints = {"constrainP" : "",
                       "constrainF" : ""}
        cpass = list(filter(lambda x :  "constrainP" in x, constrainPars))
        cfail = list(filter(lambda x :  "constrainF" in x, constrainPars))
        constraints = {"constrainP" : ",".join([x.split("::")[1].split("(")[0] for x in cpass]),
                       "constrainF" : ",".join([x.split("::")[1].split("(")[0] for x in cfail])}
        for key in constraints.keys():
            if len(constraints[key]):
                fitter.updateConstraints(key, constraints[key])

    ### set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace, sample.isMonteCarlo(), analyticPhysicsShape, modelFSR)

    title = tnpBin['title'].replace(';',' - ')
    fitter.fits(title)

    #rootfile.Close()


'''
# TO BE REWRITTEN: IT SHOULD BE A FUNCTION THAT SETS THE INITIAL PARAMETERS BASING ON A PREVIOUS FIT IN GENERAL,
# NOT ONLY FOR ALTERNATIVE SIGNAL FIT, AND DOES NOT LIMIT THEIR RANGE; A CONSTRAIN CAN INDEED BE CREATED
def createWorkspaceForAltSig( sample, tnpBin, tnpWorkspaceParam, constrainSignalFailFromMC=False):
    
    fileref = sample.mcRef.altSigFit
    filemc  = safeOpenFile(fileref, mode='READ')

    fitresP = safeGetObject(filemc, '%s_resP' % tnpBin['name'], detach=False)
    fitresF = safeGetObject(filemc, '%s_resF' % tnpBin['name'], detach=False)

    listOfParamP = ['meanP', 'sigmaP', 'nP', 'alphaP', 'sigmaP', 'sigmaP_2']
    listOfParamF = ['meanF', 'sigmaF', 'nF', 'alphaF', 'sigmaF', 'sigmaF_2', 'fsrMeanF', 'fsrSigmaF']
    #listOfParamP = ['nP', 'alphaP', 'sigmaP', 'sigmaP_2']
    #listOfParamF = ['nF', 'alphaF', 'sigmaF', 'sigmaF_2', 'fsrMeanF', 'fsrSigmaF']

    # set central value of signal parameters as in MC alt sig fit, but do not fix them
    # while those for background from the nominal fit in data (but only for failing probes, passing ones are good)
    
    # first for failing probes and signal parameters
    fitPar = fitresF.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        #print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
        x = re.compile('%s\[.*?' % pName)
        for par in listOfParamF:
            if pName == par:
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                if len(listToRM):
                    # for sigma since there is also sigma_2, but it usually picks the correct value when the first element is taken
                    #if len(listToRM) > 1:
                    #    print(f"Error: listToRM has more than 1 element: {listToRM}")
                    #    quit()
                    ir = listToRM[0] # should always be only 1 element, otherwise adding it back below becomes a problem
                    #print(f">>>>> old {ir}")
                    tnpWorkspaceParam.remove(ir)
                    if constrainSignalFailFromMC:
                        newval = fitPar[ipar].getVal()
                        newerr = fitPar[ipar].getError()
                        new_ir = "%s[%2.3f,%2.3f,%2.3f]" % (pName, newval, newval-3.0*newerr, newval+3.0*newerr)
                    else:
                        parRange = ir.split("[")[1].split("]")[0].split(",")  # get elements within square brackets, "ir" is like "name[value,low,high]"
                        parRange = ",".join(parRange[1:]) # concatenate everything except the first element                    
                        new_ir = "%s[%2.3f,%s]" % (pName, fitPar[ipar].getVal(), parRange)
                    #print(f">>>>> new {new_ir}")
                    tnpWorkspaceParam.append( new_ir )
                    
    # now for passing probes
    fitPar = fitresP.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        #print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
        x = re.compile('%s\[.*?' % pName)
        for par in listOfParamP:
            if pName == par:
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                if len(listToRM):
                    # for sigma since there is also sigma_2, but it usually picks the correct value when the first element is taken
                    #if len(listToRM) > 1:
                    #    print(f"Error: listToRM has more than 1 element: {listToRM}")
                    #    quit()
                    ir = listToRM[0] # should always be only 1 element, otherwise adding it back below becomes a problem
                    #print(f">>>>> old {ir}")
                    tnpWorkspaceParam.remove(ir)
                    parRange = ir.split("[")[1].split("]")[0].split(",")  # get elements within square brackets, "ir" is like "name[value,low,high]"
                    parRange = ",".join(parRange[1:]) # concatenate everything except the first element
                    new_ir = "%s[%2.3f,%s]" % (pName, fitPar[ipar].getVal(), parRange)
                    #print(f">>>>> new {new_ir}")
                    tnpWorkspaceParam.append( new_ir )
                    
    filemc.Close()

    return tnpWorkspaceParam
'''