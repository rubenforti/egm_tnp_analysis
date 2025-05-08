### python specific import

## example
# python libPython/checkFitStatus.py plots/results_test_globalMuons_byCharge_noMinos_RooMinimizerMinuit2//efficiencies_GtoH/mu_iso_plus/mu_RunGtoH_mu_iso_plus.nominalFit.root

import os, sys
import copy
import argparse
import ROOT

## safe batch mode
args = sys.argv[:]
sys.argv = ['-b']
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from libPython.rootUtils import safeGetObject, safeOpenFile
from libPython.plotUtils import createPlotDirAndCopyPhp, drawTH2

sys.path.append(os.getcwd() + "/libPython/")


def checkFit(infile, outdir, fitName, hbins):

    # TODO: avoid repeating code for each histogram, but fine for now
    
    hStatusPass = copy.deepcopy(hbins.Clone(f"{fitName}_pass_status"))
    hStatusPass.Reset("ICESM")
    hStatusPass.SetTitle(f"{fitName} pass")
    hStatusFail = copy.deepcopy(hbins.Clone(f"{fitName}_fail_status"))
    hStatusFail.Reset("ICESM")
    hStatusFail.SetTitle(f"{fitName} fail")

    hCovQualPass = copy.deepcopy(hbins.Clone(f"{fitName}_pass_covQual"))
    hCovQualPass.Reset("ICESM")
    hCovQualPass.SetTitle(f"{fitName} pass")
    hCovQualFail = copy.deepcopy(hbins.Clone(f"{fitName}_fail_covQual"))
    hCovQualFail.Reset("ICESM")
    hCovQualFail.SetTitle(f"{fitName} fail")

    hMeanPass = copy.deepcopy(hbins.Clone(f"{fitName}_pass_mean"))
    hMeanPass.Reset("ICESM")
    hMeanPass.SetTitle(f"{fitName} pass")
    hMeanFail = copy.deepcopy(hbins.Clone(f"{fitName}_fail_mean"))
    hMeanFail.Reset("ICESM")
    hMeanFail.SetTitle(f"{fitName} fail")

    # sigma only for nominal fit, it means how much the MC template is smeared, for the alternate fits there can be more sigma parameters instead 
    if "nominal" in fitName:
        hSigmaPass = copy.deepcopy(hbins.Clone(f"{fitName}_pass_sigma"))
        hSigmaPass.Reset("ICESM")
        hSigmaPass.SetTitle(f"{fitName} pass")
        hSigmaFail = copy.deepcopy(hbins.Clone(f"{fitName}_fail_sigma"))
        hSigmaFail.Reset("ICESM")
        hSigmaFail.SetTitle(f"{fitName} fail")
    
    nEtaBins = hStatusPass.GetNbinsX()
    nPtBins = hStatusPass.GetNbinsY()
    
    f = safeOpenFile(infile)
    for k in f.GetListOfKeys():
        name = k.GetName()
        if "_res" not in name:
            continue
        obj = safeGetObject(f, name, detach=False)
        nbin = int(name.split("_")[0].lstrip("bin"))
        neta = int((nbin % nEtaBins) + 1)
        npt = int((nbin / nEtaBins) + 1)
        if "_resP" in name:
            hStatusPass.SetBinContent(neta, npt, obj.status())
            hCovQualPass.SetBinContent(neta, npt, obj.covQual())
        else:
            hStatusFail.SetBinContent(neta, npt, obj.status())
            hCovQualFail.SetBinContent(neta, npt, obj.covQual())
        for par in obj.floatParsFinal():
            if par.GetName() == "meanP":
                hMeanPass.SetBinContent(neta, npt, par.getVal())
            elif par.GetName() == "meanF":
                hMeanFail.SetBinContent(neta, npt, par.getVal())
            if "nominal" in fitName:
                if par.GetName() == "sigmaP":
                    hSigmaPass.SetBinContent(neta, npt, par.getVal())
                elif par.GetName() == "sigmaF":
                    hSigmaFail.SetBinContent(neta, npt, par.getVal())

    f.Close()

    canvas = ROOT.TCanvas("canvas","",800,800)

    maxZval = int(hStatusPass.GetBinContent(hStatusPass.GetMaximumBin()))
    zrange = f" max={maxZval}::-0.5,4.5"    

    drawTH2(hStatusPass, "Muon #eta", "Muon p_{T} (GeV)", f"Fit status {zrange}",
            hStatusPass.GetName(), plotLabel="ForceTitle", outdir=outdir, 
            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
            palette=87, nContours=5, drawOption="colz")

    maxZval = int(hStatusFail.GetBinContent(hStatusFail.GetMaximumBin()))
    zrange = f" max={maxZval}::-0.5,4.5"    

    drawTH2(hStatusFail, "Muon #eta", "Muon p_{T} (GeV)", f"Fit status {zrange}",
            hStatusFail.GetName(), plotLabel="ForceTitle", outdir=outdir, 
            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
            palette=87, nContours=5, drawOption="colz")


    maxZval = int(hCovQualPass.GetBinContent(hCovQualPass.GetMaximumBin()))
    zrange = f" max={maxZval}::-1.5,4.5"    

    drawTH2(hCovQualPass, "Muon #eta", "Muon p_{T} (GeV)", f"Fit cov. quality {zrange}",
            hCovQualPass.GetName(), plotLabel="ForceTitle", outdir=outdir, 
            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
            palette=87, nContours=6, drawOption="colz")

    maxZval = int(hCovQualFail.GetBinContent(hCovQualFail.GetMaximumBin()))
    zrange = f" max={maxZval}::-1.5,4.5"    

    drawTH2(hCovQualFail, "Muon #eta", "Muon p_{T} (GeV)", f"Fit cov. quality {zrange}",
            hCovQualFail.GetName(), plotLabel="ForceTitle", outdir=outdir, 
            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
            palette=87, nContours=6, drawOption="colz")

    maxZval = round(float(hMeanPass.GetBinContent(hStatusPass.GetMaximumBin())), 3)
    zrange = f" max={maxZval}"    

    drawTH2(hMeanPass, "Muon #eta", "Muon p_{T} (GeV)", f"Gauss mean {zrange} (GeV)",
            hMeanPass.GetName(), plotLabel="ForceTitle", outdir=outdir, 
            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
            palette=87, nContours=51, drawOption="colz")

    maxZval = round(float(hMeanFail.GetBinContent(hMeanFail.GetMaximumBin())), 3)
    zrange = f" max={maxZval}"    

    drawTH2(hMeanFail, "Muon #eta", "Muon p_{T} (GeV)", f"Gauss mean {zrange} (GeV)",
            hMeanFail.GetName(), plotLabel="ForceTitle", outdir=outdir, 
            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
            palette=87, nContours=51, drawOption="colz")

    if "nominal" in fitName:

        maxZval = round(float(hSigmaPass.GetBinContent(hStatusPass.GetMaximumBin())), 3)
        zrange = f" max={maxZval}"    
        
        drawTH2(hSigmaPass, "Muon #eta", "Muon p_{T} (GeV)", f"Gauss sigma {zrange} (GeV)",
                hSigmaPass.GetName(), plotLabel="ForceTitle", outdir=outdir, 
                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                palette=87, nContours=51, drawOption="colz")
        
        maxZval = round(float(hSigmaFail.GetBinContent(hSigmaFail.GetMaximumBin())), 3)
        zrange = f" max={maxZval}"    
    
        drawTH2(hSigmaFail, "Muon #eta", "Muon p_{T} (GeV)", f"Gauss sigma {zrange} (GeV)",
                hSigmaFail.GetName(), plotLabel="ForceTitle", outdir=outdir, 
                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                palette=87, nContours=51, drawOption="colz")


    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Diagnostic for fit status')
    parser.add_argument("infile", type=str, nargs=1, help="Input file")
    #parser.add_argument("outdir", type=str, nargs=1, help="Output folder to save plots")
    args = parser.parse_args()

    infile = args.infile[0]
    #outdir = args.outdir[0] 
    mainPath = os.path.dirname(infile) + "/"
    outdir = mainPath + "plots/checkFitStatus/"
    createPlotDirAndCopyPhp(outdir)    

    # get histogram to read binning, might need to do it differently if this file does not exist yet
    rootfileWithEffi = safeOpenFile(mainPath + "allEfficiencies_2D.root")
    htmp = safeGetObject(rootfileWithEffi, "SF2D_nominal", detach=True)
    rootfileWithEffi.Close()
    
    tag = "MC" if "_DY_" in infile else "Data"
    fitName = f"Eff{tag}_" + infile.split("_")[-1].split(".")[0]

    checkFit(infile, outdir, fitName, htmp)

