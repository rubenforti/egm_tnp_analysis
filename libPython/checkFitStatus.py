### python specific import

## example
# python libPython/checkFitStatus.py plots/results_test_globalMuons_byCharge_noMinos_RooMinimizerMinuit2//efficiencies_GtoH/mu_iso_plus/mu_RunGtoH_mu_iso_plus.nominalFit.root

import os, sys
import argparse
import ROOT
from array import array

## safe batch mode
args = sys.argv[:]
sys.argv = ['-b']
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from libPython.rootUtils import safeGetObject, safeOpenFile
from libPython.plotUtils import createPlotDirAndCopyPhp, drawTH2

sys.path.append(os.getcwd() + "/libPython/")


def checkFit(infile, fitName, binningDef, outdir):

    bins_eta = array("d", binningDef["eta"]["bins"])
    bins_pt  = array("d", binningDef["pt"]["bins"] )

    nEtaBins, nPtBins = len(bins_eta)-1, len(bins_pt)-1
    
    toPlot = {
        "status":  (None, None),
        "covQual": (None, None),
        "mean":    (None, None),
        "sigma":   (None, None)}

    file = safeOpenFile(infile)

    for i, plot in enumerate(toPlot):

        h2_pass = ROOT.TH2D(f"{fitName}_pass_{plot}", f"{fitName} pass - {plot}", nEtaBins, bins_eta, nPtBins, bins_pt)
        h2_fail = ROOT.TH2D(f"{fitName}_fail_{plot}", f"{fitName} fail - {plot}", nEtaBins, bins_eta, nPtBins, bins_pt)

        for k in file.GetListOfKeys():
            name = k.GetName()
            if "_resP" in name:
                flag = "P"
            elif "_resF" in name:
                flag = "F"
            else:
                continue
            obj = safeGetObject(file, name, detach=False)
            nbin = int(name.split("_")[0].lstrip("bin"))
            neta = int((nbin % nEtaBins) + 1)
            npt  = int((nbin / nEtaBins) + 1)

            if plot not in [m for m in dir(obj) if callable(getattr(obj, m, None))]:
                obj = obj.floatParsFinal()
                var = obj.find(plot+flag)
                if not isinstance(var, ROOT.RooRealVar):
                    print(f"WARNING: {plot} not found in the RooFitResult, skipping")
                    continue
                obj_plot = var.getVal()
            else:
                obj_plot = getattr(obj, plot)()

            if flag == "P":
                h2_pass.SetBinContent(neta, npt, obj_plot)
            else:
                h2_fail.SetBinContent(neta, npt, obj_plot)

        toPlot[plot] = (h2_pass, h2_fail)

    createPlotDirAndCopyPhp(outdir)

    canvas = ROOT.TCanvas("canvas","",1200, 900)

    hStatusPass, hStatusFail = toPlot["status"]
    hCovQualPass, hCovQualFail = toPlot["covQual"]
    hMeanPass, hMeanFail = toPlot["mean"]
    hSigmaPass, hSigmaFail = toPlot["sigma"]

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

