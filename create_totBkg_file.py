import ROOT

lumi_data = 16.8  # fb^-1, for 2016postVFP

BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)
xsec_ZmmPostVFP = 2001.9

cross_sections_bkg = {
    # Unit = pb
    "QCD" : 238800,
    "WW" : 12.6,
    "WZ" : 27.59,  #5.4341,
    "ZZ" : 0.60,
    "TTFullyleptonic" : 88.29,
    "TTSemileptonic" : 366.34,
    "WplusJets" : 11765.9,
    "WminusJets" : 8703.87,
    "Ztautau" : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO,
    "Zjets" : xsec_ZmmPostVFP  # obtained from the signal MC file with the "reverse" gen-matching option
}

def lumi_factor(filepath, process):
    """
    Returns the lumi factor for the process in the given file
    """
    file = ROOT.TFile(filepath)
    wsum_histo = file.Get("weightSum")
    num_init = wsum_histo.Integral()

    if "mc" in filepath:
        xsection = xsec_ZmmPostVFP*1000 # has to be put in fb
    else:
        xsection = cross_sections_bkg[process]*1000
        
    lumi_process = num_init/xsection

    scale = lumi_data/lumi_process

    return scale


def checkAddConsistency(sum, h_list):
    for i in range(1, (sum.GetNbinsX()*sum.GetNbinsY()*sum.GetNbinsZ())+1):
        external_sum = 0
        for h in h_list:
            external_sum += h.GetBinContent(i)
        internal_sum = sum.GetBinContent(i)
        if external_sum != internal_sum:
            print("UNCORRECT SUM FOR BIN", i)




if __name__ == "__main__":


    old_Steve_version = True

    eff_type = "recoplus"

    base_folder = f"/home/rforti/egm_tnp_analysis/steve_histograms_2016/{eff_type.replace('plus', '').replace('minus', '')}/"

    bkg_types = ["WW", "WZ", "ZZ", "TTSemileptonic", "TTFullyleptonic", "Ztautau", "WplusJets", "WminusJets", 
                 #"QCD", 
                 "Zjets"]

    if old_Steve_version:
        bkg_types.remove("Zjets")
        bkg_types += ["mc"]

    files = [base_folder+f"tnp_{eff_type}_{bkg_type}_vertexWeights1_genMatching0_oscharge1.root" for bkg_type in bkg_types]


    if old_Steve_version:
        files = [file.replace("genMatching0", "genMatching-1") if "_mc_" in file else file for file in files]
    else:
        files = [file.replace("_genMatching0_", "") for file in files]


    file_out = ROOT.TFile(f"{base_folder}tnp_{eff_type}_bkg_vertexWeights1_genMatching0_oscharge1.root", "recreate")



    for fl in ["pass", "fail"]:
        
        rootfiles = []
        histos = []

        for file in files:
            # print("Opening file", file)
            rootfiles.append(ROOT.TFile(file, "read"))
        
        for rf in rootfiles:
            histos.append(rf.Get(f"{fl}_mu_DY_postVFP"))
            if old_Steve_version:
                # print(rf.GetName())
                bkg_process = rf.GetName().split("/")[-1].split("_")[2]
                if bkg_process == "mc": bkg_process = "Zjets"
                histos[-1].Scale(lumi_factor(rf.GetName(), bkg_process))
                # print(f"Scaling {bkg_process} by {lumi_factor(rf.GetName(), bkg_process)}")

        hist_sum = histos[0].Clone(f"{fl}_mu_DY_postVFP")
        print(hist_sum.Integral())
        hist_sum.Reset()
        print(hist_sum.Integral())

        for h in histos:
            hist_sum.Add(h)

        checkAddConsistency(hist_sum, histos)


        file_out.cd()
        hist_sum.Write()

        external_sum = 0
        for h in histos:
            external_sum += h.Integral()

        print(hist_sum.Integral(), external_sum)

        

    file_out.Close()

