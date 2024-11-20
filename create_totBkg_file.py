import ROOT
import sys
import argparse

lumi_data = {"2016" : 16.8, # only postVFP
             "2017" : 37.99, # this includes the LS where the HLT_isoMu24 was not prescaled
             "2018" : 59.81}

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

working_points_global = {
    ## for global muons
    ##
    #'mu_reco_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_reco_mc_vertexWeights1_oscharge1.root',
    #                 '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_reco_data_vertexWeights1_oscharge1.root'],
    'mu_reco_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recoplus_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recoplus_data_vertexWeights1_oscharge1.root'],
    'mu_reco_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recominus_mc_vertexWeights1_oscharge1.root',
                      '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recominus_data_vertexWeights1_oscharge1.root'],
    # 'mu_reco_minusodd': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recominusodd_mc_vertexWeights1_oscharge1.root',
    #                  '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recominusodd_data_vertexWeights1_oscharge1.root'],
    # 'mu_reco_minuseven': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recominuseven_mc_vertexWeights1_oscharge1.root',
    #                  '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_recominuseven_data_vertexWeights1_oscharge1.root'],
    'mu_tracking_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_tracking_mc_vertexWeights1_oscharge0.root',
                         '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_tracking_data_vertexWeights1_oscharge0.root'],
    'mu_tracking_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_tracking_mc_vertexWeights1_oscharge0_tagChargeMinus.root',
                             '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_tracking_data_vertexWeights1_oscharge0_tagChargeMinus.root'],
    'mu_tracking_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_tracking_mc_vertexWeights1_oscharge0_tagChargePlus.root',
                            '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_tracking_data_vertexWeights1_oscharge0_tagChargePlus.root'],
    # 'mu_tracking_odd': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_trackingodd_mc_vertexWeights1_oscharge0.root',
    #                      '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_trackingodd_data_vertexWeights1_oscharge0.root'],
    # 'mu_tracking_even': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_trackingeven_mc_vertexWeights1_oscharge0.root',
    #                      '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_trackingeven_data_vertexWeights1_oscharge0.root'],
    'mu_idip_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_idip_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_idip_data_vertexWeights1_oscharge1.root'],
    'mu_idip_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_idipplus_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_idipplus_data_vertexWeights1_oscharge1.root'],
    'mu_idip_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_idipminus_mc_vertexWeights1_oscharge1.root',
                      '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_idipminus_data_vertexWeights1_oscharge1.root'],
    'mu_trigger_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_triggerplus_mc_vertexWeights1_oscharge1.root',
                        '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_triggerplus_data_vertexWeights1_oscharge1.root'],
    'mu_trigger_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_triggerminus_mc_vertexWeights1_oscharge1.root',
                         '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_triggerminus_data_vertexWeights1_oscharge1.root'],
    'mu_iso_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_iso_mc_vertexWeights1_oscharge1.root',
                    '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_iso_data_vertexWeights1_oscharge1.root'],
    # 'mu_iso_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_isoplus_mc_vertexWeights1_oscharge1.root',
    #                 '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_isoplus_data_vertexWeights1_oscharge1.root'],
    # 'mu_iso_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_isominus_mc_vertexWeights1_oscharge1.root',
    #                 '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_isominus_data_vertexWeights1_oscharge1.root'],
    'mu_isonotrig_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_isonotrig_mc_vertexWeights1_oscharge1.root',
                          '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_isonotrig_data_vertexWeights1_oscharge1.root'],
    'mu_veto_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_veto_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_highPurity_XYZ_1orMoreValidHitsSA/tnp_veto_data_vertexWeights1_oscharge1.root'],
}

working_points_2017 = {
    'mu_reco_plus': ['/scratch/rforti/steve_histograms_2017/tnp_recoplus_mc_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2017/tnp_recoplus_data_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2017/tnp_recoplus_bkg_vertexWeights1_oscharge1.root'],
    'mu_reco_minus': ['/scratch/rforti/steve_histograms_2017/tnp_recominus_mc_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2017/tnp_recominus_data_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2017/tnp_recominus_bkg_vertexWeights1_oscharge1.root'],
    'mu_tracking_plus': ['/scratch/rforti/steve_histograms_2017/tnp_trackingplus_mc_vertexWeights1_oscharge0.root',
                         '/scratch/rforti/steve_histograms_2017/tnp_trackingplus_data_vertexWeights1_oscharge0.root',
                         '/scratch/rforti/steve_histograms_2017/tnp_trackingplus_bkg_vertexWeights1_oscharge0.root'],
    'mu_tracking_minus': ['/scratch/rforti/steve_histograms_2017/tnp_trackingminus_mc_vertexWeights1_oscharge0.root',
                          '/scratch/rforti/steve_histograms_2017/tnp_trackingminus_data_vertexWeights1_oscharge0.root',
                          '/scratch/rforti/steve_histograms_2017/tnp_trackingminus_bkg_vertexWeights1_oscharge0.root'],
    'mu_idip_plus': ['/scratch/rforti/steve_histograms_2017/tnp_idipplus_mc_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2017/tnp_idipplus_data_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2017/tnp_idipplus_bkg_vertexWeights1_oscharge1.root'],
    'mu_idip_minus': ['/scratch/rforti/steve_histograms_2017/tnp_idipminus_mc_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2017/tnp_idipminus_data_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2017/tnp_idipminus_bkg_vertexWeights1_oscharge1.root'],
    'mu_trigger_plus': ['/scratch/rforti/steve_histograms_2017/tnp_triggerplus_mc_vertexWeights1_oscharge1.root',
                        '/scratch/rforti/steve_histograms_2017/tnp_triggerplus_data_vertexWeights1_oscharge1.root',
                        '/scratch/rforti/steve_histograms_2017/tnp_triggerplus_bkg_vertexWeights1_oscharge1.root'],
    'mu_trigger_minus': ['/scratch/rforti/steve_histograms_2017/tnp_triggerminus_mc_vertexWeights1_oscharge1.root',
                         '/scratch/rforti/steve_histograms_2017/tnp_triggerminus_data_vertexWeights1_oscharge1.root',
                         '/scratch/rforti/steve_histograms_2017/tnp_triggerminus_bkg_vertexWeights1_oscharge1.root'],
    'mu_iso_both': ['/scratch/rforti/steve_histograms_2017/tnp_iso_mc_vertexWeights1_oscharge1.root',
                    '/scratch/rforti/steve_histograms_2017/tnp_iso_data_vertexWeights1_oscharge1.root',
                    '/scratch/rforti/steve_histograms_2017/tnp_iso_bkg_vertexWeights1_oscharge1.root'],
}


working_points_2018 = {

    'mu_reco_plus': ['/scratch/rforti/steve_histograms_2018/tnp_recoplus_mc_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2018/tnp_recoplus_data_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2018/tnp_recoplus_bkg_vertexWeights1_oscharge1.root'],
    'mu_reco_minus': ['/scratch/rforti/steve_histograms_2018/tnp_recominus_mc_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2018/tnp_recominus_data_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2018/tnp_recominus_bkg_vertexWeights1_oscharge1.root'],
    'mu_tracking_plus': ['/scratch/rforti/steve_histograms_2018/tnp_trackingplus_mc_vertexWeights1_oscharge0.root',
                         '/scratch/rforti/steve_histograms_2018/tnp_trackingplus_data_vertexWeights1_oscharge0.root',
                         '/scratch/rforti/steve_histograms_2018/tnp_trackingplus_bkg_vertexWeights1_oscharge0.root'],
    'mu_tracking_minus': ['/scratch/rforti/steve_histograms_2018/tnp_trackingminus_mc_vertexWeights1_oscharge0.root',
                          '/scratch/rforti/steve_histograms_2018/tnp_trackingminus_data_vertexWeights1_oscharge0.root',
                          '/scratch/rforti/steve_histograms_2018/tnp_trackingminus_bkg_vertexWeights1_oscharge0.root'],
    'mu_idip_plus': ['/scratch/rforti/steve_histograms_2018/tnp_idipplus_mc_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2018/tnp_idipplus_data_vertexWeights1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2018/tnp_idipplus_bkg_vertexWeights1_oscharge1.root'],
    'mu_idip_minus': ['/scratch/rforti/steve_histograms_2018/tnp_idipminus_mc_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2018/tnp_idipminus_data_vertexWeights1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2018/tnp_idipminus_bkg_vertexWeights1_oscharge1.root'],
    'mu_trigger_plus': ['/scratch/rforti/steve_histograms_2018/tnp_triggerplus_mc_vertexWeights1_oscharge1.root',
                        '/scratch/rforti/steve_histograms_2018/tnp_triggerplus_data_vertexWeights1_oscharge1.root',
                        '/scratch/rforti/steve_histograms_2018/tnp_triggerplus_bkg_vertexWeights1_oscharge1.root'],
    'mu_trigger_minus': ['/scratch/rforti/steve_histograms_2018/tnp_triggerminus_mc_vertexWeights1_oscharge1.root',
                         '/scratch/rforti/steve_histograms_2018/tnp_triggerminus_data_vertexWeights1_oscharge1.root',
                         '/scratch/rforti/steve_histograms_2018/tnp_triggerminus_bkg_vertexWeights1_oscharge1.root'],
    'mu_iso_both': ['/scratch/rforti/steve_histograms_2018/tnp_iso_mc_vertexWeights1_oscharge1.root',
                    '/scratch/rforti/steve_histograms_2018/tnp_iso_data_vertexWeights1_oscharge1.root',
                    '/scratch/rforti/steve_histograms_2018/tnp_iso_bkg_vertexWeights1_oscharge1.root'],
}


working_points_new = {
    'mu_reco_plus': ['/scratch/rforti/steve_histograms_2016/reco/tnp_recoplus_mc_vertexWeights1_genMatching1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2016/reco/tnp_recoplus_data_vertexWeights1_genMatching1_oscharge1.root',
                     '/scratch/rforti/steve_histograms_2016/reco/tnp_recoplus_bkg_vertexWeights1_genMatching0_oscharge1.root'],
    'mu_reco_minus': ['/scratch/rforti/steve_histograms_2016/reco/tnp_recominus_mc_vertexWeights1_genMatching1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2016/reco/tnp_recominus_data_vertexWeights1_genMatching1_oscharge1.root',
                      '/scratch/rforti/steve_histograms_2016/reco/tnp_recominus_bkg_vertexWeights1_genMatching0_oscharge1.root'],
    'mu_tracking_plus': ['/scratch/rforti/steve_histograms_2016/tracking/tnp_trackingplus_mc_vertexWeights1_genMatching1_oscharge0.root',
                         '/scratch/rforti/steve_histograms_2016/tracking/tnp_trackingplus_data_vertexWeights1_genMatching1_oscharge0.root',
                         '/scratch/rforti/steve_histograms_2016/tracking/tnp_trackingplus_bkg_vertexWeights1_genMatching0_oscharge0.root'],
    'mu_tracking_minus': ['/scratch/rforti/steve_histograms_2016/tracking/tnp_trackingminus_mc_vertexWeights1_genMatching1_oscharge0.root',
                          '/scratch/rforti/steve_histograms_2016/tracking/tnp_trackingminus_data_vertexWeights1_genMatching1_oscharge0.root',
                          '/scratch/rforti/steve_histograms_2016/tracking/tnp_trackingminus_bkg_vertexWeights1_genMatching0_oscharge0.root'],
}

bkg_types = ["WW", "WZ", "ZZ", "TTSemileptonic", "TTFullyleptonic", "Ztautau", "WplusJets", "WminusJets", 
             #"QCD", 
             "Zjets"]

def lumi_factor(filepath, process, yr="2016"):
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
    scale = lumi_data[yr]/lumi_process
    return scale


def checkAddConsistency(sum, h_list):
    for i in range(1, (sum.GetNbinsX()*sum.GetNbinsY()*sum.GetNbinsZ())+1):
        external_sum = 0
        for h in h_list:
            external_sum += h.GetBinContent(i)
        internal_sum = sum.GetBinContent(i)
        if external_sum != internal_sum:
            print("INCORRECT SUM FOR BIN", i)


parser = argparse.ArgumentParser()
parser.add_argument('-s','--steps', default=None, nargs='*', type=str, choices=list([x.split("_")[1] for x in working_points_global.keys()]),
                    help='Default runs all working points, but can choose only some if needed')
parser.add_argument('-x','--exclude', default=None, nargs='*', type=str, choices=list([x.split("_")[1] for x in working_points_global.keys()]),
                    help='Default runs all working points, but can exclude some if needed')
parser.add_argument("-y", "--year", type=str, choices=["2016", "2017", "2018"], help="Year of data taking")
parser.add_argument('--oldSteve', action='store_true',
                    help='Use input files from Steve that have the "old" naming scheme (i.e. with the "genMatching" flag)')
parser.add_argument('--scale_bkg_by_lumi', action='store_true',
                    help='Scale background datasets to match the luminosity of data')
                    
args = parser.parse_args()
                    
          
if args.oldSteve is True: 
    bkg_types.remove("Zjets") 
    bkg_types += ["mc"]
    datasets = working_points_new
    args.scale_bkg_by_lumi = True
elif args.year=="2017":
    datasets = working_points_2017
elif args.year=="2018":
    datasets = working_points_2018
else:
    datasets = working_points_global
    
 
if args.steps is None:
    args.steps = list(set([k.split("_")[1] for k in working_points_global.keys()]))
    
for eff in args.steps:
    for ch in ["plus", "minus", "both"]:
    
        if not f"mu_{eff}_{ch}" in datasets.keys(): 
            continue
        else:
            base_folder = datasets[f"mu_{eff}_{ch}"][0].replace(datasets[f"mu_{eff}_{ch}"][0].split("/")[-1], "")
            
        os = 0 if eff=="tracking" else 1
        
        print(eff, ch, base_folder, os)
        
        ch_file = ch if ch!="both" else ""
        
        fnames = [base_folder+f"tnp_{eff}{ch_file}_{bkg_type}_vertexWeights1_genMatching0_oscharge{os}.root" for bkg_type in bkg_types]
        
        if args.oldSteve:
            fnames = [fname.replace("genMatching0", "genMatching-1") if "_mc_" in fname else fname for fname in fnames]
            gm_info = "_genMatching0"
        else:
            fnames = [fname.replace("_genMatching0", "") for fname in fnames]
            gm_info = ""

        if len(datasets[f"mu_{eff}_{ch}"])==3:    
            file_out = ROOT.TFile(datasets[f"mu_{eff}_{ch}"][2], "RECREATE")
        else:
            file_out = ROOT.TFile(f"{base_folder}tnp_{eff}{ch_file}_bkg_vertexWeights1{gm_info}_oscharge{os}.root", "RECREATE")




        for fl in ["pass", "fail"]:

            rootfiles = []
            histos = []

            for fname in fnames: rootfiles.append(ROOT.TFile(fname, "READ"))       
            
            for rf in rootfiles:
            
                histos.append(rf.Get(f"{fl}_mu_DY_postVFP"))
                
                bkg_process = rf.GetName().split("/")[-1].split("_")[2]
                if args.oldSteve and bkg_process == "mc": bkg_process = "Zjets"
                
                if args.scale_bkg_by_lumi:
                    scaling = lumi_factor(rf.GetName(), bkg_process, yr=args.year)
                    print(f"Scaling {bkg_process} by {scaling}")
                    histos[-1].Scale(scaling)
                    

            hist_sum = histos[0].Clone(f"{fl}_mu_DY_postVFP")
            hist_sum.Reset()
            hist_sum.SetName(f"{fl}_mu_mcBkg_postVFP"), hist_sum.SetTitle(f"{fl}_mu_mcBkg_postVFP")
            for h in histos:
                hist_sum.Add(h)

            checkAddConsistency(hist_sum, histos)

            file_out.cd()
            hist_sum.Write()

        file_out.Close()

