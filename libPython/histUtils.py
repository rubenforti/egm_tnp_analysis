#!/usr/bin/env python3

from array import array
from libPython.rootUtils import safeOpenFile, safeGetObject


def removeNegativeBins(h):
    for i in range(1,h.GetNbinsX()+1):
        if (h.GetBinContent(i) < 0):
            h.SetBinContent(i, 0)


def makePassFailHistograms(sample, bins, bindef, var ):

    probe_binning_eta, probe_binning_pt = bindef['eta']['bins'], bindef['pt']['bins']
    probe_var_eta, probe_var_pt         = bindef['eta']['var'] , bindef['pt']['var']

    probe_binning_pt  = array('d', probe_binning_pt)
    probe_binning_eta = array('d', probe_binning_eta)

    #binning_mass = array('d', [var['min'] + i*(var['max']-var['min'])/var['nbins'] for i in range(var['nbins']+1)])

    #print("sample.getInputPath() = ",sample.getInputPath()) 
    p = sample.getInputPath() 
    #print("p = ",p)
    infile = safeOpenFile(p, mode="READ")
    #print(infile.ls())
    h_tmp_pass = safeGetObject(infile, f"pass_{sample.getName()}", detach=False)
    h_tmp_fail = safeGetObject(infile, f"fail_{sample.getName()}", detach=False)
    
    # Passing probes evaluated with standalone variables, may be needed for tracking when using all probes to form failing MC template to fit data
    altPass = "pass_" + sample.getName() + "_alt"
    keyNames = [k.GetName() for k in infile.GetListOfKeys()]    
    h_tmp_pass_alt = safeGetObject(infile, altPass, detach=False) if altPass in keyNames else None
        
    outfile = safeOpenFile(sample.getOutputPath(), mode="RECREATE")
    
    for ii, ib in enumerate(bins):
        h_name = ib['name' ]
        h_title= ib['title']

        tmp_valpt_min  = ib['vars'][probe_var_pt ]['min']
        tmp_valeta_min = ib['vars'][probe_var_eta]['min']
        tmp_valpt_max  = ib['vars'][probe_var_pt ]['max']
        tmp_valeta_max = ib['vars'][probe_var_eta]['max']

        epsilon = 0.001 # safety thing when picking the bin edges using FindFixBin
        ibin_pt_low  = h_tmp_pass.GetYaxis().FindFixBin(tmp_valpt_min  + epsilon)
        ibin_eta_low = h_tmp_pass.GetZaxis().FindFixBin(tmp_valeta_min + epsilon)
        ibin_pt_high  = h_tmp_pass.GetYaxis().FindFixBin(tmp_valpt_max  - epsilon)
        ibin_eta_high = h_tmp_pass.GetZaxis().FindFixBin(tmp_valeta_max - epsilon)

        h_pass = h_tmp_pass.ProjectionX(h_name+'_Pass', ibin_pt_low, ibin_pt_high, ibin_eta_low, ibin_eta_high)
        h_pass.SetTitle(h_title+' passing')
        removeNegativeBins(h_pass)
        h_pass.Write(h_pass.GetName())

        h_fail = h_tmp_fail.ProjectionX(h_name+'_Fail', ibin_pt_low, ibin_pt_high, ibin_eta_low, ibin_eta_high)
        h_fail.SetTitle(h_title+' failing')
        removeNegativeBins(h_fail)
        h_fail.Write(h_fail.GetName())
        
        if h_tmp_pass_alt:
            h_pass_alt = h_tmp_pass_alt.ProjectionX(h_name+'_Pass_alt', ibin_pt_low, ibin_pt_high, ibin_eta_low, ibin_eta_high)
            h_pass_alt.SetTitle(h_title+' passing alternate')
            removeNegativeBins(h_pass_alt)
            h_pass_alt.Write(h_pass_alt.GetName())
            
    outfile.Close()