#!/usr/bin/env python3 

import ROOT, math, array, ctypes, copy, sys
import numpy as np
import os.path
from . import fitUtils
from .plotUtils import safeOpenFile, safeGetObject
import functools

def removeNegativeBins(h):
    for i in range(1,h.GetNbinsX()+1):
        if (h.GetBinContent(i) < 0):
            h.SetBinContent(i, 0)


def makePassFailHistograms(sample, bins, bindef, var ):

    probe_binning_eta, probe_binning_pt = bindef['eta']['bins'], bindef['pt']['bins']
    probe_var_eta, probe_var_pt         = bindef['eta']['var'] , bindef['pt']['var']

    probe_binning_pt  = array.array('d', probe_binning_pt)
    probe_binning_eta = array.array('d', probe_binning_eta)

    #binning_mass = array.array('d', [var['min'] + i*(var['max']-var['min'])/var['nbins'] for i in range(var['nbins']+1)])

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


def histPlotter(rootfile, tnpBin, plotDir, replica=-1, verbosePlotting=True ):

    if verbosePlotting:
        binName = tnpBin['name']
        if replica<0:
            c = safeGetObject(rootfile, f"{binName}_Canv", detach=False)
            c.SaveAs(f"{plotDir}/{binName}.png")   
            # c.Print( '%s/%s.pdf' % (plotDir,tnpBin['name']))  ## this doesn't work properly, it only saves the pad with the text
        else:
            c = rootfile.Get( '%s_Canv_Stat%d' % (tnpBin['name'],replica) )
            c.Print( '%s/%s_Stat%d.png' % (plotDir,tnpBin['name'],replica))
            c.Print( '%s/%s_Stat%d.pdf' % (plotDir,tnpBin['name'],replica))
        


def computeEffi( n1,n2,e1,e2):
    effout = []
    if (n1+n2):
        eff   = n1/(n1+n2)
        nTot = n1+n2
        e_eff = math.sqrt(e1*e1*n2*n2+e2*e2*n1*n1) / (nTot*nTot)
        #if e_eff < 0.001 : e_eff = 0.001
    else:
        eff, e_eff = 1.1, 0.01

    effout.append(eff)
    effout.append(e_eff)
    
    return effout


def getAllEffi(info, bindef, outputDirectory, saveCanvas=False):
    # FIXME: this function often leads to crashes,
    #        probably because of reading canvases but it is not clear
    effis = {}
    effis_canvas = {}
    binName = bindef["name"]
    canvasesToGet = ["Data_Nominal", "MC_Nominal_fit"] #, "Data_Alt_Sig", "Data_Alt_Bkg"] 
    #canvasesToGet = ["dataNominal", "dataAltSig", "dataAltBkg", "mcNominal"]
    for key, value in info.items():
        effis[key], effis_canvas[f"canv_{key}"] = [-1, -1], None
        if value is None or not os.path.isfile(value):
            continue

        rootfile = safeOpenFile(value, mode="READ")

        if key in canvasesToGet and saveCanvas:
                effis_canvas[f"canv_{key}"] = safeGetObject(rootfile, f"{binName}_Canv", detach=False)
        
        if key == "MC_Nominal":
            hP = safeGetObject(rootfile, f"{binName}_Pass", detach=False)
            hF = safeGetObject(rootfile, f"{binName}_Fail", detach=False)
            bin1, bin2 = 1, hP.GetXaxis().GetNbins()
            ePc, eFc = ctypes.c_double(-1.0), ctypes.c_double(-1.0)
            nP = hP.IntegralAndError(bin1, bin2, ePc)
            nF = hF.IntegralAndError(bin1, bin2, eFc)
            eP, eF = float(ePc.value), float(eFc.value)
        else:
            fitresP = safeGetObject(rootfile, f"{binName}_resP", detach=False)
            fitresF = safeGetObject(rootfile, f"{binName}_resF", detach=False)
            
            fitP, fitF = fitresP.floatParsFinal().find("nSigP"), fitresF.floatParsFinal().find("nSigF")
            nP, nF = fitP.getVal(), fitF.getVal()
            eP, eF = fitP.getError(), fitF.getError()
            if "Data" in key:
                eP, eF = max(math.sqrt(nP), eP), max(math.sqrt(nF), eF)
        
        effis[key] = computeEffi(nP, nF, eP, eF) + [nP, nF, eP, eF]
        rootfile.Close()
    ##
    ##
    ## Saving canvases here, let's test if this also leads to crashes
    ##
    # canvases = [f"canv_{x}" for x in canvasesToGet]
    # padsFromCanvas = {}
    # for c in canvases:
    #     if c not in effis_canvas.keys() or effis_canvas[c] == None:                
    #         print(f"Canvas {c} not found or not available")
    #         return 0
    #     else:
    #         #padsFromCanvas[c] = list(effis[c])
    #         ## the following was for when canvases where returned
    #         if effis_canvas[c].ClassName() ==  "TCanvas":
    #             padsFromCanvas[c] = [p for p in effis_canvas[c].GetListOfPrimitives()]
    #             # print(padsFromCanvas[c])
    #             for p in padsFromCanvas[c]:
    #                 ROOT.SetOwnership(p, False)
    #             #effis_canvas[c] = None
    #             #ROOT.SetOwnership(effis_canvas[c], False)
    #         else:
    #             print(f"SOMETHING SCREWED UP WITH TCANVAS for bin {bindef['name']}")
    #             return 0

    if not saveCanvas:
        return effis

    ncols = len(canvasesToGet)

    canv_all = ROOT.TCanvas(bindef['name'], bindef['name'], 400*ncols, 1200)
    canv_all.Draw()
    pad_title = ROOT.TPad('title', 'title', 0.0, 0.95, 1.0, 1.0) 
    pad_pass = ROOT.TPad('passing', 'passing', 0.0, 0.55, 1.0, 0.95)
    pad_fail = ROOT.TPad('failing', 'failing', 0.0, 0.15, 1.0, 0.55)
    pad_eff = ROOT.TPad('efficiency', 'efficiency', 0.0, 0.0, 1.0, 0.15)

    pad_title.Draw()
    pad_pass.Divide(ncols,1), pad_pass.Draw()
    pad_fail.Divide(ncols,1), pad_fail.Draw()
    pad_eff.Divide(ncols,1), pad_eff.Draw()

    pad_title.cd()
    txt = ROOT.TLatex()
    txt.SetTextFont(62)
    txt.SetTextSize(0.5)
    txt.SetNDC()
    txt.DrawLatex(0.1, 0.5, f'{bindef['name'].replace('_',' ').replace('To', '-').replace('probe ', '').replace('m','-').replace('pt','XX').replace('p','.').replace('XX','p_{T}')}')
    txt.SetTextSize(0.15)

    for icol, canv_type in enumerate(canvasesToGet):

        for ip, p in enumerate(effis_canvas[f"canv_{canv_type}"].GetListOfPrimitives()):
            if not ip: continue
            pad_to_use = pad_pass if ip==1 else pad_fail if ip==2 else None  # TO CHANGE WHEN CREATING THE CANVASES !!!!!
            pad_to_use.cd(icol+1)
            newp = p.Clone(f"tmp_{canv_type}_{ip}")
            newp.SetPad(0.02, 0.00, 0.98, 0.98)
            newp.Draw()
        pad_eff.cd(icol+1)
        txt.SetTextFont(62)
        txt.DrawLatex(0.10, 0.77, canv_type.replace('_', " ") + " :" if not "MC_Nominal" in canv_type else "MC counting :")
        txt.SetTextFont(42)
        tmp = effis[canv_type]
        txt.DrawLatex(0.20, 0.59, f'Passing: {tmp[2]:.1f} #pm {tmp[4]:.1f}')
        txt.DrawLatex(0.20, 0.41, f'Failing: {tmp[3]:.1f} #pm {tmp[5]:.1f}')
        txt.SetTextFont(62)
        txt.DrawLatex(0.20, 0.23, f'Efficiency: {tmp[0]*100.:.2f} #pm {tmp[1]*100.:.2f} %')
        txt.SetTextFont(42)

    canv_all.cd()
    canv_all.Update()
    canv_all.SaveAs(outputDirectory+'/plots/{n}_all.pdf'.format(n=bindef['name']))
    canv_all.SaveAs(outputDirectory+'/plots/{n}_all.png'.format(n=bindef['name']))

    return effis



def getAllScales( info, bindef, refReplica ):
    scales = {}

    for key,rfile in info.iteritems():
        if not info[key] is None and os.path.isfile(rfile) :
            rootfile = ROOT.TFile( rfile, 'read' )
            replica = int(rfile.split('_Stat')[-1].split('.root')[0]) if 'dataReplica' in key else refReplica
            fitresP = rootfile.Get( '%s_resP_Stat%d' % (bindef['name'],replica)  )
            
            fitMean = fitresP.floatParsFinal().find('meanP')
            v = fitMean.getVal()
            e = fitMean.getError()
            rootfile.Close()
        
            scales[key] = [v,e]
        else:
            scales[key] = [-999,-999]
    return scales


