#!/usr/bin/env python3 

import ROOT, math, array, ctypes, copy
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

    probe_var_pass_pt, probe_var_pass_eta = probe_var_pt, probe_var_eta

    #binning_mass = array.array('d', [var['min'] + i*(var['max']-var['min'])/var['nbins'] for i in range(var['nbins']+1)])

    #print("sample.getInputPath() = ",sample.getInputPath()) 
    p = sample.getInputPath() 
    #print("p = ",p)
    infile = safeOpenFile(p, mode="read")
    #print(infile.ls())
    h_tmp_pass = safeGetObject(infile, f"pass_{sample.getName()}", detach=False)
    h_tmp_fail = safeGetObject(infile, f"fail_{sample.getName()}", detach=False)
    altPass = "pass_" + sample.getName() + "_alt"
    keyNames = [k.GetName() for k in infile.GetListOfKeys()]
    h_tmp_pass_alt = None # will be needed for tracking when using all probes to form failing probe MC template to fit data, in this case the alternate version of passing probe made with  standalone variables is needed 
    if altPass in keyNames:
        h_tmp_pass_alt = safeGetObject(infile, altPass, detach=False)
        
    outfile = safeOpenFile(sample.getOutputPath(), mode="recreate")
    
    #print('this is the integral of the passing', h_tmp_pass.Integral())
    #print('this is the integral of the failing', h_tmp_fail.Integral())
    
    for ii,ib in enumerate(bins):
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
        h_fail = h_tmp_fail.ProjectionX(h_name+'_Fail', ibin_pt_low, ibin_pt_high, ibin_eta_low, ibin_eta_high)
        h_pass.SetTitle(h_title+' passing')
        h_fail.SetTitle(h_title+' failing')
        removeNegativeBins(h_pass)
        removeNegativeBins(h_fail)
        h_pass.Write(h_pass.GetName())
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
    #canvasesToGet = ["dataNominal", "dataAltSig", "dataAltBkg", "mcAltSig"]
    canvasesToGet = ["dataNominal", "dataAltSig", "dataAltBkg"]
    for key in info.keys():
        value = info[key]
        if value is None or not os.path.isfile(value):
            effis[key] = [-1,-1]
            effis_canvas[f"canv_{key}"] = None
        else:
            rootfile = safeOpenFile(value, mode="read")
            # some customization
            if key == "mcNominal":
                # just counting using integral
                hP = safeGetObject(rootfile, f"{binName}_Pass", detach=False)
                hF = safeGetObject(rootfile, f"{binName}_Fail", detach=False)
                bin1 = 1
                bin2 = hP.GetXaxis().GetNbins()
                ePc = ctypes.c_double(-1.0)
                eFc = ctypes.c_double(-1.0)
                nP = hP.IntegralAndError(bin1,bin2,ePc)
                nF = hF.IntegralAndError(bin1,bin2,eFc)
                eP = float(ePc.value)
                eF = float(eFc.value)
            else:
                fitresP = safeGetObject(rootfile, f"{binName}_resP", detach=False)
                fitresF = safeGetObject(rootfile, f"{binName}_resF", detach=False)
                if key in canvasesToGet and saveCanvas:
                    effis_canvas[f"canv_{key}"] = safeGetObject(rootfile, f"{binName}_Canv", detach=False)
                    #effis_canvas[f"canv_{key}"] = copy.deepcopy(canv.Clone(f"{key}_{binName}"))
                #rpPass = safeGetObject(rootfile, f"{binName}_rooplotP", detach=False)
                #rpFail = safeGetObject(rootfile, f"{binName}_rooplotF", detach=False)
                #effis[f"rpPass_{key}"] = copy.deepcopy(rpPass.Clone(f"rpP_{key}_{binName}"))
                #effis[f"rpFail_{key}"] = copy.deepcopy(rpFail.Clone(f"rpF_{key}_{binName}"))
                fitP = fitresP.floatParsFinal().find("nSigP")
                fitF = fitresF.floatParsFinal().find("nSigF")
                nP = fitP.getVal()
                nF = fitF.getVal()
                eP = fitP.getError()
                eF = fitF.getError()
                if key in ["dataNominal", "dataAltSig", "dataAltBkg"]:
                    # if the fit is reliable the uncertainty in nSignal for data must be at least equal to sqrt(nSignal)
                    # it can be larger because of the presence of parameters and backgrounds
                    eP = max(math.sqrt(nP), eP)
                    eF = max(math.sqrt(nF), eF)
            # done with all cases
            effis[key] = computeEffi(nP,nF,eP,eF) +[nP, nF, eP, eF]
            rootfile.Close() # try also not closing the file
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

    canv_all = ROOT.TCanvas(bindef['name'], bindef['name'], 1200, 1200)
    canv_all.Divide(3,3)
    canv_all.Draw()
    canv_all.cd(0)
    txt = ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.03)
    txt.SetNDC()
    txt.DrawLatex(0.01, 0.97, '{n}'.format(n=bindef['name'].replace('_',' ').replace('To', '-').replace('probe ', '').replace('m','-').replace('pt','XX').replace('p','.').replace('XX','p_{T}')))
    txt.SetTextSize(0.08)
    ipad = 1
    for ip, p in enumerate(effis_canvas["canv_dataNominal"].GetListOfPrimitives()):  # Draw the data_nominal, with info on the MC counting efficiency
        if not ip: continue
        canv_all.cd(ipad)
        newp = p.Clone(f"tmp_dataNominal_{ip}")
        newp.SetPad(0.05, 0.00, 0.95, 0.90)
        newp.Draw()
        ipad += 1
    canv_all.cd(ipad)
    tmp = effis['dataNominal']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.85, 'data nominal:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.75, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.64, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.53, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    tmp = effis['mcNominal']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.35, 'MC counting efficiency:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.24, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.13, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.02, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42) 
    ipad += 1
    for ip, p in enumerate(effis_canvas["canv_dataAltSig"].GetListOfPrimitives()):   # Draw the data_alt_signal
        if not ip: continue
        canv_all.cd(ipad)
        newp = p.Clone(f"tmp_dataAltSig_{ip}")
        newp.SetPad(0.05, 0.00, 0.95, 0.90)
        newp.Draw()
        ipad+=1
    canv_all.cd(ipad)
    tmp = effis['dataAltSig']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.65, 'data alt. signal:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.54, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.43, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.32, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    ipad+=1
    for ip, p in enumerate(effis_canvas["canv_dataAltBkg"].GetListOfPrimitives()):   # Draw the data_alt_background
        if not ip: continue
        canv_all.cd(ipad)
        newp = p.Clone(f"tmp_dataAltBkg_{ip}")
        newp.SetPad(0.05, 0.00, 0.95, 0.90)
        newp.Draw()
        ipad += 1
    canv_all.cd(ipad)
    tmp = effis['dataAltBkg']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.65, 'data alt. background:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.54, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.43, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.32, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    canv_all.SaveAs(outputDirectory+'/plots/{n}_all.pdf'.format(n=bindef['name']))
    canv_all.SaveAs(outputDirectory+'/plots/{n}_all.png'.format(n=bindef['name']))
    
    '''
    ## ----------------------------------------------------------------------------
    ## OLD VERSION: the first canvases are of the fit on MC, that is no longer used
    ## ----------------------------------------------------------------------------          

    ipad = 1
    # for ip, p in enumerate(padsFromCanvas["canv_mcAltSig"]):
    for ip, p in enumerate(effis_canvas["canv_mcAltSig"].GetListOfPrimitives()):
        if not ip: continue
        canv_all.cd(ipad)
        newp = p.Clone(f"tmp_mcAltSig_{ip}")
        newp.SetPad(0.05, 0.00, 0.95, 0.90)
        newp.Draw()
        ipad += 1
    canv_all.cd(ipad)
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.85, 'MC counting efficiency:')
    txt.SetTextFont(42)
    tmp = effis['mcNominal']
    txt.DrawLatex(0.10, 0.75, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.64, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.53, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    tmp = effis['mcAltSig']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.35, 'MC fitted signal:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.24, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.13, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.02, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    ipad+=1
    #for ip, p in enumerate(padsFromCanvas["canv_dataNominal"]):
    for ip, p in enumerate(effis_canvas["canv_dataNominal"].GetListOfPrimitives()):
        if not ip: continue
        canv_all.cd(ipad)
        newp = p.Clone(f"tmp_dataNominal_{ip}")
        newp.SetPad(0.05, 0.00, 0.95, 0.90)
        newp.Draw()
        ipad += 1
    canv_all.cd(ipad)
    tmp = effis['dataNominal']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.65, 'data nominal:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.54, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.43, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.32, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    ipad += 1
    #for ip, p in enumerate(padsFromCanvas["canv_dataAltSig"]):
    for ip, p in enumerate(effis_canvas["canv_dataAltSig"].GetListOfPrimitives()):
        if not ip: continue
        canv_all.cd(ipad)
        newp = p.Clone(f"tmp_dataAltSig_{ip}")
        newp.SetPad(0.05, 0.00, 0.95, 0.90)
        newp.Draw()
        ipad+=1
    canv_all.cd(ipad)
    tmp = effis['dataAltSig']
    txt.SetTextFont(62)
    txt.DrawLatex(0.00, 0.65, 'data alternative:')
    txt.SetTextFont(42)
    txt.DrawLatex(0.10, 0.54, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
    txt.DrawLatex(0.10, 0.43, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
    txt.SetTextFont(62)
    txt.DrawLatex(0.10, 0.32, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
    txt.SetTextFont(42)
    canv_all.SaveAs(outputDirectory+'/plots/{n}_all.pdf'.format(n=bindef['name']))
    canv_all.SaveAs(outputDirectory+'/plots/{n}_all.png'.format(n=bindef['name']))
    '''
    # not sure if the following helps removing crashes
    #canv_all = None
    #canv_all.Clear()
    #for c in effis_canvas:
    #    c.Clear()
    #effis_canvas = {}
    ###
    ###
    # only return numbers and not canvases
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


