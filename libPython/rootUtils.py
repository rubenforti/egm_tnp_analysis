#!/usr/bin/env python3 

import os
import ctypes
import math
import ROOT


def addStringToEnd(name, matchToAdd, notAddIfEndswithMatch=False):
    if notAddIfEndswithMatch and name.endswith(matchToAdd):
        return name
    elif not name.endswith(matchToAdd):
        return name + matchToAdd


def compileMacro(x): #, basedir=os.environ['PWD']):
    #ROOT.gROOT.ProcessLine(".L %s/%s+" % (os.environ['CMSSW_BASE'],x));
    success = ROOT.gSystem.CompileMacro("%s" % (x), "k")
    if not success:
        print("Loading and compiling %s failed! Exit" % x)
        quit()


def compileFileMerger(x):
    y=x.strip(".C")
    print(f"Compiling {x} into {y}")
    res = os.system(f"g++ `root-config --libs --cflags --glibs` -O3 {x} -o {y}")
    if res:
        print("Compiling %s failed! Exit" % x)
        quit()


def safeGetObject(fileObject, objectName, quitOnFail=True, silent=False, detach=True):
    obj = fileObject.Get(objectName)
    if obj == None:
        if not silent:
            print(f"Error getting {objectName} from file {fileObject.GetName()}")
        if quitOnFail:
            quit()
        return None
    else:
        if detach:
            obj.SetDirectory(0)
        return obj


def safeOpenFile(fileName, quitOnFail=True, silent=False, mode="READ"):
    fileObject = ROOT.TFile.Open(fileName, mode)
    if not fileObject or fileObject.IsZombie():
        if not silent:
            print(f"Error when opening file {fileName}")
        if quitOnFail:
            quit()
        else:
            return None
    elif not fileObject.IsOpen():
        if not silent:
            print(f"File {fileName} was not opened")
        if quitOnFail:
            quit()
        else:
            return None
    else:
        return fileObject


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


def computeEffi(n1, n2, e1, e2):
    if (n1+n2):
        eff  = n1/(n1+n2)
        nTot = n1+n2
        e_eff = math.sqrt(e1*e1*n2*n2+e2*e2*n1*n1) / (nTot*nTot)
    else:
        eff, e_eff = 1.1, 0.01

    effout = [eff, e_eff]
 
    return effout


def getAllEffi(info, bindef):

    effis = {}
    binName = bindef["name"]
    for key, value in info.items():
        effis[key] = [-1, -1]       

        if not value or not os.path.isfile(value):
            continue
        
        rootfile = safeOpenFile(value, mode="READ")

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
        
        effis[key] = computeEffi(nP, nF, eP, eF) + [nP, eP, nF, eF]
        rootfile.Close()

    return effis

def plotAllEffi(info, bindef, outputDirectory, effis):
    effis_canvas = {}
    binName = bindef["name"]
    canvasesToGet = ["MC_Nominal_fit", "Data_Nominal", "Data_Alt_Sig", "Data_Alt_Bkg"] 
    
    for key, value in info.items():
        if not (value and os.path.isfile(value)):
            if key in canvasesToGet:
                print(f"Warning: {key} fit not found, skipping...")
                canvasesToGet.remove(key)
            continue

        rootfile = safeOpenFile(value, mode="READ")
        effis_canvas[f"canv_{key}"] = safeGetObject(rootfile, f"{binName}_Canv", detach=False) if key in canvasesToGet else None
        rootfile.Close()


    ncols = len(canvasesToGet)

    canv_all = ROOT.TCanvas(bindef['name'], bindef['name'], 400*ncols, 1200)
    canv_all.Draw()
    pad_title = ROOT.TPad('title', 'title', 0.0, 0.95, 1.0, 1.0) 
    pad_pass = ROOT.TPad('passing', 'passing', 0.0, 0.55, 1.0, 0.95)
    pad_fail = ROOT.TPad('failing', 'failing', 0.0, 0.15, 1.0, 0.55)
    pad_eff = ROOT.TPad('efficiency', 'efficiency', 0.0, 0.0, 1.0, 0.15)

    pad_title.Draw()
    pad_pass.Divide(ncols, 1), pad_pass.Draw()
    pad_fail.Divide(ncols, 1), pad_fail.Draw()
    pad_eff.Divide(ncols, 1),  pad_eff.Draw()

    pad_title.cd()
    txt = ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.5)
    txt.SetNDC()
    txt.DrawLatex(0.4, 0.5, f'{bindef['name'].replace('_','  ').replace('To', ' - ').replace('probe ', '').replace('m',' - ').replace('pt','XX').replace('p','.').replace('XX','p_{T}')}')
    txt.SetTextSize(0.15)

    for icol, canv_type in enumerate(canvasesToGet):

        for ip, p in enumerate(effis_canvas[f"canv_{canv_type}"].GetListOfPrimitives()):
            if not ip: continue
            pad_to_use = pad_pass if ip==1 else pad_fail if ip==2 else None
            pad_to_use.cd(icol+1)
            newp = p.Clone(f"tmp_{canv_type}_{ip}")
            newp.SetPad(0.02, 0.00, 0.98, 0.98)
            newp.Draw()
        pad_eff.cd(icol+1)
        txt.SetTextFont(62)
        txt.DrawLatex(0.10, 0.77, canv_type.replace('_', " ") + " :" if not "MC_Nominal" in canv_type else "MC counting :")
        txt.SetTextFont(42)
        tmp = effis[canv_type]
        txt.DrawLatex(0.20, 0.59, f'Passing: {tmp[2]:.1f} #pm {tmp[3]:.1f}')
        txt.DrawLatex(0.20, 0.41, f'Failing: {tmp[4]:.1f} #pm {tmp[5]:.1f}')
        txt.SetTextFont(62)
        txt.DrawLatex(0.20, 0.23, f'Efficiency: {tmp[0]*100.:.2f} #pm {tmp[1]*100.:.2f} %')
        txt.SetTextFont(42)

    canv_all.cd()
    canv_all.Update()
    canv_all.SaveAs(f"{outputDirectory}/{bindef['name']}_all.pdf")
    canv_all.SaveAs(f"{outputDirectory}/{bindef['name']}_all.png")


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


