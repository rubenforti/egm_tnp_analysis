#!/usr/bin/env python

import sys,os
import ROOT

from libPython import CMS_lumi, tdrstyle
from libPython.efficiencyUtils import efficiency, efficiencyManager, makeTGraphFromListEff

#tdrstyle.setTDRStyle()


effiMin = 0.68
effiMax = 1.07

sfMin = 0.78
sfMax = 1.12


def isFloat(myFloat):
    try:
        float(myFloat)
        return True
    except:
        return False

graphColors = [ROOT.kBlack,    ROOT.kGray+1,   ROOT.kBlue-3,    ROOT.kBlue-9,    ROOT.kAzure-4, #ROOT.kCyan-3, ROOT.kCyan-7,ROOT.kTeal-5, ROOT.kTeal+8,
               ROOT.kAzure+8,  ROOT.kGreen-3,  ROOT.kSpring+10, ROOT.kOrange-2,  ROOT.kOrange+1,
               ROOT.kRed-3,    ROOT.kRed-9,    ROOT.kPink-2,    ROOT.kMagenta-3, ROOT.kViolet, 
               ROOT.kCyan-7,   ROOT.kTeal-5,   ROOT.kYellow+1,  ROOT.kYellow-4,  ROOT.kOrange-8, 
               ROOT.kPink-4,   ROOT.kGray+3,   ROOT.kMagenta+3, ROOT.kCyan+7,    ROOT.kTeal+5, 
               ROOT.kYellow+4, ROOT.kOrange+8, ROOT.kPink+4,
               33, 36, 38, 40, 41, 42, 43, 45, 46, 48, 30, 32, 20, 25, 27, 28, 29, 9, 8, 7, 6, 2, 3, 4]


def findMinMax(effis):
    """
    """
    mini = 999
    maxi = -999

    for eff_key, eff_bin in effis.items():
        for eff in eff_bin:
            mini = min(mini, eff['val']-eff['err'])
            maxi = max(maxi, eff['val']+eff['err'])

    mini -= 0.01
    maxi += 0.01

    return mini, maxi

    

def EffiGraph1D(effDataList, effMCList, sfList, nameout, xAxis='pT', yAxis='eta'):
            
    W = 800
    H = 800
    yUp = 0.45
    canName = 'toto' + xAxis

    c = ROOT.TCanvas(canName, canName, 50, 50, H, W)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.10)
    c.SetLeftMargin(0.1)
    c.SetRightMargin(0.1)
        
    p1 = ROOT.TPad(f"{canName}_up", f"{canName}_up", 0, yUp, 1, 1,   0, 0, 0)
    p2 = ROOT.TPad(f"{canName}_do", f"{canName}_do", 0, 0,   1, yUp, 0, 0, 0)
    p1.SetBottomMargin(0.0075)
    p1.SetTopMargin(c.GetTopMargin() * 1 / (1 - yUp))
    p2.SetTopMargin(0.0075)
    p2.SetBottomMargin(c.GetBottomMargin() * 1 / yUp)
    p1.SetLeftMargin(c.GetLeftMargin())
    p2.SetLeftMargin(c.GetLeftMargin())
    nGraphs = len(effDataList.keys())
    nCol = 2 if nGraphs < 10 else 3 if nGraphs < 20 else 4
    nRow = 1 + int((nGraphs-1) / nCol)
    ymaxLeg = 0.945
    yminLeg = ymaxLeg - nRow*0.025
    leg = ROOT.TLegend(0.15, yminLeg, 0.85, ymaxLeg)
    leg.SetNColumns(nCol)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)


    xMin = 10
    xMax = 200
    if 'pT' in xAxis or 'pt' in xAxis:
        xMin = 20
        xMax = 70
    elif 'eta' in xAxis or 'Eta' in xAxis:
        xMin = -2.50
        xMax = +2.50
    else:
        print("xAxis not recognized")
        sys.exit(1)
    
    if 'abs' in xAxis or 'Abs' in xAxis:
        xMin = 0.0

    effiMin, effiMax = findMinMax(effDataList)
    sfMin,   sfMax   = findMinMax(sfList)

    print(f"effiMin, effiMax = {effiMin}, {effiMax}")
    print(f"sfMin, sfMax = {sfMin}, {sfMax}")

    listTGraphsData, listTGraphsMC, listTGraphsSF = [], [], []

    for igr, key in enumerate(sorted(effDataList.keys())):

        desc = 'To'.join([str(i) for i in key])
        desc = desc.replace(' ','').replace('.','p').replace('-','m')
        
        grBinsEffData = makeTGraphFromListEff(effDataList[key], 'min', 'max')
        grBinsEffData.SetName(f"{grBinsEffData.GetName()}_effDATA_{xAxis}_{yAxis}_{desc}")
        grBinsEffData.GetHistogram().GetXaxis().SetLimits(xMin, xMax)
        grBinsEffData.GetHistogram().GetYaxis().SetTitleOffset(1)
        grBinsEffData.GetHistogram().GetYaxis().SetTitle("Efficiency Data")
        grBinsEffData.GetHistogram().GetYaxis().SetRangeUser(effiMin, effiMax+(nRow*0.006))
        
        grBinsSF = makeTGraphFromListEff(sfList[key], 'min', 'max')
        grBinsSF.SetName(f"{grBinsSF.GetName()}_SF_{xAxis}_{yAxis}_{desc}")
        

        grBinsSF.GetHistogram().GetXaxis().SetLimits(xMin, xMax)
        if 'eta' in xAxis or 'Eta' in xAxis:
            grBinsSF.GetHistogram().GetXaxis().SetTitle("#eta^{#mu}")
        elif 'pt' in xAxis or 'pT' in xAxis:
            grBinsSF.GetHistogram().GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        grBinsSF.GetHistogram().GetXaxis().SetTitleOffset(1)
        grBinsSF.GetHistogram().GetYaxis().SetTitle("Data / MC")
        grBinsSF.GetHistogram().GetYaxis().SetTitleOffset(1)
        grBinsSF.GetHistogram().GetYaxis().SetRangeUser(sfMin, sfMax)

        grBinsEffMC = makeTGraphFromListEff(effMCList[key], 'min', 'max')
        grBinsEffMC.SetName(f"{grBinsEffMC.GetName()}_effMC_{xAxis}_{yAxis}_{desc}")
        grBinsEffMC.SetLineStyle(ROOT.kDashed)
        grBinsEffMC.SetLineColor(graphColors[igr])
        grBinsEffMC.SetMarkerSize(0)
        grBinsEffMC.SetLineWidth(2)

        option = "P" if igr != 0 else "AP"

        p1.cd()
        grBinsEffData.SetMarkerColor(graphColors[igr])
        grBinsEffData.SetLineColor(graphColors[igr])
        grBinsEffData.SetLineWidth(2)
        grBinsEffData.Draw(option)
        if grBinsEffData.GetTitle() == 'Graph':
            grBinsEffData.SetTitle('')
        if grBinsEffMC is not None:
            grBinsEffMC.Draw("ez")

        p2.cd()
        grBinsSF.SetMarkerColor(graphColors[igr])
        grBinsSF.SetLineColor(graphColors[igr])
        grBinsSF.SetLineWidth(2)
        if grBinsSF.GetTitle() == 'Graph':
            grBinsSF.SetTitle('')
        if 'pT' in xAxis or 'pt' in xAxis:
            grBinsSF.GetHistogram().GetXaxis().SetMoreLogLabels()
        grBinsSF.GetHistogram().GetXaxis().SetNoExponent()
        grBinsSF.Draw(option)

        if 'eta' in yAxis or 'Eta' in yAxis:
            leg.AddEntry(grBinsEffData, f"{float(key[0]):1.1f} #leq |#eta| #leq {float(key[1]):1.1f}", "PL")
        elif 'pt' in yAxis or 'pT' in yAxis:
            leg.AddEntry(grBinsEffData, f"{float(key[0]):3.0f} #leq p_{{T}} #leq {float(key[1]):3.0f} GeV", "PL")

        listTGraphsData.append(grBinsEffData)
        listTGraphsMC.append(grBinsEffMC)
        listTGraphsSF.append(grBinsSF)
        

    lineAtOne = ROOT.TLine(xMin, 1, xMax,1)
    lineAtOne.SetLineStyle(ROOT.kDashed)
    lineAtOne.SetLineWidth(2)
    
    p2.cd()
    lineAtOne.Draw()

    c.cd()
    p2.Draw()
    p1.Draw()
    CMS_lumi.CMS_lumi(c, 4, 0)
    

    leg.Draw()    
    c.RedrawAxis("sameaxis")
    
    c.Print(nameout)

    return listTGraphsData+listTGraphsSF+listTGraphsMC

    #################################################    


def diagnosticErrorPlot( effgr, ierror, nameout ):
    errorNames = efficiency.getSystematicNames()
    c2D_Err = ROOT.TCanvas('canScaleFactor_%s' % errorNames[ierror] ,'canScaleFactor: %s' % errorNames[ierror],1000,600)    
    c2D_Err.Divide(2,1)
    c2D_Err.GetPad(1).SetLogy()
    c2D_Err.GetPad(2).SetLogy()
    c2D_Err.GetPad(1).SetRightMargin(0.15)
    c2D_Err.GetPad(1).SetLeftMargin( 0.15)
    c2D_Err.GetPad(1).SetTopMargin(  0.10)
    c2D_Err.GetPad(2).SetRightMargin(0.15)
    c2D_Err.GetPad(2).SetLeftMargin( 0.15)
    c2D_Err.GetPad(2).SetTopMargin(  0.10)

    h2_sfErrorAbs = effgr.ptEtaScaleFactor_2DHisto(ierror+1, False )
    h2_sfErrorRel = effgr.ptEtaScaleFactor_2DHisto(ierror+1, True  )
    h2_sfErrorAbs.SetMinimum(0)
    h2_sfErrorAbs.SetMaximum(min(h2_sfErrorAbs.GetMaximum(),0.2))
    h2_sfErrorRel.SetMinimum(0)
    h2_sfErrorRel.SetMaximum(1)
    h2_sfErrorAbs.SetTitle('lepton absolute SF syst: %s ' % errorNames[ierror])
    h2_sfErrorRel.SetTitle('lepton relative SF syst: %s ' % errorNames[ierror])
    c2D_Err.cd(1)
    h2_sfErrorAbs.DrawCopy("colz TEXT45")
    c2D_Err.cd(2)
    h2_sfErrorRel.DrawCopy("colz TEXT45")
    
    c2D_Err.Print(nameout)



def doSFs(filein, lumi, axis=['pT','eta'], plotdir='' ):
    """
    """
    print(f"Opening file: {filein} (plot lumi: {lumi:.1f})")
    CMS_lumi.lumi_13TeV = "%3.1f fb^{-1}" % lumi 

    nameOutBase = filein.replace('.txt','')

    if not os.path.exists( filein ) :
        print(f'File {filein} does not exist')
        sys.exit(1)

    fileWithEff = open(filein, 'r')
    effGraph = efficiencyManager()
    
    for line in fileWithEff :
        modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
        numbers = modifiedLine.split('\t')

        if len(numbers) > 0 and isFloat(numbers[0]):
            etaKey = ( float(numbers[0]), float(numbers[1]) )
            ptKey  = ( float(numbers[2]), float(numbers[3]) ) 
        
            myeff = efficiency(ptKey, etaKey,
                               float(numbers[4]), float(numbers[5]), ## data eff and error
                               float(numbers[6]), float(numbers[7]), ## mc eff and error
                               float(numbers[8]), float(numbers[9]),   ## data_alt_sig eff and error
                               float(numbers[10]), float(numbers[11]), ## data_alt_bkg eff and error
                               float(numbers[12]), float(numbers[13]), ## eff_alt_sig eff and error
                               float(numbers[14]) )
            effGraph.addEfficiency(myeff)

    fileWithEff.close()

    pdfout = nameOutBase + '_efficiencyPlots.pdf'
    cDummy = ROOT.TCanvas()
    cDummy.Print( pdfout + "[" )
    
        
    listOfSF1D = EffiGraph1D(effGraph.get1DGraphList("pt"),  # eff Data
                             effGraph.get1DGraphList("pt", effMC=True),  # eff MC
                             effGraph.get1DGraphList("pt", doScaleFactor=True),  # SF
                             pdfout, xAxis=axis[0], yAxis=axis[1])

    listOfSF1D.extend(
        EffiGraph1D(effGraph.get1DGraphList("eta"),  # eff Data
                    effGraph.get1DGraphList("eta", effMC=True),  # eff MC
                    effGraph.get1DGraphList("eta", doScaleFactor=True),  # SF
                    pdfout, xAxis = axis[1], yAxis = axis[0]) )


    listOf2DKeys = [
        #typePlot, isMC, typeEff
        ("eff", False, "nominal"),
        ("err", False, "nominal"),
        ("eff", True,  "nominal"),
        ("err", True,  "nominal"),
        ("eff", False, "altSig"),
        ("err", False, "altSig"),
        ("eff", False, "altBkg"),
        ("err", False, "altBkg"),
        ("syst", False, "altSig"),
        ("syst", False, "altBkg"),
        ("syst", False, "allSyst"),
        ("sf", False, "nominal"),
        ("sf", False, "altSig"),
        ("sf", False, "altBkg"),
    ]

    h2D_dict = {}
    for key in listOf2DKeys:
        h2 = effGraph.ptEtaScaleFactor_2DHisto(typePlot=key[0], isMC=key[1], typeEff=key[2])
        h2D_dict[h2.GetName()] = h2


    ROOT.gStyle.SetPalette(87)
    ROOT.gStyle.SetPaintTextFormat('1.3f')
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat(0)

    c2D = ROOT.TCanvas('canScaleFactor','canScaleFactor',900,600)
    c2D.Divide(2,1)
    c2D.GetPad(1).SetRightMargin(0.15)
    c2D.GetPad(1).SetLeftMargin( 0.15)
    c2D.GetPad(1).SetTopMargin(  0.10)
    c2D.GetPad(2).SetRightMargin(0.15)
    c2D.GetPad(2).SetLeftMargin( 0.15)
    c2D.GetPad(2).SetTopMargin(  0.10)
    c2D.GetPad(1).SetLogy()
    c2D.GetPad(2).SetLogy()
    

    c2D.cd(1)
    h2D_dict['h2_sfDataNominal'].DrawCopy("colz TEXT45")
    
    c2D.cd(2)
    h2Error = h2D_dict["h2_statUncDataNominal"].Clone('h2Error')
    h2Error.SetMinimum(0)
    h2Error.SetMaximum(min(h2Error.GetMaximum(), 0.2))    
    h2Error.DrawCopy("colz TEXT45")

    c2D.Print(pdfout)

    print(f"NameOutBase = {nameOutBase}_2D.root")
    rootout = ROOT.TFile(f"{nameOutBase}_2D.root", "RECREATE")
    rootout.cd()

    for h2D in h2D_dict.values():
        h2D.Write(h2D.GetName(), ROOT.TObject.kOverwrite)
    for igr in listOfSF1D:
        igr.Write(igr.GetName(), ROOT.TObject.kOverwrite) #'grSF1D_{ib}'.format(ib=igr), ROOT.TObject.kOverwrite)
    rootout.Close()

    #for isyst in range(len(efficiency.getSystematicNames())):
    #    diagnosticErrorPlot( effGraph, isyst, pdfout )

    cDummy.Print( pdfout + "]" )

    canv = ROOT.TCanvas('c','c', 1200, 900)
    canv.SetRightMargin(0.20)
    ROOT.gStyle.SetPalette(55)
    ROOT.gStyle.SetNumberContours(51)
    for hist in h2D_dict.values():
        canv.Clear()
        hist.Draw('colz')
        canv.SaveAs(plotdir+'/'+hist.GetName()+'.png')
        canv.SaveAs(plotdir+'/'+hist.GetName()+'.pdf')


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='tnp EGM scale factors')
    parser.add_argument('--lumi'  , type = float, default = -1, help = 'Lumi (just for plotting purpose)')
    parser.add_argument('txtFile' , default = None, help = 'EGM formatted txt file')
    parser.add_argument('--PV'    , action  = 'store_true', help = 'plot 1 vs nVtx instead of pT' )
    args = parser.parse_args()

    if args.txtFile is None:
        print(' - Needs EGM txt file as input')
        sys.exit(1)
    

    CMS_lumi.lumi_13TeV = "5.5 fb^{-1}"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.lumi_S = "13 TeV"
    
    axis = ['pT','eta']
    if args.PV:
        axis = ['nVtx','eta']

    doSFs(args.txtFile, args.lumi,axis,'')
