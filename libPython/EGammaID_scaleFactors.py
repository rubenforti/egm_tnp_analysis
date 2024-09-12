#!/usr/bin/env python

import sys,os
from math import sqrt
import ROOT as rt
from . import CMS_lumi, tdrstyle

from .efficiencyUtils import efficiency
from .efficiencyUtils import efficiencyList
from . import efficiencyUtils as effUtil

tdrstyle.setTDRStyle()


effiMin = 0.68
effiMax = 1.07

sfMin = 0.78
sfMax = 1.12


def isFloat( myFloat ):
    try:
        float(myFloat)
        return True
    except:
        return False



## graphColors = [rt.kBlack, rt.kGray+1, rt.kRed +1, rt.kRed-2, rt.kAzure+2, rt.kAzure-1, 
##                rt.kSpring-1, rt.kYellow -2 , rt.kYellow+1,
##                rt.kBlack, rt.kBlack, rt.kBlack, 
##                rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack ]

graphColors = [rt.kBlack, rt.kGray+1, rt.kBlue-3, rt.kBlue-9, rt.kAzure-4, rt.kAzure+8, #rt.kCyan-3, rt.kCyan-7,rt.kTeal-5, rt.kTeal+8,
               rt.kGreen-3, rt.kSpring+10, 
               rt.kOrange-2, rt.kOrange+1, rt.kRed-3, rt.kRed-9, rt.kPink-2,
               rt.kMagenta-3, rt.kViolet, rt.kCyan-7, rt.kTeal-5, rt.kYellow+1, rt.kYellow-4, rt.kOrange-8, rt.kPink-4, rt.kGray+3,
               rt.kMagenta+3, rt.kCyan+7, rt.kTeal+5, rt.kYellow+4, rt.kOrange+8, rt.kPink+4,
               33, 36, 38, 40, 41, 42, 43, 45, 46, 48, 30, 32, 20, 25, 27, 28, 29, 9, 8, 7, 6, 2, 3, 4, ]




def findMinMax( effis ):
    mini = +999
    maxi = -999

    for key in effis.keys():
        for eff in effis[key]:
            if eff['val'] - eff['err'] < mini:
                mini = eff['val'] - eff['err']
            if eff['val'] + eff['err'] > maxi:
                maxi = eff['val'] + eff['err']

    if mini > 0.18 and mini < 0.28:
        mini = 0.18
    if mini > 0.28 and mini < 0.38:
        mini = 0.28
    if mini > 0.38 and mini < 0.48:
        mini = 0.38
    if mini > 0.48 and mini < 0.58:
        mini = 0.48
    if mini > 0.58 and mini < 0.68:
        mini = 0.58
    if mini > 0.68 and mini < 0.78:
        mini = 0.68
    if mini > 0.78 and mini < 0.88:
        mini = 0.78
    if mini > 0.88:
        mini = 0.88
    if mini > 0.92:
        mini = 0.92

        
    if  maxi > 0.95:
        maxi = 1.17        
    elif maxi < 0.87:
        maxi = 0.87
    else:
        maxi = 1.07

    if maxi-mini > 0.5:
        maxi = maxi + 0.2
        
    return (mini,maxi)

    

def EffiGraph1D(effDataList, effMCList, sfList ,nameout, xAxis = 'pT', yAxis = 'eta'):
            
    W = 800
    H = 800
    yUp = 0.45
    canName = 'toto' + xAxis

    c = rt.TCanvas(canName,canName,50,50,H,W)
    c.SetTopMargin(0.055)
    c.SetBottomMargin(0.10)
    c.SetLeftMargin(0.12)
        
    p1 = rt.TPad( canName + '_up', canName + '_up', 0, yUp, 1,   1, 0,0,0)
    p2 = rt.TPad( canName + '_do', canName + '_do', 0,   0, 1, yUp, 0,0,0)
    p1.SetBottomMargin(0.0075)
    p1.SetTopMargin(   c.GetTopMargin()*1/(1-yUp))
    p2.SetTopMargin(   0.0075)
    p2.SetBottomMargin( c.GetBottomMargin()*1/yUp)
    p1.SetLeftMargin( c.GetLeftMargin() )
    p2.SetLeftMargin( c.GetLeftMargin() )
    firstGraph = True
    nGraph = len(effDataList.keys())
    nCol = 2 if nGraph < 10 else 3 if nGraph < 20 else 4
    nRow = 1 + int((nGraph-1) / nCol)
    ymaxLeg = 0.92
    yminLeg = max(0.72, ymaxLeg - nRow * 0.03)
    leg = rt.TLegend(0.15, yminLeg, 0.9, ymaxLeg)
    leg.SetNColumns(nCol)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    igr = 0
    listOfTGraph1 = []
    listOfTGraph2 = []
    listOfMC      = []

    xMin = 10
    xMax = 200
    if 'pT' in xAxis or 'pt' in xAxis:
        ## p1.SetLogx()
        ## p2.SetLogx()    
        xMin = 20
        xMax = 80
    elif 'vtx' in xAxis or 'Vtx' in xAxis or 'PV' in xAxis:
        xMin =  3
        xMax = 42
    elif 'eta' in xAxis or 'Eta' in xAxis:
        xMin = -2.60
        xMax = +2.60
    
    if 'abs' in xAxis or 'Abs' in xAxis:
        xMin = 0.0

    effminmax =  findMinMax( effDataList )
    effiMin = effminmax[0]
    effiMax = effminmax[1]

    sfminmax =  findMinMax( sfList )
    sfMin = sfminmax[0]
#    sfMin = 0.94
#    sfMax = 1.02

    for key in sorted(effDataList.keys()):
        #print('this is key', key)
        desc = 'To'.join([str(i) for i in key])#str(key[0])+str(key[1])
        desc = desc.replace(' ','').replace('.','p').replace('-','m')
        grBinsEffData = effUtil.makeTGraphFromList(effDataList[key], 'min', 'max')
        grBinsSF      = effUtil.makeTGraphFromList(sfList[key]     , 'min', 'max')
        grBinsEffData.SetName(grBinsEffData.GetName()+'_effDATA_{x}_{y}_{d}'.format(x=xAxis,y=yAxis,d=desc))
        grBinsSF.SetName(grBinsSF.GetName()+'_SF_{x}_{y}_{d}'.format(x=xAxis,y=yAxis,d=desc))
        grBinsEffMC = None
        if not effMCList is None:
            grBinsEffMC = effUtil.makeTGraphFromList(effMCList[key], 'min', 'max')
            grBinsEffMC.SetName(grBinsEffMC.GetName()+'_effMC_{x}_{y}_{d}'.format(x=xAxis,y=yAxis,d=desc))
            grBinsEffMC.SetLineStyle( rt.kDashed )
            grBinsEffMC.SetLineColor( graphColors[igr] )
            grBinsEffMC.SetMarkerSize( 0 )
            grBinsEffMC.SetLineWidth( 2 )

        grBinsSF     .SetMarkerColor( graphColors[igr] )
        grBinsSF     .SetLineColor(   graphColors[igr] )
        grBinsSF     .SetLineWidth(2)
        grBinsEffData.SetMarkerColor( graphColors[igr] )
        grBinsEffData.SetLineColor(   graphColors[igr] )
        grBinsEffData.SetLineWidth(2) 
                
        grBinsEffData.GetHistogram().SetMinimum(effiMin)
        grBinsEffData.GetHistogram().SetMaximum(effiMax)

        grBinsEffData.GetHistogram().GetXaxis().SetLimits(xMin,xMax)
        grBinsSF.GetHistogram()     .GetXaxis().SetLimits(xMin,xMax)
        grBinsSF.GetHistogram().SetMinimum(sfMin)
        grBinsSF.GetHistogram().SetMaximum(sfMax)
        
        grBinsSF.GetHistogram().GetXaxis().SetTitleOffset(1)
        if 'eta' in xAxis or 'Eta' in xAxis:
            grBinsSF.GetHistogram().GetXaxis().SetTitle("lepton #eta")
        elif 'pt' in xAxis or 'pT' in xAxis:
            grBinsSF.GetHistogram().GetXaxis().SetTitle("p_{T} (GeV)")  
        elif 'vtx' in xAxis or 'Vtx' in xAxis or 'PV' in xAxis:
            grBinsSF.GetHistogram().GetXaxis().SetTitle("n_{vtx}")  
            
        grBinsSF.GetHistogram().GetYaxis().SetTitle("data / MC " )
        grBinsSF.GetHistogram().GetYaxis().SetTitleOffset(1)
            
        grBinsEffData.GetHistogram().GetYaxis().SetTitleOffset(1)
        grBinsEffData.GetHistogram().GetYaxis().SetTitle("efficiency" )
        grBinsEffData.GetHistogram().GetYaxis().SetRangeUser( effiMin, effiMax )

            
        ### to avoid loosing the TGraph keep it in memory by adding it to a list
        listOfTGraph1.append( grBinsEffData )
        listOfTGraph2.append( grBinsSF ) 
        listOfMC.append( grBinsEffMC   )
        if 'eta' in yAxis or 'Eta' in yAxis:
            leg.AddEntry( grBinsEffData, '%1.1f #leq |#eta| #leq  %1.1f' % (float(key[0]),float(key[1])), "PL")        
        elif 'pt' in yAxis or 'pT' in yAxis:
            leg.AddEntry( grBinsEffData, '%3.0f #leq p_{T} #leq  %3.0f GeV' % (float(key[0]),float(key[1])), "PL")        
        elif 'vtx' in yAxis or 'Vtx' in yAxis or 'PV' in yAxis:
            leg.AddEntry( grBinsEffData, '%3.0f #leq nVtx #leq  %3.0f'      % (float(key[0]),float(key[1])), "PL")        

        
    for igr in range(len(listOfTGraph1)+1):

        option = "P"
        if igr == 1:
            option = "AP"

        use_igr = igr
        if use_igr == len(listOfTGraph1):
            use_igr = 0
            
        #print('use_igr and len(graphColors()', use_igr, len(graphColors))
        listOfTGraph1[use_igr].SetLineColor(graphColors[use_igr])
        listOfTGraph1[use_igr].SetMarkerColor(graphColors[use_igr])
        if not listOfMC[use_igr] is None:
            listOfMC[use_igr].SetLineColor(graphColors[use_igr])

        listOfTGraph1[use_igr].GetHistogram().SetMinimum(effiMin)
        listOfTGraph1[use_igr].GetHistogram().SetMaximum(effiMax)
        p1.cd()
        listOfTGraph1[use_igr].Draw(option)
        if listOfTGraph1[use_igr].GetTitle() == 'Graph':
            listOfTGraph1[use_igr].SetTitle('')
        if not listOfMC[use_igr] is None:
            listOfMC[use_igr].Draw("ez")

        p2.cd()            
        listOfTGraph2[use_igr].SetLineColor(graphColors[use_igr])
        listOfTGraph2[use_igr].SetMarkerColor(graphColors[use_igr])
        listOfTGraph2[use_igr].GetHistogram().SetMinimum(sfMin)
        listOfTGraph2[use_igr].GetHistogram().SetMaximum(sfMax)
        if listOfTGraph2[use_igr].GetTitle() == 'Graph':
            listOfTGraph2[use_igr].SetTitle('')
        if 'pT' in xAxis or 'pt' in xAxis :
            listOfTGraph2[use_igr].GetHistogram().GetXaxis().SetMoreLogLabels()
        listOfTGraph2[use_igr].GetHistogram().GetXaxis().SetNoExponent()
        listOfTGraph2[use_igr].Draw(option)
        

    lineAtOne = rt.TLine(xMin,1,xMax,1)
    lineAtOne.SetLineStyle(rt.kDashed)
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

    return listOfTGraph1+listOfTGraph2+listOfMC

    #################################################    


def diagnosticErrorPlot( effgr, ierror, nameout ):
    errorNames = efficiency.getSystematicNames()
    c2D_Err = rt.TCanvas('canScaleFactor_%s' % errorNames[ierror] ,'canScaleFactor: %s' % errorNames[ierror],1000,600)    
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



def doSFs(filein, lumi, axis = ['pT','eta'], plotdir='' ):
    print(" Opening file: {f} (plot lumi: {l:.1f})".format(f=filein, l=lumi ))
    CMS_lumi.lumi_13TeV = "%3.1f fb^{-1}" % lumi 

    nameOutBase = filein.replace('.txt','')
    if not os.path.exists( filein ) :
        print('file {f} does not exist'.format(f=filein))
        sys.exit(1)


    fileWithEff = open(filein, 'r')
    effGraph = efficiencyList()
    
    for line in fileWithEff :
        modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
        numbers = modifiedLine.split('\t')

        if len(numbers) > 0 and isFloat(numbers[0]):
            etaKey = ( float(numbers[0]), float(numbers[1]) )
            ptKey  = ( float(numbers[2]), min(500,float(numbers[3])) )
        
            myeff = efficiency( ptKey, etaKey,
                                float(numbers[4]), float(numbers[5]), ## data eff and error
                                float(numbers[6]), float(numbers[7]), ## mc eff and error
                                float(numbers[8]), float(numbers[9]),    ## data_alt_sig eff and error
                                float(numbers[10]), float(numbers[11]),  ## data_alt_bkg eff and error
                                float(numbers[12]), float(numbers[13]),  ## eff_alt_sig eff and error
                                float(numbers[14]) )
            effGraph.addEfficiency(myeff)

    fileWithEff.close()

### massage the numbers a bit
    #marc test effGraph.symmetrizeSystVsEta()
    #marc test effGraph.combineSyst()

    print(" ------------------------------- ")

    ## marc customEtaBining = []
    ## marc customEtaBining.append( (0.000,0.800))
    ## marc customEtaBining.append( (0.800,1.444))
    ## marc customEtaBining.append( (1.444,1.566))
    ## marc customEtaBining.append( (1.566,2.000))
    ## marc customEtaBining.append( (2.000,2.500))


    pdfout = nameOutBase + '_efficiencyPlots.pdf'
    cDummy = rt.TCanvas()
    cDummy.Print( pdfout + "[" )

    #print('this is effgraph', effGraph)
    #print('this is effGraph.pt_1DGraph_list)', effGraph.pt_1DGraph_list)

    listOfSF1D = EffiGraph1D( effGraph.pt_1DGraph_list( False ) , #eff Data
                 #None, 
                 effGraph.pt_1DGraph_list( False, typeGR = -1 ),
                 effGraph.pt_1DGraph_list( True ) , #SF
                 pdfout,
                 xAxis = axis[0], yAxis = axis[1] )

    #print('this is length of listOfSF1D', len(listOfSF1D))
#EffiGraph1D( effGraph.pt_1DGraph_list_customEtaBining(customEtaBining,False) , 
#             effGraph.pt_1DGraph_list_customEtaBining(customEtaBining,True)   , False, pdfout )
#    EffiGraph1D( effGraph.eta_1DGraph_list(False), effGraph.eta_1DGraph_list(True), True , pdfout )
    listOfSF1D.extend( EffiGraph1D( effGraph.eta_1DGraph_list( typeGR =  0 ) , # eff Data
                              effGraph.eta_1DGraph_list( typeGR = -1 ) , # eff MC
                              effGraph.eta_1DGraph_list( typeGR = +1 ) , # SF
                              pdfout, 
                              xAxis = axis[1], yAxis = axis[0] ) )

    h2EffData       = effGraph.ptEtaScaleFactor_2DHisto(40)
    h2EffMC         = effGraph.ptEtaScaleFactor_2DHisto(41)
    h2EffDataAltSig = effGraph.ptEtaScaleFactor_2DHisto(42)
    h2EffDataAltBkg   = effGraph.ptEtaScaleFactor_2DHisto(43)
    h2EffMCAltSig   = effGraph.ptEtaScaleFactor_2DHisto(44)
    h2StatUncEffData       = effGraph.ptEtaScaleFactor_2DHisto(50)
    h2StatUncEffMC         = effGraph.ptEtaScaleFactor_2DHisto(51)
    h2StatUncEffDataAltSig = effGraph.ptEtaScaleFactor_2DHisto(52)
    h2StatUncEffDataAltBkg   = effGraph.ptEtaScaleFactor_2DHisto(53)
    h2StatUncEffMCAltSig   = effGraph.ptEtaScaleFactor_2DHisto(54)
    #h2SF      = effGraph.ptEtaScaleFactor_2DHisto(-1)
    h2SF = h2EffData.Clone('h2_SF_nominal');
    h2SF.SetTitle('SF nominal data / nominal MC')
    h2SF.Divide(h2EffMC)
    h2SFDataAltSig = h2EffDataAltSig.Clone('h2_SF_dataAltSig');
    h2SFDataAltSig.SetTitle('SF altSig data / nominal MC')
    h2SFDataAltSig.Divide(h2EffMC)
    h2SFDataAltBkg = h2EffDataAltBkg.Clone('h2_SF_dataAltBkg');
    h2SFDataAltBkg.SetTitle('SF altBkg data / nominal MC')
    h2SFDataAltBkg.Divide(h2EffMC)
    h2SFMCAltSig = h2EffData.Clone('h2_SF_MCAltSig');
    h2SFMCAltSig.SetTitle('SF nominal data / altSig MC')
    h2SFMCAltSig.Divide(h2EffMCAltSig)
    h2SFDataMCAltSig = h2EffDataAltSig.Clone('h2_SF_dataMCAltSig');
    h2SFDataMCAltSig.SetTitle('SF altSig data / altSig MC')
    h2SFDataMCAltSig.Divide(h2EffMCAltSig)
    h2Error   = effGraph.ptEtaScaleFactor_2DHisto( 0)  ## only error bars

    rt.gStyle.SetPalette(87)
    rt.gStyle.SetPaintTextFormat('1.3f');
    rt.gStyle.SetOptTitle(1)

    c2D = rt.TCanvas('canScaleFactor','canScaleFactor',900,600)
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
    dmin = 1.0 - h2SF.GetMinimum()
    dmax = h2SF.GetMaximum() - 1.0
    #dall = max(dmin,dmax)
    #h2SF.SetMinimum(1-dall)
    #h2SF.SetMaximum(1+dall)
    h2SF.DrawCopy("colz TEXT45")
    
    c2D.cd(2)
    h2Error.SetMinimum(0)
    h2Error.SetMaximum(min(h2Error.GetMaximum(),0.2))    
    h2Error.DrawCopy("colz TEXT45")

    c2D.Print( pdfout )

    #print('this is listOfSF1D', listOfSF1D)
    print("NameOutBase = ",nameOutBase + '_2D.root')
    rootout = rt.TFile(nameOutBase + '_2D.root','recreate')
    rootout.cd()
    h2SF.Write('SF2D_nominal',rt.TObject.kOverwrite)
    h2EffData.Write('EffData2D',rt.TObject.kOverwrite)
    h2EffMC  .Write('EffMC2D'  ,rt.TObject.kOverwrite)
    h2EffDataAltSig.Write('EffDataAltSig2D', rt.TObject.kOverwrite)
    h2EffMCAltSig  .Write('EffMCAltSig2D'  , rt.TObject.kOverwrite)
    h2SFDataAltSig.Write('SF2D_dataAltSig', rt.TObject.kOverwrite)
    h2SFDataAltBkg.Write('SF2D_dataAltBkg', rt.TObject.kOverwrite)
    h2SFMCAltSig.Write('SF2D_MCAltSig', rt.TObject.kOverwrite)
    h2SFDataMCAltSig.Write('SF2D_dataMCAltSig', rt.TObject.kOverwrite)
    for igr in listOfSF1D:
        igr.Write( igr.GetName(), rt.TObject.kOverwrite) #'grSF1D_{ib}'.format(ib=igr), rt.TObject.kOverwrite)
    rootout.Close()

    #for isyst in range(len(efficiency.getSystematicNames())):
    #    diagnosticErrorPlot( effGraph, isyst, pdfout )

    cDummy.Print( pdfout + "]" )

    allEffsAndSFs= [h2SF, h2EffData, h2EffMC, h2EffDataAltSig, h2EffDataAltBkg, h2EffMCAltSig, h2SFDataAltSig, h2SFDataAltBkg, h2SFMCAltSig, h2SFDataMCAltSig, h2StatUncEffData, h2StatUncEffMC, h2StatUncEffDataAltSig, h2StatUncEffDataAltBkg, h2StatUncEffMCAltSig]
    canv = rt.TCanvas('c','c', 800, 800)
    canv.SetRightMargin(0.20)
    rt.gStyle.SetPalette(55)
    rt.gStyle.SetNumberContours(51)
    for hist in allEffsAndSFs:
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
    CMS_lumi.lumi_sqrtS = "13 TeV"
    
    axis = ['pT','eta']
    if args.PV:
        axis = ['nVtx','eta']

    doSFs(args.txtFile, args.lumi,axis,'')
