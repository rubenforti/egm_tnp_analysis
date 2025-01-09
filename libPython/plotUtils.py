#!/usr/bin/env python3

import os
import copy
import shutil
from array import array
import numpy as np
import re
import math

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from libPython.CMS_lumi_new import *

def printLine(marker='-', repeat=30):
    print(marker*repeat)

    
def addStringToEnd(name, matchToAdd, notAddIfEndswithMatch=False):
    if notAddIfEndswithMatch and name.endswith(matchToAdd):
        return name
    elif not name.endswith(matchToAdd):
        return name + matchToAdd


def createPlotDirAndCopyPhp(outdir):
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    htmlpath = "etc/inputs/index.php"
    shutil.copy(htmlpath, outdir)


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


def testBinning(bins, testbins, var="var", flag="workingPoint", allowRebin=False):
    """
    """
    if bins!=testbins:
        if bins[0] in testbins:
            firstTestIdx = testbins.index(bins[0])
            if all(bins[i]==testbins[i+firstTestIdx] for i in range(len(bins))):
                print(f"\nWarning: {var} binning not consistent with the one in histograms for {flag}")
                print(f"{bins}")
                print(f"{testbins}")
                print(f"However it seems to be a slice of it, so I will continue assuming it is intentional. Proceed with caution!\n\n")
                return 0
            elif allowRebin and all(bins[i] in testbins for i in range(len(bins))):
                print(f"\nWarning: {var} binning not consistent with the one in histograms for {flag}")
                print(f"{bins}")
                print(f"{testbins}")
                print(f"However it seems to be a subset of it, so I will continue assuming you wanted to rebin. Proceed with caution!\n\n")
                return 0
        print(f"\nError: {var} binning not consistent with the one in histograms for {flag}")
        print(f"{bins}")
        print(f"{testbins}")
        print("Please check!\n\n")
        return -1
    else:
        return 0
    
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


def adjustSettings_CMS_lumi():

    ## dummy function to be called before using any other fucntion calling CMS_lumi
    ## for some reason, the settings of the very first plot are screwed up.
    ## To fix this issue, it is enough to call it to a dummy plot
    
    dummy = ROOT.TH1D("dummy","",10,0,10)
    for i in range(1,1+dummy.GetNbinsX()):
        dummy.SetBinContent(i,i)
    dummy.GetXaxis().SetTitle("x axis")
    dummy.GetYaxis().SetTitle("y axis")
    cdummy = ROOT.TCanvas("cdummy","",600,600)
    dummy.Draw("HE")
    CMS_lumi_new(cdummy,"",True,False)
    setTDRStyle()
    ## no need to save the canvas

def getAxisRangeFromUser(axisNameTmp="",
                         separator="::",
                         rangeSeparator=","
                         ):

    setXAxisRangeFromUser = False;
    fields = axisNameTmp.split(separator)
    axisName = fields[0]

    if len(fields) > 1:
        setXAxisRangeFromUser = True;
        xmin = float(fields[1].split(rangeSeparator)[0])
        xmax = float(fields[1].split(rangeSeparator)[1])
    else:
        xmin = 0
        xmax = 0

    return axisName,setXAxisRangeFromUser,xmin,xmax

    
def drawTH2(h2D_tmp,
            labelXtmp="xaxis", labelYtmp="yaxis", labelZtmp="zaxis",
            canvasName="default", plotLabel="", outdir="./",
            rebinFactorX=0,
            rebinFactorY=0,
            smoothPlot=False,
            drawProfileX=False,
            scaleToUnitArea=False,
            draw_both0_noLog1_onlyLog2=1,
            leftMargin=0.16,
            rightMargin=0.20,
            nContours=51,
            palette=55,
            invertePalette=False,
            canvasSize="700,625",
            passCanvas=None,
            bottomMargin=0.1,
            plotError=False,
            plotRelativeError=False,
            lumi=None,
            drawOption = "colz",
            skipLumi=False
):


    ROOT.TH1.SetDefaultSumw2()
    adjustSettings_CMS_lumi()

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h2D_tmp.RebinX(rebinFactorX)
        else:                             h2D_tmp.RebinX(len(rebinFactorX)-1,"",array('d',rebinFactorX)) # case in which rebinFactorX is a list of bin edges

    if (rebinFactorY): 
        if isinstance(rebinFactorY, int): h2D_tmp.RebinY(rebinFactorY)
        else:                             h2D_tmp.RebinY(len(rebinFactorY)-1,"",array('d',rebinFactorY)) # case in which rebinFactorX is a list of bin edges

    if plotError or plotRelativeError:
        herr = h2D_tmp.Clone(h2D_tmp.GetName()+"_err")
        herr.Reset("ICESM")
        for i in range(1,herr.GetNbinsX()+1):
            for j in range(1,herr.GetNbinsY()+1):
                errval = h2D_tmp.GetBinError(i,j)
                if plotRelativeError:
                    if h2D_tmp.GetBinContent(i,j) != 0.0:
                        errval = errval/h2D_tmp.GetBinContent(i,j)
                    else:
                        errval = 1.0 if errval == 0 else 0.0
                herr.SetBinContent(i,j,errval)
        h2D = herr
    else:
        h2D = h2D_tmp

    # dark blue to red
    ROOT.TColor.CreateGradientColorTable(4,
                                         array ("d", [0.00, 0.45, 0.55, 1.00]),
                                         array ("d", [0.00, 1.00, 1.00, 1.00]),
                                         array ("d", [0.00, 1.00, 1.00, 0.00]),
                                         array ("d", [1.00, 1.00, 1.00, 0.00]),
                                         255,  1.0)
                                         # array ("d", [0.00, 0.50, 1.00]),
                                         # array ("d", [0.00, 1.00, 1.00]),
                                         # array ("d", [0.00, 1.00, 0.00]),
                                         # array ("d", [1.00, 1.00, 0.00]),
                                         # 255,  0.95)

    if palette > 0:
        ROOT.gStyle.SetPalette(palette)  # 55:raibow palette ; 57: kBird (blue to yellow, default) ; 107 kVisibleSpectrum ; 77 kDarkRainBow 
    ROOT.gStyle.SetNumberContours(nContours) # default is 20 
    if invertePalette:
        ROOT.TColor.InvertPalette()

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    labelZ,setZAxisRangeFromUser,zmin,zmax = getAxisRangeFromUser(labelZtmp)
        
    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)    
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.SetBottomMargin(bottomMargin)
    canvas.cd()

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)
    # normalize to 1
    if (scaleToUnitArea): h2D.Scale(1./h2D.Integral())

    h2DGraph = 0

    h2DPlot = 0
    if (not smoothPlot): h2DPlot = h2D
    else:
        h2DGraph = ROOT.TGraph2D()
        h2DGraph.SetNpx(300)
        h2DGraph.SetNpy(300)
        nPoint = 0
        for iBinX in range (1,1+h2D.GetNbinsX()):
            for iBinY in range(1,1+h2D.GetNbinsY()):
                h2DGraph.SetPoint(nPoint,h2D.GetXaxis().GetBinCenter(iBinX),h2D.GetYaxis().GetBinCenter(iBinY),h2D.GetBinContent(iBinX,iBinY))
                nPoint += 1
            

        h2DPlot = h2DGraph.GetHistogram()

    if plotLabel == "ForceTitle":
        h2DPlot.SetTitle(h2D_tmp.GetTitle())
  
    h2DPlot.GetXaxis().SetTitle(labelX)
    h2DPlot.GetYaxis().SetTitle(labelY)
    h2DPlot.GetXaxis().SetTitleSize(0.05)
    h2DPlot.GetXaxis().SetLabelSize(0.04)
    h2DPlot.GetXaxis().SetTitleOffset(0.95) # 1.1 goes outside sometimes, maybe depends on root version or canvas width
    h2DPlot.GetYaxis().SetTitleSize(0.05)
    h2DPlot.GetYaxis().SetLabelSize(0.04)
    h2DPlot.GetYaxis().SetTitleOffset(1.1)
    h2DPlot.GetZaxis().SetTitleSize(0.05)
    h2DPlot.GetZaxis().SetLabelSize(0.04)
    h2DPlot.GetZaxis().SetTitleOffset(1.2)

    h2DPlot.GetZaxis().SetTitle(labelZ) 
    h2DPlot.Draw(drawOption)
    
    if (setXAxisRangeFromUser): h2DPlot.GetXaxis().SetRangeUser(xmin,xmax)
    if (setYAxisRangeFromUser): h2DPlot.GetYaxis().SetRangeUser(ymin,ymax)
    if (setZAxisRangeFromUser): h2DPlot.GetZaxis().SetRangeUser(zmin,zmax)


    h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.4)


    h2DProfile = 0
    if drawProfileX:
        h2DProfile = h2D.ProfileX("%s_pfx" %h2D.GetName())
        h2DProfile.SetMarkerColor(ROOT.kBlack)
        h2DProfile.SetMarkerStyle(20)
        h2DProfile.SetMarkerSize(1)
        h2DProfile.Draw("EPsame")
        
    # not yet implemented
    setTDRStyle()
    if not skipLumi and not plotLabel == "ForceTitle": 
        if lumi != None: CMS_lumi_new(canvas,lumi,True,False)
        else:            CMS_lumi_new(canvas,"",True,False)

    if plotLabel == "ForceTitle":
        ROOT.gStyle.SetOptTitle(1)        

    #h2DPlot.GetZaxis().SetMaxDigits(1)  #for N>99, should use scientific notation, I'd like to make it work only with negative exponential but haven't succeeded yet
    # canvas.Modified()
    # canvas.Update()

    leg = ROOT.TLegend(0.39,0.75,0.89,0.95)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)
    if plotLabel not in ["", "ForceTitle"]: leg.AddEntry(0,plotLabel,"")
    if drawProfileX: leg.AddEntry(0,"Correlation = %.2f" % h2DPlot.GetCorrelationFactor(),"")
    leg.Draw("same")

    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 1):
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        
    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 2):
        canvas.SetLogz()
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_logZ.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        canvas.SetLogz(0)
