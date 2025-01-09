#!/usr/bin/env python3 

### python specific import
import os
import pickle
import shutil
from multiprocessing import Pool
import datetime
import copy
import argparse
import re
from array import array

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gInterpreter.ProcessLine(".O3")

from libPython.tnpClassUtils import tnpSample
from libPython.plotUtils import compileMacro, testBinning, safeGetObject, safeOpenFile, createPlotDirAndCopyPhp, compileFileMerger
from libPython.checkFitStatus import checkFit

compileMacro("libCpp/RooCBExGaussShape.cc")
compileMacro("libCpp/RooCMSShape.cc")
compileMacro("libCpp/histFitter.C")
#compileFileMerger("libCpp/FileMerger.C")
compileMacro("libCpp/FileMerger.C")

### tnp library
import libPython.binUtils  as tnpBiner
import libPython.rootUtils as tnpRoot
import libPython.fitUtils as fitUtils
        
parser = argparse.ArgumentParser()
parser.add_argument('--flag'       , default = None       , help ='WP to test')
parser.add_argument('--inputMC'    , type=str, default = '', help = "MC input file which contains 3d histograms")
parser.add_argument('--inputData'  , type=str, default = '', help = "Data input file which contains 3d histograms")
parser.add_argument('--inputBkg'   , type=str, default = '', help = "Background input file which contains 3d histograms")
parser.add_argument('--era'        , type=str, default = '', choices=["BtoF", "GtoH"], help ='era to perform tnp fits for')
parser.add_argument('--checkBins'  , action='store_true'  ,  help = 'check  bining definition')
parser.add_argument('--createBins' , action='store_true'  ,  help = 'create bining definition')
parser.add_argument('--createHists', action='store_true'  ,  help = 'create histograms')
parser.add_argument('--sample'     , default='all'        ,  help = 'create histograms (per sample, expert only)')
parser.add_argument('--altSig'     , action='store_true'  ,  help = 'alternate signal model fit')
parser.add_argument('--altBkg'     , action='store_true'  ,  help = 'alternate background model fit')
parser.add_argument('--doFit'      , action='store_true'  ,  help = 'fit sample (sample should be defined in settings.py)')
parser.add_argument('--mcSig'      , action='store_true'  ,  help = 'fit MC nom [to init fit params]')
parser.add_argument('--mergeFiles'     , action='store_true'  ,  help = 'plotting')
parser.add_argument('--sumUp'      , action='store_true'  ,  help = 'sum up efficiencies')
parser.add_argument('--iBin'       , dest = 'binNumber'   , type=int,  default=-1, help='bin number (to refit individual bin)')
parser.add_argument('--outdir'     , type=str, default=None,
                    help="name of the output folder (if not passed, a default one is used, which has the time stamp in it)")
parser.add_argument('--useTrackerMuons', action='store_true'  , help = 'Measuring efficiencies specific for tracker muons (different tunings needed')

args = parser.parse_args()

if args.flag is None:
    print('[tnpEGM_fitter] flag is MANDATORY, this is the working point as defined in the settings.py')
    sys.exit(0)

## put most of the stuff here and make it configurable...
## AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
## ===========================================================================


massbins, massmin, massmax = 60, 60, 120

binning_mass = [round(massmin + (massmax - massmin)/massbins*i, 1) for i in range(massbins+1) ]

## define the binning here, much easier...
binning_eta = [round(-2.4+0.1*i,2) for i in range(49) ]
#binning_eta = [round(-2.4+0.4*i,2) for i in range(13) ]

binning_pt  = [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.]
#binning_pt  = [24., 28., 32., 36., 40., 47., 55., 65.]

typeflag = args.flag.split('_')[1]

print("typeflag = ",typeflag)

if typeflag == 'tracking':
    #binning_pt  = [15., 25.,35.,45.,55.,65.,80.]
    #massbins, massmin, massmax = 100, 40, 140
    #binning_pt  = [55., 65.]
    binning_pt  = [10., 15., 24., 35., 45., 55., 65.]  # [24., 65.]
    #massbins, massmin, massmax = 100, 50, 150
    massbins, massmin, massmax = 80, 50, 130
    binningDef = {
        'eta' : {'var' : 'eta', 'type': 'float', 'bins': binning_eta},
        'pt'  : {'var' : 'pt' , 'type': 'float', 'bins': binning_pt }
    }

elif typeflag == 'reco':
    #binning_pt   = [24., 65.]
    #massbins, massmin, massmax = 52, 68, 120
    massbins, massmin, massmax = 60, 60, 120
    if args.useTrackerMuons:
        binning_pt  = [24., 26., 30., 34., 38., 42., 46., 50., 55., 65.]
    else:
        #binning_pt  = [10., 15., 20., 24., 26., 30., 34., 38., 42., 46., 50., 55., 60., 65.]
        binning_pt  = [24., 26., 30., 34., 38., 42., 46., 50., 55., 60., 65.]
    binningDef = {
        'eta' : {'var' : 'eta', 'type': 'float', 'bins': binning_eta},
        'pt'  : {'var' : 'pt' , 'type': 'float', 'bins': binning_pt }
    }

elif typeflag == 'veto':
    binning_pt  = [10., 15., 20., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.]
    #binning_pt = [(15. + 5.*i) for i in range(11)]
    binningDef = {
        'eta' : {'var' : 'eta', 'type': 'float', 'bins': binning_eta},
        'pt'  : {'var' : 'pt' , 'type': 'float', 'bins': binning_pt }
    }

else:
    binningDef = {
        'eta' : {'var' : 'eta', 'type': 'float', 'bins': binning_eta},
        'pt'  : {'var' : 'pt' , 'type': 'float', 'bins': binning_pt }
    }


#############################################################
########## fitting params to tune fit by hand if necessary
#############################################################
if typeflag == 'tracking':

    bkgParFit = [
        "expalphaP[0.,-5.,5.]",
        "expalphaF[0.,-5.,5.]",
        "acmsF[60.,40.,130.]","betaF[5.,0.1,40.]","gammaF[0.1, 0, 1]","peakF[90.0]",
        "c1F[0.0,-1.0,1.0]","c2F[-0.5,-1.0,1.0]","c3F[0.0,-1.0,1.0]","c4F[-0.5,-1.0,1.0]",
        #"argMassF[100.,80.,120.]","argSlopeF[5.,0.,20.]","argPowF[1.0,0.0,3.0]"
    ]
    bkgShapes = [
        "Exponential::bkgPass(x, expalphaP)",
        #"RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        #"ArgusBG::bkgFail(x,argMassF,argSlopeF,argPowF)",
        "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F,c4F})",
        #"Bernstein::bkgFailBackup(x,{b0F[0.5,0,1.0],b1F[0.5,0,1.0],b2F[0.5,0,1.0],b3F[0.5,0,1.0],b4F[0.5,0,1.0]})",
        "Exponential::bkgFailMC(x, expalphaF)"
    ]
    
    tnpParNomFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.02,3.0]",
    ]

    # these might be partially overridden when running the fit to data by taking the values from the MC fit and narrowing the range in which they can float to help convergence
    tnpParAltSigFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]",'nP[3,0.01,5]',"sigmaP_2[1.5,0.5,6.0]",
        "meanF[-0.0,-12.0,12.0]","sigmaF[2,0.7,12.0]","alphaF[2.0,1.2,3.5]",'nF[3,0.01,5]',"sigmaF_2[2.0,0.5,6.0]",
    ]
    
    tnpParAltBkgFit = [
        "meanP[-0.0,-5.0,5.0]", "sigmaP[0.5,0.1,5.0]", "expalphaP[0.,-5.,5.]",
        "meanF[-0.0,-5.0,5.0]", "sigmaF[0.5,0.02,3.0]",
    ]

    # for pt >= 55 and tracking (se also note above)
    tnpParAltSigFitTrackingHighPt = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]",'nP[3,0,5]',"sigmaP_2[1.5,0.5,6.0]",
        "meanF[4.0,-1.0,15.0]","sigmaF[2,0.7,15.0]","alphaF[2.0,1.2,3.5]",'nF[3,0,5]',"sigmaF_2[2.0,0.5,6.0]",
    ]

    tnpParNomFit.extend(bkgParFit)
    tnpParAltSigFit.extend(bkgParFit)
    tnpParAltSigFitTrackingHighPt.extend(bkgParFit)
    
    if not args.mcSig:
        if args.useTrackerMuons:
            # for tracker muons the fraction is much larger
            tnpParNomFit.extend(["maxFracSigF[0.5]"])
            tnpParAltSigFit.extend(["maxFracSigF[0.5]"])
            tnpParAltSigFitTrackingHighPt.extend(["maxFracSigF[0.5]"])
        else:
            # use 0.2 if not using MergedStandAlone_nValidHits > 0
            tnpParNomFit.extend(["maxFracSigF[0.5]"])
            tnpParAltSigFit.extend(["maxFracSigF[0.5]"])
            tnpParAltSigFitTrackingHighPt.extend(["maxFracSigF[0.5]"])

    # ## Try to constrain some background parameters (for tracking might need to do it for signal instead, since S/B is small)
    parConstraints = [
        # Passing
        #"Gaussian::constrainP_acmsP(acmsP,90,50)",
        #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
        #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
        # failing
        "Gaussian::constrainF_acmsF(acmsF,90,50)",
        #"Gaussian::constrainF_betaF(betaF,0.05,0.25)",
        "Gaussian::constrainF_betaF(betaF,5.0,25.0)",
        "Gaussian::constrainF_gammaF(gammaF,0.5,0.8)",
    ]

            
elif typeflag == 'reco':

    bkgParFit = [
        "expalphaP[0.,-5.,5.]",
        "expalphaF[0.,-5.,5.]",
        "acmsF[60.,40.,130.]","betaF[5.,0.1,40.]","gammaF[0.1, 0, 1]","peakF[90.0]",
        "c1F[0.0,-1.0,1.0]","c2F[-0.5,-1.0,1.0]","c3F[0.0,-1.0,1.0]",
        #"b0F[0.1,0,100]","b1F[0.1,0,100]","b2F[0.1,0,100]","b3F[0.1,0,100]"
    ]
    bkgShapes = [
        "Exponential::bkgPass(x, expalphaP)",
        #"RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        #"Bernstein::bkgFailBackup(x,{b0F,b1F,b2F,b3F})",
        #"Exponential::bkgFailBackup(x, expalphaF)"
        "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F})",
        #"Chebychev::bkgFailBackup(x,{c1F[0.0,-1.0,1.0],c2F[-0.5,-1.0,1.0],c3F[0.0,-1.0,1.0],c4F[-0.5,-1.0,1.0]})",
        "Exponential::bkgFailMC(x, expalphaF)"
    ]

    tnpParNomFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,3.0]",
        "meanF[-0.0,-3.0,3.0]","sigmaF[0.5,0.01,2.0]",
    ]

    # was to tune few bins for reco, but currently used everywhere
    tnpParAltSigFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]",'nP[3,0.01,5]',"sigmaP_2[1.5,0.5,6.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[2,0.7,5.0]","alphaF[2.0,1.2,3.5]",'nF[3,0.1,5]',"sigmaF_2[2.0,0.5,6.0]",
    ]
    
    tnpParAltBkgFit = [
        "meanP[-0.0,-5.0,5.0]", "sigmaP[0.5,0.1,3.0]", "expalphaP[0.,-5.,5.]",
        "meanF[-0.0,-3.0,3.0]", "sigmaF[0.5,0.01,2.0]",
    ]

    tnpParNomFit.extend(bkgParFit)
    tnpParAltSigFit.extend(bkgParFit)

    if not args.mcSig and args.useTrackerMuons:
        # for tracker muons
        tnpParNomFit.extend(["maxFracSigF[0.05]"] if args.binNumber in [24] else ["maxFracSigF[0.05]"] if args.binNumber in [50, 69, 79, 133, 420] else ["maxFracSigF[0.3]"])
        tnpParAltSigFit.extend(["maxFracSigF[0.1]"] if args.binNumber in [29] else ["maxFracSigF[0.3]"])

    # ## Try to constrain some background parameters (for tracking might need to do it for signal instead, since S/B is small)
    parConstraints = [
        # Passing
        #"Gaussian::constrainP_acmsP(acmsP,90,50)",
        #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
        #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
        # failing
        "Gaussian::constrainF_acmsF(acmsF,90,50)",
        "Gaussian::constrainF_betaF(betaF,5.0,25.0)",
        "Gaussian::constrainF_gammaF(gammaF,0.5,0.8)",
    ]

else:

    bkgParFit = [
        "expalphaP[0.,-5.,5.]",
        "expalphaF[0.,-5.,5.]",
        "acmsF[60.,40.,130.]","betaF[5.,0.1,40.]","gammaF[0.1, 0, 1]","peakF[90.0]",
        "c1F[0.0,-1.0,1.0]","c2F[-0.5,-1.0,1.0]","c3F[0.0,-1.0,1.0]",
        #"argMassF[90.,60.,120.]","argSlopeF[5.,0.,20.]","argPowF[1.0,0.0,3.0]"
    ]
    bkgShapes = [
        "Exponential::bkgPass(x, expalphaP)",
        "Exponential::bkgFail(x, expalphaF)",
        "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F})",
        #"ArgusBG::bkgFailBackup(x,argMassF,argSlopeF,argPowF)",
     ]

    tnpParNomFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
    ]
    
    tnpParAltBkgFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]", "expalphaP[0.,-5.,5.]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
    ]
    
    tnpParAltSigFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]",'nP[3,0.01,5]',"sigmaP_2[1.5,0.5,6.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[2,0.7,15.0]","alphaF[2.0,1.2,3.5]",'nF[3,0.01,5]',"sigmaF_2[2.0,0.5,6.0]",
    ]
    tnpParNomFit.extend(bkgParFit)
    tnpParAltSigFit.extend(bkgParFit)
    
    parConstraints = []
#####

# add second gaussian at low mass around 70 to model FSR bump for working points with isolation
flagsWithFSR = ["iso", "trigger", "isonotrig"]
if any(x in typeflag for x in flagsWithFSR):
    fsrGauss = ["fsrMeanF[70.0,65.0,80.0]", "fsrSigmaF[1.0,1.2,5.0]"]
    tnpParAltSigFit.extend(fsrGauss)

########################



if args.outdir:
    baseOutDir = f'{args.outdir}/efficiencies_{args.era}/'
else:
    baseOutDir = f'plots/results_{datetime.date.isoformat(datetime.date.today())}/efficiencies_{args.era}/'

outputDirectory = f'{baseOutDir}/{args.flag}'

print('===>  Output directory: ', outputDirectory)
print()

luminosity = 16.8 if args.era == "GtoH" else 19.5
eraMC = "postVFP" if args.era == "GtoH" else "preVFP"
dataName = f"mu_Run{args.era}"
mcName = f"mu_DY_{eraMC}"
bkg_name = f"mu_mcBkg_{eraMC}"

samples_data = tnpSample(dataName, args.inputData, f"{outputDirectory}/{dataName}_{args.flag}.root", False)
samples_dy   = tnpSample(mcName,   args.inputMC,   f"{outputDirectory}/{mcName}_{args.flag}.root",   True)
samples_bkg  = tnpSample(bkg_name, args.inputBkg,  f"{outputDirectory}/{bkg_name}_{args.flag}.root", True)


## check binning in histogram and consistency with settings above
## FIXME: should be done for each step, but histograms are not always passed
if args.createHists:
    ftest = ROOT.TFile(args.inputData, "READ")
    htest = ftest.Get(f"pass_{dataName}")
    this_binning_mass = [round(htest.GetXaxis().GetBinLowEdge(i), 1) for i in range(1, htest.GetNbinsX()+2) ]
    resTestMass = testBinning(binning_mass, this_binning_mass, "mass", typeflag, allowRebin=True)
    this_binning_pt = [round(htest.GetYaxis().GetBinLowEdge(i), 1) for i in range(1, htest.GetNbinsY()+2) ]
    resTestPt = testBinning(binning_pt, this_binning_pt, "pt", typeflag, allowRebin=True)
    this_binning_eta = [round(htest.GetZaxis().GetBinLowEdge(i), 1) for i in range(1, htest.GetNbinsZ()+2) ]
    resTestEta = testBinning(binning_eta, this_binning_eta, "eta", typeflag, allowRebin=True)
    ftest.Close()
    if resTestMass or resTestPt or resTestEta:
        exit(-1)


####################################################################
##### Create (check) Bins
####################################################################
if args.checkBins:
    print(">>> Check bins")
    tnpBins = tnpBiner.createBins(binningDef, None)
    for ib in range(len(tnpBins['bins'])):
        print(tnpBins['bins'][ib]['name'])
        print('  - cut: ',tnpBins['bins'][ib]['cut'])
    sys.exit(0)

if args.createBins:
    print(">>> Create bins")
    if os.path.exists( outputDirectory ):
        shutil.rmtree( outputDirectory )
    createPlotDirAndCopyPhp(outputDirectory)
    tnpBins = tnpBiner.createBins(binningDef, None)
    pickle.dump( tnpBins, open( f'{outputDirectory}/bining.pkl', 'wb') )
    print(f'Created dir: {outputDirectory} ')
    print(f'Bining created successfully... ')
    print(f'Note than any additional call to createBins will overwrite directory {outputDirectory}')
    sys.exit(0)

tnpBins = pickle.load( open( f'{outputDirectory}/bining.pkl', 'rb') )


####################################################################
##### Create Histograms
####################################################################

samplesDef = {
    'data'   : samples_data,
    'mcNom'  : samples_dy,
    'mcBkg'  : samples_bkg,
    #'tagSel' : None,
}

samplesDef["data"].printConfig()

if args.createHists:
    print()
    print(">>> Create histograms")
    def parallel_hists(sampleType):
        sample = samplesDef[sampleType]
        if sample is not None and (sampleType == args.sample or args.sample == 'all'):
            print('Creating histogram for sample', sample.getName())
            sample.printConfig()
            if typeflag == 'tracking':
                var = { 'namePassing' : 'pair_mass', 'nameFailing' : 'pair_massStandalone', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            else:
                var = { 'name' : 'pair_mass', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            tnpRoot.makePassFailHistograms(sample, tnpBins['bins'], binningDef, var)

    pool = Pool()
    pool.map(parallel_hists, samplesDef.keys())
    sys.exit(0)


####################################################################
##### Actual Fitter
####################################################################


sampleMC = samplesDef['mcNom']
if sampleMC is None:
    print('[tnpEGM_fitter, prelim checks]: MC sample not available... check your settings')
    sys.exit(1)

for s in samplesDef.keys():
    sample =  samplesDef[s]
    if sample is None: continue
    setattr( sample, 'mcRef'     , sampleMC )
    setattr( sample, 'nominalFit', '%s/%s_%s_nominalFit.root' % ( outputDirectory , sample.getName(), args.flag ) )
    setattr( sample, 'altSigFit' , '%s/%s_%s_altSigFit.root'  % ( outputDirectory , sample.getName(), args.flag ) )
    setattr( sample, 'altBkgFit' , '%s/%s_%s_altBkgFit.root'  % ( outputDirectory , sample.getName(), args.flag ) )

sampleToFit = samplesDef['data'] if not args.mcSig else sampleMC

if sampleToFit is None:
    print('[tnpEGM_fitter, prelim checks]: data sample not available... check your settings')
    sys.exit(1)

if args.altSig:
    fileName = sampleToFit.altSigFit
    fitType  = 'altSigFit'
elif args.altBkg:
    fileName = sampleToFit.altBkgFit
    fitType  = 'altBkgFit'
else:
    fileName = sampleToFit.nominalFit
    fitType  = 'nominalFit'


plottingDir = '%s/plots/%s/%s' % (outputDirectory,sampleToFit.getName(),fitType)
createPlotDirAndCopyPhp(plottingDir)

    
if args.doFit:
    print()
    print(">>> Running fits")

    #useAllTemplateForFail = True if typeflag not in flagsWithFSR else False # use all probes to build MC template for failing probes when fitting data nominal
    #maxFailIntegralToUseAllProbe = 300 if typeflag not in ["tracking"] else -1 # use all probes for the failing template only when stat is very small, otherwise sometimes the fit doesn't work well
    useAllTemplateForFail = False
    maxFailIntegralToUseAllProbe = -1
    altSignalFail = True if typeflag in ["tracking", "reco", "veto"] else False # use Gaussian as resolution function for altSig model
    modelFSR = True if typeflag in flagsWithFSR else False # add Gaussian to model low mass bump from FSR, in altSig fit
    def parallel_fit(ib): ## parallel
        #print("tnpBins['bins'][ib] = ",tnpBins['bins'][ib])
        if (args.binNumber >= 0 and ib == args.binNumber) or args.binNumber < 0:
            
            if not (args.altSig or args.altBkg):
                fitUtils.histFitterNominal(sampleToFit, tnpBins['bins'][ib], tnpParNomFit, massbins, massmin, massmax,
                                            useAllTemplateForFail, maxFailIntegralToUseAllProbe, 
                                            constrainPars=parConstraints, bkgShapes=bkgShapes)
            elif args.altSig:
                if typeflag == 'tracking':
                    # constrainSignalFailFromMC sets the data fit parameters to MC value +/- 3*uncertainty
                    if fitUtils.ptMin(tnpBins['bins'][ib]) > 54.0: # force peak mean more on the right for high pt bins and tracking efficiency
                        fitUtils.histFitterAltSig(sampleToFit, tnpBins['bins'][ib], tnpParAltSigFitTrackingHighPt, massbins, massmin, massmax,
                                                  altSignalFail=altSignalFail, modelFSR=False, constrainSignalFailFromMC=False, constrainPars=parConstraints, bkgShapes=bkgShapes)
                    else:
                        fitUtils.histFitterAltSig(sampleToFit, tnpBins['bins'][ib], tnpParAltSigFit, massbins, massmin, massmax,
                                                  altSignalFail=altSignalFail, modelFSR=False, constrainSignalFailFromMC=False, constrainPars=parConstraints, bkgShapes=bkgShapes)
                elif typeflag == 'reco': 
                        fitUtils.histFitterAltSig(sampleToFit, tnpBins['bins'][ib], tnpParAltSigFit, massbins, massmin, massmax,
                                                  altSignalFail=altSignalFail, modelFSR=False, constrainPars=parConstraints, bkgShapes=bkgShapes)
                else:
                    #fitUtils.histFitterAltSig(sampleToFit, tnpBins['bins'][ib], tnpParAltSigFit, massbins, massmin, massmax,
                    #                          altSignalFail=altSignalFail, modelFSR=modelFSR, zeroBackground=True)
                    fitUtils.histFitterAltSig(sampleToFit, tnpBins['bins'][ib], tnpParAltSigFit, massbins, massmin, massmax,
                                              altSignalFail=altSignalFail, modelFSR=modelFSR, constrainPars=parConstraints, bkgShapes=bkgShapes)
            elif args.altBkg and not args.mcSig:
                fitUtils.histFitterAltBkg(sampleToFit, tnpBins['bins'][ib], tnpParAltBkgFit, massbins, massmin, massmax,
                                                      useAllTemplateForFail, maxFailIntegralToUseAllProbe, constrainPars=parConstraints, bkgShapes=bkgShapes)
            else:
                # This is the case for fitting on MC with altSig/altBkg models. Not used for now
                pass

    pool = Pool() ## parallel
    pool.map(parallel_fit, range(len(tnpBins['bins']))) ## parallel
    #parallel_fit(0)
    args.mergeFiles = True

####################################################################
##### Merging the files with the results into a single one
####################################################################
if args.mergeFiles:

    print()
    print(">>> Merging root files...")
    print(f"Output: {fileName}")
    fileNameNoExt = fileName.rstrip(".root")
    if args.binNumber >= 0:
        thisbin = tnpBins['bins'][args.binNumber]['name']
        rootfileBin = safeOpenFile(f"{fileNameNoExt}_bin_{thisbin}.root")
        rootfile = safeOpenFile(f"{fileName}", mode="UPDATE")
        for k in rootfileBin.GetListOfKeys():
            obj = safeGetObject(rootfileBin, k.GetName(), detach=False)
            rootfile.cd()
            obj.Write(k.GetName(), ROOT.TObject.kOverwrite) # write in merged root file overwriting keys if they already existed
        rootfile.Close()
        rootfileBin.Close()
        os.system(f"rm {fileNameNoExt}_bin_bin*.root")
    else:
        numberOfBins = len(tnpBins['bins'])
        regexp = re.compile(f"{os.path.basename(fileNameNoExt)}_bin_bin.*.root")
        outpath = os.path.dirname(fileName) + "/"
        vec = ROOT.std.vector["std::string"]()
        for fname in os.listdir(outpath):
            if os.path.isfile(os.path.join(outpath, fname)) and regexp.match(fname):
                vec.push_back(os.path.join(outpath, fname))
        if vec.size() != numberOfBins:
            print(f"Error: number of files ({vec.size()}) does not coincide with number of bins ({numberOfBins})")
            vec.clear()
            os.system("sleep 3")
            print("Deleting temporary files and exiting")
            os.system(f"rm {fileNameNoExt}_bin_bin*.root")
            quit()
        ROOT.FileMerger(numberOfBins, fileName, vec)
        vec.clear()
        os.system("sleep 3")
        os.system(f"rm {fileNameNoExt}_bin_bin*.root")
    print("Done with merging :-)")
    print()


####################################################################
##### Producing the sumUp .txt files and plots
####################################################################
if args.sumUp:
    from pprint import pprint

    info = {
        'Data_Nominal'  : sampleToFit.nominalFit,
        'Data_Alt_Sig'  : sampleToFit.altSigFit ,
        'Data_Alt_Bkg'  : sampleToFit.altBkgFit ,
        'MC_Nominal'    : sampleToFit.mcRef.getOutputPath(),
        'MC_Nominal_fit': sampleToFit.mcRef.nominalFit,
        'MC_Alt_Sig'    : sampleToFit.mcRef.altSigFit, # not used
        'MC_Alt_Bkg'    : sampleToFit.mcRef.altBkgFit, # not used
        'tagSel'      : None
        }

    if 'MC_Alt_Sig' in samplesDef.keys() and not (samplesDef['MC_Alt_Sig'] is None):

        info['MC_Alt_Sig'] = samplesDef["MC_Alt_Sig"].getOutputPath()
    if 'tagSel' in samplesDef.keys() and not (samplesDef['tagSel'] is None):
        info['tagSel'] = samplesDef['tagSel'].getOutputPath()

    effFileName = outputDirectory+'/allEfficiencies.txt'

    # security check, if the code crashes the temporary files are still present, let's remove them before executing parallel_sumUp
    if any ("_tmpTMP_" in f for f in os.listdir(outputDirectory)):
        os.system(f"rm {effFileName}_tmpTMP_*")
            

    def parallel_sumUp(_bin):
        effis = tnpRoot.getAllEffi(info, _bin, outputDirectory, saveCanvas=False)
        v1Range = _bin['title'].split(';')[1].split('<')
        v2Range = _bin['title'].split(';')[2].split('<')

        ib = int(_bin['name'].split('_')[0].replace('bin',''))

        fOut = open(effFileName+'_tmpTMP_'+str(ib), 'w')

        if not ib:
            fOut.write( '### var1 : %s\n' % v1Range[1] )
            fOut.write( '### var2 : %s\n' % v2Range[1] )
            exp = ''
            for v in ['var1min', 'var1max', 'var2min', 'var2max']:
                exp += f'{v:8s}\t'
            for v in ['eff data', 'err data', 'eff mc', 'err mc']:
                exp += f'{v:10s}\t'
            for v in ['effD altS', 'errD altS', 'effD altB', 'errD altB', 'eff mc altS', 'err mc altS', 'eff tag sel']:
                exp += f'{v:12s}\t'

            fOut.write(exp + '\n')

        vals = ''
        for r in [v1Range, v2Range]:
            vals += f'{float(r[0]):<+8.3f}\t{float(r[2]):<+8.3f}\t'
        for v in ['Data_Nominal', 'MC_Nominal']:
            vals += f'{effis[v][0]:<+10.6f}\t{effis[v][1]:<+10.6f}\t'
        for v in ['Data_Alt_Sig', 'Data_Alt_Bkg', 'MC_Alt_Sig', 'tagSel']:
            vals += f'{effis[v][0]:<+12.6f}\t{effis[v][1]:<+12.6f}\t'

        fOut.write( vals + '\n' )
        fOut.close()
        effis = {}
        
    #odllevel = ROOT.gErrorIgnoreLevel
    #ROOT.gErrorIgnoreLevel = ROOT.kWarning
    pool = Pool()
    pool.map(parallel_sumUp, tnpBins['bins'])

    lsfiles = []
    alltmpfiles = os.listdir(outputDirectory)
    for ifile in alltmpfiles:
        if not 'tmpTMP' in ifile: continue
        lsfiles.append(outputDirectory+'/'+ifile)

    lsfiles = sorted(lsfiles, key = lambda x: int(x.split('_')[-1]))

    os.system('cat '+' '.join(lsfiles)+' > '+effFileName)
    os.system('rm  '+' '.join(lsfiles))

    print('Efficiencies saved in file : ',  effFileName)

    outputDirectoryPlots = f"{outputDirectory}/plots/"
    createPlotDirAndCopyPhp(outputDirectoryPlots)
    #fOut.close()
    outputDirectoryTH2 = f"{outputDirectoryPlots}/histo2D/"
    createPlotDirAndCopyPhp(outputDirectoryTH2)

    # Everything works, but the syntax of the various functions is awful
    # Will change it in a second time
    import libPython.EGammaID_scaleFactors as makesf
    makesf.doSFs(effFileName,luminosity,['pt', 'eta'], outputDirectoryTH2)

    # plotting sanity-check plots on fitresults
    print()
    print("Plotting sanity-check histograms from fit results in this file")
    print(f">>> {fileName}")
    basePath = os.path.dirname(fileName)
    outdirCheckPlots = basePath + "/plots/checkFitStatus/"
    createPlotDirAndCopyPhp(outdirCheckPlots)
    # get histogram to read binning, it is passed to checkFit
    #rootfileWithEffi = safeOpenFile(basePath + "/allEfficiencies_2D.root")
    #th2ForBinning = safeGetObject(rootfileWithEffi, "SF2D_nominal", detach=True)
    #rootfileWithEffi.Close()
    th2ForBinning = ROOT.TH2D("th2ForBinning", "", len(binning_eta)-1, array("d", binning_eta), len(binning_pt)-1, array("d", binning_pt))
    print(f"Eta binning: {binning_eta}")
    print(f"Pt binning : {binning_pt}")
    tag = "MC" if "_DY_" in fileName else "Data"
    fitName = f"Eff{tag}_" + fitType
    checkFit(fileName, outdirCheckPlots, fitName, th2ForBinning)
    print("================================================")
    print("THE END!")
    print("================================================")
