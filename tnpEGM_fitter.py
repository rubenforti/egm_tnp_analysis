#!/usr/bin/env python3 

import os
import sys
import re
import argparse
import datetime
import pickle
import shutil
from array import array
from multiprocessing import Pool
import ROOT

## Safe batch mode
args = sys.argv[:]
sys.argv = ['-b']
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gInterpreter.ProcessLine(".O3")


## TnP libraries and functions
import libPython.fitUtils as tnpFit
import libPython.MuonWmass_scaleFactors as makesf
from libPython.binUtils import createBins, testBinning, binMinPt
from libPython.checkFitStatus import checkFit
from libPython.histUtils import makePassFailHistograms
from libPython.plotUtils import createPlotDirAndCopyPhp
from libPython.rootUtils import compileMacro, safeGetObject, safeOpenFile, getAllEffi, plotAllEffi
from libPython.tnpClassUtils import tnpSample


## Import fitting settings and strategies
from config.fit_settings import fitBinning, fitParsAndShapes
from config.fitting_strategies import fitStrategies

## Compile C++ libraries
compileMacro("libCpp/RooCBExGaussShape.cc")
compileMacro("libCpp/RooCMSShape.cc")
compileMacro("libCpp/histFitter.C")
compileMacro("libCpp/FileMerger.C")

#################################################
#### General settings for the fit, hardcoded here



#################################################

## Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('--flag'       , type=str, default='', help = 'WP to test')
parser.add_argument('--inputMC'    , type=str, default='', help = 'MC input file which contains 3d histograms')
parser.add_argument('--inputData'  , type=str, default='', help = 'Data input file which contains 3d histograms')
parser.add_argument('--inputBkg'   , type=str, default='', help = 'Background input file which contains 3d histograms')
parser.add_argument('--era'        , type=str, default='', choices=["BtoF", "GtoH"], 
                    help = 'Era to perform tnp fits for')
parser.add_argument('--year'       , type=str, default='2016', choices=["2016", "2017", "2018"],
                    help = 'Year of data taking')
parser.add_argument('--createBins' , action='store_true' ,  help = 'Create binning definition')
parser.add_argument('--checkBins'  , action='store_true' ,  help = 'Check binning definition')
parser.add_argument('--createHists', action='store_true' ,  help = 'Create histograms')
parser.add_argument('--doFit'      , action='store_true' ,  help = 'Fit sample')
parser.add_argument('--sample'     , default='all'       ,  help = 'Create histograms (per sample, expert only)')
parser.add_argument('--altSig'     , action='store_true' ,  help = 'Perform "alternate-signal" fit')
parser.add_argument('--altBkg'     , action='store_true' ,  help = 'Perform "alternate-background" fit')
parser.add_argument('--mcSig'      , action='store_true' ,  help = 'Fit MC sample')
parser.add_argument('--mergeFiles' , action='store_true' ,  help = 'Merge files of individual fitted bins into one')
parser.add_argument('--sumUp'      , action='store_true' ,  help = 'Produce the sum-up .txt file and plots')
parser.add_argument('--outdir'     , type=str, default=None,
                    help="name of the output folder (if not passed, a default one is used, which has the time stamp in it)")
parser.add_argument('--iBin'       , dest = 'binNumber' , type=int,  default=-1, help='bin number (to refit individual bin)')
parser.add_argument('--useTrackerMuons', action='store_true', help = 'Measuring efficiencies specific for tracker muons (different tunings needed)')
parser.add_argument('--useVetoBins', action='store_true', help = 'Use veto bins for tracking and reco efficiencies')
args = parser.parse_args()

if args.flag=="":
    print('[tnpEGM_fitter] flag argument is MANDATORY')
    sys.exit(0)

typeflag = args.flag.split('_')[1]
print("typeflag = ", typeflag)


###############################################################################
##### General configuration steps
###############################################################################

## Binning
binningDef_mass, binningDef = fitBinning(typeflag)
massbins, massmin, massmax = binningDef_mass['nbins'], binningDef_mass['min'], binningDef_mass['max']

if args.useVetoBins: 
    if typeflag == "tracking":
        binningDef["pt"]["bins"] = [10., 15.] + binningDef["pt"]["bins"] 
    else:
        binningDef["pt"]["bins"] = [10., 15., 20.] + binningDef["pt"]["bins"]

if args.useTrackerMuons:
    pass

binning_mass = [round(massmin + i*(massmax-massmin)/massbins, 1) for i in range(massbins+1)]
binning_eta = binningDef["eta"]["bins"]
binning_pt = binningDef["pt"]["bins"]


## Output directory
outputDirectory = args.outdir if args.outdir else f"plots/results_{datetime.date.isoformat(datetime.date.today())}/"
outputDirectory += f"efficiencies_{args.era}/{args.flag}/"
print('===>  Output directory: ', outputDirectory)
print()

## Names and sample definitions
if args.year == "2016":
    luminosity = 16.8 if args.era == "GtoH" else 19.5
    eraMC = "postVFP" if args.era == "GtoH" else "preVFP"
    dataName = f"mu_Run{args.era}"
else:
    luminosity = 59.8 if args.year=="2018" else 38.0  # 2017 is not used
    eraMC = args.year
    dataName = f"mu_{args.year}"
mcName = f"mu_DY_{eraMC}"
bkg_name = f"mu_mcBkg_{eraMC}"

samples_data = tnpSample(dataName, args.inputData, f"{outputDirectory}/{dataName}_{args.flag}.root", False)
samples_dy   = tnpSample(mcName,   args.inputMC,   f"{outputDirectory}/{mcName}_{args.flag}.root",   True)
samples_bkg  = tnpSample(bkg_name, args.inputBkg,  f"{outputDirectory}/{bkg_name}_{args.flag}.root", True)

samplesDef = {
    'data'   : samples_data,
    'mcNom'  : samples_dy,
    'mcAltSig' : None,
    'mcBkg'  : samples_bkg if fitStrategies[typeflag]['AltBkg'] else None
    #'tagSel' : None,
}

if samplesDef['mcNom'] is None:
    print('[tnpEGM_fitter, prelim checks]: MC sample not available... check your settings')
    sys.exit(1)

for sample in samplesDef.values():
    if sample is None:
        continue
    setattr(sample, 'mcRef',  samplesDef['mcNom'] )
    setattr(sample, 'bkgRef', samplesDef['mcBkg'] if fitStrategies[typeflag]['AltBkg'] else None)
    setattr(sample, 'nominalFit', f"{outputDirectory}/{sample.getName()}_nominalFit.root")
    setattr(sample, 'altSigFit',  f"{outputDirectory}/{sample.getName()}_altSigFit.root")
    setattr(sample, 'altBkgFit',  f"{outputDirectory}/{sample.getName()}_altBkgFit.root")

sampleToFit = samplesDef['data'] if not args.mcSig else samplesDef['mcNom']

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


###############################################################################
##### Create (check) Bins
###############################################################################
if args.checkBins:
    print()
    print(">>> Check bins")
    tnpBins = createBins(binningDef, None)
    for ib in range(len(tnpBins['bins'])):
        print(tnpBins['bins'][ib]['name'])
        print('  - cut: ', tnpBins['bins'][ib]['cut'])
        print('')
    sys.exit(0)

## Check if the binning is consistent with the one in the histograms
if args.createHists:
    print()
    print(">>> Check bin consistency with histograms")
    ftest = ROOT.TFile(args.inputData, "READ")
    htest = ftest.Get(f"pass_{dataName}")
    this_binning_mass = [round(htest.GetXaxis().GetBinLowEdge(i), 1) for i in range(1, htest.GetNbinsX()+2) ]
    resTestMass = testBinning(binning_mass, this_binning_mass, "mass", typeflag, allowRebin=True)
    this_binning_pt = [round(htest.GetYaxis().GetBinLowEdge(i), 1) for i in range(1, htest.GetNbinsY()+2) ]
    resTestPt = testBinning(binning_pt, this_binning_pt, "pt", typeflag, allowRebin=True)
    this_binning_eta = [round(htest.GetZaxis().GetBinLowEdge(i), 1) for i in range(1, htest.GetNbinsZ()+2) ]
    resTestEta = testBinning(binning_eta, this_binning_eta, "eta", typeflag, allowRebin=True)
    ftest.Close()
    if any(res<0 for  res in [resTestMass, resTestPt, resTestEta]):
        sys.exit(0)

if args.createBins:
    print(">>> Create bins")
    if os.path.exists( outputDirectory ):
        shutil.rmtree( outputDirectory )
    createPlotDirAndCopyPhp(outputDirectory)
    tnpBins = createBins(binningDef, None)
    pickle.dump(tnpBins, open(f'{outputDirectory}/binning.pkl', 'wb') )
    print(f'Created dir: {outputDirectory} ')
    print(f'Bining created successfully... ')
    print(f'Note than any additional call to "createBins" will overwrite directory {outputDirectory}')
    sys.exit(0)

with open(f'{outputDirectory}/binning.pkl', 'rb') as f:
    tnpBins = pickle.load(f)


###############################################################################
##### Create Histograms
###############################################################################
if args.createHists:
    print()
    print(">>> Create histograms")
    def parallel_hists(sampleType):
        sample = samplesDef[sampleType]
        if sample is not None and (sampleType==args.sample or args.sample=='all'):
            print('Creating histogram for sample', sample.getName())
            sample.printConfig()
            if typeflag == 'tracking':
                var = { 'namePassing' : 'pair_mass', 'nameFailing' : 'pair_massStandalone', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            else:
                var = { 'name' : 'pair_mass', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            makePassFailHistograms(sample, tnpBins['bins'], binningDef, var)

    pool = Pool()
    pool.map(parallel_hists, samplesDef.keys())
    sys.exit(0)


####################################################################
##### Actual Fitter
####################################################################

## Create plot directory
plottingDir = f"{outputDirectory}/plots/{sampleToFit.getName()}/{fitType}"
createPlotDirAndCopyPhp(plottingDir)

## Load the shapes and parameters
ps = fitParsAndShapes(typeflag)

## General fit settings
symmConvSigFail = False  # use Gaussian as resolution function for altSig model
modelFSR = True if typeflag in ["iso", "trigger", "isonotrig"] else False # add Gaussian to model low mass bump from FSR (for fits with analytic signal model)
useBBFail = True  # use Barlow-Beeston method for the bkg template
doPrefit = "sig-both"  # do prefit on sig/bkg, when it is analytic
constrainMode = ""  # constrain the parameters of the analytic model to the prefit values (see "_applyPrefitOnPar" in the fitterHandler class for details)

fitterNominal = getattr(tnpFit, fitStrategies[typeflag]['Nominal'])
fitterAltSig  = getattr(tnpFit, fitStrategies[typeflag]['AltSig'])
fitterAltBkg  = getattr(tnpFit, fitStrategies[typeflag]['AltBkg']) if fitStrategies[typeflag]['AltBkg'] else None

if fitterNominal is None:
    sys.exit("[tnpEGM_fitter]: No fitting strategy defined for nominal fit")

if fitterAltSig is None and fitterAltBkg is None:
    sys.exit("[tnpEGM_fitter]: No fitting strategy defined for alternative fit")

if args.doFit:
    print()
    print(">>> Running fits")

    def parallel_fit(ib):
        
        if not ((args.binNumber>=0 and ib==args.binNumber) or (args.binNumber<0)): return
            
        if not (args.altSig or args.altBkg):
            fitterNominal(sampleToFit, tnpBins['bins'][ib], typeflag, "Nominal", massbins, massmin, massmax,
                          ps[f"tnpParNominal"], ps[f"tnpShapesNominal"], ps[f"parConstraints"],
                          isMC=args.mcSig, useBBFail=useBBFail,
                          modelFSR=modelFSR, symmConvSigFail=symmConvSigFail,
                          doPrefit=doPrefit, constrainMode=constrainMode)
        
        elif args.altSig and fitterAltSig:
            key_pars_altSig = f"tnpParAltSig_trackingHighPt" if binMinPt(tnpBins['bins'][ib])>54.0 and typeflag=="tracking" \
                              else f"tnpParAltSig" # force peak mean more on the right for high pt bins and tracking efficiency
            fitterAltSig(sampleToFit, tnpBins['bins'][ib], typeflag, "AltSig", massbins, massmin, massmax,
                         ps[key_pars_altSig], ps[f"tnpShapesAltSig"], ps[f"parConstraints"], 
                         isMC=args.mcSig, useBBFail=useBBFail, 
                         modelFSR=modelFSR, symmConvSigFail=symmConvSigFail,
                         doPrefit=doPrefit, constrainMode=constrainMode)
        
        elif args.altBkg and fitterAltBkg:
            fitterAltBkg(sampleToFit, tnpBins['bins'][ib], typeflag, "AltBkg", massbins, massmin, massmax,
                         ps[f"tnpParAltBkg"], ps[f"tnpShapesAltBkg"], ps["parConstraints"],
                         isMC=args.mcSig, useBBFail=useBBFail,
                         modelFSR=modelFSR, symmConvSigFail=symmConvSigFail,
                         doPrefit=doPrefit, constrainMode=constrainMode)
        
        else:
            print(f"[tnpEGM_fitter]: No fitting strategy defined for alternative fit")
            return
    
    #parallel_fit(166)
    
    pool = Pool() ## parallel
    pool.map(parallel_fit, range(len(tnpBins['bins']))) ## parallel
    args.mergeFiles = True
    

####################################################################
##### Merging the files with the results into a single one
####################################################################
if args.mergeFiles:

    print()
    print(">>> Merging root files...")
    print(f"Output: {fileName}")
    fileNameNoExt = fileName.replace(".root", "")
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
    print()
    print(">>> Summing up...")

    info = {
        'Data_Nominal'  : sampleToFit.nominalFit,
        'Data_Alt_Sig'  : sampleToFit.altSigFit ,
        'Data_Alt_Bkg'  : sampleToFit.altBkgFit ,
        'MC_Nominal'    : sampleToFit.mcRef.getOutputPath(),
        'MC_Nominal_fit': sampleToFit.mcRef.nominalFit,
        'MC_Alt_Sig'    : sampleToFit.mcRef.altSigFit,
        'MC_Alt_Bkg'    : sampleToFit.mcRef.altBkgFit,
        'tagSel'        : None
        }

    effFileName = outputDirectory+'/allEfficiencies.txt'

    outputDirectoryPlots = f"{outputDirectory}/plots/"
    createPlotDirAndCopyPhp(outputDirectoryPlots)

    # security check, if the code crashes the temporary files are still present, let's remove them before executing parallel_sumUp
    if any("_tmpEFF_" in f for f in os.listdir(outputDirectory)):
        os.system(f"rm {effFileName.rstrip(".txt")}_tmpEFF_*")

    
    def parallel_sumUp(_bin):

        effis = getAllEffi(info, _bin)
        plotAllEffi(info, _bin, outputDirectoryPlots, effis)
        v1Range = _bin['title'].split(';')[0].split('<')
        v2Range = _bin['title'].split(';')[1].split('<')

        ib = int(_bin['name'].split('_')[0].replace('bin',''))

        fOut = open(f'{effFileName.rstrip(".txt")}_tmpEFF_{ib}.txt', 'w')

        if not ib:
            fOut.write('### var1 : %s\n' % v1Range[1])
            fOut.write('### var2 : %s\n' % v2Range[1])
            exp = ''
            for v in ['var1min', 'var1max', 'var2min', 'var2max']:
                exp += f'{v:8s}\t'
            for v in ['eff data', 'err data', 'eff mc', 'err mc']:
                exp += f'{v:10s}\t'
            for v in ['effD altS', 'errD altS', 'effD altB', 'errD altB', 'effMC altS', 'errMC altS', 'effMC altB', 'errMC altB','eff tagSel']:
                exp += f'{v:12s}\t'

            fOut.write(exp + '\n')

        vals = ''
        for r in [v1Range, v2Range]:
            vals += f'{float(r[0]):<+8.3f}\t{float(r[2]):+8.3f}\t'
        for v in ['Data_Nominal', 'MC_Nominal']:
            vals += f'{effis[v][0]:<10.6f}\t{effis[v][1]:<10.6f}\t'
        for v in ['Data_Alt_Sig', 'Data_Alt_Bkg', 'MC_Alt_Sig', 'MC_Alt_Bkg']:
            vals += f'{effis[v][0]:<12.6f}\t{effis[v][1]:<12.6f}\t'
        vals += f'{effis["tagSel"][0]:<12.6f}\t'

        fOut.write( vals + '\n' )
        fOut.close()
        effis = {}
        
    pool = Pool()
    pool.map(parallel_sumUp, tnpBins['bins'])

    lsfiles = [outputDirectory+"/"+f for f in os.listdir(outputDirectory) if "_tmpEFF_" in f]
    lsfiles = sorted(lsfiles, key = lambda x: int(x.rstrip(".txt").split("_")[-1]))

    os.system('cat '+' '.join(lsfiles)+' > ' + effFileName)
    os.system('rm  '+' '.join(lsfiles))
    print('Efficiencies saved in file : ', effFileName)
    
    outputDirectoryTH2 = f"{outputDirectoryPlots}/histo2D/"
    createPlotDirAndCopyPhp(outputDirectoryTH2)
    makesf.doSFs(effFileName, luminosity, ['pt','eta'], outputDirectoryTH2)

    # plotting sanity-check plots on fitresults
    outputDirectoryCheckPlots = f"{outputDirectoryPlots}/checkFitStatus/"
    createPlotDirAndCopyPhp(outputDirectoryCheckPlots)
    for typeFit in [st_name for st_name, st in fitStrategies[typeflag].items() if st is not None]:
        print()
        print("Plotting sanity-check histograms from fit results in this file")
        print(f">>> {fileName}")
        outDirChecksFit = f"{outputDirectoryCheckPlots}/{typeFit}/"
        createPlotDirAndCopyPhp(outDirChecksFit)
        print(f"Output: {outDirChecksFit}")
        
        checkFit(fileName, typeFit, binningDef, f"{outDirChecksFit}/")

    print("================================================")
    print("THE END!")
    print("================================================")
