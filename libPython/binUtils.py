#!/usr/bin/env python3

import copy


## whoever wrote this code originally should take a few lessons in basic python operations
## this is, by all means, terrible....
## i'm not proud that i didn't fix it, i'm just lazy

## LETSGOSKI

def createBins(bining, cut=None):

    nbin = 1
    index = [-1 for i in bining.keys()]
    listOfIndex = []    
    listOfIndex.append( index )
    print('this is listOfIndex', listOfIndex)
    ### first map nD bins in a single list
    for iiv, (var, bin_dict) in enumerate(bining.items()):
        print(iiv, var, bin_dict)
        if not any([x in list(bin_dict.keys()) for x in ["type", "bins"]]):
            print(f'Bining is not complete for var {var}')
            return listOfIndex
        nb1D = 1
        if   bin_dict['type']=='float':
            nb1D = len(bin_dict['bins'])-1
        elif bin_dict['type']=='int':
            nb1D = len(bin_dict['bins'])
        nbin = nbin * nb1D

        print(f"Variable {var} -> {nb1D} bins")
        listOfIndexInit = copy.deepcopy(listOfIndex)
        for ib_var in range(nb1D):
            if ib_var == 0 :
                for ib in range(len(listOfIndex)):
                    listOfIndex[ib][iiv] = ib_var            
            else: 
                for ib in range(len(listOfIndexInit)):
                    listOfIndexInit[ib][iiv] = ib_var
                           
                listOfIndex.extend(copy.deepcopy(listOfIndexInit))

    listOfBins = []
    nbins = len(listOfIndex)
    print('this is listOfIndex', listOfIndex)
    for ibin, ix in enumerate(listOfIndex):
        if   nbins <= 100 :   binName  = f'bin{ibin:02d}_'
        elif nbins <= 1000 :  binName  = f'bin{ibin:03d}_'
        elif nbins <= 10000 : binName  = f'bin{ibin:04d}_'
        else:                 binName  = f'bin{ibin}_'

        binTitle = ''
        binCut   = f"{cut} && " if cut is not None else ""
        binVars  = {}

        for iv, (var, bin_dict) in enumerate(bining.items()):
            varType, bins1D = bin_dict['type'], bin_dict['bins']
            lowEdge, highEdge = bins1D[ix[iv]], bins1D[ix[iv]+1]

            if iv!=0:
                binCut   += ' && '
                binTitle += '; '
                binName  += '_'
            
            if varType == 'float' :
                binCut   += f'{var} >= {lowEdge} && {var} < {highEdge}'
                binTitle += f'{lowEdge:1.3f} < {var} < {highEdge:1.3f}'
                binName  += f'{var}_{lowEdge:1.2f}To{highEdge:1.2f}'
                binVars[var] = {'min': lowEdge, 'max': highEdge}

            if varType == 'int' :
                binCut   += f'{var} == {lowEdge}'
                binTitle += f'{var} = {lowEdge}'
                binName  += f'{var}Eq{lowEdge}'
                binVars[var] = {'min': lowEdge, 'max': lowEdge}

            binName = binName.replace('-', 'm')
            binName = binName.replace('.', 'p')

        listOfBins.append({'cut' : binCut, 'title': binTitle, 'name' : binName, 'vars' : binVars })
        
    binDefinition = {
        'vars' : list(bining.keys()),
        'bins' : listOfBins
    }

    return binDefinition


def binMinPt(tnpBin):
    ptmin = 1
    if tnpBin['name'].find('pt_') >= 0:
        ptmin = float(tnpBin['name'].split('pt_')[1].split('p')[0])
    elif tnpBin['name'].find('et_') >= 0:
        ptmin = float(tnpBin['name'].split('et_')[1].split('p')[0])
    return ptmin


def testBinning(bins, testbins, var, flag, allowRebin=False):
    """
    """
    
    if bins == testbins:
        return 0

    if bins[0] in testbins:
        first_test_idx = testbins.index(bins[0])
        is_slice = all(bins[i] == testbins[i + first_test_idx] for i in range(len(bins)))
        is_subset = allowRebin and all(bin_edge in testbins for bin_edge in bins)

        if is_slice:
            print(f"\nWarning: {var} binning not consistent with the one in histograms for {flag}")
            print(f"Bins: {bins}")
            print(f"Test bins: {testbins}")
            print("However, it seems to be a slice of it. Proceeding with caution!\n")
            return 0

        if is_subset:
            print(f"\nWarning: {var} binning not consistent with the one in histograms for {flag}")
            print(f"Bins: {bins}")
            print(f"Test bins: {testbins}")
            print("However, it seems to be a subset of it. Proceeding with caution!\n")
            return 0

    print(f"\nError: {var} binning not consistent with the one in histograms for {flag}")
    print(f"Bins: {bins}")
    print(f"Test bins: {testbins}")
    print("Please check!\n")
    return -1


def tuneCuts( bindef, cuts ) :
    if cuts is None:
        return
    
    for ibin in cuts.keys():
        cut0 = bindef['bins'][ibin]['cut']
        cut1 = cuts[ibin]
        bindef['bins'][ibin]['cut'] = '%s && %s ' % (cut0,cut1)
    

