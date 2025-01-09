#!/usr/bin/env python3 

import subprocess
#from multiprocessing import Process
import time
import argparse
from libPython import working_points

def runCommands(wp, era, options, inputMC, inputData, inputBkg):
    '''
    '''
    print()
    print("-"*30)
    print("Working point = ", wp)
    print("-"*30)
    opt_e = '--era=' + era
    opt_f = '--flag=' + wp
    opt_iMC = '--inputMC=' + inputMC
    opt_iData = '--inputData=' + inputData
    opt_iBkg = '--inputBkg=' + inputBkg if inputBkg else ''

    cmds = []
    ex = 'tnpEGM_fitter.py'
    cmds.append(['python', ex, opt_e, opt_f, '--createBins'])
    cmds.append([el for el in ['python', ex, opt_e, opt_f, opt_iMC , opt_iData, opt_iBkg, '--createHists'] if el!=""])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit'                        ])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig'             ])
    #cmds.append(['python', ex, opt_e, opt_f, '--doFit',             '--altSig'])
    #cmds.append(['python', ex, opt_e, opt_f, '--doFit',             '--altBkg'])
    
    
    # cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig',  '--altSig'])
    # cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig',  '--altBkg'])

    cmds.append(['python', ex, opt_e, opt_f, '--sumUp'                        ])


    for cmd in cmds:
        if args.onlySumUp and "--sumUp" not in cmd:
            continue
        
        if options.outdir:          cmd.append(f"--outdir={options.outdir}")
        if options.useTrackerMuons: cmd.append("--useTrackerMuons")
        
        if options.dryRun:
            print(' '.join(cmd))
        else:
            subprocess.run(cmd, check=True)




all_wp = working_points.wp_list["working_points_global"]

parser = argparse.ArgumentParser()
parser.add_argument('-o',  '--outdir', default=None, type=str,
                    help='name of the output folder (if not passed, a default one is used, which has the time stamp in it)')
# parser.add_argument('-e',  '--era', default=['GtoH'], nargs='+', type=str, choices=['GtoH', 'BtoF'],
#                     help='Choose the era')
parser.add_argument("-y", "--year", type=str, default="2016", choices=["2016", "2017", "2018"], 
                    help="Year of data taking")
parser.add_argument('-d',  '--dryRun', action='store_true',
                    help='Do not execute commands, just print them')
parser.add_argument('-s','--steps', default=None, nargs='*', type=str, choices=list([x.split("_")[1] for x in all_wp.keys()]),
                    help='Default runs all working points, but can choose only some if needed')
parser.add_argument('-x','--exclude', default=None, nargs='*', type=str, choices=list([x.split("_")[1] for x in all_wp.keys()]),
                    help='Default runs all working points, but can exclude some if needed')
parser.add_argument('--forVeto', action='store_true', help='Use the working points related to veto')               
parser.add_argument('--useTrackerMuons', action='store_true'  , help = 'Measuring efficiencies specific for tracker muons (different tunings needed')
parser.add_argument('--onlySumUp', action='store_true',
                    help='Execute only the final command with --sumUp')
args = parser.parse_args()

tstart = time.time()
cpustrat = time.process_time()

if args.exclude and args.steps:
    print("Warning: --exclude and --steps are not supposed to be used together. Try again")
    quit()

era = "GtoH"  # fixed to this value since we will not run on 2016BtoF. Eventually it will be "removed" when all the code will be adapted to run on 2017-18

wp_sel = "working_points_"
if args.useTrackerMuons:
    wp_sel += "tracker"  
else:
    wp_sel += args.year
    if args.forVeto:
        wp_sel += "_forVeto"

    wp_sel.replace("2016", "global")  # maintained this convention for backward compatibility    

wp_dict = working_points.wp_list[wp_sel]

if args.steps:
    stepsToRun = [x for x in wp_dict.keys() if x.split("_")[1] in args.steps]
elif args.exclude:
    stepsToRun = [x for x in wp_dict.keys() if x.split("_")[1] not in args.exclude]
else:
    stepsToRun = wp_dict.keys()
    
for wp_name in wp_dict:
    if wp_name not in stepsToRun:
        continue
    inputMC = wp_dict[wp_name][0]
    inputData = wp_dict[wp_name][1]
    inputBkg = wp_dict[wp_name][2] if wp_name.split("_")[1] in working_points.wp_with_bkg else None
    runCommands( wp_name, era, args, inputMC, inputData, inputBkg)


print('\n\n')
print('Execution time:', time.time()-tstart, 'seconds')
print('CPU Execution time:', time.process_time()-cpustrat , 'seconds')
print('\n')
