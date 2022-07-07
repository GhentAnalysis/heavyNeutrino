#!/bin/env python

import json, os, subprocess, sys

# suffix for the job names of the recovery tasks
SUFFIX = 'recoveryTask'

# recover MC job for a single sample
# this assumes that all jobs have been killed
def recover_sample(dirname):
    # process directory name
    dirname = os.path.realpath(dirname)
    print 'Recovering jobs for directory:', dirname
    # extract DAS data set name from directory name
    parts = dirname.split('/')
    samplename = '/'+'/'.join([
        parts[-2],                              # process name
        parts[-1][5:].split('_'+parts[-3])[0],  # campaign name
        'MINIAODSIM',                           # datatier name
    ])
    print 'Identified sample name:', samplename

    if any(SUFFIX in prt for prt in parts):
        print("Note: you are tryinh to recover a recovery task. Are you sure?")
        answ = input("y/n? ")
        if (answ == "n") return
    recoveryExistsName = dirname + '_' + SUFFIX
    
    if os.path.exists(recoveryExistsName):
        print("Error: {} is already recovered. Kill previous recovery and move/delete the crab folder before retrying.".format(dirname))
        return

    if (len(recoveryExistsName.split('/')[-1]) > 100):
        shortened = recoveryExistsName.split('/')[-1][:100]
        dirnameTempToCheck = '/'.join(parts[:-1]) + shortened
            
        if os.path.exists(dirnameTempToCheck):
            print("Error: {} is already recovered. Kill previous recovery and move/delete the crab folder before retrying.".format(dirname))
            return

    # get DAS information about all lumisections in sample
    lines = subprocess.check_output(
        'dasgoclient -query="lumi dataset={}"'.format(samplename),
        shell=True,
    )
    all_lumis = sorted(map(int, lines.split()))
    # run crab to get list of processed lumisections
    subprocess.call('crab report -d {}'.format(dirname), shell=True)

    if (not os.path.exists(os.path.join(dirname, 'results/processedLumis.json'))):
        print("Error: file {} does not exist. Crab report likely executed > 30 days after initial submission. Skipping dataset {}".format(os.path.join(dirname, 'results/processedLumis.json'), dirname))
        return

    with open(os.path.join(dirname, 'results/processedLumis.json')) as f:
        processed_lumis = sum(map(lambda (a,b): range(a,b+1), json.load(f)["1"]), [])
    # evaluate missing lumisections
    missing_lumis = sorted(set(all_lumis)-set(processed_lumis))
    if len(missing_lumis)==0:
        return
    missing_lumis_json = []
    first, last = missing_lumis[0], missing_lumis[0]
    for l in missing_lumis[1:]:
        if l-1==last:
            last = l
            continue
        missing_lumis_json.append([first, last])
        first, last = l, l
    missing_lumis_file = os.path.join(dirname, 'results/notFinishedLumis.json')
    with open(missing_lumis_file, 'w') as f:
        json.dump({"1": missing_lumis_json}, f)
    # create crab config file for recovery job
    with open('crab_temp.py', 'w') as config, open(os.path.join(dirname, 'crab.log')) as f:
        foundConfigLines = False
        for line in f:
            if foundConfigLines:
                if line.count('requestName'):
                    newname = '_'.join([parts[-1].split('crab_')[-1], SUFFIX])
                    if len(newname)>100:
                        newname = newname[:100]
                    config.write('config.General.requestName="{}"\n'.format(newname))
                elif line.count('lumiMask'):
                    config.write('config.Data.lumiMask="{}"\n'.format(missing_lumis_file))
                elif line.count("config.section_('Site')"):
                    config.write("config.section_('Site')\n")
                    config.write("config.Site.blacklist = ['T2_IN_TIFR', 'T2_US_Florida', 'T2_US_UCSD', 'T2_IT_Pisa']\n")
                elif 'DEBUG' in line:
                    break
                else:
                    config.write(line)
            elif line.count('from WMCore.Configuration import Configuration'):
                foundConfigLines = True
                config.write('from WMCore.Configuration import Configuration\n')
    # submit recovery jobs on crab
    subprocess.call('crab submit -c crab_temp.py', shell=True)
    os.remove('crab_temp.py')
    os.remove('crab_temp.pyc')

# Run from test/ directory of heavyNeutrino
# Example usages:
# ./scripts/recoverMonteCarlo.py crab/singlelepton_MC_2018_ULv5/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/crab_RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2_singlelepton_MC_2018_ULv5
# ./scripts/recoverMonteCarlo.py crab/singlelepton_MC_2018_ULv5/*/crab_*/
# ./scripts/recoverMonteCarlo.py `cat samples.txt`
#     (if samples.txt is a file with one directory name per line)
if __name__ == '__main__':
    for dirname in sys.argv[1:]:
        recover_sample(dirname)
