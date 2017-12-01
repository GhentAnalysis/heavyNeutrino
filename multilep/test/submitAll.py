#!/usr/bin/env python
import os, glob, sys
#import function returning JSON from multilep.py
from multilep import getJSON

datasetsFile    = sys.argv[1]                                                                             # Input file with datasets
productionLabel = os.path.basename(datasetsFile.split('/')[-1].split('.')[0])                             # Label to keep track of the tuple version (is taken from the name of the above input file)
outDir          = '/user/' + os.environ['USER'] + '/public/heavyNeutrino'                                 # Output directory in case of local submission
datasets        = [dataset.strip() for dataset in open(datasetsFile)]                                     # Get list of datasets from file given as first argument
datasets        = [dataset.split()[0] for dataset in datasets if dataset and not dataset.startswith('#')] # Clean empty and comment lines
#check if call asked for local submission
submitLocal     = ""
#Use third argument to specify the number of jobs per file
filesPerJob     = 10

if len(sys.argv) > 2:
    submitLocal = sys.argv[2]
    if len(sys.argv) > 3:
        filesPerJob = sys.argv[3]

for dataset in datasets:
    skim, dataset = dataset.split(':')

    if 'pnfs' in dataset or 'user' in dataset or submitLocal == "local":
        dir        = os.getcwd()
        outputDir = outDir + "/"
        if 'pnfs' in dataset or 'user' in dataset: outputDir = outputDir + dataset.split('/')[-1]      
        else: outputDir = outputDir + dataset.split('/')[1]
        if 'ext' in dataset:
            outputDir = outputDir + "_"
            index = dataset.find('ext')
            for i in range (index, len(dataset)):
                if dataset[i] == "/": break
                else : outputDir = outputDir + dataset[i]

        print outputDir
        #cut out the first part of /pnfs path for official sample if needed
        if 'pnfs' in dataset and 'user' not in dataset:
            dataset = dataset.replace("/pnfs/iihe/cms/ph/sc4/store/mc", "")
            #naming in pnfs directories is slightly altered compared to the CMSDAS name, rever this to the CMSDAS version:
            period = dataset[1:].split('/')[0]
            name   = dataset[1:].split('/')[1]
            form   = dataset[1:].split('/')[2]
            puScen = dataset[1:].split('/')[3]
            dataset = '/' + name + '/' + period + '-' + puScen + '/' + form

        os.system('bash runLocal.sh ' + dataset + ' ' + outputDir + ' ' + skim + ' ' + str(filesPerJob) )

    else: # use crab
        print 'Submitting ' + dataset + ' using crab:'
        os.environ['CRAB_PRODUCTIONLABEL'] = productionLabel
        os.environ['CRAB_DATASET']         = dataset
        os.environ['CRAB_OUTPUTFILE']      = skim + '.root'
        if 'Run2017' in dataset : os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/PromptReco/" + getJSON( ('Run2017' in dataset) )
        else : os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/" + getJSON( ('Run2017' in dataset) )
        os.system('crab submit -c crab.py')
