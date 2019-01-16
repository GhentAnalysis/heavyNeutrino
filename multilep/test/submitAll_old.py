#!/usr/bin/env python
import os, glob, sys

datasetsFile    = sys.argv[1]                                                                             # Input file with datasets
productionLabel = os.path.basename(datasetsFile.split('.')[0].split('/')[-1])                             # Label to keep track of the tuple version (is taken from the name of the above input file)
outDir          = '/user/' + os.environ['USER'] + '/public/heavyNeutrino/taaaalk'                                 # Output directory in case of local submission
datasets        = [dataset.strip() for dataset in open(datasetsFile)]                                     # Get list of datasets from file given as first argument
datasets        = [dataset.split()[0] for dataset in datasets if dataset and not dataset.startswith('#')] # Clean empty and comment lines
groupFiles      = 1                                                                                       # Group files together when running locally


for dataset in datasets:
  outputName, dataset = dataset.split(':')

  if 'pnfs' in dataset or 'user' in dataset:
    if 'user' in dataset: datasetName = dataset.split('/')[-1]
    else:        datasetName = dataset.split('/MINIAOD')[0].split('/')[-1]
    print dataset

    i = 0
    j = 0
    inputFiles = []
    for file in glob.glob(dataset + ('/*.root' if 'user' in dataset else '/*/*.root')):
      j          += 1
      inputFiles += [('dcap://maite.iihe.ac.be' if not 'user' in dataset else 'file://') + file]
      if j%groupFiles!=0: continue

      dir        = os.getcwd()
      wallTime   = '05:00:00'
      inputFile  = ','.join(inputFiles)
      outputFile = os.path.join(outDir, datasetName, 'local_' + productionLabel, outputName + '_' + str(i) + '.root')
      logFile    = os.path.join(outDir, datasetName, 'local_' + productionLabel, 'log', outputName + '_' + str(i) + '.log')

      for filename in [outputFile, logFile]:
        try:    os.makedirs(os.path.dirname(filename))
        except: pass
      
      print 'Submitting ' + inputFile + ' to cream02:'
      args  = 'dir=' + dir + ',inputFile=\"' + inputFile + '\",outputFile=' + outputFile + ',events=-1'
      args += ',isData='+ ('False' if ('SIM' in dataset or 'HeavyNeutrino' in dataset) else 'True')
      os.system('qsub -v ' + args + ' -q localgrid@cream02 -o ' + logFile + ' -e ' + logFile + ' -l walltime=' + wallTime + ' runOnCream02.sh')
      i += 1
      inputFiles = []

  else: # use crab
    print 'Submitting ' + dataset + ' using crab:'
    os.environ['CRAB_PRODUCTIONLABEL'] = productionLabel
    os.environ['CRAB_DATASET']         = dataset
    os.environ['CRAB_OUTPUTFILE']      = outputName + '.root'
    os.system('crab submit -c crab.py')

