#!/usr/bin/env python
import os, glob, sys, subprocess

#
# Check if we have already done cmsenv
#
if os.environ['CMSSW_BASE'].replace('/storage_mnt/storage','') not in os.getcwd():
  print '\033[1m\033[91mPlease do cmsenv first!'
  exit(0)

#
# Special git diff check for special people like Willem
#
def confirm(prompt, resp=False):
  prompt = '%s %s|%s: ' % (prompt, 'y', 'n')
  while True:
    ans = raw_input(prompt)
    if not ans: return resp
    if ans not in ['y', 'Y', 'n', 'N']:
      print 'please enter y or n.'
      continue
    if ans == 'y' or ans == 'Y': return True
    if ans == 'n' or ans == 'N': return False

def system(command):
    return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

if len(system('git diff')):
  print '\033[1m\033[91mYou have local changes not yet committed to git!'
  print
  print os.system('git diff') # directly using os system to keep colorized output
  print
  if not confirm('\033[1m\033[94mAre you sure you would want to continue?'):
    exit(0)
#
# Proceed with the submission script
#
from multilep import getJSON

datasetsFile    = sys.argv[1]                                                                             # Input file with datasets
submitLocal     = (sys.argv[2]=='local') if len(sys.argv) > 2 else False                                  # Check if call asked for local submission
filesPerJob     = sys.argv[3]            if len(sys.argv) > 3 else 10                                     # Use third argument to specify the number of jobs per file
productionLabel = os.path.basename(datasetsFile.split('/')[-1].split('.')[0])                             # Label to keep track of the tuple version (is taken from the name of the above input file)
datasets        = [dataset.strip() for dataset in open(datasetsFile)]                                     # Get list of datasets from file given as first argument
datasets        = [dataset.split()[0] for dataset in datasets if dataset and not dataset.startswith('#')] # Clean empty and comment lines
extraContent    = []

for dataset in datasets:
  if dataset.startswith('+'):
    extraContent += dataset.replace('+','').split(',')
  elif dataset.startswith('-'):
    for toRemove in dataset.replace('-','').split(','):
      extraContent.remove(toRemove)
  elif dataset.startswith('%'):
    continue
  else:
    skim, dataset = dataset.split(':')

    if submitLocal or 'user' in dataset:
      print 'Submitting ' + dataset + ' on the local resources:'
      if 'user' in dataset:   outputDir = dataset.split('/')[-1]
      else:                   outputDir = dataset.split('/')[1]

      if 'Run201' in dataset: outputDir = os.path.join(outputDir, dataset.split('-')[0].split('/')[-1])

      if 'ext' in dataset:    outputDir = os.path.join(outputDir, 'localSubmission_ext' + dataset.split('ext')[-1].split('/')[0] + '_' + productionLabel)
      else:                   outputDir = os.path.join(outputDir, 'localSubmission_' + productionLabel)

      extra  = ('extraContent=' + ','.join(extraContent)) if len(extraContent) else ''
      os.system('bash runLocal.sh ' + dataset + ' ' + outputDir + ' ' + skim + ' ' + str(filesPerJob) + ' ' + extra)

    else:
      print 'Submitting ' + dataset + ' using crab:'
      os.environ['CRAB_PRODUCTIONLABEL'] = productionLabel
      os.environ['CRAB_DATASET']         = dataset
      os.environ['CRAB_OUTPUTFILE']      = skim + '.root'
      os.environ['CRAB_EXTRACONTENT']    = ','.join(extraContent)
      if 'Run2017' in dataset :   os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/" +       getJSON(True, False)
      elif 'Run2018' in dataset : os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/" +   getJSON(False, True)
      else :                      os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/" + getJSON(False, False)
      os.system('crab submit -c crab.py')
