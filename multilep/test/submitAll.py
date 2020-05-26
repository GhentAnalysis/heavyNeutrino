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
  else:
    skim, dataset = dataset.split(':')


    if 'heavyN' in dataset: # just make some similar requestName as standard samples
      if 'Moriond17_aug2018_miniAODv3' in dataset: miniAODver  = 'MiniAOD2016v3'
      if 'Fall17' in dataset:                      miniAODver  = 'RunIIFall17MiniAODv2'
      if 'Autumn18' in dataset:                    miniAODver  = 'RunIIAutumn18MiniAOD'
      requestName = '%s_%s' % (miniAODver, productionLabel)
    else:
      requestName = dataset.split('/')[2] + '_' + productionLabel
      requestName = requestName.replace('RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3', 'MiniAOD2016v3')
      requestName = requestName.replace('RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14', 'MiniAOD2017v2')
      requestName = requestName.replace('RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14', 'MiniAOD2017v2NewPMX')
      requestName = requestName.replace('RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15', 'MiniAOD2018')

    if submitLocal or 'user' in dataset:
      print 'Submitting %s (requestname %s) on the local resources:' % (dataset, requestName)
      if 'user' in dataset:   outputDir = dataset.split('/')[-1]
      else:                   outputDir = dataset.split('/')[1]

      if 'Run201' in dataset: outputDir = os.path.join(outputDir, dataset.split('-')[0].split('/')[-1])

      outputDir = os.path.join(outputDir, 'localSubmission_%s' % requestName)

      extra  = ('extraContent=' + ','.join(extraContent)) if len(extraContent) else ''
      os.system('bash runLocal.sh ' + dataset + ' ' + outputDir + ' ' + skim + ' ' + str(filesPerJob) + ' ' + extra)

    else:
      print 'Submitting %s (requestname %s) to crab:' % (dataset, requestName)
      os.environ['CRAB_PRODUCTIONLABEL'] = productionLabel
      os.environ['CRAB_DATASET']         = dataset
      os.environ['CRAB_OUTPUTFILE']      = skim + '.root'
      os.environ['CRAB_EXTRACONTENT']    = ','.join(extraContent)
      os.environ['CRAB_REQUESTNAME']     = requestName
      if 'Run2017' in dataset :   os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/" +       getJSON(True, False)
      elif 'Run2018' in dataset : os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/" +       getJSON(False, True)
      else :                      os.environ['CRAB_LUMIMASK'] = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/" + getJSON(False, False)
      os.system('crab submit -c crab.py')
