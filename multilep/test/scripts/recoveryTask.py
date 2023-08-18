#!/usr/bin/env python3

# Wrapper around crab status which allows to use filters on output and automatic resubmissions
# crabStatus.py <options> DIRECTORY
# DIRECTORY is also optional, it will check everything in the crab directory

import os, subprocess
from optparse import OptionParser

#
# Check if we have already done cmsenv
#
if os.environ['CMSSW_BASE'].replace('/ada_mnt/ada','') not in os.getcwd():
  print('\033[1m\033[91mPlease do cmsenv first!')
  exit(0)

def system(command):
  try: 
    return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT).decode()
  except Exception as e:
    return e.output


#Option parser
parser = OptionParser()
parser.add_option("-d", "--dirs",         dest="dirs",         default=None,           help="Directories for which to make a recovery task, could take wildcards")
parser.add_option("-r", "--recoveryTask", dest="recoveryTask", default='recoveryTask', help="Set recoveryTask=<name> to launch recovery tasks, former tasks will be killed")
parser.add_option("-s", "--splitting",    dest="splitting",    default=None,           help="Change the splitting type with respect to the original task")
parser.add_option("-u", "--units",        dest="units",        default=None,           help="Change the splitting unit with respect to the orginal task")
(options, args) = parser.parse_args()

# Get the crab config from the log file
def recoveryJob(dir):
  print('Making a recovery job for %s' % dir)
  if not 'KILLED' in system('crab status -d ' + dir):
    system('crab kill -d ' + dir)
  print(system('crab report -d ' + dir))
  missingLumis = os.path.join(dir, 'results/notFinishedLumis.json')
  newName      = dir.split('/crab_')[-1] + '_%s' % options.recoveryTask

  cwd = os.getcwd()
  os.chdir(os.path.join(os.environ['CMSSW_BASE'], 'src/heavyNeutrino/multilep/test/'))
  with open('crab_temp.py', 'w') as config:
    with open(os.path.join(dir, 'crab.log')) as f:
      foundConfigLines = False
      for line in f:
        if foundConfigLines:
          if newName and line.count('requestName'):
            config.write('config.General.requestName="%s"\n' % newName)
          elif line.count('lumiMask'):
            config.write('config.Data.lumiMask="%s"\n' % missingLumis)
          elif options.splitting and line.count('splitting'):
            config.write('config.Data.splitting="%s"\n' % options.splitting)
          elif options.units and line.count('unitsPerJob'):
            config.write('config.Data.unitsPerJob=%s\n' % options.units)
          elif 'DEBUG' in line:
            break
          else:
            config.write(line)
        if line.count('from WMCore.Configuration import Configuration'):
          foundConfigLines = True
          config.write('from WMCore.Configuration import Configuration\n')
  print(system('crab submit -c crab_temp.py'))
  os.remove('crab_temp.py')
  os.remove('crab_temp.pyc')
  os.chdir(cwd)


import glob
if len(args): topDir = os.path.join(os.getcwd(), args[0])
else:         topDir = os.path.join(os.environ['CMSSW_BASE'], 'src/heavyNeutrino/multilep/test/crab/')
for requestCache in glob.glob(os.path.join(topDir, options.dirs, '.requestcache'), recursive=True):
  recoveryJob(os.path.dirname(requestCache))
