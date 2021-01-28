#!/usr/bin/env python3

# Wrapper around crab status which allows to use filters on output and automatic resubmissions
# crabStatus.py <options> DIRECTORY
# DIRECTORY is also optional, it will check everything in the crab directory

import os, shutil, subprocess
from optparse import OptionParser

#
# Check if we have already done cmsenv
#
if os.environ['CMSSW_BASE'].replace('/storage_mnt/storage','') not in os.getcwd():
  print('\033[1m\033[91mPlease do cmsenv first!')
  exit(0)


def system(command):
  return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT).decode()


#Option parser
parser = OptionParser()
parser.add_option("-n", "--noStatusCheck", action="store_false", dest="checkCrab", default = True, help="Crab status is not re-checked (fast debug option)")
parser.add_option("-j", "--job", dest="jobs", default = "All", help="filter on job numbers")
parser.add_option("-s", "--state", dest="state", default="All", help="filter on state")
parser.add_option("-e", "--exitCode", dest="exitCode", default="All", help="filter on exitCode")
#parser.add_option("-b", "--blacklist", dest="blacklist", default="", help="blacklist E_HOST") # OPTIONAL, not yet implemented
#parser.add_option("-k", "--kill", action="store_true", dest="kill", default = False, help="Kill the (selected) jobs")
#parser.add_option("-r", "--resubmit", action="store_true", dest="resubmit", default = False, help="Resubmit the (selected) jobs")
#parser.add_option("-f", "--forceResubmit", action="store_true", dest="forceResubmit", default = False, help="Force resubmit the (selected) jobs")
(options, args) = parser.parse_args()

# changing python list to crab list
def jobsToCrabList(jobs):
  crabList = ""
  for i in jobs:
    if crabList == "": crabList = format(i)
    else:
      pythonRange = crabList.split(',')[-1]
      last = pythonRange.split('-')[-1]
      if i == int(last) + 1:
        if pythonRange.split('-')[0] == last: crabList += '-' + format(i)
        else: crabList = crabList.replace(last, format(i))
      else: crabList += "," + format(i)
  return crabList

# changing crab list to python list
def jobsToPythonList(jobs):
  pythonList = []
  for crabRange in jobs.split(','):
    i = int(crabRange.split('-')[0])
    j = int(crabRange.split('-')[-1])
    while i <= j:
      pythonList.append(i)
      i+=1
  return pythonList

# Get the crab config from the log file
def resubmitCrabConfig(dir):
  cwd = os.getcwd()
  os.chdir(os.path.join(os.environ['CMSSW_BASE'], 'src/heavyNeutrino/multilep/test/'))
  with open('crab_temp.py', 'w') as config:
    with open(os.path.join(dir, 'crab.log')) as f:
      foundConfigLines = False
      for line in f:
        if foundConfigLines:
          if line.count('DEBUG'): break
          config.write(line)
        if line.count('from WMCore.Configuration import Configuration'):
          foundConfigLines = True
          config.write('from WMCore.Configuration import Configuration\n')
  shutil.rmtree(dir)
  print(system('crab submit -c crab_temp.py'))
  os.remove('crab_temp.py')
  os.remove('crab_temp.pyc')
  os.chdir(cwd)

def checkCrabDir(dir):
  if options.checkCrab:
    print("Getting the crab status for %s..." % dir)
    if len(dir) > 0: system("crab status --long " + dir + "> .status.txt")
    else:            system("crab status --long > .status.txt")


  jobs = []
  outputIsGood = False
  submitFailed = False
  with open(".status.txt") as statusFile:
    for str in statusFile:
      if 'SUBMITFAILED' in str:
        submitFailed = True
      if 'Jobs status' in str:
        outputIsGood = True
        if 'finished' in str and '100.0%' in str:
          print('All jobs finished with exit code 0')
          return
      if 'jobs' in str: continue
      if 'finished' in str: continue
      if len(str.split()) and str.split()[0].isdigit():
         jobs.append({ 'Job' : int(str.split()[0]), 'State' : str.split()[1], 'Exit code' : str.split()[-1], 'Str' : str})

  if not outputIsGood:
    print('Could not get crab status:')
    if submitFailed:
      with open(".status.txt") as statusFile:
        for str in statusFile:
          if "is not 'VALID' but 'DELETED'" in str or "is not 'VALID' but 'INVALID'" in str:
            print("    SUBMITFAILED --> Dataset is not 'VALID' but 'DELETED' or 'INVALID")
            for str in statusFile:
              if 'dataset=' in str: print('                     %s' % str.split('dataset=')[-1].replace('.',''))
            shutil.rmtree(dir)
            break
          if 'splitting on your task generated' in str and 'maximum number of jobs' in str:
            print("    SUBMITFAILED --> Too much jobs created, try to raise the splitting unit!")
            shutil.rmtree(dir)
            break
        else:
          print('   SUBMITFAILED --> resubmiting!')
          resubmitCrabConfig(dir)
    else:
      with open(".status.txt") as statusFile:
        for str in statusFile: print(str, end='')

  # Filter on job numbers
  if options.jobs != 'All':
    jobsToKeep = jobsToPythonList(options.jobs)
    jobs = [job for job in jobs if job['Job'] in jobsToKeep]

  # Filter on state
  if options.state != 'All':
    print(options.state)
    jobs = [job for job in jobs if job['State'] in options.state]

  # Filter on exit code
  if options.exitCode != 'All':
    jobs = [job for job in jobs if job['Exit code'] in options.exitCode]

  # Print out
  print(" Job State        Most Recent Site        Runtime   Mem (MB)      CPU %    Retries   Restarts      Waste       Exit Code")
  for job in jobs: print(job['Str'], end='')

  # Automatic resubmissions (if not killed)
  if not 'KILLED' in system('crab status -d ' + dir):
    jobsToResubmit = []
    raiseMemoryLimit = False
    raiseWallTime    = False
    for job in jobs:
      if job['Exit code'] in ['50660']: raiseMemoryLimit = True
      if job['Exit code'] in ['50664']: raiseWallTime    = True
      if job['Exit code'] in ['-1','65','83', '86', '97', '60311','60318','8001','8002','8012','8022','8021','8020','8028','134','135','8004','-15','139','60317','60307','60302','10030','10031','10034','10040','50115','50664','50660','50662','8010','7002','50665', '60324', '47']: jobsToResubmit.append(job['Job'])
      if job['State'] == 'failed' and job['Exit code'] == '0':                                  jobsToResubmit.append(job['Job'])
      if job['State'] == 'failed' and job['Exit code'] == '-1':                                 jobsToResubmit.append(job['Job'])
      if 'failed' in job['Exit code']:                                                          jobsToResubmit.append(job['Job'])
      if job['State'] == 'failed' and job['Exit code'] == 'Unknown':                            jobsToResubmit.append(job['Job'])

    if len(jobsToResubmit) > 0:
      print("Resubmitting " + jobsToCrabList(jobsToResubmit))
      if len(args) > 1: os.system("crab resubmit --jobids=" + jobsToCrabList(jobsToResubmit) + (' --maxjobruntime=2800' if raiseWallTime else '') + (' --maxmemory=2800' if raiseMemoryLimit else '') + " " +  args[1])
      else:             os.system("crab resubmit --jobids=" + jobsToCrabList(jobsToResubmit) + (' --maxjobruntime=2800' if raiseWallTime else '') + (' --maxmemory=2800' if raiseMemoryLimit else '') + " --siteblacklist=T2_US_UCSD")


# Run this script for all found crab directories
import glob
if len(args): topDir = os.path.join(os.getcwd(), args[0])
else:         topDir = os.path.join(os.environ['CMSSW_BASE'], 'src/heavyNeutrino/multilep/test/crab/')
for requestCache in glob.glob(os.path.join(topDir, '**/.requestcache'), recursive=True):
  checkCrabDir(os.path.dirname(requestCache))
