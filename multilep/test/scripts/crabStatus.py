#!/usr/bin/env python2.6

# crabStatus.py <options> DIRECTORY
# Wrapper around crab status which allows to use filters on output and automatic resubmissions

import os
from optparse import OptionParser

#Option parser
parser = OptionParser()
parser.add_option("-n", "--noStatusCheck", action="store_false", dest="checkCrab", default = True, help="Crab status is not re-checked (fast debug option)")
parser.add_option("-j", "--job", dest="jobs", default = "All", help="filter on job numbers")
parser.add_option("-s", "--state", dest="state", default="All", help="filter on state")
parser.add_option("-e", "--exitCode", dest="exitCode", default="All", help="filter on exitCode")
#parser.add_option("-b", "--blacklist", dest="blacklist", default="", help="blacklist E_HOST")
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

if options.checkCrab:
  print "Getting the crab status..."
  if len(args) > 0: os.system("crab status --long " + args[0] + "> .status.txt")
  else:             os.system("crab status --long > .status.txt")


jobs = []
outputIsGood = False
with open(".status.txt") as statusFile:
  for str in statusFile:
    if 'Jobs status' in str:
      outputIsGood = True
      if 'finished' in str and '100.0%' in str:
	print 'All jobs finished with exit code 0'
	exit(0)
    if 'jobs' in str: continue
    if 'finished' in str: continue
    if len(str.split()) and str.split()[0].isdigit():
       jobs.append({ 'Job' : int(str.split()[0]), 'State' : str.split()[1], 'Exit code' : str.split()[-1], 'Str' : str})

if not outputIsGood:
  print 'Could not get crab status:'
  with open(".status.txt") as statusFile:
    for str in statusFile: print str,

# Filter on job numbers
if options.jobs != 'All':
  jobsToKeep = jobsToPythonList(options.jobs)
  jobs = [job for job in jobs if job['Job'] in jobsToKeep]

# Filter on state
if options.state != 'All':
  print options.state
  jobs = [job for job in jobs if job['State'] in options.state]

# Filter on exit code
if options.exitCode != 'All':
  jobs = [job for job in jobs if job['Exit code'] in options.exitCode]

# Print out
print " Job State        Most Recent Site        Runtime   Mem (MB)      CPU %    Retries   Restarts      Waste       Exit Code"
for job in jobs: print job['Str'],

# Automatic resubmissions
jobsToResubmit = []
raiseMemoryLimit = False
for job in jobs:
  if job['Exit code'] in ['50660']: raiseMemoryLimit = True
  if job['Exit code'] in ['-1','83','60311','60318','8001','8002','8022','8021','8020','8028','134','135','8004','-15','139','60317','60307','60302','10030','10031','10034','10040','50115','50664','50660','50662', '8010']: jobsToResubmit.append(job['Job'])
  if job['State'] == 'failed' and job['Exit code'] == '0':                                  jobsToResubmit.append(job['Job'])
  if job['State'] == 'failed' and job['Exit code'] == '-1':                                 jobsToResubmit.append(job['Job'])
  if 'failed' in job['Exit code']:                                                          jobsToResubmit.append(job['Job'])
  if job['State'] == 'failed' and job['Exit code'] == 'Unknown':                            jobsToResubmit.append(job['Job'])

if len(jobsToResubmit) > 0:
  print "Resubmitting " + jobsToCrabList(jobsToResubmit)
  if len(args) > 0: os.system("crab resubmit --jobids=" + jobsToCrabList(jobsToResubmit) + (' --maxmemory=4000' if raiseMemoryLimit else '') + " " +  args[0])
  else:             os.system("crab resubmit --jobids=" + jobsToCrabList(jobsToResubmit) + (' --maxmemory=4000' if raiseMemoryLimit else '') + " --siteblacklist=T2_US_Wisconsin,T2_US_Vanderbilt,T2_CH_CERN,T2_BE_IIHE,T2_IT_Rome")
