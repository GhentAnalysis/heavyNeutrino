#!/bin/env python

##################################################################################
# extension to monitorCrabJobs.py which allows continuous automatic resubmission #
##################################################################################

# how to use:
# - make sure to have done cmsenv, created a proxy, etc.
#   (it might be good to check that the crab status command is working as expected 
#    by running it on a single sample).
# - change the settings below to your liking.
# - run with "nohup python monitorCrabJobsAuto.py &> [name of a log file] &".
#   some explanation:
#   - "python monitorCrabJobsAuto.py" trivially runs this script.
#   - "&> [name of a log file]" redirects both stdout and stderr to the log file you chose.
#   - "nohup [...] &" runs the command in the background, 
#     and keeps it running even after you close the terminal.
# - alternatively, you could open a screen session and simply run 
#   "python monitorCrabJobsAuto.py &> [name of a log file]".
#   this gives you more control over terminating the command,
#   since with "nohup &", the process ID seems to be irretrievable 
#   after logging out of the m-machine (?).
# notes:
# - you will need a valid proxy for the entire duration of this script,
#   so create one with a long enough lifetime. 
#   if the proxy becomes invalid during the runtime of this script,
#   you risk overwriting the information with 'finished 0%' for each sample,
#   since the crab status command does not work anymore.
#   TO DO: implement a safety against this overwriting, see also monitorCrabJobs.py.


import sys
import os
import time


if __name__=='__main__':

    # settings
    niterations = 1
    # (number of iterations to do)
    tsleep = 4*3600
    # (time to sleep in seconds)
    crabdir = '../crab'
    # (folder to monitor, passed down to monitorCrabJobs.py)
    resubmit = True
    # (whether or not to resubmit, passed down to monitorCrabJobs.py)
    webpage = 'crab_status_auto'
    # (name of the webpage where to put the result)

    # make the command
    cmd = 'python monitorCrabJobs.py'
    cmd += ' --crabdir {}'.format(crabdir)
    cmd += ' --resubmit {}'.format(resubmit)
    cmd += ' --webpage {}'.format(webpage)
    # temp for testing
    cmd += ' --istest True'

    # loop
    idx = 0
    while idx < niterations:
        os.system(cmd)
        sys.stdout.flush()
        sys.stderr.flush()
        time.sleep(tsleep)
        idx += 1
