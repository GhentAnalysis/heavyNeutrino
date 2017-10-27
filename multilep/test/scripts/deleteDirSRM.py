#!/usr/bin/env python
import os
import re
import lsFilesSRM
import getpass
import sys

userdir = "srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/" + getpass.getuser() + "/"
directory = sys.argv[1]
print "List of files to delete:"
list = lsFilesSRM.ls(directory, '.*', getpass.getuser(), True, onlyFailed=('onlyFailed' in sys.argv))
directory = directory.replace("/pnfs/iihe/cms/store/user/" + getpass.getuser() + "/", "")
if not directory.endswith('/'): directory += '/'

# Function for confirmation
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

if(confirm('Delete these files?')):
  lsfile = open(".ls.txt")
  for file in list:
    print "Delete srm://maite.iihe.ac.be:8443" + file.split()[-1] 
    os.system("srmrm srm://maite.iihe.ac.be:8443" + file.split()[-1])
#  if os.listdir(userdir + directory) == []: os.system("srmrmdir -recursive=true " + userdir + directory) try to put some check in here if everything is empty
  if True: os.system("srmrmdir -recursive=true " + userdir + directory)
print "DONE"

