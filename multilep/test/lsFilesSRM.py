#!/usr/bin/env python
#Tom Cornelis, 16/01/2012
import os
import re
import sys
from optparse import OptionParser

def ls(directory, search = ".*", user = "tomc", doPrint = True, onlyFailed=True):
  files = []
  # old command "srmls   -recursion_depth=10 -offset=" + format(offset) + " -count=1000 " + userdir + directory + " >& .ls.txt"
  os.system("find " + directory + " -type f &> .ls.txt")
  lsfile = open(".ls.txt")
  for line in lsfile:
    if line.strip('\n').endswith('/') or len(line.split()) == 0: continue
    if onlyFailed and not line.count('failed'): continue
    files.append(line.split()[-1])
  lsfile.close()

  def sortWithNumbers(str):
    pieces = re.split(r'(\d+)', str)
    pieces[1::2] = map(int, pieces[1::2])
    return pieces 

  if doPrint: 
    for file in sorted(files, key=sortWithNumbers): print file
  return sorted(files, key=sortWithNumbers) 


if __name__ == "__main__": 
  parser = OptionParser()
  parser.add_option("-d", "--directory", dest="directory", default="", help="directory")
  parser.add_option("-s", "--search", dest="search", default="Tuple.*root", help="search")
  parser.add_option("-u", "--user", dest="user", default="tomc")
  parser.add_option("-n", action="store_false", dest="SRMCheck", default=True)
  (options, args) = parser.parse_args()

  ls(options.directory, options.search)
