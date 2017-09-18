#!/bin/bash

#read command-line arguments
input=$1
output=$2
jobsPerFile=$3

#list all files present in sample in a txt file
python ~/das_client.py --query="file dataset=$input" --limit=0 > fileNames.txt

#location of files on local cluster
location=/pnfs/iihe/cms/ph/sc4
dcap=dcap://maite.iihe.ac.be

#xrootd redirector for files not present on local cluster

#make temporary txt file to store final list of files

#loop over files in sample
    
    #check if file is locally available

    #if not available add xrootd redirector to remote location





#remove temporary files
rm fileNames.txt
rm fileList.txt
