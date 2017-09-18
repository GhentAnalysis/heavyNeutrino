#!/bin/bash

input=$1
output=$2
filesPerJob=$3
#find all files contained in sample and write them to "filenames.txt"
python ~/das_client.py --query="file dataset=$input" --limit=0 > filenames.txt

#location of samples on local cluster
location=/pnfs/iihe/cms/ph/sc4
dcap=dcap://maite.iihe.ac.be
#name of xrootd redictector for files that are not present
redirector=root://cmsxrootd.fnal.gov//
#new temporary txt file containing files to run over
touch fileList.txt
#loop over all files in the sample
while read f;
    #check if file is present on the local cluster
    do if [ -e ${location}${f} ] 
        #add path to local file to list
        then echo "${dcap}${location}${f}" >> fileList.txt
    else 
        #add remote path to list
        echo "${redirector}${f}" >> fileList.txt
    fi
done < filenames.txt

#loop over new file list to submit jobs
while read f
    do echo "$f"
done < fileList.txt

#clean up temporary files
rm filenames.txt
rm fileList.txt
