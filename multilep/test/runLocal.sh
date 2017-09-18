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
redirector=root://cmsxrootd.fnal.gov//

#make temporary txt file to store final list of files
touch fileList.txt

while read f
    #check if file is locally available
    do if [ -e ${location}${f} ]
        then echo "${dcap}${location}${f}" >> fileList.txt
    #if the file is not available add an xrootd redirector to a remote location
    else 
        echo "${redirector}${f}" >> fileList.txt
    fi
done < fileNames.txt

cat fileNames.txt



#remove temporary files
rm fileNames.txt
rm fileList.txt
