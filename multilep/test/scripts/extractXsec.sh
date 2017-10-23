#!/bin/bash

#include function to make list of files in sample
source makeFileList.sh

#give input sample and output txt file to hold the cross section info
input=$1
output=$2

#make list of files
fileList $input

#make one big string as an argument to the cross section extraction script
total=""
while read f; 
    do total="$total,$f"
done  < fileList.txt
#cut out final characters
total=${total#*,}

#download xsec extraction script
if ! [ -f ana.py ]
    then curl https://raw.githubusercontent.com/syuvivida/generator/master/cross_section/runJob/ana.py  -o ana.py
fi
cmsRun ana.py inputFiles="$total" maxEvents=-1 > $output 2>>$output

#clean up temporary files
rm fileList.txt
