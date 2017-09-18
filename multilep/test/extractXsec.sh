#!/bin/bash

input=$1
output=$2
python ./das_client.py --query="file dataset=$input" --limit=0 > filenames.txt
total=""
while read f; 
    do f="root://cmsxrootd.fnal.gov//$f"
    total="$total,$f"
done  < filenames.txt
total=${total#*,}
if ! [ -f ana.py ]
    then curl https://raw.githubusercontent.com/syuvivida/generator/master/cross_section/runJob/ana.py  -o ana.py
fi
cmsRun ana.py inputFiles="$total" maxEvents=-1 > $output 2>>$output

