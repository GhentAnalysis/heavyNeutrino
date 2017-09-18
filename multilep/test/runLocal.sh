#!/bin/bash

#read command-line arguments
input=$1
output=$3
filesPerJob=$2

#check if output exists

#if no output directory given, automatically initialize one

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

#loop over new list of files and submit jobs
count=0
submit=submit.sh
while read f
    do echo "$f"
    #submit a job for every few files, as specified in the input
    if (( $count % $filesPerJob == 0 ))
        then if (( $count != 0)) 
            #then qsub $submit -l walltime=40:00:00;
            then cat $submit
        fi
        #initialize temporary submission script
        if [ -e $submit ]; then rm $submit; fi
        touch $submit
        #initialize CMSSW environment in submission script
        echo "cd ${CMSSW_BASE}/src" >> $submit
        echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> $submit
        echo "eval \`scram runtime -sh\`" >> $submit
    fi
    echo "cmsRun ./heavyNeutrino/multilep/test/multilep.py $f ${output}/Job_${count}.root > ${output}/logs/Job_${count}.txt 2> ${output}/errs/Job_${count}.txt" >> $submit
    count=$((count + 1))
done < fileList.txt
qsub $submit -l walltime=40:00:00;
#remove temporary files
rm $submit
rm fileNames.txt
rm fileList.txt
