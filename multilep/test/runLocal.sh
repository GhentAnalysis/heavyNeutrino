#!/bin/bash

#include function to make list of all files in given sample
source makeFileList.sh

#function to set up CMSSW in a job
setCMSSW(){
    echo "cd ${CMSSW_BASE}/src" >> $1
    echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> $1
    echo "eval \`scram runtime -sh\`" >> $1
}

#read command-line arguments
input=$1
output=$2
skim=$3             #skim condition for the sample
                    #This will determine the name of the output files.
                    #The names of the output files will determine the
                    #skim, via multilep.py
filesPerJob=$4
                    

#if no output directory given, automatically initialize one
if [[ -z "$output" ]]; then
    #strip sample name from input
    output=${input///}    
    #set output directory to default 
    output=~/heavyNeutrino/${output}    
fi

#make output directory structure if needed
mkdir -p $output

#initialize filesPerJob to a default value if the argument is unspecified
if [[ -z "$filesPerJob" ]]; then
    filesPerJob=10
fi


#check if output exists, if not make the directory
mkdir -p $output

#make list of all files in input sample
fileList $input

#loop over new list of files and submit jobs
count=0
submit=submit.sh
while read f
    do echo "$f"
    #submit a job for every few files, as specified in the input
    if (( $count % $filesPerJob == 0 ))
        then if (( $count != 0)) 
            #then qsub $submit -l walltime=40:00:00;
            then cat $submit #temporary for testing
        fi
        #initialize temporary submission script
        if [ -e $submit ]; then rm $submit; fi
        touch $submit
        #initialize CMSSW environment in submission script
        setCMSSW $1
    fi
    echo "cmsRun ./heavyNeutrino/multilep/test/multilep.py dir=$output , inputFile=$f, outputFile=${output}/Job_${count}_${skim}.root, events=-1 > ${output}/logs/Job_${count}.txt 2> ${output}/errs/Job_${count}.txt" >> $submitS
    count=$((count + 1))
done < fileList.txt
qsub $submit -l walltime=40:00:00;
#remove temporary files
rm fileList.txt
