#!/bin/bash

#include function to make list of all files in given sample
source scripts/makeFileList.sh

#function to set up CMSSW in a job
setCMSSW(){
    echo "cd ${CMSSW_BASE}/src" >> $1
    echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> $1
    echo "eval \`scram runtime -sh\`" >> $1
    echo "cd heavyNeutrino/multilep/test/" >> $1
}

#function to submit a job and catch invalid credentials
submitJob(){
    qsub $1 -l walltime=40:00:00 > outputCheck.txt 2>> outputCheck.txt
    while grep "Invalid credential" outputCheck.txt; do
        echo "Invalid credential caught, resubmitting"
        sleep 2  #sleep 2 seconds before attemtping resubmission
        qsub $1 -l walltime=40:00:00 > outputCheck.txt 2>> outputCheck.txt
    done
    cat outputCheck.txt
    rm outputCheck.txt
}

##################################
#Change to your proxy!
proxy=/user/wverbeke/x509up_u20640
##################################

if ! [[ $proxy = *"$USER"* ]]; then
    echo "Change path to proxy in RunLocal.sh before submitting jobs!"
    exit 1
fi

exportProxy(){
   echo "export X509_USER_PROXY=$proxy" >> $1
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
    output=$input
    if [[ $input == *"/user/"* ]] || [[ $input == *"/pnfs/"* ]]; then
        if [[  -z "${a##/}" ]]; then
            output=${output%:*}
        fi
        output=${input##/}
    else 
        #strip sample name from input
        output=${input:1}    
        output=${output%%/}
        #set output directory to default 
    fi
    echo "OUTPUT = $output"
    output=~/public/heavyNeutrino/${output}    
fi

#make output directory structure if needed
mkdir -p $output

#initialize filesPerJob to a default value if the argument is unspecified
if [[ -z "$filesPerJob" ]]; then
    filesPerJob=10
fi


#check if output exists, if not make the directory
mkdir -p $output
mkdir -p ${output}/errs
mkdir -p ${output}/logs

#make list of all files in input sample
fileList $input

#loop over new list of files and submit jobs
fileCount=0
#jobCount=0
submit=crabSucks.sh
fileList=""
while read f; do
    #fileList="${fileList}${f},"
    fileCount=$((fileCount + 1))
    #submit a job for every few files, as specified in the input
    if (( $fileCount % $filesPerJob == 0 )) || (( $fileCount == 1 ))
        then if (( $fileCount % $filesPerJob == 0 )); then
            #then fileList="${fileList%,}" #remove trailing comma from fileList
            #echo "cmsRun ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/multilep.py inputFile=$fileList outputFile=${output}/Job_${jobCount}_${skim}.root events=-1 > ${output}/logs/Job_${jobCount}.txt 2> ${output}/errs/Job_${jobCount}.txt" >> $submit
            #submit job
            submitJob $submit
            #cat $submit
            jobCount=$((jobCount + 1))
            fileList=""
        fi
        #initialize temporary submission script
        if [ -e $submit ]; then rm $submit; fi
        touch $submit
        #initialize CMSSW environment in submission script
        setCMSSW $submit
        exportProxy $submit
    fi
    echo "cmsRun ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/multilep.py inputFile=$f outputFile=${output}/${skim}_Job_${fileCount}.root events=-1 > ${output}/logs/Job_${fileCount}.txt 2> ${output}/errs/Job_${fileCount}.txt" >> $submit
done < fileList.txt
if (( $fileCount % $filesPerJob != 0 )); then
    #fileList="${fileList%,}" #remove trailing comma from fileList
    #echo "cmsRun ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/multilep.py inputFile=$fileList outputFile=${output}/${skim}_Job_${jobCount}.root events=-1 > ${output}/logs/Job_${jobCount}.txt 2> ${output}/errs/Job_${jobCount}.txt" >> $submit
    #submit job
    submitJob $submit
    #cat $submit
fi
#remove temporary files
rm $submit
rm fileList.txt
