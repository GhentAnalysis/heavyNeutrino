#!/bin/bash

#
# The script needs a proxy, so make sure to have your preferred way of accessing the proxy here!
#
if [[ $USER == "tutran" ]]; then     proxy=/user/tutran/private/x509up_u23068
elif [[ $USER == "wverbeke" ]]; then proxy=/user/wverbeke/x509up_u20640
elif [[ $USER == "gmestdac" ]]; then proxy=/user/gmestdac/x509up_u20676
elif [[ $USER == "lwezenbe" ]]; then proxy=/user/lwezenbe/proxylwezenbe
elif [[ $USER == "tomc" ]]; then     proxy=/user/$USER/production/proxyExpect.sh
else
  echo "Add your proxy in RunLocal.sh before submitting jobs!"
  exit 1
fi

exportProxy(){
  if [[ $proxy == *.sh ]]; then
    echo "$proxy" >> $1
  else
    echo "export X509_USER_PROXY=$proxy" >> $1
  fi
}

DATE=`date '+%Y%m%d_%H%M%S'`
#
# Read command-line arguments, typically given by submitAll.py
#
input=$1
output="$2/$DATE/0000" # similar structure as crab
skim=$3
filesPerJob=$4
extraContent=$5


#include function to make list of all files in given sample
source scripts/makeFileList.sh

#function to set up CMSSW in a job
setCMSSW(){
    echo "cd ${CMSSW_BASE}/src" >> $1
    echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> $1
    echo "eval \`scram runtime -sh\`" >> $1
    echo "cd heavyNeutrino/multilep/test/" >> $1
}

transfer(){
    echo "gfal-mkdir -p srm://maite.iihe.ac.be:8443/$2/$(dirname $3)" >> $4
    echo "gfal-copy -f file://$1/$3 srm://maite.iihe.ac.be:8443/$2/$3" >> $4
    if [[ $proxy == *.sh ]]; then # Only clean-up the local files when a script is used for the proxy, then we are sure the proxy is valid and files should be copied correctly
      echo "rm $1/$3" >> $4
    fi
}

#function to submit a job and catch invalid credentials
submitJob(){
    qsub $1 -l walltime=40:00:00 > outputCheck.txt 2>> outputCheck.txt
    while !(grep ".cream02.iihe.ac.be" outputCheck.txt); do
        echo "Submit failed, resubmitting"
        sleep 2  #sleep 2 seconds before attemtping resubmission
        qsub $1 -l walltime=40:00:00 > outputCheck.txt 2>> outputCheck.txt
    done
    rm outputCheck.txt
}

#make list of all files in input sample
fileList $input


#loop over new list of files and submit jobs
fileCount=0
submit=localSubmission.sh
fileList=""
while read f; do
    fileCount=$((fileCount + 1))
    #submit a job for every few files, as specified in the input
    if (( $fileCount % $filesPerJob == 0 )) || (( $fileCount == 1 ))
        then if (( $fileCount % $filesPerJob == 0 )); then
            submitJob $submit
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
    userDir="/user/$USER/public/heavyNeutrino/$output"
    pnfsDir="/pnfs/iihe/cms/store/user/$USER/heavyNeutrino/$output"
    outputFile="${skim}_${fileCount}.root"
    logFile="logs/${skim}_${fileCount}.txt"
    errFile="errs/${skim}_${fileCount}.txt"
    mkdir -p ${userDir}/errs
    mkdir -p ${userDir}/logs
    echo "cmsRun ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/multilep.py inputFile=$f outputFile=$userDir/$outputFile events=-1 ${extraContent} > $userDir/$logFile 2> $userDir/$errFile" >> $submit
    exportProxy $submit
    transfer $userDir $pnfsDir $outputFile $submit
    transfer $userDir $pnfsDir $logFile $submit
    transfer $userDir $pnfsDir $errFile $submit
done < fileList.txt
if (( $fileCount % $filesPerJob != 0 )); then
    submitJob $submit
fi

#remove temporary files
rm $submit
rm fileList.txt
