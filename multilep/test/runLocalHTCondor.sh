#!/bin/bash

#
# The script needs a proxy, so make sure to have your preferred way of accessing the proxy here!
#
if [[ $USER == "tutran" ]]; then     proxy=/user/tutran/private/x509up_u23068
elif [[ $USER == "wverbeke" ]]; then proxy=/user/wverbeke/x509up_u20640
elif [[ $USER == "gmestdac" ]]; then proxy=/user/$USER/production/proxyExpect.sh
elif [[ $USER == "kskovpen" ]]; then proxy=/tmp/x509up_u20657
elif [[ $USER == "lwezenbe" ]]; then proxy=/user/lwezenbe/x509up_u20675
elif [[ $USER == "tomc" ]]; then     proxy=/user/$USER/production/proxyExpect.sh
elif [[ $USER == "nivanden" ]]; then proxy=/user/nivanden/x509up_u23145
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
datasetName=$(echo "$2" | sed 's/\/.*//')
skim=$3
filesPerJob=$4
extraContent=$5


#include function to make list of all files in given sample
source scripts/makeFileList.sh

#function to set up CMSSW in a job
setCMSSW(){
    echo "#!/bin/bash" >> $1
    echo "cd ${CMSSW_BASE}/src" >> $1
    echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> $1
    echo "eval \`scram runtime -sh\`" >> $1
    echo "cd heavyNeutrino/multilep/test/" >> $1
    echo "export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/python/2.7.15-pafccj" >> $1
}



#function to submit a job and catch invalid credentials
submitJob(){
    if [ -e HTCondorSubmission/${datasetName}_localSubmission.sub ]; then rm HTCondorSubmission/${datasetName}_localSubmission.sub; fi

    touch HTCondorSubmission/${datasetName}_localSubmission.sub
    if [[ ! -d "HTCondorSubmission/logs" ]]; then
      mkdir "HTCondorSubmission/logs"
      mkdir "HTCondorSubmission/err"
      mkdir "HTCondorSubmission/out"
    fi   

    echo "universe    =   vanilla" >> HTCondorSubmission/${datasetName}_localSubmission.sub
    echo "executable  =   HTCondorSubmission/${datasetName}_localSubmission\$(Process).sh" >> HTCondorSubmission/${datasetName}_localSubmission.sub

    dir=$(pwd)
    echo "log         =   $dir/HTCondorSubmission/logs/${datasetName}_localSubmission_\$(ClusterId)_\$(Process).log" >> HTCondorSubmission/${datasetName}_localSubmission.sub
    echo "error       =   $dir/HTCondorSubmission/err/${datasetName}_localSubmission_\$(ClusterId)_\$(Process).err" >> HTCondorSubmission/${datasetName}_localSubmission.sub
    echo "output      =   $dir/HTCondorSubmission/out/${datasetName}_localSubmission_\$(ClusterId)_\$(Process).out" >> HTCondorSubmission/${datasetName}_localSubmission.sub

    echo "queue $1" >> HTCondorSubmission/${datasetName}_localSubmission.sub


    condor_submit HTCondorSubmission/${datasetName}_localSubmission.sub
}

#make list of all files in input sample
fileList $input

if [[ ! -d "HTCondorSubmission" ]]; then
    mkdir "HTCondorSubmission"
fi        

#loop over new list of files and submit jobs
fileCount=0
submitCount=-1
jobCount=0
submit=""

while read f; do
    fileCount=$((fileCount + 1))
    #submit a job for every few files, as specified in the input
    if (( $fileCount % $filesPerJob == 0 )) || (( $fileCount == 1 ))
        then if (( $fileCount % $filesPerJob == 0 )); then
            jobCount=$((jobCount + 1))
        fi
        #initialize temporary submission script
        submitCount=$((submitCount + 1))

        submit="HTCondorSubmission/${datasetName}_localSubmission${submitCount}.sh"
        
        if [ -e $submit ]; then rm $submit; fi
        
        touch $submit
        chmod +x $submit
        #initialize CMSSW environment in submission script
        setCMSSW $submit
        exportProxy $submit
    fi
    pnfsDir="/pnfs/iihe/cms/store/user/$USER/heavyNeutrino/$output"
    outputFile="${skim}_${fileCount}.root"
    logFile="logs/${skim}_${fileCount}.txt"
    errFile="errs/${skim}_${fileCount}.txt"

    mkdir -p ${pnfsDir}
    mkdir -p ${pnfsDir}/errs
    mkdir -p ${pnfsDir}/logs
    echo "cmsRun ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/multilep.py inputFile=file://$f outputFile=$pnfsDir/$outputFile events=-1 ${extraContent} > $pnfsDir/$logFile 2> $pnfsDir/$errFile" >> $submit

done < fileList.txt

submitCount=$((submitCount + 1))
submitJob $submitCount

rm fileList.txt
