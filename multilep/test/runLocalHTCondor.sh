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

#transfer(){
#    #echo "gfal-mkdir -p srm://maite.iihe.ac.be:8443/$2/$(dirname $3)" >> $4
#    #echo "cp -f file://$1/$3 srm://maite.iihe.ac.be:8443/$2/$3" >> $4
#    #if [[ $proxy == *.sh ]]; then # Only clean-up the local files when a script is used for the proxy, then we are sure the proxy is valid and files should be copied correctly
#    #  echo "rm $1/$3" >> $4
#    #fi
#}

#function to submit a job and catch invalid credentials
submitJob(){
    if [ -e localSubmission.sub ]; then rm localSubmission.sub; fi

    touch localSubmission.sub

    echo "universe    =   vanilla" >> localSubmission.sub
    echo "executable  =   HTCondorSubmission/localSubmission\$(Process).sh" >> localSubmission.sub

    echo "log         =   /user/nivanden/condor/logs/localSubmission_\$(ClusterId)_\$(Process).log" >> localSubmission.sub
    echo "error       =   /user/nivanden/condor/error/localSubmission_\$(ClusterId)_\$(Process).err" >> localSubmission.sub
    echo "output      =   /user/nivanden/condor/output/localSubmission_\$(ClusterId)_\$(Process).out" >> localSubmission.sub

    # echo "x509userproxy = /user/nivanden/x509up_u23145" >> localSubmission.sub
    # echo "use_x509userproxy = True" >> localSubmission.sub

    echo "queue $1" >> localSubmission.sub


    condor_submit localSubmission.sub
    # quick change to make it run. Right now no real error checking
}

#make list of all files in input sample
fileList $input

mkdir "HTCondorSubmission"

#loop over new list of files and submit jobs
fileCount=0
submitCount=-1
submit=localSubmission0.sh
fileList=""
while read f; do
    fileCount=$((fileCount + 1))
    #submit a job for every few files, as specified in the input
    if (( $fileCount % $filesPerJob == 0 )) || (( $fileCount == 1 ))
        then if (( $fileCount % $filesPerJob == 0 )); then
            #submitJob $submit
            jobCount=$((jobCount + 1))
            fileList=""
            submitCount=$((submitCount + 1))
        fi
        #initialize temporary submission script
        submitCount=$((submitCount + 1))
        submit="HTCondorSubmission/localSubmission$submitCount.sh"
        
        if [ -e $submit ]; then rm $submit; fi
        
        touch $submit
        chmod +x $submit
        #initialize CMSSW environment in submission script
        setCMSSW $submit
        exportProxy $submit
    fi
    userDir="/user/$USER/public/heavyNeutrino/$output"
    pnfsDir="/pnfs/iihe/cms/store/user/$USER/heavyNeutrino/$output"
    outputFile="${skim}_${fileCount}.root"
    logFile="logs/${skim}_${fileCount}.txt"
    errFile="errs/${skim}_${fileCount}.txt"

    mkdir -p ${pnfsDir}
    mkdir -p ${pnfsDir}/errs
    mkdir -p ${pnfsDir}/logs
    echo "cmsRun ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/multilep.py inputFile=$f outputFile=$pnfsDir/$outputFile events=-1 ${extraContent} > $pnfsDir/$logFile 2> $pnfsDir/$errFile" >> $submit

done < fileList.txt

submitCount=$((submitCount + 1))
submitJob $submitCount

rm fileList.txt
