#!/bin/bash

#This function will make a list of all files present in a given sample
fileList(){
    #sample path given as input
    input=$1
    
    #initialize file which will contain all files in the sample
    if [[ -e fileList.txt ]]; then
        rm fileList.txt
    fi
    touch fileList.txt
    
    #check if sample is private or official production
    if [[ input == *"/user/"* ]]; then  #private sample

        #add all files to list 
        for file in $input/*
            do echo "$file" >> fileList.txt
        done

    else                                #official sample

        #use CMSDAS query to list all files present in sample in a txt file
        python ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/scripts/das_client.py --query="file dataset=$input" --limit=0 > fileNames.txt
        
        #location of files on local cluster
        location=/pnfs/iihe/cms/ph/sc4
        #T2-BE-IIHE redirector
        dcap=dcap://maite.iihe.ac.be
        
        #xrootd redirector for files not present on local cluster
        redirector=root://cmsxrootd.fnal.gov//
        
        while read f
            #check if file is locally available
            do if [[ -e ${location}${f} ]]
                then echo "${dcap}${location}${f}" >> fileList.txt
            #if the file is not available add an xrootd redirector to a remote location
            else
                echo "${redirector}${f}" >> fileList.txt
            fi
        done < fileNames.txt

        #clean up temporary file
        rm fileNames.txt

    fi
}
