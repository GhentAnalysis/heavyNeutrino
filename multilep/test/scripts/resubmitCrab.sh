#!/bin/bash
checkCrab(){
    for dir in ${1}/*/*;
        do if [[ $dir != *"done"* ]]
            then echo $dir 
            ${CMSSW_BASE}/src/heavyNeutrino/multilep/test/crabStatus.py $dir
        fi
    done
}
if [[ -z $1 ]] 
    then echo "Empty string given, will resubmit ALL failed crab jobs"
    for d in ../crab/*
        do if [[ -d $d ]]
            then checkCrab $d
        fi
    done
else 
    checkCrab $1
fi
