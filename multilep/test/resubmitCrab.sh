#!/bin/bash
checkCrab(){
    for dir in ./*/*;
        do if [[ $dir != *"done"* ]]
            then (echo $dir; ./checkCrab.py $dir)
        fi
    done
}
if [[ -z $1 ]] 
    echo "Empty string given, will resubmit ALL failed crab jobs"
    then for d in ./crab
        do checkCrab $d
    done
else 
    checkCrab $1
fi
