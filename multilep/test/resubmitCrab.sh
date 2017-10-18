#!/bin/bash
checkCrab(){
    for dir in ./*/*;
        do if [[ $dir != *"done"* ]]
            then (echo $dir; ./checkCrab.py $dir)
        fi
    done
}
checkCrab $1
