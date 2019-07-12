#!/bin/bash
echo "Getting status of all jobs in ./crab directory"
if [[ -d $1 ]]; then
    cd $1
    for dir2 in *
        do if [[ -d ${dir2} ]]; then
            cd ${dir2}
            for dir3 in *
                do if [[ -d ${dir3} ]]; then
                    cd ${dir3}
                    crab status --dir=./
                    echo -e "\n\n\n"
                    cd ..
                fi
            done
            cd ..
        fi
    done
    cd ..
fi
