#!/bin/bash
echo "Killing all running crab jobs listed in ./crab directory"
cd ./crab
for dir in *
    do if [[ -d $dir ]]; then
        cd $dir
        for dir2 in *
            do if [[ -d ${dir2} ]]; then
                cd ${dir2}
                for dir3 in *
                    do if [[ -d ${dir3} ]]; then
                        cd ${dir3}
                        crab kill --dir=./
                        cd ..
                    fi
                done
                cd ..
            fi
        done
        cd ..
    fi
done
