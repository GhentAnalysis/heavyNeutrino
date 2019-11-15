#!/bin/bash
#See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial#Re_analyzing_the_missing_luminos for more info on the method to resubmit notFinishedLumis.json
#Method:
#Change campaign name to xx_missingLumis
#Change lumiMask to notFinishedLumis.json from the unfinished campaign
#submitAll for same data files
echo "Making reports for a specific crab campaign"
if [[ -d $1 ]]; then
    cd $1
    for dir in *
        do if [[ -d $dir ]]; then
            cd $dir
            for dir2 in *
                do if [[ -d $dir2 ]]; then
                    cd $dir2
                    crab report -d .
                    echo -e "\n\n\n"
                    cd ..
                fi
            done
            cd ..
        fi
    done
    cd ..
fi
