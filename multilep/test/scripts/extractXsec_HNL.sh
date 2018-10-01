#!/bin/bash

#calculate cross sections for all HNL samples in a directory at once
#normal use: "./extractXsec.sh HNL_list.txt"
outputpath="/user/bvermass/heavyNeutrino/Dileptonprompt/CMSSW_9_4_0/src/samesign_Analysis/xsecs_HNL/" #output in samesign_analysis directory

#where all HNL cross sections will be put
xsec_list="${outputpath}HNL_Xsec_List.txt"
> $xsec_list
#sample path of collection of HNL samples (use HNL_list.txt)
#line="/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/prompt/"
while IFS='' read -r line || [[ -n "$line" ]]; do #read line per line from input text file
    for D in ${line}HeavyNeutrino_lljj*; do
        if [ -d "${D}" ]; then
            filename=$(echo ${D}| cut -d'/' -f 11) #get name of HNL sample
            filename=${filename}"_"$(echo ${D}| cut -d'/' -f 10)"_"$(echo ${D}| cut -d'/' -f 9) #attach prompt/displaced and production period to name
            output=${outputpath}"log_xsec_"${filename}".txt"
            echo "output can be found in :${output}"
            if [[ -e ${output} ]]; then
                rm ${output}
            fi
            ./extractXsec.sh ${D} ${output}
            echo -e "${filename} \t \t \t \t \t \t $(grep -F "Before matching" ${output} | cut -d' ' -f 7)" >> ${xsec_list}
        fi
    done
done < "$1"
