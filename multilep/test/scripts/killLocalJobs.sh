#!/bin/bash
temp="jobs.txt"
if [[ -e $temp ]]; then
    rm $temp
fi
touch $temp
qstat -u bvermass > $temp
while IFS='' read -r line || [[ -n "$line" ]]; do
    #echo $line 
    if echo "$line" | grep -q "bvermass"; then
        jobnumber=$(echo ${line}| cut -d'.' -f 1)
        qdel "$jobnumber" 
    fi
done < "$temp"
rm $temp
