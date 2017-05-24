#PBS -l nodes=1:ppn=1
#!/bin/zsh

echo "Creating tuple for $inputFile" >&2
echo "Id: $PBS_JOBID" >&2

cd $dir 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cmsRun multilep.py inputFile=$inputFile outputFile=$outputFile events=-1 isData=$isData
