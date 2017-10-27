# Testing and submitting jobs

### multilep.py
The main python config of the module is multilep.py. The file contains of following subsections:
* default arguments: the ones which are taken when a simply test job is run (i.e. just "cmsRun multilep.py")
* argument evaluation: replaces the default arguments which are given on the fly (i.e. "cmsRun multilep.py isData=True outputFile=dilep.root events=10000")
* standard cms loads and options like the messagelogger, global tag etc..
* load of goodOfflinePrimaryVertex module
* loading jet and egm sequences
* loading additional MET filters
* initializing the multilep plugin with all of its configuration parameters
* lumiList for data
* the path defining the order of all CMS modules to be run

Note that the outputFile argument is transformed into the skim argument of the multilep plugin, and skims are therefore done based on the filename given in the outputFile argument.

### submitting jobs (both on local T2 as with crab)
Run
```
./submitAll.py <datasetsFile> (<runLocal>) (<files per job>)
```
where <datasetsFile> is a file structures like test-v1.txt. The submitAll.py will submit jobs to the local T2 when a local dataset path is given (i.e. starting with /pnfs) and to crab when a DAS-style path is given. The identifier before ":" is the outputFileName without .root, and hence the skim. The name of the datasetsfile is taken as the productionlabel. 

By giving an extra argument "True" to the submission script, all files will be submitted to the local T2 cluster, and files missing locally will use xrootd redirectors to other sites so that every file in the sample is guaranteed to be processed. Finally one can give the number of files to be processed in each job, for which the default is 10.

  The outputfiles will be found at /user/$USER/public/heavyNeutrino/\<datasetName\>/local\_\<productionLabel\>/\<skim\>\_i.root. Logfiles in /user/$USER/public/heavyNeutrino/\<datasetName\>/local\_\<productionLabel\>/log/\<skim\>\_i.log.

* When running on crab, jobs are submitted using the default crab.py in this directory, filling in the productionlabel, dataset and outputfile.
  The outputfiles will be found at /pnfs/iihe/cms/store/user/$USER/heavyNeutrino/\<datasetName\>/crab\_<campaign\>\_\<productionLabel\>/\*/\*/\<skim\>\_i.root

* When jobs are submitted with crab, one can check the crab status on the subdirectories in the newly created ./crab/\<productionLabel\>/. A tool (on T2\_BE\_IIHE) to check the crab status and automatically resubmit failed jobs is /user/tomc/public/scripts/crabStatus.py <crabDirectory>. A few problematic error codes are which typically not go away by simply resubmitting:
  - *50664* (walltime exceeded) --> To avoid this one can lower the unitsPerJob in crab.py. A lower value of unitsPerJob is used for data because processing time per lumi is much larger there
  - *50660* (memory exceeded) --> Even though we have (at moment of writing) no memory leaks in the heavyNeutrino/multilep module, inefficiencies in CMSSW or other modules could raise the needed memory for some events.
                                 Could happen randomly.
                                 Best solved by resubmitting with the --maxmemory option as is done automatically by the crabStatus.py script.
  - *60318* (stageout) --> Probably a problem with the T2 storage, try again
  - *10034* (release not available) --> Check for more recent releases in the same cycle
