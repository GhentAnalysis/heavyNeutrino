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
./submitAll.py <datasetsFile>
```
where <datasetsFile> is a file structures like test-v1.txt. The submitAll.py will submit jobs to the local T2 when a local dataset path is given (i.e. starting with /pnfs) and to crab when a DAS-style path is given. The identifier before ":" is the outputFileName without .root, and hence the skim. The name of the datasetsfile is taken as the productionlabel.

* When running on the local Belgium T2, the submitAll.py script will look for all available .root files within the given path (i.e. make sure the sample is complete). For each .root file a job is sent to the T2, simply by using runOnCream02.sh as wrapper script which simply provides the inputFile, outputFile, isData, events=-1 arguments to the cmsRun multilep.py command. One can group multiple .root files by changing the groupFiles parameter in submitAll.py.
  The outputfiles will be found at /user/$USER/public/heavyNeutrino/\<datasetName\>/local\_\<productionLabel\>/\<skim\>\_i.root. Logfiles in /user/$USER/public/heavyNeutrino/\<datasetName\>/local\_\<productionLabel\>/log/\<skim\>\_i.log.

* When running on crab, jobs are submitted using the default crab.py in this directory, filling in the productionlabel, dataset and outputfile.
  The outputfiles will be found at /pnfs/iihe/cms/store/user/$USER/heavyNeutrino/\<datasetName\>/crab\_<campaign\>\_\<productionLabel\>/\*/\*/\<skim\>\_i.root

* When jobs are submitted with crab, one can check the crab status on the subdirectories in the newly created ./crab/\<productionLabel\>/
