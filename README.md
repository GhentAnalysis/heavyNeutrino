# heavyNeutrino
New repository for Ghent CMS analysis tuplizer framework (initially created for heavy neutrino analysis).
In order to have a working copy, one needs to set up the CMSSW release as instructed below.
The tuplizer itself is entirely contained in the multilep directory.

# Set-up instructions for branch *displaced*
```
cmsrel CMSSW_9_4_3
cd CMSSW_9_4_3/src
cmsenv
git cms-init
git clone https://github.com/GhentAnalysis/heavyNeutrino
cd heavyNeutrino
git checkout displaced
./setup.sh
```

# Running a test job
```
cd heavyNeutrino/multilep/test
cmsRun multilep.py
```

# Mass production of tuples
```
cd heavyNeutrino/multilep/test
./submitAll.py <datasetsFile>  (<run locally>) (<files per job>)
```
