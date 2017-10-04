# heavyNeutrino
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
This is a temporary version trying to implement all necessary changes for running on 2017 data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
New repository for Ghent CMS analysis tuplizer framework. (Initially created for heavy neutrino analysis)
In order to have a working copy, one needs to set up the CMSSW release as instructed below and run the ./setup.sh script in order to have additional packages downloaded.
The tuplizer itself is entirely contained in the multilep directory.
# Set-up instructions
```
cmsrel CMSSW_9_2_4 #Warning: CMSSW_9_2_3 and lower will not work due to bug in packedCandidate::hasTrackDetails() function (is non-const)
#Warning: framework does not work in CMSSW_9_2_12_patch1 for unknown reasons
cd CMSSW_9_2_4/src
cmsenv
git cms-init
git clone https://github.com/GhentAnalysis/heavyNeutrino
cd heavyNeutrino
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
./submitAll.py <datasetsFile>
```
