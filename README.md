# heavyNeutrino
New repository for Ghent CMS analysis tuplizer framework. (Initially created for heavy neutrino analysis)
In order to have a working copy, one needs to set up the CMSSW release as instructed below and run the ./setup.sh script in order to have additional packages downloaded.
The tuplizer itself is entirely contained in the multilep directory.

# Set-up instructions
```
cmsrel CMSSW_8_0_30
cd CMSSW_8_0_30/src
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
