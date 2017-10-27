# heavyNeutrino
New repository for Ghent CMS analysis tuplizer framework (initially created for heavy neutrino analysis).
In order to have a working copy, one needs to set up the CMSSW release as instructed below. For CMSSW\_8, the ./setup.sh script downloads some needed additional packages.
The tuplizer itself is entirely contained in the multilep directory.

# Set-up instructions for 2016 data
```
cmsrel CMSSW_8_0_30
cd CMSSW_8_0_30/src
cmsenv
git cms-init
git clone https://github.com/GhentAnalysis/heavyNeutrino
cd heavyNeutrino
./setup.sh
```

# Set-up instructions for 2017 data
```
cmsrel CMSSW_9_2_4
cd CMSSW_9_2_4/src
cmsenv
git cms-init
git clone https://github.com/GhentAnalysis/heavyNeutrino
cd heavyNeutrino
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
