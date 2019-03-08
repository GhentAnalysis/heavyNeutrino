# heavyNeutrino
New repository for Ghent CMS analysis tuplizer framework (initially created for heavy neutrino analysis).
In order to have a working copy, following the instructions below which sets up the needed CMSSW release and dependency packages.
The tuplizer itself is entirely contained in the multilep directory.
**Please take a moment to read the [wiki](https://github.com/GhentAnalysis/heavyNeutrino/wiki) with useful guidelines**

# Set-up instructions (displacedLegacy branch)
```
wget https://raw.githubusercontent.com/GhentAnalysis/heavyNeutrino/displacedLegacy/setup.sh
source setup.sh
```

# The cmsRun executable
```
cd heavyNeutrino/multilep/test
cmsRun multilep.py
```

# Running test jobs
Preferably ran before an important git push:
```
cd heavyNeutrino/multilep/test/testing
./runTests.py
```

# Mass production of tuples
```
cd heavyNeutrino/multilep/test
./submitAll.py <datasetsFile>  (<run locally>) (<files per job>)
```
