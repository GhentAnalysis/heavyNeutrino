# heavyNeutrino
New repository for Ghent CMS analysis tuplizer framework (initially created for heavy neutrino analysis).
In order to have a working copy, following the instructions below which sets up the needed CMSSW release and dependency packages.
The tuplizer itself is entirely contained in the multilep directory.

# Set-up instructions (master branch)
```
wget https://raw.githubusercontent.com/GhentAnalysis/heavyNeutrino/master/setup.sh
chmod +x setup.sh
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
