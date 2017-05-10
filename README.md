# heavyNeutrino
New repository for the heavyNeutrino analysis

# Set-up instructions
```
cmsrel CMSSW_8_0_27
cd CMSSW_8_0_27/src
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
./submitAll.py <datasetsFile>
```
