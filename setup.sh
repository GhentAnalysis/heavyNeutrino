
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init

#for EE noise fix of 2017 MET 
git cms-merge-topic cms-met:METFixEE2017_949
scram b -j 10

git clone https://github.com/GhentAnalysis/heavyNeutrino
cd heavyNeutrino
scram b -j 10
