
cmsrel CMSSW_9_4_12
cd CMSSW_9_4_12/src
cmsenv
git cms-init

#for EE noise fix of 2017 MET 
git cms-merge-topic cms-met:METFixEE2017_949_v2
scram b -j 10

cd heavyNeutrino
scram b -j 10
