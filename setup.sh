# Setup for new EGamma ID code
eval `scram runtime -sh`
cd $CMSSW_BASE/src
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
scram b -j 10
