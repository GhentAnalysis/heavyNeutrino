# Setup for new EGamma ID code
eval `scram runtime -sh`
cd $CMSSW_BASE/src
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
scram b -j 10

# Setting up new photon ID
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/src
scram b -j10
