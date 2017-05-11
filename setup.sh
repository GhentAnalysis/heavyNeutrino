# Script to setup for new electron and photon MVA
eval `scram runtime -sh`
cd $CMSSW_BASE/src
git cms-merge-topic ikrav:egm_id_80X_v3_photons
scram b -j 10

# Not needed anymore?
#cd $CMSSW_BASE/external
#cd $SCRAM_ARCH
#git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
#cd data/RecoEgamma/ElectronIdentification/data
#cmsenv
#git checkout egm_id_80X_v1
#cd $CMSSW_BASE/src/heavyNeutrino
