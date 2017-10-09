if [[ $CMSSW_BASE == *CMSSW_8* ]]
then
  # Setup for new electron and photon MVA
  eval `scram runtime -sh`
  cd $CMSSW_BASE/src
  git cms-merge-topic ikrav:egm_id_80X_v3_photons
  scram b -j 10
fi
# nothing to do for CMSSW_9_X
