if [[ $CMSSW_BASE == *CMSSW_8* ]]
then
  eval `scram runtime -sh`
  cd $CMSSW_BASE/src

  # First EGM smearing
  git cms-merge-topic cms-egamma:EGM_gain_v1
  cd EgammaAnalysis/ElectronTools/data
  git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git 
  cd -

  # Setup for new electron and photon MVA
  git cms-merge-topic ikrav:egm_id_80X_v3_photons
  scram b -j 10
else
  echo "Wrong branch! Setup script is made for CMSSW_8_X"
fi
