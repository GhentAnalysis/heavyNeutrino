# Setup script for branch: master
RELEASE=CMSSW_9_4_12

# If the release is already available using cmsenv, use it, otherwise set up a new one
if [[ $CMSSW_BASE == *$RELEASE ]]; then
  echo "Setting up heavyNeutrino package in current release: $CMSSW_BASE"
  cd $CMSSW_BASE/src
else
  scram project CMSSW $RELEASE
  cd $RELEASE/src
  eval `scram runtime -sh`
  echo "Creating release for heavyNeutrino package: $CMSSW_BASE"
fi

# The git commands
git cms-init
git clone https://github.com/GhentAnalysis/heavyNeutrino
git cms-merge-topic cms-met:METFixEE2017_949_v2 #for EE noise fix of 2017 MET

# Compile and move into package
scram b -j 10
cd $CMSSW_BASE/src/heavyNeutrino
echo "Setup finished"
