# Setup script for branch: master
BRANCH=master

# Release: take CMSSW_10_6_X as default (works for both UL and older reprocessings), fall back to CMSSW_10_2_X on T2_BE_IIHE
if [[ $SCRAM_ARCH == *slc6* ]]; then
  RELEASE=CMSSW_10_2_25
else
  RELEASE=CMSSW_10_6_19_patch2
fi

# If the release is already available using cmsenv, use it, otherwise set up a new one
if [[ $CMSSW_BASE == *$RELEASE ]] && [[ -d $CMSSW_BASE ]]; then
  echo "Setting up heavyNeutrino package in current release: $CMSSW_BASE"
  cd $CMSSW_BASE/src
else
  scram project CMSSW $RELEASE
  cd $RELEASE/src
  eval `scram runtime -sh`
  echo "Creating release for heavyNeutrino package: $CMSSW_BASE"
fi

# The git commands
gitUser=$(git config --get user.github)
remoteLs=$(git ls-remote https://github.com/$gitUser/heavyNeutrino 2>&1)
if [[ "$remoteLs" = *'Repository not found'* ]]; then
  echo "WARNING: You do not have a heavyNeutrino fork yet, checking out the GhentAnalysis remote instead"
  gitUser="GhentAnalysis"
fi
git cms-init
git clone https://github.com/$gitUser/heavyNeutrino
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git EgammaUser/EgammaPostRecoTools
git clone https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer.git TopQuarkAnalysis/BFragmentationAnalyzer
cd $CMSSW_BASE/src/heavyNeutrino
git checkout --track origin/$BRANCH
if [[ "$gitUser" != "GhentAnalysis" ]]; then
  git remote add ghent https://github.com/GhentAnalysis/heavyNeutrino
  git fetch ghent
  git rebase ghent/$BRANCH
fi

# Compile and move into package
cd $CMSSW_BASE
scram b -j 10
cd $CMSSW_BASE/src/heavyNeutrino
echo "Setup finished"
