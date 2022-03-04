# Setup script for branch: master
BRANCH=UL_master

# Release: take CMSSW_10_6_X as default (works for both UL and older reprocessings), fall back to CMSSW_10_2_X on T2_BE_IIHE
# if [[ $SCRAM_ARCH == *slc6* ]]; then
#   RELEASE=CMSSW_10_2_25
# else
#   RELEASE=CMSSW_10_6_20
# fi

export SCRAM_ARCH=slc7_amd64_gcc820
RELEASE=CMSSW_10_6_27

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

git cms-addpkg RecoEgamma/EgammaTools  ### essentially just checkout the package from CMSSW
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools

cd $CMSSW_BASE/src/heavyNeutrino
git checkout --track origin/$BRANCH
if [[ "$gitUser" != "GhentAnalysis" ]]; then
  git remote add ghent https://github.com/GhentAnalysis/heavyNeutrino
  git fetch ghent
  git rebase ghent/$BRANCH
fi

# Install additional packages
cd $CMSSW_BASE
cp heavyNeutrino/multilep/data/tools/rabit.xml ../config/toolbox/slc7_amd64_gcc820/tools/selected/
cp heavyNeutrino/multilep/data/tools/xgboost.xml ../config/toolbox/slc7_amd64_gcc820/tools/selected/
scram setup rabit
scram setup xgboost
cmsenv

# Compile and move into package
scram b -j 10
cd $CMSSW_BASE/src/heavyNeutrino
echo "Setup finished"
