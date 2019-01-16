import sys, os
import FWCore.ParameterSet.Config as cms

# Default input file (could be overwritten by parameters given on the command line and by crab), some examples:
#inputFile      = 'file:///pnfs/iihe/cms/ph/sc4/store/data/Run2017F/DoubleMuon/MINIAOD/17Nov2017-v1/70000/E4B6F7A1-7BDE-E711-8C42-02163E019DE8.root'
#inputFile      = "root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAOD/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/02041699-0BFB-E711-AAD4-FA163E965751.root"
#inputFile      = '/store/data/Run2017F/SingleMuon/MINIAOD/17Nov2017-v1/00000/3E7C07F9-E6F1-E711-841A-0CC47A4C8E46.root'
#inputFile      = '/store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v3/000/316/569/00000/0085320B-4E64-E811-A2D3-FA163E2A55D6.root'
#inputFile      = '/store/data/Run2018A/MET/MINIAOD/PromptReco-v3/000/316/666/00000/0CC8EDCD-FD64-E811-BCA8-02163E01A020.root'
inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root'

# Other default arguments
nEvents         = 1000
extraContent    = ''
outputFile      = 'noskim.root' # trilep    --> skim three leptons (basic pt/eta criteria)
                                # displtrilep --> skim three leptons (basic pt/eta criteria)
                                # dilep     --> skim two leptons
                                # singlelep --> skim one lepton
                                # ttg       --> skim two light leptons + one photon
                                # singlejet --> one jet
                                # FR        --> one jet and one light lepton

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "outputFile"     in sys.argv[i]: outputFile   = getVal(sys.argv[i])
    elif "inputFile"    in sys.argv[i]: inputFile    = getVal(sys.argv[i])
    elif "extraContent" in sys.argv[i]: extraContent = getVal(sys.argv[i])
    elif "events"       in sys.argv[i]: nEvents      = int(getVal(sys.argv[i]))

isData = not ('SIM' in inputFile or 'HeavyNeutrino' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile
is2018 = "Run2018" in inputFile or "18MiniAOD" in inputFile
isSUSY = "SMS-T" in inputFile

process = cms.Process("BlackJackAndHookers")

# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")), duplicateCheckMode = cms.untracked.string('noDuplicateCheck'))
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if is2018 and 'PromptReco' in inputFile: process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v11'
elif is2018:                             process.GlobalTag.globaltag = '102X_dataRun2_Sep2018Rereco_v1' if isData else '102X_upgrade2018_realistic_v12'
elif is2017:                             process.GlobalTag.globaltag = '94X_dataRun2_v11'               if isData else '94X_mc2017_realistic_v17'
else:                                    process.GlobalTag.globaltag = '94X_dataRun2_v10'               if isData else '94X_mcRun2_asymptotic_v3'

#
# TrackingComponentsRecord 
#
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi')

#
# Vertex collection
#
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')
process.goodOfflinePrimaryVertices.src    = cms.InputTag('offlineSlimmedPrimaryVertices')
process.goodOfflinePrimaryVertices.filter = cms.bool(False)                          #Don't use any EDFilter when relying on hCounter!

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#
# Import some objectsequences sequence (details in cff files)
#
from heavyNeutrino.multilep.jetSequence_cff import addJetSequence
addJetSequence(process, isData, is2017, is2018)

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if is2018:   setupEgammaPostRecoSeq(process, runEnergyCorrections=False, era='2018-Prompt')      # No scale and smearings available yet
elif is2017: setupEgammaPostRecoSeq(process, runEnergyCorrections=True,  era='2017-Nov17ReReco') # Rerun scale and smearings for shiftscale bug
else:        setupEgammaPostRecoSeq(process, runEnergyCorrections=False, era='2016-Legacy')      # Default scale and smearings are ok

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  allowMatchingToAllIds         = cms.bool(False),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffectiveAreas           = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'), # TODO: check if muon POG has updates on effective areas
  muonsEffectiveAreasFall17     = cms.FileInPath('heavyNeutrino/multilep/data/effAreas_cone03_Muons_Fall17.txt'), # TODO
  electrons                     = cms.InputTag("slimmedElectrons"),
  electronsEffectiveAreas       = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'), # Recommended, used by standard IDs (the difference with the outdated effective areas is typically small)
  leptonMvaWeightsMuSUSY16      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_SUSY16_BDTG.weights.xml"), # TODO: clean-up old trainings here?
  leptonMvaWeightsEleSUSY16     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_SUSY16_BDTG.weights.xml"),
  leptonMvaWeightsMuttH16       = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH16_BDTG.weights.xml"),
  leptonMvaWeightsElettH16      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH16_BDTG.weights.xml"),
  leptonMvaWeightsMuSUSY17      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_SUSY17_BDTG.weights.xml"),
  leptonMvaWeightsEleSUSY17     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_SUSY17_BDTG.weights.xml"),
  leptonMvaWeightsMuttH17       = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH17_BDTG.weights.xml"),
  leptonMvaWeightsElettH17      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH17_BDTG.weights.xml"),
  leptonMvaWeightsEletZqTTV16   = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV16_BDTG.weights.xml"),
  leptonMvaWeightsMutZqTTV16    = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV16_BDTG.weights.xml"),
  leptonMvaWeightsEletZqTTV17   = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV17_BDTG.weights.xml"),
  leptonMvaWeightsMutZqTTV17    = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV17_BDTG.weights.xml"),
  photons                       = cms.InputTag("slimmedPhotons"),
  photonsChargedEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt'),
  photonsNeutralEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt'),
  photonsPhotonsEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt'),
  taus                          = cms.InputTag("slimmedTaus"),
  packedCandidates              = cms.InputTag("packedPFCandidates"),
  rho                           = cms.InputTag("fixedGridRhoFastjetAll"),
  met                           = cms.InputTag("slimmedMETs"),
  jets                          = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jecUncertaintyFile16          = cms.FileInPath("heavyNeutrino/multilep/data/JEC/Summer16_07Aug2017_V9_MC_Uncertainty_AK4PFchs.txt"),
  jecUncertaintyFile17          = cms.FileInPath("heavyNeutrino/multilep/data/JEC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt"), # TODO: add 2018
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  recoResultsPrimary            = cms.InputTag("TriggerResults::PAT"),
  recoResultsSecondary          = cms.InputTag("TriggerResults::RECO"),
  triggerObjects                = cms.InputTag("slimmedPatTrigger"),   # displaced specific --> this was "selectedPatTrigger" in the displaced branch but does not work for the tests TODO please check
  SingleEleTriggers             = cms.vstring(),                       # displaced specific
  SingleMuoTriggers             = cms.vstring(),                       # displaced specific
  SingleEleTriggers2017         = cms.vstring(),                       # displaced specific
  SingleMuoTriggers2017         = cms.vstring(),                       # displaced specific
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  is2018                        = cms.untracked.bool(is2018),
  isSUSY                        = cms.untracked.bool(isSUSY),
  storeLheParticles             = cms.untracked.bool('storeLheParticles' in extraContent),
)

def getJSON(is2017, is2018):
    if is2018:   return "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
    elif is2017: return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else:        return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"


## Single triggers for matching
if 'FR' in outputFile:
    process.blackJackAndHookers.SingleEleTriggers.extend(["HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*"])
    process.blackJackAndHookers.SingleMuoTriggers.extend(["HLT_Mu3_PFJet40_v*", "HLT_Mu8_v*", "HLT_Mu17_v*"])
    process.blackJackAndHookers.SingleEleTriggers2017.extend(["HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*"])
    process.blackJackAndHookers.SingleMuoTriggers2017.extend(["HLT_Mu3_PFJet40_v*", "HLT_Mu8_v*", "HLT_Mu17_v*"])
else:
    process.blackJackAndHookers.SingleEleTriggers.extend(["HLT_Ele27_WPTight_Gsf_v*"])
    process.blackJackAndHookers.SingleMuoTriggers.extend(["HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*"])
    process.blackJackAndHookers.SingleEleTriggers2017.extend(["HLT_Ele32_WPTight_Gsf_v*"])
    process.blackJackAndHookers.SingleMuoTriggers2017.extend(["HLT_IsoMu27_v*"])


if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  jsonDir = os.path.expandvars('$CMSSW_BASE/src/heavyNeutrino/multilep/data/JSON')
  process.source.lumisToProcess = LumiList.LumiList(filename = os.path.join(jsonDir, getJSON(is2017, is2018))).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.egammaPostRecoSeq *
                     process.jetSequence *
                     process.fullPatMetSequence *
                     process.blackJackAndHookers)
