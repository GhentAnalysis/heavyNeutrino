import sys, os
import FWCore.ParameterSet.Config as cms

# Default input file (could be overwritten by parameters given on the command line and by crab), some examples:
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-105To160_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/2E242480-5C0D-E911-B9A6-90E2BACBAA90.root"
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017RECOPF_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0A1754A2-256F-E811-AD07-6CC2173CAAE0.root"
inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root'

# Other default arguments
nEvents         = 1000
extraContent    = ''
outputFile      = 'noskim.root' # trilep    --> skim three leptons (basic pt/eta criteria)
                                # dilep     --> skim two leptons
                                # singlelep --> skim one lepton
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

# Print a warning if a different release is used as the one in the setup script (i.e. probably something will be broken)
os.system('if [[ "$(grep RELEASE= $CMSSW_BASE/src/heavyNeutrino/setup.sh)" != *"$CMSSW_VERSION"* ]];then echo ">>> WARNING: you are using a different release as the one specified in the setup.sh script! <<< ";fi')

# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")))
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if is2018 and 'PromptReco' in inputFile: process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v13'
elif is2018:                             process.GlobalTag.globaltag = '102X_dataRun2_Sep2018ABC_v2' if isData else '102X_upgrade2018_realistic_v18'
elif is2017:                             process.GlobalTag.globaltag = '94X_dataRun2_v11'            if isData else '94X_mc2017_realistic_v17'
else:                                    process.GlobalTag.globaltag = '94X_dataRun2_v10'            if isData else '94X_mcRun2_asymptotic_v3'

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
if is2018:   jecUncertaintyFile = 'Autumn18_V8_MC_Uncertainty_AK4PFchs.txt'
elif is2017: jecUncertaintyFile = 'Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt'
else:        jecUncertaintyFile = 'Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt'

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if is2018:   setupEgammaPostRecoSeq(process, runEnergyCorrections=True,  era='2018-Prompt')      # Updated scale and smearings
elif is2017: setupEgammaPostRecoSeq(process, runEnergyCorrections=True,  era='2017-Nov17ReReco') # Rerun scale and smearings for shiftscale bug
else:        setupEgammaPostRecoSeq(process, runEnergyCorrections=False, era='2016-Legacy')      # Default scale and smearings are ok

#
# L1 prefiring (only needed for 2016/2017, use empty sequence for 2018)
#
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
if not is2018:
  process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
      DataEra                      = cms.string("2017BtoF" if is2017 else "2016BtoH"),
      UseJetEMPt                   = cms.bool(False),
      PrefiringRateSystematicUncty = cms.double(0.2),
      SkipWarnings                 = False
  )
else:
  process.prefiringweight = cms.Sequence()

#
# For the particleLevelProducer (useful for rivet implementation and/or unfolding)
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/ParticleLevelProducer
#
if 'storeParticleLevel' in extraContent and not isData:
  process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
  process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
  process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
  process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
  process.genParticles2HepMC.signalParticlePdgIds = cms.vint32(6,-6) # for top analyses, though not yet sure what it exactlye does, I think it is only relevant to find the signal vertex which we currently do not save
  process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
  process.particleLevelSequence = cms.Sequence(process.mergedGenParticles * process.genParticles2HepMC * process.particleLevel)
else:
  process.particleLevelSequence = cms.Sequence()

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  particleLevelLeptons          = cms.InputTag("particleLevel:leptons"),
  particleLevelPhotons          = cms.InputTag("particleLevel:photons"),
  particleLevelJets             = cms.InputTag("particleLevel:jets"),
  particleLevelMets             = cms.InputTag("particleLevel:mets"),
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
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintyFile),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  recoResultsPrimary            = cms.InputTag("TriggerResults::PAT"),
  recoResultsSecondary          = cms.InputTag("TriggerResults::RECO"),
  triggerObjects                = cms.InputTag("slimmedPatTrigger"),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  is2018                        = cms.untracked.bool(is2018),
  isSUSY                        = cms.untracked.bool(isSUSY),
  storeLheParticles             = cms.untracked.bool('storeLheParticles' in extraContent),
  storeParticleLevel            = cms.untracked.bool('storeParticleLevel' in extraContent),
)

def getJSON(is2017, is2018):
    if is2018:   return "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
    elif is2017: return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else:        return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  jsonDir = os.path.expandvars('$CMSSW_BASE/src/heavyNeutrino/multilep/data/JSON')
  process.source.lumisToProcess = LumiList.LumiList(filename = os.path.join(jsonDir, getJSON(is2017, is2018))).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.egammaPostRecoSeq *
                     process.jetSequence *
                     process.fullPatMetSequence *
                     process.prefiringweight *
                     process.particleLevelSequence *
                     process.blackJackAndHookers)
