import sys
import FWCore.ParameterSet.Config as cms

# Default arguments
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/tZq_W_lept_Z_hadron_4f_ckm_NLO_13TeV_amcatnlo_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0CAF6064-664E-E711-ADA0-1CB72C1B63C2.root'
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/TTGamma_Dilept_TuneCUETP8M2T4_13TeV-madgraph-fsrup-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/02748F9D-F640-E711-BF1F-D067E5F91695.root'
inputFile       = '/store/mc/RunIISummer16MiniAODv2/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/A27FACBE-7659-E711-9ABB-14DDA90900BB.root'
isData          = not 'SIM' in inputFile
nEvents         = 1000
outputFile      = 'ttg.root'     # trilep    --> skim three leptons (basic pt/eta criteria)
                                 # dilep     --> skim two leptons
                                 # singlelep --> skim one lepton
                                 # ttg       --> skim two leptons + one photon
                                 # fakerate  --> not implemented

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "isData"            in sys.argv[i]: isData          = (getVal(sys.argv[i]) == "True")
    elif "outputFile"      in sys.argv[i]: outputFile      = getVal(sys.argv[i])
    elif "inputFile"       in sys.argv[i]: inputFile       = getVal(sys.argv[i])
    elif "events"          in sys.argv[i]: nEvents         = int(getVal(sys.argv[i]))


process = cms.Process("BlackJackAndHookers")

# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Unscheduled mode makes no difference for us
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")))
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if isData: process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
else:      process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'


#
# Vertex collection
#
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')
process.goodOfflinePrimaryVertices.src    = cms.InputTag('offlineSlimmedPrimaryVertices')
process.goodOfflinePrimaryVertices.filter = cms.bool(True) # This was false in the old Majorana code, why?


#
# Import some objectsequences sequence (details in cff files)
#
from heavyNeutrino.multilep.jetSequence_cff import addJetSequence
from heavyNeutrino.multilep.egmSequence_cff import addElectronAndPhotonSequence
addJetSequence(process, isData)
addElectronAndPhotonSequence(process)



#
# Read additional MET filters not stored in miniAOD
#
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
for module in [process.BadPFMuonFilter, process.BadChargedCandidateFilter]:
  module.muons        = cms.InputTag("slimmedMuons")
  module.PFCandidates = cms.InputTag("packedPFCandidates")

# Propagate JEC to MET (need to add fullPatMetSequence to path) (maybe good to add here link to a twiki page, if it exist)
# from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD   # currently broken
from heavyNeutrino.multilep.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, isData=isData)


# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffectiveAreas           = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
  electrons                     = cms.InputTag("slimmedElectrons"),
  electronsEffectiveAreas       = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'), # WARNING this is spring 15, following SUSY-standard, i.e. not the most up-to-date values
  electronsMva                  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
  electronsMvaHZZ               = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
  electronsCutBasedVeto         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  electronsCutBasedLoose        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  electronsCutBasedMedium       = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  electronsCutBasedTight        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  photons                       = cms.InputTag("slimmedPhotons"),
  photonsChargedEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt'),
  photonsNeutralEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt'),
  photonsPhotonsEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
  photonsCutBasedLoose          = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
  photonsCutBasedMedium         = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
  photonsCutBasedTight          = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),
  photonsMva                    = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
  photonsChargedIsolation       = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
  photonsNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
  photonsPhotonIsolation        = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
  taus                          = cms.InputTag("slimmedTaus"),
  packedCandidates              = cms.InputTag("packedPFCandidates"),
  rhoCentralNeutral             = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
  rhoAll                        = cms.InputTag("fixedGridRhoFastjetAll"),
  met                           = cms.InputTag("slimmedMETs"),
  jets                          = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt"),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  recoResults                   = cms.InputTag("TriggerResults::RECO"),
  badPFMuonFilter               = cms.InputTag("BadPFMuonFilter"),
  badChargedCandFilter          = cms.InputTag("BadChargedCandidateFilter"),
  skim                          = cms.untracked.string(outputFile.split('.')[0]),
  isData                        = cms.untracked.bool(isData)
)


if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(filename = '../data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()


process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.BadPFMuonFilter *
                     process.BadChargedCandidateFilter *
                     process.egmSequence *
                     process.jetSequence *
                     process.fullPatMetSequence *
                     process.blackJackAndHookers)

