import sys
import FWCore.ParameterSet.Config as cms

# Default arguments
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00A9113F-15D6-E611-9142-047D7B881D3A.root'
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/003658EE-77E6-E611-ACB1-7CD30ABD295A.root'
#inputFile       = '/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/100000/52779EE0-F4ED-E611-BF87-70106F49CD3C.root'
inputFile       = "root://cmsxrootd.fnal.gov///store/data/Run2017C/MuonEG/MINIAOD/PromptReco-v3/000/300/780/00000/86494C82-EA7E-E711-ACCC-02163E01441B.root"
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/prompt/HeavyNeutrino_trilepton_M-100_V-0.01_2l_NLO/heavyNeutrino_1.root'
isData          = not ('SIM' in inputFile or 'HeavyNeutrino' in inputFile)
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
  photonsFull5x5SigmaIEtaIPhi   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIPhi"),
  taus                          = cms.InputTag("slimmedTaus"),
  packedCandidates              = cms.InputTag("packedPFCandidates"),
  rho                           = cms.InputTag("fixedGridRhoFastjetAll"),
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
  data2017 = "Run2017" in inputFile
  if data2017: JSON = "../data/JSON/Cert_294927-302663_13TeV_PromptReco_Collisions17_JSON.txt"
  else:        JSON = "../data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
  process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()


process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.BadPFMuonFilter *
                     process.BadChargedCandidateFilter *
                     process.egmSequence *
                     process.jetSequence *
                     process.fullPatMetSequence *
                     process.blackJackAndHookers)

