import sys, os
import FWCore.ParameterSet.Config as cms

#function to return JSON file
def getJSON(is2017):
    if is2017: return "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
    else: return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

# Default arguments
outputFile      = 'trilepTree_M_2_V_0p01_m_10x.root'
inputFile       = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_100.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_103.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_105.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_101.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_102.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_104.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_107.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_108.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_109.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_10.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_106.root"

nEvents         = -1
#outputFile      = 'trilep.root'     # trilep    --> skim three leptons (basic pt/eta criteria)
                                 # dilep     --> skim two leptons
                                 # singlelep --> skim one lepton
                                 # ttg       --> skim two leptons + one photon
                                 # fakerate  --> not implemented

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "outputFile"  in sys.argv[i]: outputFile = getVal(sys.argv[i])
    elif "inputFile" in sys.argv[i]: inputFile  = getVal(sys.argv[i])
    elif "events"    in sys.argv[i]: nEvents    = int(getVal(sys.argv[i]))

isData = not ('SIM' in inputFile or 'HeavyNeutrino' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile
isSUSY = "SMS-T" in inputFile

process = cms.Process("BlackJackAndHookers")

# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Unscheduled mode makes no difference for us
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")), duplicateCheckMode = cms.untracked.string('noDuplicateCheck'))
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if   isData and is2017: process.GlobalTag.globaltag = '92X_dataRun2_2017Repro_v4'   
elif is2017:            process.GlobalTag.globaltag = '93X_upgrade2023_realistic_v3'
elif isData:            process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
else:                   process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

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

#clean 2016 data met from spurious muons and ECAL slew rate
metCollection = "slimmedMETs"
if (not is2017) and isData : metCollection = "slimmedMETsMuEGClean"

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
  leptonMvaWeightsMu            = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_BDTG.weights.xml"),
  leptonMvaWeightsEle           = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_BDTG.weights.xml"),
  leptonMvaWeightsMuttH         = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH_BDTG.weights.xml"),
  leptonMvaWeightsElettH        = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH_BDTG.weights.xml"),   
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
  met                           = cms.InputTag(metCollection),
  jets                          = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt"),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  triggerObjets                 = cms.InputTag("selectedPatTrigger"),
  SingleEleTriggers             = cms.untracked.vstring(),
  SingleMuoTriggers             = cms.untracked.vstring(),
  SingleEleTriggers2017         = cms.untracked.vstring(),
  SingleMuoTriggers2017         = cms.untracked.vstring(),
  recoResults                   = cms.InputTag("TriggerResults::RECO"),
  badPFMuonFilter               = cms.InputTag("BadPFMuonFilter"),
  badChargedCandFilter          = cms.InputTag("BadChargedCandidateFilter"),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  isSUSY                        = cms.untracked.bool(isSUSY)
)

## Single triggers for matching
process.blackJackAndHookers.SingleEleTriggers.extend([
        "HLT_Ele27_WPTight_Gsf_v*"
        ])
process.blackJackAndHookers.SingleMuoTriggers.extend([
        "HLT_IsoMu24_v*", 
        "HLT_IsoTkMu24_v*" 
        ])
process.blackJackAndHookers.SingleEleTriggers2017.extend([
        "HLT_Ele32_WPTight_Gsf_v*"
        ])
process.blackJackAndHookers.SingleMuoTriggers2017.extend([
        "HLT_IsoMu27_v*"
        ])


if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(filename = "../data/JSON/" + getJSON(is2017)).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.BadPFMuonFilter *
                     process.BadChargedCandidateFilter *
                     process.egmSequence *
                     process.jetSequence *
                     process.fullPatMetSequence *
                     process.blackJackAndHookers)
