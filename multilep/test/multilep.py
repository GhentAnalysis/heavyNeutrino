import sys, os
import FWCore.ParameterSet.Config as cms

#function to return JSON file
def getJSON(is2017):
    if is2017: return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else: return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

# Default arguments
nEvents         = -1
## TTbar
inputFile       = "root://xrootd-cms.infn.it///store/mc/RunIISummer16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/70000/A03D5CD4-04C8-E611-8199-02163E019B3D.root"
## DY
#inputFile       = "root://xrootd-cms.infn.it///store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/CA5B0759-7EE5-E611-A6AB-0025905A6126.root"
#inputFile       = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_12.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_13.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_14.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_18.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_2.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_29.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_31.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_35.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_38.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_4.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_40.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_41.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_43.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_48.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_49.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_50.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_55.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_57.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_58.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_60.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_61.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_68.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_69.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_73.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_78.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_8.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_86.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_87.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_88.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_9.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_92.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_94.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_95.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_97.root,file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_trilepton_M-5_V-0.00707106781187_e_massiveAndCKM_LO/heavyNeutrino_98.root"

###                            ###
### Skim from output file name ###
###                            ###
# trilep      --> skim three leptons (basic pt/eta criteria)
# displtrilep --> skim three leptons (basic pt/eta criteria)
# dilep       --> skim two leptons
# singlelep   --> skim one lepton
# ttg         --> skim two leptons + one photon
# fakerate    --> not implemented
outputFile      = 'trilep_00.root'

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

process.source       = cms.Source("PoolSource", 
                                  fileNames = cms.untracked.vstring(inputFile.split(",")), 
                                  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                                  )
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if   isData and is2017: process.GlobalTag.globaltag = '94X_dataRun2_v6'   
elif is2017:            process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'
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
addJetSequence(process, isData, is2017)
addElectronAndPhotonSequence(process)

#
# Read additional MET filters not stored in miniAOD
#
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
for module in [process.BadPFMuonFilter, process.BadChargedCandidateFilter]:
  module.muons        = cms.InputTag("slimmedMuons")
  module.PFCandidates = cms.InputTag("packedPFCandidates")

process.BadPFMuonFilter.filter = cms.bool(False)
process.BadChargedCandidateFilter.filter = cms.bool(False)

#clean 2016 data met from spurious muons and ECAL slew rate
metCollection = "slimmedMETs"
#if (not is2017) and isData : metCollection = "slimmedMETsMuEGClean"

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  allowMatchingToAllIds         = cms.bool(False),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffectiveAreas           = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
  muonsEffectiveAreasFall17     = cms.FileInPath('heavyNeutrino/multilep/data/effAreas_cone03_Muons_Fall17.txt'),
  electrons                     = cms.InputTag("slimmedElectrons"),
  # WARNING this is spring 15, following SUSY-standard, i.e. not the most up-to-date values
  electronsEffectiveAreas       = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'),
  electronsEffectiveAreasFall17 = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt'),

  electronsMva                  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
  electronsMvaHZZ               = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
  # electronMvaFall17Iso          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
  # electronMvaFall17NoIso        = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
  electronsCutBasedVeto         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  electronsCutBasedLoose        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  electronsCutBasedMedium       = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  electronsCutBasedTight        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  # leptonMvaWeightsMuttH         = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH_BDTG.weights.xml"),
  # leptonMvaWeightsElettH        = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH_BDTG.weights.xml"),   
  # leptonMvaWeightsMutZqTTV      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_leptonMva_tZqTTV.weights.xml"),
  # leptonMvaWeightsEletZqTTV     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/ele_leptonMva_tZqTTV.weights.xml"),
  # leptonMvaWeightsMuSUSY16      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_SUSY16_BDTG.weights.xml"),
  # leptonMvaWeightsEleSUSY16     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_SUSY16_BDTG.weights.xml"),
  # leptonMvaWeightsMuttH16       = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH16_BDTG.weights.xml"),
  # leptonMvaWeightsElettH16      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH16_BDTG.weights.xml"),   
  # leptonMvaWeightsMuSUSY17      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_SUSY17_BDTG.weights.xml"),
  # leptonMvaWeightsEleSUSY17     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_SUSY17_BDTG.weights.xml"),
  # leptonMvaWeightsMuttH17       = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH17_BDTG.weights.xml"),
  # leptonMvaWeightsElettH17      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH17_BDTG.weights.xml"),   
  # leptonMvaWeightsEletZqTTV16   = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV16_BDTG.weights.xml"),
  # leptonMvaWeightsMutZqTTV16    = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV16_BDTG.weights.xml"),
  # leptonMvaWeightsEletZqTTV17   = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV17_BDTG.weights.xml"),
  # leptonMvaWeightsMutZqTTV17    = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV17_BDTG.weights.xml"),
  JECtxtPath                    = cms.FileInPath("heavyNeutrino/multilep/data/JEC/dummy.txt"), 
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
  #jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  #jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  #jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt"),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  triggerObjects                = cms.InputTag("selectedPatTrigger"),
  #triggerObjects                = cms.InputTag("slimmedPatTrigger"),
  SingleEleTriggers             = cms.vstring(),
  SingleMuoTriggers             = cms.vstring(),
  SingleEleTriggers2017         = cms.vstring(),
  SingleMuoTriggers2017         = cms.vstring(),
  #recoResultsPrimary            = cms.InputTag("TriggerResults::PAT"),
  #recoResultsSecondary          = cms.InputTag("TriggerResults::RECO"),
  badPFMuonFilter               = cms.InputTag("BadPFMuonFilter"),
  badChargedCandFilter          = cms.InputTag("BadChargedCandidateFilter"),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  isSUSY                        = cms.untracked.bool(isSUSY)
)

## Change trigger-object label for signal
#if 'HeavyNeutrino' in inputFile: process.blackJackAndHookers.triggerObjects = cms.InputTag("slimmedPatTrigger")

## Single triggers for matching
if 'FR' in outputFile:
    process.blackJackAndHookers.SingleEleTriggers.extend([
            "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*",
            "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*"
            ])
    process.blackJackAndHookers.SingleMuoTriggers.extend([
            "HLT_Mu3_PFJet40_v*",
            "HLT_Mu8_v*",
            "HLT_Mu17_v*"
            ])
    process.blackJackAndHookers.SingleEleTriggers2017.extend([
            "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*",
            "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*"
            ])
    process.blackJackAndHookers.SingleMuoTriggers2017.extend([
            "HLT_Mu3_PFJet40_v*",
            "HLT_Mu8_v*",
            "HLT_Mu17_v*"
            ])
else:
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
