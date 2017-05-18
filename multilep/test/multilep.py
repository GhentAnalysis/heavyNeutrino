import sys
import FWCore.ParameterSet.Config as cms

# Default arguments
inputFile       = '/store/mc/RunIISummer16MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/50D92A94-D1D0-E611-BEA6-D4AE526A023A.root'
isData          = False
nEvents         = 1000
outputFile      = 'trilepton.root'

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "isData"            in sys.argv[i]: isData          = getVal(sys.argv[i])
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

#
# Vertex collection
#
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')
process.goodOfflinePrimaryVertices.src    = cms.InputTag('offlineSlimmedPrimaryVertices')
process.goodOfflinePrimaryVertices.filter = cms.bool(True) # This was false in the old Majorana code, why?


#
# Latest JEC through globaltag, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
# Has currently no effect (Moriond2017 miniAOD contains latest JEC)
#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' if isData else '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
if isData:
  jetCorrectorLevels          = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
  process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
else:
  jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(jetCorrectorLevels), 'None')
)
process.jetSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.ak4PFCHSL1FastL2L3CorrectorChain)



#
# JER, see https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
#
for (i, j) in [(0, ''), (-1, 'Down'), (1, 'Up')]:
  jetSmearing = cms.EDProducer('SmearedPATJetProducer',
        src          = cms.InputTag('updatedPatJetsUpdatedJEC'),
        enabled      = cms.bool(True),
        rho          = cms.InputTag("fixedGridRhoFastjetAll"),
        algo         = cms.string('AK4PFchs'),
        algopt       = cms.string('AK4PFchs_pt'),
        genJets      = cms.InputTag('slimmedGenJets'),
        dRMax        = cms.double(0.2),
        dPtMaxFactor = cms.double(3),
        debug        = cms.untracked.bool(False),
        variation    = cms.int32(i),
  )
  setattr(process, 'jetSmearing'+j, jetSmearing)
  process.jetSequence *= jetSmearing


#
# Set up electron and photon identifications
#
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
switchOnVIDPhotonIdProducer(  process, DataFormat.MiniAOD)
electronModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                   'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                   'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
photonModules   = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff',
                   'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
for idmod in electronModules: setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in photonModules:   setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)



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
# fakeRateTree                  = cms.untracked.bool(outputFile.count('fakeRate')), # TO BE IMPLEMENTED
# dileptonTree                  = cms.untracked.bool(outputFile.count('dilepton')), # TO BE IMPLEMENTED
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffectiveAreas           = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
  electrons                     = cms.InputTag("slimmedElectrons"),
  electronsEffectiveAreas       = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'), # WARNING this is spring 15, following SUSY-standard, i.e. not the most up-to-date values
  electronsMva                  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
  electronsMvaHZZ               = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
  electronsCutBasedLoose        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  electronsCutBasedMedium       = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  electronsCutBasedTight        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  photons                       = cms.InputTag("slimmedPhotons"),
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
  jets                          = cms.InputTag("updatedPatJetsUpdatedJEC"),
  jetsSmeared                   = cms.InputTag("jetSmearing"),
  jetsSmearedUp                 = cms.InputTag("jetSmearingUp"),
  jetsSmearedDown               = cms.InputTag("jetSmearingDown"),
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt"),
  triggers                      = cms.InputTag("TriggerResults","","HLT"),
  recoResults                   = cms.InputTag("TriggerResults", "", "RECO"),
  badPFMuonFilter               = cms.InputTag("BadPFMuonFilter"),
  badChargedCandFilter          = cms.InputTag("BadChargedCandidateFilter"),
)


if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(filename = 'TO_BE_ADDED').getVLuminosityBlockRange()


process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.BadPFMuonFilter *
                     process.BadChargedCandidateFilter *
                     process.egmGsfElectronIDSequence *
                     process.egmPhotonIDSequence *
                     process.jetSequence *
                     process.blackJackAndHookers)

