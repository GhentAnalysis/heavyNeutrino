import sys
import FWCore.ParameterSet.Config as cms

# Default arguments
inputFile       = '/store/mc/RunIISummer16MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/50D92A94-D1D0-E611-BEA6-D4AE526A023A.root'
isData          = False
nEvents         = -1
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

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO' # Options: INFO, WARNING, ERROR

#allow unscheduled mode  # aargh I hate this unscheduled mode, it is evil! Why do we need it?
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.source    = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nEvents))

#define globaltag for JEC
#process.load('Configuration.StandardSequences.Services_cff') # do we need this? I don't think so
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '' if isData else '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

#load JEC
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')


# Not working yet
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
)

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))


#Read additional MET filters not stored in miniAOD
#Bad muon filter
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
for module in [process.BadPFMuonFilter, process.BadChargedCandidateFilter]:
  module.muons        = cms.InputTag("slimmedMuons")
  module.PFCandidates = cms.InputTag("packedPFCandidates")


# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
# fakeRateTree         = cms.untracked.bool(outputFile.count('fakeRate')), # TO BE IMPLEMENTED
# dileptonTree         = cms.untracked.bool(outputFile.count('dilepton')), # TO BE IMPLEMENTED
  vertices             = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons                = cms.InputTag("slimmedMuons"),
  electrons            = cms.InputTag("slimmedElectrons"),
  taus                 = cms.InputTag("slimmedTaus"),
  packedCandidates     = cms.InputTag("packedPFCandidates"),
  rhoCentralNeutral    = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
  rhoAll               = cms.InputTag("fixedGridRhoFastjetAll"),
  met                  = cms.InputTag("slimmedMETs"),
  jets                 = cms.InputTag("slimmedJets"),
 #jets                 = cms.InputTag("updatedPatJetsUpdatedJEC"),
  triggers             = cms.InputTag("TriggerResults","","HLT"),
  recoResults          = cms.InputTag("TriggerResults", "", "RECO"),
  badPFMuonFilter      = cms.InputTag("BadPFMuonFilter"),
  badChargedCandFilter = cms.InputTag("BadChargedCandidateFilter")
)


if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(filename = 'TO_BE_ADDED').getVLuminosityBlockRange()


process.p = cms.Path(process.BadPFMuonFilter * 
                     process.BadChargedCandidateFilter *
                     process.ak4PFCHSL1FastL2L3CorrectorChain *
                     process.blackJackAndHookers)

