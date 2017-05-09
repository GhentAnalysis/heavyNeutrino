import sys
import FWCore.ParameterSet.Config as cms

isData          = False
treeForFakeRate = False # To be implemented
isDiLep         = False # To be implemented (for 1l+2l trigger efficiencies in MET, 3l trigger efficiencies could be done with trilep samples in DoubleEG and DoubleMu)
inputFile       = '/store/mc/RunIISummer16MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/50D92A94-D1D0-E611-BEA6-D4AE526A023A.root'
nEvents         = -1
outputFile      = None

def getVal(arg):
    return arg.split('=')[-1]

## loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if   "isData"          in sys.argv[i]: isData          = (getVal(sys.argv[i]) == "True")
    elif "treeForFakeRate" in sys.argv[i]: treeForFakeRate = (getVal(sys.argv[i]) == "True")
    elif "isDiLep"         in sys.argv[i]: isDiLep         = (getVal(sys.argv[i]) == "True")
    elif "output"          in sys.argv[i]: outputFile      = getVal(sys.argv[i])
    elif "input"           in sys.argv[i]: inputFile       = getVal(sys.argv[i])
    elif "events"          in sys.argv[i]: nEvents         = int(getVal(sys.argv[i]))


process = cms.Process("BlackJackAndHookers")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO' # Options: INFO, WARNING, ERROR
process.load("FWCore.MessageService.MessageLogger_cfi")

#allow unscheduled mode  # aargh I hate this unscheduled mode, it is evil! Why do we need it?
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.source    = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nEvents))

#define globaltag for JEC
process.load('Configuration.StandardSequences.Services_cff') # do we need this?
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '???' if isData else '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

#load JEC
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')



from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
)


if not outputFile:
  if treeForFakeRate: outputFile = 'fakeRate.root'
  elif singleLep:     outputFile = 'singleLep.root'
  else:               outputFile = 'trilepton.root'
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

#Read additional MET filters not stored in miniAOD
#Bad muon filter
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#Bad chargedCandidate filter
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.blackJackAndHookers = cms.EDAnalyzer('multilep',
	vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
	muons = cms.InputTag("slimmedMuons"),
	electrons = cms.InputTag("slimmedElectrons"),
	taus = cms.InputTag("slimmedTaus"),
	packedCandidates = cms.InputTag("packedPFCandidates"),
	rhoCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
	rhoAll = cms.InputTag("fixedGridRhoFastjetAll"),
	met = cms.InputTag("slimmedMETs"),
	#jets = cms.InputTag("slimmedJets")
	jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
	triggers = cms.InputTag("TriggerResults","","HLT"),
	recoResults = cms.InputTag("TriggerResults", "", "RECO"),
	badPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
	badChargedCandFilter = cms.InputTag("BadChargedCandidateFilter")
)


process.p = cms.Path(process.BadPFMuonFilter * 
					process.BadChargedCandidateFilter *
					process.ak4PFCHSL1FastL2L3CorrectorChain *
					process.blackJackAndHookers)

