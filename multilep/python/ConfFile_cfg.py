import FWCore.ParameterSet.Config as cms #load all CMS specific python modules

process = cms.Process("BlackJackAndHookers")

process.load("FWCore.MessageService.MessageLogger_cfi")
#allow unscheduled mode
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
#change number of events back to -1 (100 used for testing the code)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#make an output file
process.TFileService = cms.Service("TFileService", fileName = cms.string("../test/output.root") )
#define globaltag for JEC
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

#load JEC
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')



from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
)


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring('dcap://maite.iihe.ac.be:/pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A4140169-0FC2-E611-9809-002590FD030A.root')
	fileNames = cms.untracked.vstring('dcap://maite.iihe.ac.be:/pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer16MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/50D92A94-D1D0-E611-BEA6-D4AE526A023A.root')
)

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

