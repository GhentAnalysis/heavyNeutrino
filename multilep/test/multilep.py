import sys, os
import FWCore.ParameterSet.Config as cms

# Default input file (could be overwritten by parameters given on the command line and by crab), some examples:
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-105To160_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/2E242480-5C0D-E911-B9A6-90E2BACBAA90.root"
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017RECOPF_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0A1754A2-256F-E811-AD07-6CC2173CAAE0.root"
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root'
#inputFile        = '/store/mc/RunIIAutumn18MiniAOD/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/70000/F447BDAD-6642-BD46-B8E9-750F7F961BA7.root'
inputFile       = '/store/mc/RunIISummer16MiniAODv3/SMS-TChiWZ_ZToLL_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/100000/502F9078-3296-E911-BFB6-0025905B85EC.root'
#inputFile       = '/store/mc/RunIIFall17MiniAODv2/SMS-TChiWZ_ZToLL_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/10000/00071D11-4C7A-E911-8E48-0CC47A1E0484.root'
#inputFile       = '/store/mc/RunIIAutumn18MiniAOD/SMS-TChiWZ_ZToLL_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall18Fast_102X_upgrade2018_realistic_v15-v1/50000/FBA243F4-DA64-3A40-8DD8-58A72064AD86.root'
# Other default arguments

nEvents         =-1
extraContent    = 'storeAllTauID'
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

isData = not ('SIM' in inputFile or '/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile or 'Fall17' in inputFile
is2018 = "Run2018" in inputFile or "18MiniAOD" in inputFile or 'Autumn18' in inputFile
isSUSY = "SMS-T" in inputFile
isFastSim = ( 'PUSummer16v3Fast' in inputFile ) or ( 'PUFall17Fast' in inputFile ) or ( 'PUFall18Fast' in inputFile )

process = cms.Process("BlackJackAndHookers")

# Print a warning if a different release is used as the one in the setup script (i.e. probably something will be broken)
os.system('if [[ "$(grep RELEASE= $CMSSW_BASE/src/heavyNeutrino/setup.sh)" != *"$CMSSW_VERSION"* ]];then echo ">>> WARNING: you are using a different release as the one specified in the setup.sh script! <<< ";fi')

# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")), duplicateCheckMode = cms.untracked.string('noDuplicateCheck'))
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

# Latest recommended global tags can always be checked here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#TO FIX RUN DEPENDENT JEC IN 2018!!! 102X_dataRun2_v12 (ABC)/ 102X_dataRun2_Prompt_v15 (D)
if is2018 and 'PromptReco' in inputFile: process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v15'
elif is2018:                             process.GlobalTag.globaltag = '102X_dataRun2_v12' if isData else '102X_upgrade2018_realistic_v20'
elif is2017:                             process.GlobalTag.globaltag = '94X_dataRun2_v11'  if isData else '94X_mc2017_realistic_v17'
else:                                    process.GlobalTag.globaltag = '94X_dataRun2_v10'  if isData else '94X_mcRun2_asymptotic_v3'

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
addJetSequence( process, isData, is2017, is2018, isFastSim )
if is2018 and ( not isFastSim ):
  jecUncertaintyFile = 'Autumn18_V19_MC_Uncertainty_AK4PFchs.txt'
elif is2018 and isFastSim:
  jecUncertaintyFile = 'Autumn18_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
elif is2017 and ( not isFastSim ):
  jecUncertaintyFile = 'Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt'
elif is2017 and isFastSim:
 jecUncertaintyFile = 'Fall17_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
elif not isFastSim:
  jecUncertaintyFile = 'Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt'
else:
  jecUncertaintyFile = 'Summer16_FastSimV1_MC_Uncertainty_AK4PFchs.txt'


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

yy = '17' if is2017 or is2018 else '16'

#
#Latest tau ID
#
updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = False,
                    updatedTauName = updatedTauName,
                    toKeep = [ "2017v2", "newDM2017v2", #classic MVAIso tau-Ids
                                "deepTau2017v2p1" #latest deepTau2017v2p1
                               ])
tauIdEmbedder.runTauID()

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  allowMatchingToAllIds         = cms.bool(False),
  particleLevelLeptons          = cms.InputTag("particleLevel:leptons"),
  particleLevelPhotons          = cms.InputTag("particleLevel:photons"),
  particleLevelJets             = cms.InputTag("particleLevel:jets"),
  particleLevelMets             = cms.InputTag("particleLevel:mets"),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffectiveAreas           = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
  muonsEffectiveAreasFall17     = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt'),
  electrons                     = cms.InputTag("slimmedElectrons"),
  # FIXME: currently the next two lines are used for all electron values, not only the ttH leptonMva!
  electronsEffectiveAreasMiniIso2016  = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'),
  electronsEffectiveAreasRelIso2016  = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'), #TODO: fix this is the master branch, for some reason (probably because of ttH syncing) the old effective areas slipped back into here; it seems we need a separate ttH
  electronsEffectiveAreasFall17 = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'), # Recommended, used by standard IDs (the difference with the outdated effective areas is typically small)
  leptonMvaWeightsMuttH         = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH"+yy+"_BDTG.weights.xml"),
  leptonMvaWeightsElettH        = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH"+yy+"_BDTG.weights.xml"),
  leptonMvaWeightsEletZqTTV     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV"+yy+"_BDTG.weights.xml"),
  leptonMvaWeightsMutZqTTV      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV"+yy+"_BDTG.weights.xml"),
  photons                       = cms.InputTag("slimmedPhotons"),
  photonsChargedEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt'),
  photonsNeutralEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt'),
  photonsPhotonsEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt'),
#  taus                          = cms.InputTag("slimmedTaus"),
  taus                          = cms.InputTag("slimmedTausNewID"),
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
  SingleEleTriggers             = cms.vstring(),                       # displaced specific
  SingleMuoTriggers             = cms.vstring(),                       # displaced specific
  SingleEleTriggers2017         = cms.vstring(),                       # displaced specific
  SingleMuoTriggers2017         = cms.vstring(),                       # displaced specific
  SingleEleTriggers2018         = cms.vstring(),                       # displaced specific
  SingleMuoTriggers2018         = cms.vstring(),                       # displaced specific
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  is2018                        = cms.untracked.bool(is2018),
  isFastSim                     = cms.untracked.bool(isFastSim),
  isSUSY                        = cms.untracked.bool(isSUSY),
  storeLheParticles             = cms.untracked.bool('storeLheParticles' in extraContent),
  storeParticleLevel            = cms.untracked.bool('storeParticleLevel' in extraContent),
  storeAllTauID                 = cms.untracked.bool('storeAllTauID' in extraContent),
)

def getJSON(is2017, is2018):
    if is2018:   return "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
    elif is2017: return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else:        return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"


## Single triggers for matching


process.blackJackAndHookers.SingleEleTriggers.extend(["HLT_Ele27_WPTight_Gsf_v*"])
process.blackJackAndHookers.SingleMuoTriggers.extend(["HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*"])
process.blackJackAndHookers.SingleEleTriggers2017.extend(["HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*","HLT_Ele35_WPTight_Gsf_v*","HLT_Ele32_WPTight_Gsf"])
process.blackJackAndHookers.SingleMuoTriggers2017.extend(["HLT_IsoMu24_*","HLT_IsoMu27_v*"])
process.blackJackAndHookers.SingleEleTriggers2018.extend(["HLT_Ele32_WPTight_Gsf_v*"])
process.blackJackAndHookers.SingleMuoTriggers2018.extend(["HLT_IsoMu24_*","HLT_IsoMu27_v*"])

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
                     process.rerunMvaIsolationSequence *
                     getattr(process,updatedTauName) *                    
                     process.blackJackAndHookers)
