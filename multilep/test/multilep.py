import sys, os
import FWCore.ParameterSet.Config as cms

# Default input file (could be overwritten by parameters given on the command line and by crab), some examples:
#inputFile      = 'file:///pnfs/iihe/cms/ph/sc4/store/data/Run2017F/DoubleMuon/MINIAOD/17Nov2017-v1/70000/E4B6F7A1-7BDE-E711-8C42-02163E019DE8.root'
#inputFile      = "root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAOD/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/02041699-0BFB-E711-AAD4-FA163E965751.root"
#inputFile      = '/store/data/Run2017F/SingleMuon/MINIAOD/17Nov2017-v1/00000/3E7C07F9-E6F1-E711-841A-0CC47A4C8E46.root'
#inputFile      = '/store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v3/000/316/569/00000/0085320B-4E64-E811-A2D3-FA163E2A55D6.root'
#inputFile      = '/store/data/Run2018A/MET/MINIAOD/PromptReco-v3/000/316/666/00000/0CC8EDCD-FD64-E811-BCA8-02163E01A020.root'
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root'
inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-5_V-0.00126095202129_mu_massiveAndCKM_LO/heavyNeutrino_45.root'
#inputFile       = '/store/data/Run2016C/SingleMuon/MINIAOD/17Jul2018-v1/20000/FEC97F81-0097-E811-A7B9-90E2BACC5EEC.root'
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-2_V-0.0141421356237_e_massiveAndCKM_LO/heavyNeutrino_108.root'
#inputFile       = 'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/20000/76EBF3B8-67EF-E811-95B0-0CC47AC52AFC.root'
#inputFile       = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/00A25ADE-DFD4-E611-8EAC-0025905A48B2.root"
 
#2018D dataset files for unskimmed synchronisation with Mohamed
#inputFile       = '/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/325/159/00000/0731BA23-725E-8540-BA91-D6B06D8D7826.root'#60723 events
#inputFile       = '/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/320/569/00000/3C8C28E7-1A96-E811-BA8D-02163E012DD8.root'#999 events
#inputFile       = '/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/320/570/00000/34C040A0-4096-E811-A1B7-FA163EEC3A41.root'#1426 events

# Other default arguments

nEvents         = 1000
extraContent    = ''
outputFile      = 'dilep.root'  # trilep    --> skim three leptons (basic pt/eta criteria)
                                # dilep     --> skim two leptons
                                # singlelep --> skim one lepton
                                # singlejet --> one jet
                                # FR        --> one jet and one light lepton
                                # noskim    --> no skim

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "outputFile"     in sys.argv[i]: outputFile   = getVal(sys.argv[i])
    elif "inputFile"    in sys.argv[i]: inputFile    = getVal(sys.argv[i])
    elif "extraContent" in sys.argv[i]: extraContent = getVal(sys.argv[i])
    elif "events"       in sys.argv[i]: nEvents      = int(getVal(sys.argv[i]))

isData = not ('SIM' in inputFile or 'heavyNeutrinoMiniAOD' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile or 'Fall17' in inputFile
is2018 = "Run2018" in inputFile or "18MiniAOD" in inputFile or 'Autumn18' in inputFile
isSUSY = "SMS-T" in inputFile
isFastSim = 'Fast' in inputFile

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

# Latest recommended global tags can always be checked here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if is2018 and 'PromptReco' in inputFile: process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v16'

#TO FIX RUN DEPENDENT JEC IN 2018!!! 102X_dataRun2_v12 (ABC)/ 102X_dataRun2_Prompt_v15 (D)
elif is2018:                             process.GlobalTag.globaltag = '102X_dataRun2_v14' if isData else '102X_upgrade2018_realistic_v21'
elif is2017:                             process.GlobalTag.globaltag = '102X_dataRun2_v14'  if isData else '102X_mc2017_realistic_v8'
else:                                    process.GlobalTag.globaltag = '102X_dataRun2_v14'  if isData else '102X_mcRun2_asymptotic_v8'

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
if isFastSim:
  if is2018:   jecUncertaintyFile = 'Autumn18_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
  elif is2017: jecUncertaintyFile = 'Fall17_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
  else:        jecUncertaintyFile = 'Summer16_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
else:
  if is2018:   jecUncertaintyFile = 'Autumn18_V19_MC_Uncertainty_AK4PFchs.txt'
  elif is2017: jecUncertaintyFile = 'Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt'
  else:        jecUncertaintyFile = 'Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt'


from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
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

yy = '16'
if is2017: yy = '17'
elif is2018: yy = '18'

yyy = '16'
if is2017 or is2018: yyy = '17'

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

process.load('heavyNeutrino.multilep.displacedInclusiveVertexing_cff')

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  offlineBeamSpot		= cms.InputTag("offlineBeamSpot"),
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  packedGenParticles            = cms.InputTag("packedGenParticles"),
  particleLevelLeptons          = cms.InputTag("particleLevel:leptons"),
  particleLevelPhotons          = cms.InputTag("particleLevel:photons"),
  particleLevelJets             = cms.InputTag("particleLevel:jets"),
  particleLevelMets             = cms.InputTag("particleLevel:mets"),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffAreas                 = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt'),
  muonsEffAreas_80X             = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
  electrons                     = cms.InputTag("slimmedElectrons"),
  electronsEffAreas             = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'),    # Recommended, used by standard IDs (the difference with the outdated effective areas is typically small)
  electronsEffAreas_Summer16    = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),  # old effective aras are used in 2016 computation of ttH MVA
  electronsEffAreas_Spring15    = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'), # prehistoric effective aras are used in 2016 computation of ttH MVA
  leptonMvaWeightsMuttH         = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH"+yyy+"_BDTG.weights.xml"),
  leptonMvaWeightsElettH        = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH"+yyy+"_BDTG.weights.xml"),
  leptonMvaWeightsEletZqTTV     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV"+yyy+"_BDTG.weights.xml"),
  leptonMvaWeightsMutZqTTV      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV"+yyy+"_BDTG.weights.xml"),
  leptonMvaWeightsEleTOP        = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_TOP"+yy+"_BDTG.weights.xml"),
  leptonMvaWeightsMuTOP         = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_TOP"+yy+"_BDTG.weights.xml"),
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
  secondaryVertices             = cms.InputTag("displacedInclusiveSecondaryVertices"),
  triggerObjects                = cms.InputTag("slimmedPatTrigger"),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  is2018                        = cms.untracked.bool(is2018),
  isFastSim                     = cms.untracked.bool(isFastSim),
  isSUSY                        = cms.untracked.bool(isSUSY),
  storeLheParticles             = cms.untracked.bool('storeLheParticles' in extraContent),
  storeGenParticles             = cms.untracked.bool('storeGenParticles' in extraContent),
  storeParticleLevel            = cms.untracked.bool('storeParticleLevel' in extraContent),
  storeAllTauID                 = cms.untracked.bool('storeAllTauID' in extraContent),
  headerPart1                   = cms.FileInPath("heavyNeutrino/multilep/data/header/soviet.txt"),
  headerPart2                   = cms.FileInPath("heavyNeutrino/multilep/data/header/text.txt")
)

#
# Additional dilep filter and trigger filter to be run before the IVF in order to avoid running out of walltime
#
process.dilepFilter = cms.EDFilter('multilepFilter',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  muons                         = cms.InputTag("slimmedMuons"),
  electrons                     = cms.InputTag("slimmedElectrons"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  triggerObjects                = cms.InputTag("slimmedPatTrigger"),
  is2017                        = cms.untracked.bool(is2017),
  is2018                        = cms.untracked.bool(is2018),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0])
)

#
# if using a filter before running blackJackAndHookers module, run this module first to calculate hCounter and other 'global sample' variables before filtering events
#
process.blackJackAndHookersGlobal = cms.EDAnalyzer('multilepGlobal',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  isData                        = cms.untracked.bool(isData),
)

def getJSON(is2017, is2018):
    if is2018:   return "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
    elif is2017: return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else:        return "Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt"

if isData:
  print('Sample is found to be 20%s data, will process using %s' % (yy, getJSON(is2017, is2018)))
  import FWCore.PythonUtilities.LumiList as LumiList
  jsonDir = os.path.expandvars('$CMSSW_BASE/src/heavyNeutrino/multilep/data/JSON')
  process.source.lumisToProcess = LumiList.LumiList(filename = os.path.join(jsonDir, getJSON(is2017, is2018))).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.blackJackAndHookersGlobal *
                     process.dilepFilter *
                     process.egammaPostRecoSeq *
                     process.jetSequence *
                     process.fullPatMetSequence *
                     process.prefiringweight *
                     process.particleLevelSequence *
                     process.rerunMvaIsolationSequence *
                     process.displacedInclusiveVertexing *
                     getattr(process,updatedTauName) *                    
                     process.blackJackAndHookers)
