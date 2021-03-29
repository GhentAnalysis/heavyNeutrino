import sys, os
import FWCore.ParameterSet.Config as cms

# Default input file (could be overwritten by parameters given on the command line and by crab), some examples:
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-105To160_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/2E242480-5C0D-E911-B9A6-90E2BACBAA90.root"
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017RECOPF_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0A1754A2-256F-E811-AD07-6CC2173CAAE0.root"
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root'
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIAutumn18MiniAOD/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v2/110000/00707922-8E6F-3042-A709-2DD4DB9AEDED.root'
#inputFile       = '/store/mc/RunIISummer19UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FBC4E43E-2C35-974A-84C9-29AD0430BD39.root'
#inputFile       = '/store/data/Run2017E/SingleElectron/MINIAOD/09Aug2019_UL2017-v1/130000/1428BFE4-BCE6-EA4D-9E49-CE5DC19AFFCE.root'
inputFile       = '/store/mc/RunIISummer19UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/70000/D01D206D-918C-4242-B9C3-07C3D8C106D3.root'
#inputFile       = '/store/data/Run2018B/DoubleMuon/MINIAOD/12Nov2019_UL2018-v2/70000/F5018BF9-184A-524E-A209-BC7A6BB7D8D4.root'
#inputFile       = '/store/mc/RunIISummer19UL18MiniAOD/DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/110000/5AA65A39-AEA4-6444-8B0E-1EC4846CCB3B.root'
#inputFile       = 'root://ccxrootdcms.in2p3.fr:1094/pnfs/in2p3.fr/data/cms/t2data/store/data/Run2018C/EGamma/MINIAOD/12Nov2019_UL2018-v2/270000/CD1A6E73-92DF-454F-BC1C-31199CB2E02D.root'
#inputFile        = '/store/mc/RunIIAutumn18MiniAOD/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/70000/F447BDAD-6642-BD46-B8E9-750F7F961BA7.root'
#inputFile       = '/store/mc/RunIIFall17MiniAODv2/SMS-TChiWZ_ZToLL_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/10000/00071D11-4C7A-E911-8E48-0CC47A1E0484.root'
#inputFile       = '/store/mc/RunIIAutumn18MiniAOD/SMS-TChiWZ_ZToLL_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall18Fast_102X_upgrade2018_realistic_v15-v1/50000/FBA243F4-DA64-3A40-8DD8-58A72064AD86.root'

# Other default arguments

nEvents         = 1000
extraContent    = ''
outputFile      = 'dilep.root' # trilep    --> skim three leptons (basic pt/eta criteria)
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

isData = not ('SIM' in inputFile or 'heavyNeutrinoMiniAOD' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile or 'Fall17' in inputFile
is2018 = "Run2018" in inputFile or "18MiniAOD" in inputFile or 'Autumn18' in inputFile
isSUSY = "SMS-T" in inputFile
isFastSim = 'Fast' in inputFile
isUL   = ("Summer19UL" in inputFile) or (isData and ("21Feb2020_UL2016" in inputFile or "09Aug2019_UL2017" in inputFile or "12Nov2019_UL2018" in inputFile))

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
# or if that doesn't give the necessary info, try here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
# or on conddb: https://cms-conddb.cern.ch/cmsDbBrowser/index/Prod
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if isUL:
    if is2017:                               process.GlobalTag.globaltag = '106X_dataRun2_v28' if isData else '106X_mc2017_realistic_v7'
    elif is2018:                             process.GlobalTag.globaltag = '106X_dataRun2_v28' if isData else '106X_upgrade2018_realistic_v11_L1v1'
else:
    #TO FIX RUN DEPENDENT JEC IN 2018!!! 102X_dataRun2_v12 (ABC)/ 102X_dataRun2_Prompt_v15 (D)
    if is2018 and 'PromptReco' in inputFile: process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v15'
    elif is2018:                             process.GlobalTag.globaltag = '102X_dataRun2_v12' if isData else '102X_upgrade2018_realistic_v20'
    elif is2017:                             process.GlobalTag.globaltag = '94X_dataRun2_v11'  if isData else '94X_mc2017_realistic_v17'
    else:                                    process.GlobalTag.globaltag = '94X_dataRun2_v10'  if isData else '94X_mcRun2_asymptotic_v3'

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
addJetSequence( process, inputFile, isData, is2017, is2018, isFastSim, isUL )
if isFastSim:
  if is2018:   jecUncertaintyFile = 'Autumn18_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
  elif is2017: jecUncertaintyFile = 'Fall17_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
  else:        jecUncertaintyFile = 'Summer16_FastSimV1_MC_Uncertainty_AK4PFchs.txt'
else:
  if is2018:   jecUncertaintyFile = 'Autumn18_V19_MC_Uncertainty_AK4PFchs.txt'
  elif is2017: jecUncertaintyFile = 'Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt'
  else:        jecUncertaintyFile = 'Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt'
#JEC uncertainty for Puppi:
jecUncertaintyFilePuppi = jecUncertaintyFile.replace('PFchs', 'PFPuppi')

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

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("offlineSlimmedPrimaryVertices"), #MET XY corr: make sure to use offlineSlimmedPrimaryVertices, not goodOfflinePrimaryVertices, since those are only a subset of vertices
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
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
  metPuppi                      = cms.InputTag("slimmedMETsPuppi"),
  jets                          = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jetsPuppi                     = cms.InputTag("slimmedJetsPuppi"),
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintyFile),
  jecUncertaintyFilePuppi       = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintyFilePuppi),
  JECtxtPath                    = cms.FileInPath("heavyNeutrino/multilep/data/JEC/dummy.txt"),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  recoResultsPrimary            = cms.InputTag("TriggerResults::PAT"),
  recoResultsSecondary          = cms.InputTag("TriggerResults::RECO"),
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

def getJSON(isUL, is2017, is2018):
    if is2018:
        if isUL: return "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        else:    return "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
    elif is2017:
        if isUL: return "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        else:    return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else:
        if isUL: return ""#not available yet
        else:    return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

if isData:
  print('Sample is found to be 20%s data, will process using %s' % (yy, getJSON(isUL, is2017, is2018)))
  import FWCore.PythonUtilities.LumiList as LumiList
  jsonDir = os.path.expandvars('$CMSSW_BASE/src/heavyNeutrino/multilep/data/JSON')
  process.source.lumisToProcess = LumiList.LumiList(filename = os.path.join(jsonDir, getJSON(isUL, is2017, is2018))).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.egammaPostRecoSeq *
                     process.jetSequence *
                     process.prefiringweight *
                     process.particleLevelSequence *
                     process.rerunMvaIsolationSequence *
                     getattr(process,updatedTauName) *                    
                     process.blackJackAndHookers)
