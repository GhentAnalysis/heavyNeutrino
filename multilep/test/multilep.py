import sys, os
import FWCore.ParameterSet.Config as cms

# Default input file (could be overwritten by parameters given on the command line and by crab), some examples:
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-105To160_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/2E242480-5C0D-E911-B9A6-90E2BACBAA90.root"
#inputFile      = "file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017RECOPF_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0A1754A2-256F-E811-AD07-6CC2173CAAE0.root"
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrino/testFiles/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/100000/42EFAC9D-DC91-DB47-B931-B6B816C60C21.root'
# inputFile        = '/store/mc/RunIISummer20UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260000/0071F930-6376-7A48-89F1-74E189BD3BFC.root'
# inputFile        = '/store/mc/RunIISummer20UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v1/280000/00085E9D-5FBB-B840-82A6-D5CE33C0F4C4.root'
inputFile        = '/store/mc/RunIISummer20UL16MiniAODAPV/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v8-v2/00000/01181D3D-D52C-4E45-9A25-2F5B5CC381D1.root'
#inputFile = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/prompt/ttGamma_Dilept_5f_ckm_LO_1line/heavyNeutrino_429.root'
#inputFile = '/store/mc/RunIISummer16MiniAODv3/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/50000/FEB45979-7B72-E911-9CD8-0242AC1C0505.root'
#inputFile       = '/store/mc/RunIISummer16MiniAODv3/SMS-TChiWZ_ZToLL_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/100000/502F9078-3296-E911-BFB6-0025905B85EC.root'
#inputFile       = '/store/mc/RunIIFall17MiniAODv2/SMS-TChiWZ_ZToLL_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/10000/00071D11-4C7A-E911-8E48-0CC47A1E0484.root'
#inputFile       = '/store/mc/RunIIAutumn18MiniAOD/SMS-TChiWZ_ZToLL_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall18Fast_102X_upgrade2018_realistic_v15-v1/50000/FBA243F4-DA64-3A40-8DD8-58A72064AD86.root'

# Other default arguments

nEvents         = 1000
extraContent    = 'storeAllTauID'
# extraContent    = 'storeJecSources'
# extraContent    = ''
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

isData = not ('SIM' in inputFile or 'heavyNeutrinoMiniAOD' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile or 'Fall17' in inputFile
is2018 = "Run2018" in inputFile or "18MiniAOD" in inputFile or 'Autumn18' in inputFile
is2016preVFP = "preVFP" in inputFile or "_HIPM" in inputFile
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
if is2018:                             process.GlobalTag.globaltag = '106X_dataRun2_v32' if isData else '106X_upgrade2018_realistic_v15_L1v1'
elif is2017:                             process.GlobalTag.globaltag = '106X_dataRun2_v32' if isData else '106X_mc2017_realistic_v8'
elif is2016preVFP:                       process.GlobalTag.globaltag = '106X_dataRun2_v32' if isData else '106X_mcRun2_asymptotic_preVFP_v9'
else:                                    process.GlobalTag.globaltag = '106X_dataRun2_v32' if isData else '106X_mcRun2_asymptotic_v15'

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
addJetSequence( process, inputFile, isData, is2017, is2018, is2016preVFP, isFastSim)
unc_prefix = ''
if isFastSim:
  if is2018:   unc_prefix = 'Autumn18_FastSimV1_MC'
  elif is2017: unc_prefix = 'Fall17_FastSimV1_MC'
  else:        unc_prefix = 'Summer16_FastSimV1_MC'
else:
  if is2018:
    if isData:
      if 'Run2018A' in inputFile: unc_prefix = 'Summer19UL18_RunA_V5_DATA'
      if 'Run2018B' in inputFile: unc_prefix = 'Summer19UL18_RunB_V5_DATA'
      if 'Run2018C' in inputFile: unc_prefix = 'Summer19UL18_RunC_V5_DATA'
      if 'Run2018D' in inputFile: unc_prefix = 'Summer19UL18_RunD_V5_DATA'
    else:
      unc_prefix = 'Summer19UL18_V5_MC'
  elif is2017: 
    if isData:
      if 'Run2017B' in inputFile: unc_prefix = 'Summer19UL17_RunB_V5_DATA'
      if 'Run2017C' in inputFile: unc_prefix = 'Summer19UL17_RunC_V5_DATA'
      if 'Run2017D' in inputFile: unc_prefix = 'Summer19UL17_RunD_V5_DATA'
      if 'Run2017E' in inputFile: unc_prefix = 'Summer19UL17_RunE_V5_DATA'
      if 'Run2017F' in inputFile: unc_prefix = 'Summer19UL17_RunF_V5_DATA'
    else: 
      unc_prefix = 'Summer19UL17_V5_MC'
  elif is2016preVFP: 
    if isData:
      if 'Run2017B' in inputFile or 'Run2017C' in inputFile or 'Run2017D' in inputFile: unc_prefix = 'Summer19UL16APV_RunBCD_V7_DATA'
      if 'Run2017E' in inputFile or 'Run2017E' in inputFile: unc_prefix = 'Summer19UL16APV_RunEF_V7_DATA'
    else:
      unc_prefix = 'Summer19UL16APV_V7_MC'
  else:
    if isData:
      unc_prefix = 'Summer19UL16_RunFGH_V7_DATA'
    else:
      unc_prefix = 'Summer19UL16_V7_MC'

  jecUncertaintyFile = '{0}/{0}_Uncertainty_AK4PFchs.txt'.format(unc_prefix)
  jecUncertaintySourcesFile =  '{0}/{0}_UncertaintySources_AK4PFchs.txt'.format(unc_prefix)
  jecUncertaintyRegroupedFile = '{0}/RegroupedV2_{0}_UncertaintySources_AK4PFchs.txt'.format(unc_prefix)
  jecL1FastJetFile =  "{0}/{0}_L1FastJet_AK4PFchs.txt".format(unc_prefix)
  jecL2RelativeFile =  "{0}/{0}_L2Relative_AK4PFchs.txt".format(unc_prefix)
  jecL3AbsoluteFile =  "{0}/{0}_L3Absolute_AK4PFchs.txt".format(unc_prefix)
  jecL2L3ResidualFile =  "{0}/{0}_L2L3Residual_AK4PFchs.txt".format(unc_prefix)

  jecUncertaintySourcesFilePuppi = jecUncertaintySourcesFile.replace('PFchs', 'PFPuppi')
  jecUncertaintyRegroupedFilePuppi = jecUncertaintyRegroupedFile.replace('PFchs', 'PFPuppi')
  jecL1FastJetFilePuppi = jecL1FastJetFile.replace('PFchs', 'PFPuppi')
  jecL2RelativeFilePuppi = jecL2RelativeFile.replace('PFchs', 'PFPuppi')
  jecL3AbsoluteFilePuppi = jecL3AbsoluteFile.replace('PFchs', 'PFPuppi')
  jecL2L3ResidualFilePuppi = jecL2L3ResidualFile.replace('PFchs', 'PFPuppi')


#JEC uncertainty for Puppi:
jecUncertaintyFilePuppi = jecUncertaintyFile.replace('PFchs', 'PFPuppi')

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

if is2018: setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2018-UL')  
if is2017: setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2017-UL')  
if is2016preVFP: setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2016preVFP-UL')  
else: setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2016postVFP-UL')  

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
  process.load("heavyNeutrino.multilep.particleLevelTTG_cfi")
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
                    toKeep = [ "2017v2", #classic MVAIso tau-Ids
                                "deepTau2017v2p1" #latest deepTau2017v2p1
                               ])
tauIdEmbedder.runTauID()

#rochester correction file to use
if is2017:
    rochesterCorrectionFile = 'UL/RoccoR2017UL.txt'
elif is2018:
    rochesterCorrectionFile = 'UL/RoccoR2018UL.txt'
elif is2016preVFP:
    rochesterCorrectionFile = 'UL/RoccoR2016aUL.txt'
else:
    rochesterCorrectionFile = 'UL/RoccoR2016bUL.txt'


# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
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
  taus                          = cms.InputTag("slimmedTausNewID"),
  packedCandidates              = cms.InputTag("packedPFCandidates"),
  rho                           = cms.InputTag("fixedGridRhoFastjetAll"),
  met                           = cms.InputTag("slimmedMETs"),
  metPuppi                      = cms.InputTag("slimmedMETsPuppi"),
  jets                          = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  jetsPuppi                     = cms.InputTag("slimmedJetsPuppi"),
  jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jecUncertaintyFile            = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintyFile),
  jecUncertaintyFilePuppi       = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintyFilePuppi),
  jecUncertaintySourcesFile     = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintySourcesFile),
  jecUncertaintyRegroupedFile   = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecUncertaintyRegroupedFile),
  jecL1FastJetFile              = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecL1FastJetFile),
  jecL2RelativeFile             = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecL2RelativeFile),
  jecL3AbsoluteFile             = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecL3AbsoluteFile),
  jecL2L3ResidualFile           = cms.FileInPath("heavyNeutrino/multilep/data/JEC/" + jecL2L3ResidualFile),
  rochesterCorrectionFile       = cms.FileInPath("heavyNeutrino/multilep/data/RochesterCorrections/" + rochesterCorrectionFile ),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  recoResultsPrimary            = cms.InputTag("TriggerResults::PAT"),
  recoResultsSecondary          = cms.InputTag("TriggerResults::RECO"),
  triggerObjects                = cms.InputTag("slimmedPatTrigger"),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  is2018                        = cms.untracked.bool(is2018),
  is2016preVFP                  = cms.untracked.bool(is2016preVFP),
  isFastSim                     = cms.untracked.bool(isFastSim),
  isSUSY                        = cms.untracked.bool(isSUSY),
  storeLheParticles             = cms.untracked.bool('storeLheParticles' in extraContent),
  storeGenParticles             = cms.untracked.bool('storeGenParticles' in extraContent),
  storeParticleLevel            = cms.untracked.bool('storeParticleLevel' in extraContent),
  storeJecSources               = cms.untracked.bool('storeJecSources' in extraContent),
  storeAllTauID                 = cms.untracked.bool('storeAllTauID' in extraContent),
  headerPart1                   = cms.FileInPath("heavyNeutrino/multilep/data/header/soviet.txt"),
  headerPart2                   = cms.FileInPath("heavyNeutrino/multilep/data/header/text.txt"),
)

def getJSON(is2017, is2018):
    if is2018:   return "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
    elif is2017: return "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
    else:        return "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

if isData:
  print('Sample is found to be 20%s data, will process using %s' % (yy, getJSON(is2017, is2018)))
  import FWCore.PythonUtilities.LumiList as LumiList
  jsonDir = os.path.expandvars('$CMSSW_BASE/src/heavyNeutrino/multilep/data/JSON')
  process.source.lumisToProcess = LumiList.LumiList(filename = os.path.join(jsonDir, getJSON(is2017, is2018))).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.egammaPostRecoSeq *
                     process.pileupJetIdUpdated *
                     process.jetSequence *
                     process.prefiringweight *
                     process.particleLevelSequence *        
                     process.rerunMvaIsolationSequence *
                     getattr(process,updatedTauName) *           
                     process.blackJackAndHookers)
