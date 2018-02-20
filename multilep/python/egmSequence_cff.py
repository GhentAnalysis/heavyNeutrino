import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

def addElectronAndPhotonSequence(process, isData):
  #
  # EGM Regression
  #
  from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
  process = regressionWeights(process)
  process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

  #
  # EGM scale corrections (on data) and smearing (on MC)
  #
  from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import files
  process.load('Configuration.StandardSequences.Services_cff')
  process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81), engineName = cms.untracked.string('TRandom3'),),
    calibratedPatPhotons    = cms.PSet( initialSeed = cms.untracked.uint32(81), engineName = cms.untracked.string('TRandom3'),),
  )
  process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
  process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
  for calmod in [process.calibratedPatElectrons, process.calibratedPatPhotons]:
    calmod.correctionFile = cms.string(files['Moriond17_23Jan'])
    calmod.isMC           = cms.bool(not isData)

  #
  # Set up electron and photon identifications
  #
  process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
  switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
  switchOnVIDPhotonIdProducer(  process, DataFormat.MiniAOD)
  electronModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff']
  photonModules   = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff',
                     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
  for idmod in electronModules: setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
  for idmod in photonModules:   setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

  process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
  process.electronIDValueMapProducer.srcMiniAOD  = cms.InputTag('slimmedElectrons')
  process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
  process.photonIDValueMapProducer.srcMiniAOD    = cms.InputTag('slimmedPhotons')
  process.photonMVAValueMapProducer.srcMiniAOD   = cms.InputTag('slimmedPhotons')

  process.egmSequence = cms.Sequence(process.regressionApplication * process.calibratedPatElectrons * process.calibratedPatPhotons * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence)
