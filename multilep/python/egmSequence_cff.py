import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

def addElectronAndPhotonSequence(process):
  from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import files
  process.load('Configuration.StandardSequences.Services_cff')
  process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81), engineName = cms.untracked.string('TRandom3'),),
    calibratedPatPhotons    = cms.PSet( initialSeed = cms.untracked.uint32(81), engineName = cms.untracked.string('TRandom3'),),
  )
  process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
  process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
  process.calibratedPatElectrons.correctionFile = cms.string(files['Moriond17_23Jan'])
  process.calibratedPatPhotons.correctionFile   = cms.string(files['Moriond17_23Jan'])

  #
  # Set up electron and photon identifications
  #
  process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
  switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
  switchOnVIDPhotonIdProducer(  process, DataFormat.MiniAOD)
  electronModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
  photonModules   = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff',
                     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
  for idmod in electronModules: setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
  for idmod in photonModules:   setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

  process.egmSequence = cms.Sequence(process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.calibratedPatElectrons * process.calibratedPatPhotons)
