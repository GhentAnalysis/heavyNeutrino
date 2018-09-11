import FWCore.ParameterSet.Config as cms
import os

def addJetSequence(process, isData, is2017):

  #
  # Latest JEC through globaltag, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  # Has (at time of writing) no effect (Moriond2017 miniAOD contains latest JEC)
  #
  process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
  process.load('Configuration.StandardSequences.MagneticField_cff')  # needed for pfImpactParameterTagInfos
  if isData: jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
  else:      jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']

  from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

  if not is2017 :
    updateJetCollection(
      process,
      jetSource = cms.InputTag('slimmedJets'),
      labelName = 'UpdatedJEC',
      jetCorrections = ('AK4PFchs', cms.vstring(jetCorrectorLevels), 'None'),
      # DeepCSV twiki: https://twiki.cern.ch/twiki/bin/view/CMS/DeepFlavour
      btagDiscriminators = [
        'pfCombinedSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probudsg',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb',
       #'pfDeepCSVJetTags:probcc', # not available in CMSSW_9_X_Y, also not really needed for us
      ]
     )
  else: 
     updateJetCollection(
       process,
       jetSource = cms.InputTag('slimmedJets'),
       labelName = 'UpdatedJEC',
       jetCorrections = ('AK4PFchs', cms.vstring(jetCorrectorLevels), 'None')
    ) 

  if os.environ['CMSSW_BASE'].count('CMSSW_9'):                                                                     # From CMSSW_9_X, the patAlgosToolsTask contains the whole updateJetCollection sequence
    process.jetSequence = cms.Sequence(process.patAlgosToolsTask)
  else:
    process.jetSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC *
                                       process.pfImpactParameterTagInfosUpdatedJEC *
                                       process.pfSecondaryVertexTagInfosUpdatedJEC *
                                       process.pfCombinedSecondaryVertexV2BJetTagsUpdatedJEC *
                                       process.patJetCorrFactorsTransientCorrectedUpdatedJEC *
                                       process.pfInclusiveSecondaryVertexFinderTagInfosUpdatedJEC *
                                       process.pfDeepCSVTagInfosUpdatedJEC *
                                       process.pfDeepCSVJetTagsUpdatedJEC *
                                       process.updatedPatJetsTransientCorrectedUpdatedJEC *
                                       process.selectedUpdatedPatJetsUpdatedJEC)

  #
  # Jet energy resolution, see https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
  # Run three times the SmeredPATJetProducer for nominal, up and down variations
  #
  if not isData:
    for (i, j) in [(0, ''), (-1, 'Down'), (1, 'Up')]:
      jetSmearing = cms.EDProducer('SmearedPATJetProducer',
            src          = cms.InputTag('selectedUpdatedPatJetsUpdatedJEC'),
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
      setattr(process, 'slimmedJetsCorrectedAndSmeared'+j, jetSmearing)
      process.jetSequence *= jetSmearing

  # Propagate JEC to MET (need to add fullPatMetSequence to path) (maybe good to add here link to a twiki page, if it exist)
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD   # currently broken

  print( "isData = {}".format( isData) )
  print( "is2017 = {}".format( is2017) )
  runMetCorAndUncFromMiniAOD (
    process,
    isData = isData, # false for MC
    fixEE2017 = is2017
    #postfix = "ModifiedMET" if is2017 else ""
  )
