import FWCore.ParameterSet.Config as cms

def addJetSequence(process, isData):

  #
  # Latest JEC through globaltag, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  # Has (at time of writing) no effect (Moriond2017 miniAOD contains latest JEC)
  #
  process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
  if isData: jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
  else:      jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']

  from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
  updateJetCollection(
     process,
     jetSource = cms.InputTag('slimmedJets'),
     labelName = 'UpdatedJEC',
     jetCorrections = ('AK4PFchs', cms.vstring(jetCorrectorLevels), 'None')
  )

  # TODO: try add deepCSV, though no twiki existing for it?
  process.jetSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

  #
  # Jet energy resolution, see https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
  # Run three times the SmeredPATJetProducer for nominal, up and down variations
  #
  if not isData:
    for (i, j) in [(0, ''), (-1, 'Down'), (1, 'Up')]:
      jetSmearing = cms.EDProducer('SmearedPATJetProducer',
            src          = cms.InputTag('updatedPatJetsUpdatedJEC'),
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
