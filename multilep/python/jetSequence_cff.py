import FWCore.ParameterSet.Config as cms
import os

def addJetSequence( process, inputFile, isData, is2017, is2018, isFastSim, isUL ):
  #
  # Latest JEC through globaltag, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  #
  process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
  process.load('Configuration.StandardSequences.MagneticField_cff')  # needed for pfImpactParameterTagInfos
  if isData: jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
  else:      jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']

  from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
  from CondCore.CondDB.CondDB_cfi import CondDB
 
  #
  # Load specific JEC through sqlite file
  #
  if isData:
    if isUL:
        if is2017:
            if 'Run2017B' in inputFile: JECVersion = 'Summer19UL17_RunB_V5_DATA'
            if 'Run2017C' in inputFile: JECVersion = 'Summer19UL17_RunC_V5_DATA'
            if 'Run2017D' in inputFile: JECVersion = 'Summer19UL17_RunD_V5_DATA'
            if 'Run2017E' in inputFile: JECVersion = 'Summer19UL17_RunE_V5_DATA'
            if 'Run2017F' in inputFile: JECVersion = 'Summer19UL17_RunF_V5_DATA'
        elif is2018:
            if 'Run2018A' in inputFile: JECVersion = 'Summer19UL18_RunA_V4_DATA'
            if 'Run2018B' in inputFile: JECVersion = 'Summer19UL18_RunB_V4_DATA'
            if 'Run2018C' in inputFile: JECVersion = 'Summer19UL18_RunC_V4_DATA'
            if 'Run2018D' in inputFile: JECVersion = 'Summer19UL18_RunD_V4_DATA'

        else: JECVersion = ''
    else: JECVersion = ''
  else:
    if isUL:
        if is2017:   JECVersion = 'Summer19UL17_V5_MC'
        #if is2017:   JECVersion = 'Summer19UL17_V1_ComplexL1_MC'
        elif is2018: JECVersion = 'Summer19UL18_V2_MC'
        else:        JECVersion = ''
    else:
      if isFastSim : 
        if is2018:   JECVersion = 'Autumn18_FastSimV1_MC'
        elif is2017: JECVersion = 'Fall17_FastSimV1_MC'
        else :       JECVersion = 'Summer16_FastSimV1_MC'
      else:
        if is2018:   JECVersion = 'Autumn18_V19_MC'
        elif is2017: JECVersion = 'Fall17_17Nov2017_V32_102X_MC'
        else:        JECVersion = 'Summer16_07Aug2017_V11_MC'


  if not JECVersion == '':
    print 'using JEC: {}'.format(JECVersion)
    CondDBJECFile = CondDB.clone( connect = cms.string('sqlite_fip:heavyNeutrino/multilep/data/JEC/{}.db'.format( JECVersion ) ) )
    process.jec = cms.ESSource('PoolDBESSource',
      CondDBJECFile,
      toGet = cms.VPSet(
        cms.PSet(
          record = cms.string('JetCorrectionsRecord'),
          tag    = cms.string('JetCorrectorParametersCollection_{}_AK4PFchs'.format( JECVersion ) ),
          label  = cms.untracked.string('AK4PFchs')
        ),
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_{}_AK4PFPuppi'.format( JECVersion ) ),
            label  = cms.untracked.string('AK4PFPuppi')
        )
      )
    )
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', cms.vstring(jetCorrectorLevels), 'None'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = [
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg'
    ],
  )

  process.jetSequence = cms.Sequence(process.patAlgosToolsTask)

  #
  # Recompute MET with preferred JEC (need to add fullPatMetSequence to path)
  # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10
  # Make sure the correct JEC has been loaded
  #
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  #Type 1 PFMET
  runMetCorAndUncFromMiniAOD (
      process,
      isData = isData
  )
  
  #PUPPI MET
  from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
  makePuppiesFromMiniAOD( process, True );
  runMetCorAndUncFromMiniAOD(process,
                             isData= isData,
                             metType="Puppi",
                             postfix="Puppi",
                             jetFlavor="AK4PFPuppi",
                             )
  process.puppiNoLep.useExistingWeights = False
  process.puppi.useExistingWeights = False
  
  from CommonTools.PileupAlgos.customizePuppiTune_cff import UpdatePuppiTuneV15
  UpdatePuppiTuneV15(process, not isData)

  #Add MET sequences to path
  process.jetSequence *= process.puppiMETSequence
  process.jetSequence *= process.fullPatMetSequence
  process.jetSequence *= process.fullPatMetSequencePuppi
  
  #
  # Jet energy resolution, see https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
  # Run three times the SmearedPATJetProducer for nominal, up and down variations
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

  #
  # To get updated ecalBadCalibReducedMINIAODFilter
  # See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
  # Recipe is preliminary, i.e. recommended to check for updates
  # UltraLegacy doesn't have working filters yet! add them when available
  #
  if(is2017 or is2018):
    process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

    baddetEcallist = cms.vuint32(
        [872439604,872422825,872420274,872423218,
         872423215,872416066,872435036,872439336,
         872420273,872436907,872420147,872439731,
         872436657,872420397,872439732,872439339,
         872439603,872422436,872439861,872437051,
         872437052,872420649,872422436,872421950,
         872437185,872422564,872421566,872421695,
         872421955,872421567,872437184,872421951,
         872421694,872437056,872437057,872437313])


    process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
        ecalMinEt        = cms.double(50.),
        baddetEcal       = baddetEcallist,
        taggingMode      = cms.bool(True),
        debug            = cms.bool(False)
        )

    process.jetSequence *= process.ecalBadCalibReducedMINIAODFilter
