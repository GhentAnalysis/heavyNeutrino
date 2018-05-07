#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TLorentzVector.h"

LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer), 
    electronsEffectiveAreas(multilepAnalyzer->is2017 ? (iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas")).fullPath() ),
    muonsEffectiveAreas    (multilepAnalyzer->is2017 ? (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath() ),
    singleEleTrigs(multilepAnalyzer->is2017 ? iConfig.getParameter<std::vector<std::string> >("SingleEleTriggers2017") : iConfig.getParameter<std::vector<std::string> >("SingleEleTriggers")),
    singleMuoTrigs(multilepAnalyzer->is2017 ? iConfig.getParameter<std::vector<std::string> >("SingleMuoTriggers2017") : iConfig.getParameter<std::vector<std::string> >("SingleMuoTriggers"))
{
    if(!multilepAnalyzer->isData) genMatcher = new GenMatching(iConfig, multilepAnalyzer);
    if(multilepAnalyzer->isData){
        jecLevel = "L2L3Residual";
    } else {
	jecLevel = "L3Absolute";
    }
};

LeptonAnalyzer::~LeptonAnalyzer(){
  if(!multilepAnalyzer->isData) delete genMatcher;
}

void LeptonAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nL",                           &_nL,                           "_nL/i");
  outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/i");
  outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/i");
  outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/i");
  outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/i");
  outputTree->Branch("_nVFit",                        &_nVFit,                        "_nVFit/i");
  outputTree->Branch("_nGoodLeading",                 &_nGoodLeading,                 "_nGoodLeading/i");
  outputTree->Branch("_nGoodDisplaced",               &_nGoodDisplaced,               "_nGoodDisplaced/i");
  outputTree->Branch("_lIndex",                       &_lIndex,                       "_lIndex[_nL]/i");
  outputTree->Branch("_vertices",                     &_vertices,                     "_vertices[_nVFit][12]/D");
  outputTree->Branch("_lDisplaced",                   &_lDisplaced,                   "_lDisplaced[_nVFit][24]/D");
  outputTree->Branch("_lHasTrigger",                  &_lHasTrigger,                  "_lHasTrigger[_nL]/O");
  outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
  outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
  outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
  outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
  outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
  outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/i");
  outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
  outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
  outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
  outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
  outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
  outputTree->Branch("_2dIP",                         &_2dIP,                         "_2dIP[_nL]/D");
  outputTree->Branch("_2dIPSig",                      &_2dIPSig,                      "_2dIPSig[_nL]/D");
  outputTree->Branch("_lSimType",                     &_lSimType,                     "_lSimType[_nL]/I");
  outputTree->Branch("_lSimExtType",                  &_lSimExtType,                  "_lSimExtType[_nL]/I");
  outputTree->Branch("_lSimFlavour",                  &_lSimFlavour,                  "_lSimFlavour[_nL]/I");
  outputTree->Branch("_lElectronMva",                 &_lElectronMva,                 "_lElectronMva[_nLight]/F");
  outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
  outputTree->Branch("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", &_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, "_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL]/O");
  outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
  outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
  outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
  outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
  outputTree->Branch("_lpassConversionVeto",          &_lpassConversionVeto,          "_lpassConversionVeto[_nL]/O");
  outputTree->Branch("_eleNumberInnerHitsMissing",    &_eleNumberInnerHitsMissing,    "_eleNumberInnerHitsMissing[_nL]/D");
  outputTree->Branch("_lGlobalMuon",                  &_lGlobalMuon,                  "_lGlobalMuon[_nL]/O");
  outputTree->Branch("_lTrackerMuon",                 &_lTrackerMuon,                 "_lTrackerMuon[_nL]/O");
  outputTree->Branch("_lInnerTrackValidFraction",     &_lInnerTrackValidFraction,     "_lInnerTrackValidFraction[_nL]/D");
  outputTree->Branch("_lGlobalTrackNormalizeChi2",    &_lGlobalTrackNormalizeChi2,    "_lGlobalTrackNormalizeChi2[_nL]/D");
  outputTree->Branch("_lCQChi2Position",              &_lCQChi2Position,              "_lCQChi2Position[_nL]/D");
  outputTree->Branch("_lCQTrackKink",                 &_lCQTrackKink,                 "_lCQTrackKink[_nL]/D");
  outputTree->Branch("_muonSegComp",                  &_muonSegComp,                  "_muonSegComp[_nMu]/D");
  outputTree->Branch("_lNumberOfMatchedStation",      &_lNumberOfMatchedStation,      "_lNumberOfMatchedStation[_nL]/i");
  outputTree->Branch("_lNumberOfValidPixelHits",      &_lNumberOfValidPixelHits,      "_lNumberOfValidPixelHits[_nL]/i");
  outputTree->Branch("_muNumberInnerHits",            &_muNumberInnerHits,            "_muNumberInnerHits[_nL]/i");
  outputTree->Branch("_lTrackerLayersWithMeasurement",&_lTrackerLayersWithMeasurement,"_lTrackerLayersWithMeasurement[_nL]/i");

 	
  outputTree->Branch("_lMuTime",                    &_lMuTime,         "_lMuTime[_nL]/D");
  outputTree->Branch("_lMuTimeErr",                 &_lMuTimeErr,         "_lMuTimeErr[_nL]/D");
  outputTree->Branch("_lMuRPCTime",                 &_lMuRPCTime,         "_lMuRPCTime[_nL]/D");
  outputTree->Branch("_lMuRPCTimeErr",              &_lMuRPCTimeErr,         "_lMuRPCTimeErr[_nL]/D");	
  outputTree->Branch("_lMuTimenDof",                 &_lMuTimenDof,                 "_lMuTimenDof[_nL]/I");		
  outputTree->Branch("_lMuRPCTimenDof",                 &_lMuRPCTimenDof,                 "_lMuRPCTimenDof[_nL]/I");		

  outputTree->Branch("_lEleIsEB",                     &_lEleIsEB ,                    "_lEleIsEB[_nL]/O");
  outputTree->Branch("_lEleIsEE",                     &_lEleIsEE ,                    "_lEleIsEE[_nL]/O");
  outputTree->Branch("_lEleSuperClusterOverP",        &_lEleSuperClusterOverP ,       "_lEleSuperClusterOverP[_nL]/D");
  outputTree->Branch("_lEleEcalEnergy",               &_lEleEcalEnergy ,              "_lEleEcalEnergy[_nL]/D");
  outputTree->Branch("_lElefull5x5SigmaIetaIeta",     &_lElefull5x5SigmaIetaIeta ,    "_lElefull5x5SigmaIetaIeta[_nL]/D");
  outputTree->Branch("_lEleDEtaInSeed",               &_lEleDEtaInSeed ,              "_lEleDEtaInSeed[_nL]/D");
  outputTree->Branch("_lEleDeltaPhiSuperClusterTrackAtVtx", &_lEleDeltaPhiSuperClusterTrackAtVtx , "_lEleDeltaPhiSuperClusterTrackAtVtx[_nL]/D");
  outputTree->Branch("_lElehadronicOverEm",           &_lElehadronicOverEm ,          "_lElehadronicOverEm[_nL]/D");
  outputTree->Branch("_lEleInvMinusPInv",             &_lEleInvMinusPInv ,            "_lEleInvMinusPInv[_nL]/D");
  outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
  outputTree->Branch("_puCorr",                       &_puCorr,                       "_puCorr[_nL]/D");
  outputTree->Branch("_absIso03",                     &_absIso03,                     "_absIso03[_nL]/D");
  outputTree->Branch("_absIso04",                     &_absIso04,                     "_absIso04[_nL]/D");
  outputTree->Branch("_sumNeutralHadronEt04",         &_sumNeutralHadronEt04,         "_sumNeutralHadronEt04[_nL]/D");
  outputTree->Branch("_sumChargedHadronPt04",         &_sumChargedHadronPt04,         "_sumChargedHadronPt04[_nL]/D");
  outputTree->Branch("_sumPhotonEt04",                &_sumPhotonEt04,                "_sumPhotonEt04[_nL]/D");
  outputTree->Branch("_sumNeutralHadronEt03",         &_sumNeutralHadronEt03,         "_sumNeutralHadronEt03[_nL]/D");
  outputTree->Branch("_sumChargedHadronPt03",         &_sumChargedHadronPt03,         "_sumChargedHadronPt03[_nL]/D");
  outputTree->Branch("_sumPhotonEt03",                &_sumPhotonEt03,                "_sumPhotonEt03[_nL]/D");
  outputTree->Branch("_trackIso",                     &_trackIso ,                    "_trackIso[_nL]/D");
  outputTree->Branch("_ecalIso",                      &_ecalIso ,                     "_ecalIso[_nL]/D");
  outputTree->Branch("_hcalIso",                      &_hcalIso ,                     "_hcalIso[_nL]/D");
  outputTree->Branch("_deltaBIso",                    &_deltaBIso,                    "_deltaBIso[_nL]/D");
  outputTree->Branch("_ecalPFClusterIso",             &_ecalPFClusterIso ,            "_ecalPFClusterIso[_nL]/D");
  outputTree->Branch("_hcalPFClusterIso",             &_hcalPFClusterIso ,            "_hcalPFClusterIso[_nL]/D");
  outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
  outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
  outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
  outputTree->Branch("_tauMuonVeto",                  &_tauMuonVeto,                  "_tauMuonVeto[_nL]/O");
  outputTree->Branch("_tauEleVeto",                   &_tauEleVeto,                   "_tauEleVeto[_nL]/O");
  outputTree->Branch("_decayModeFindingNew",          &_decayModeFindingNew,          "_decayModeFindingNew[_nL]/O");
  outputTree->Branch("_tauVLooseMvaNew",              &_tauVLooseMvaNew,              "_tauVLooseMvaNew[_nL]/O");
  outputTree->Branch("_tauLooseMvaNew",               &_tauLooseMvaNew,               "_tauLooseMvaNew[_nL]/O");
  outputTree->Branch("_tauMediumMvaNew",              &_tauMediumMvaNew,              "_tauMediumMvaNew[_nL]/O");
  outputTree->Branch("_tauTightMvaNew",               &_tauTightMvaNew,               "_tauTightMvaNew[_nL]/O");
  outputTree->Branch("_tauVTightMvaNew",              &_tauVTightMvaNew,              "_tauVTightMvaNew[_nL]/O");
  outputTree->Branch("_tauVTightMvaOld",              &_tauVTightMvaOld,              "_tauVTightMvaOld[_nL]/O");
  outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
  outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
  outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
  
 
  if(!multilepAnalyzer->isData){
    outputTree->Branch("_lGenIndex",                  &_lGenIndex,                    "_lGenIndex[_nL]/i");
    outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
    outputTree->Branch("_lIsPromptFinalState",        &_lIsPromptFinalState,          "_lIsPromptFinalState[_nL]/O");
    outputTree->Branch("_lIsPromptDecayed",           &_lIsPromptDecayed,             "_lIsPromptDecayed[_nL]/O");
    outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
    outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
    outputTree->Branch("_lProvenanceCompressed",      &_lProvenanceCompressed,        "_lProvenanceCompressed[_nL]/i");
    outputTree->Branch("_lMatchPt",                   &_lMatchPt,                     "_lMatchPt[_nL]/D");
    outputTree->Branch("_lMatchEta",                  &_lMatchEta,                    "_lMatchEta[_nL]/D");
    outputTree->Branch("_lMatchPhi",                  &_lMatchPhi,                    "_lMatchPhi[_nL]/D");
    outputTree->Branch("_lMatchVertexX",              &_lMatchVertexX,                "_lMatchVertexX[_nL]/D");
    outputTree->Branch("_lMatchVertexY",              &_lMatchVertexY,                "_lMatchVertexY[_nL]/D");
    outputTree->Branch("_lMatchVertexZ",              &_lMatchVertexZ,                "_lMatchVertexZ[_nL]/D");
  }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& primaryVertex){
  edm::Handle<std::vector<pat::Electron>> electrons;            iEvent.getByToken(multilepAnalyzer->eleToken,               electrons);
  edm::Handle<edm::ValueMap<float>> electronsMva;               iEvent.getByToken(multilepAnalyzer->eleMvaToken,            electronsMva);
  edm::Handle<edm::ValueMap<float>> electronsMvaHZZ;            iEvent.getByToken(multilepAnalyzer->eleMvaHZZToken,         electronsMvaHZZ);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedVeto;       iEvent.getByToken(multilepAnalyzer->eleCutBasedVetoToken,   electronsCutBasedVeto);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedLoose;      iEvent.getByToken(multilepAnalyzer->eleCutBasedLooseToken,  electronsCutBasedLoose);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedMedium;     iEvent.getByToken(multilepAnalyzer->eleCutBasedMediumToken, electronsCutBasedMedium);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedTight;      iEvent.getByToken(multilepAnalyzer->eleCutBasedTightToken,  electronsCutBasedTight);
  edm::Handle<std::vector<pat::Muon>> muons;                    iEvent.getByToken(multilepAnalyzer->muonToken,              muons);
  edm::Handle<std::vector<pat::Tau>> taus;                      iEvent.getByToken(multilepAnalyzer->tauToken,               taus);
  edm::Handle<std::vector<pat::PackedCandidate>> packedCands;   iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,  packedCands);
  edm::Handle<double> rho;                                      iEvent.getByToken(multilepAnalyzer->rhoToken,               rho);
  edm::Handle<std::vector<pat::Jet>> jets;                      iEvent.getByToken(multilepAnalyzer->jetToken,               jets);
  edm::Handle<edm::TriggerResults> trigBits;                    iEvent.getByToken(multilepAnalyzer->triggerToken,           trigBits);
  edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjs; iEvent.getByToken(multilepAnalyzer->trigObjToken,           trigObjs);
  iSetup.get<IdealMagneticFieldRecord>().get(_bField);
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", _shProp);

  _nL     = 0;
  _nLight = 0;
  _nMu    = 0;
  _nEle   = 0;
  _nTau   = 0;
  _nVFit  = 0;
  _nGoodLeading = 0;
  _nGoodDisplaced = 0;

  //bool good_leading=false; // to check 1 leading-well_isolated lepton
  unsigned counter_index_leptons   = 0;
  //int counter_number_vertices = 0;

  //set up generator matching
  if(!multilepAnalyzer->isData) genMatcher->setGenParticles(iEvent);

  //loop over muons
  for(const pat::Muon& mu : *muons){
    // Check if muon passes minimum criteria
    if(passMuonPreselection(mu)==false) continue;

    counter_index_leptons++  ;                               // unique index to identify the 2 tracks for each vertex
    _lIndex[_nL] = counter_index_leptons;
    _lPFMuon[_nL]=  mu.isPFMuon();
	  
     const reco::MuonTime cmb = mu.time();
     const reco::MuonTime rpc = mu.rpcTime();  
		//csc + dt
	  	_lMuTimenDof[_nL] = cmb.nDof;
		if (cmb. Direction () == -1) {
		    _lMuTime[_nL] = cmb.timeAtIpOutIn;
		    _lMuTimeErr[_nL] = cmb.timeAtIpOutInErr;
		}
		if (cmb. Direction () == InsideOut) {
		    _lMuTime[_nL] = cmb.timeAtIpInOut;
		    _lMuTimeErr[_nL] = cmb.timeAtIpInOutErr;
		}
		if (cmb. Direction () == Undefined) {
		    _lMuTime[_nL] = 999999;
		    _lMuTimeErr[_nL] = 999999;
		}
		//RPC
		_lMuRPCTimenDof[_nL] = rpc.nDof;
		if (rpc. Direction () == OutsideIn) {
		    _lMuRPCTime[_nL] = rpc.timeAtIpOutIn;
		    _lMuRPCTimeErr[_nL] = rpc.timeAtIpOutInErr;
		}
		if (rpc. Direction () == InsideOut) {
		    _lMuRPCTime[_nL] = rpc.timeAtIpInOut;
		    _lMuRPCTimeErr[_nL] = rpc.timeAtIpInOutErr;
		}
		if (rpc. Direction () == Undefined) {
		    _lMuRPCTime[_nL] = 999999;
		    _lMuRPCTimeErr[_nL] = 999999;
		}
	  
	
	
	  
	  
    

    fillLeptonImpactParameters(mu, primaryVertex);

    fillLeptonKinVars(mu);
    fillLeptonIsoVars(mu, *rho);
    if(!multilepAnalyzer->isData) fillLeptonGenVars(mu, genMatcher);

        fillLeptonJetVariables(mu, jets, primaryVertex, *rho);
	

    _lGlobalMuon[_nL] = mu.isGlobalMuon();
    _lTrackerMuon[_nL]= mu.isTrackerMuon();
    _lInnerTrackValidFraction[_nL] = (!mu.innerTrack().isNull()) ?   mu.innerTrack()->validFraction()  : -1;
    _lGlobalTrackNormalizeChi2[_nL]= (!mu.globalTrack().isNull()) ?   mu.globalTrack()->normalizedChi2()  : -1;
    _lCQChi2Position[_nL] = mu.combinedQuality().chi2LocalPosition;
    _lCQTrackKink[_nL] = mu.combinedQuality().trkKink;
    _lNumberOfMatchedStation[_nL] = mu.numberOfMatchedStations();
    _lNumberOfValidPixelHits[_nL] = (!mu.innerTrack().isNull()) ?   mu.innerTrack()->hitPattern().numberOfValidPixelHits()  : 0; // cannot be -1 !!
    _lTrackerLayersWithMeasurement[_nL] = (!mu.innerTrack().isNull()) ?   mu.innerTrack()->hitPattern().trackerLayersWithMeasurement()  : 0; // cannot be -1 !!

    _lFlavor[_nL]        = 1;
    _muonSegComp[_nL]    = mu.segmentCompatibility();
    _relIso[_nL]         = getRelIso03(mu, *rho);   // Isolation variables

    _lSimType[_nL]       = mu.simType();
    _lSimExtType[_nL]    = mu.simExtType();
    _lSimFlavour[_nL]    = mu.simFlavour();

    // TODO: this is a possible solution to the missing trackRef, but maybe not what you want 
    _muNumberInnerHits[_nL]= (!mu.globalTrack().isNull()) ?   mu.globalTrack()->hitPattern().numberOfValidMuonHits() : (!mu.outerTrack().isNull() ? mu.outerTrack()->hitPattern().numberOfValidMuonHits() : 0); // cannot be -1 !!!
    _lPOGVeto[_nL]   = mu.isLooseMuon();
    _lPOGLoose[_nL]  = mu.isLooseMuon();
    if ( mu.isLooseMuon()) _lPOGMedium[_nL] = mu.isMediumMuon();
    if ( mu.isLooseMuon()) _lPOGTight[_nL]  = mu.isTightMuon(primaryVertex);
 
    _eleNumberInnerHitsMissing[_nL] =-1;
    _lpassConversionVeto[_nL] = false;
    _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL] = false;
    
    _lEleIsEB [_nL] = false;
    _lEleIsEE[_nL] = false;
    _lEleSuperClusterOverP[_nL] =  -1;
    _lEleEcalEnergy[_nL]=   -1;
    _lElefull5x5SigmaIetaIeta[_nL] =   -1;
    _lEleDEtaInSeed[_nL] =   -1;
    _lEleDeltaPhiSuperClusterTrackAtVtx[_nL] =   -1;
    _lElehadronicOverEm[_nL] =   -1;
    _lEleInvMinusPInv[_nL] =   -1;
	  
	  
    if (mu.pt() >  3 && std::abs(_dxy[_nL]) > 0.02) ++_nGoodDisplaced; 
    if (mu.pt() > 22 && std::abs(_dxy[_nL]) < 0.05 && std::abs(_dz[_nL])< 0.1 && getRelIso03(mu, *rho) < 0.3 && !mu.innerTrack().isNull() && (mu.isTrackerMuon() || mu.isGlobalMuon()) ) {
      ++_nGoodLeading;
      _lHasTrigger[_nL] = matchSingleTrigger(false, _lEta[_nL], _lPhi[_nL], iEvent.triggerNames(*trigBits), trigObjs);
    }
    else {
      _lHasTrigger[_nL] = false;
    }

    ++_nMu;
    ++_nL;
    ++_nLight;
  }

  // Loop over electrons (note: using iterator we can easily get the ref too)
  for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
    // Check if electron passes minimum criteria
    if(passElectronPreselection(*ele)==false) continue;

    auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele - electrons->begin()));
    counter_index_leptons++  ;                               // unique index to identify the 2 tracks for each vertex
    _lIndex[_nL] = counter_index_leptons;
    
	  
	  
	  _lMuRPCTimenDof[_nL] = -1;
	  _lMuTimenDof[_nL] = -1;
	  _lMuRPCTime[_nL] = -1;
	  _lMuRPCTimeErr[_nL] = -1;
	  _lMuTime[_nL] = -1;
	  _lMuTimeErr[_nL] = -1;
	
	  
	  
	  
    _eleNumberInnerHitsMissing[_nL]=ele->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    _muNumberInnerHits[_nL] = 0; // cannot be -1 !!
    fillLeptonImpactParameters(*ele, primaryVertex);
    fillLeptonIsoVars(*ele, *rho);

    _lpassConversionVeto[_nL] = ele->passConversionVeto();

    // ID ele variables	  
    _lEleIsEB [_nL] = ele->isEB();
    _lEleIsEE[_nL] = ele->isEE();
    _lEleSuperClusterOverP[_nL] = ele->eSuperClusterOverP();
    _lEleEcalEnergy[_nL]= ele->ecalEnergy();
    _lElefull5x5SigmaIetaIeta[_nL] = ele->full5x5_sigmaIetaIeta();
    _lEleDEtaInSeed[_nL] = std::abs(dEtaInSeed(&*ele));
    _lEleDeltaPhiSuperClusterTrackAtVtx[_nL] = std::abs(ele->deltaPhiSuperClusterTrackAtVtx());
    _lElehadronicOverEm[_nL] = ele->hadronicOverEm();
    _lEleInvMinusPInv[_nL] = std::abs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();

    fillLeptonKinVars(*ele);
    //fillLeptonGenVars(ele->genParticle());
    if(!multilepAnalyzer->isData) fillLeptonGenVars(*ele, genMatcher);

    fillLeptonJetVariables(*ele, jets, primaryVertex, *rho);
    _lFlavor[_nL]      = 0;
    _lEtaSC[_nL]       = ele->superCluster()->eta();
    _relIso[_nL]       = getRelIso03(*ele, *rho);

    _lElectronMva[_nL] = (*electronsMva)[electronRef];
    _lElectronPassEmu[_nL] = passTriggerEmulationDoubleEG(&*ele);
    _lElectronMvaHZZ[_nL]       = (*electronsMvaHZZ)[electronRef];

   
    _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL] = isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele);
    _lPOGVeto[_nL]     = (*electronsCutBasedVeto)[electronRef];
    _lPOGLoose[_nL]    = (*electronsCutBasedLoose)[electronRef]; // cut-based cuts
    _lPOGMedium[_nL]   = (*electronsCutBasedMedium)[electronRef];
    _lPOGTight[_nL]    = (*electronsCutBasedTight)[electronRef];             // Actually in SUS-17-001 we applied addtionaly lostHists==0, probably not a big impact

    //muon variables
    _lGlobalMuon[_nL] = false;
    _lTrackerMuon[_nL]= false;
    _lInnerTrackValidFraction[_nL] = -1;
    _lGlobalTrackNormalizeChi2[_nL]=  -1;
    _lCQChi2Position[_nL] =  -1;
    _lCQTrackKink[_nL] =  -1;
    _lNumberOfMatchedStation[_nL]       = 0; // cannot be -1 !!!
    _lNumberOfValidPixelHits[_nL]       = 0; // cannot be -1 !!!
    _lTrackerLayersWithMeasurement[_nL] = 0; // cannot be -1 !!!
    _lSimType[_nL]       = -1; 
    _lSimExtType[_nL]    = -1; 
    _lSimFlavour[_nL]    = -1; 
	  
    if(ele->pt() >  7 && std::abs(_dxy[_nL]) > 0.02) ++_nGoodDisplaced; 
    if(ele->pt() > 22 && std::abs(_dxy[_nL]) < 0.05 && std::abs(_dz[_nL])< 0.1 && _relIso[_nL] < 0.3 && !ele->gsfTrack().isNull() && _eleNumberInnerHitsMissing[_nL] <=2 && ele->passConversionVeto()) {
      ++_nGoodLeading;
      _lHasTrigger[_nL] = matchSingleTrigger(true, _lEta[_nL], _lPhi[_nL], iEvent.triggerNames(*trigBits), trigObjs);
    }
    else {
      _lHasTrigger[_nL] = false;
    }

    ++_nEle;
    ++_nL;
    ++_nLight;
  }

  //loop over taus
  for(const pat::Tau& tau : *taus){
    // Check if tau passes minimum criteria
    if(passTauPreselection(tau, _nL)==false) continue;

    fillLeptonKinVars(tau);
    //fillLeptonGenVars(tau.genParticle());
    if(!multilepAnalyzer->isData) fillLeptonGenVars(tau, genMatcher);
    fillLeptonImpactParameters(tau, primaryVertex);

    _lFlavor[_nL]  = 2;
    _tauMuonVeto[_nL] = tau.tauID("againstMuonLoose3");                        //Light lepton vetos
    _tauEleVeto[_nL] = tau.tauID("againstElectronLooseMVA6");

    _lPOGVeto[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");        //old tau ID
    _lPOGLoose[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    _lPOGMedium[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    _lPOGTight[_nL] = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    _tauVTightMvaOld[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

    _decayModeFindingNew[_nL] = tau.tauID("decayModeFindingNewDMs");           //new Tau ID 
    _tauVLooseMvaNew[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
    _tauLooseMvaNew[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
    _tauMediumMvaNew[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
    _tauTightMvaNew[_nL] = tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT");
    _tauVTightMvaNew[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT");
    _lSimType[_nL]       = -1; 
    _lSimExtType[_nL]    = -1; 
    _lSimFlavour[_nL]    = -1; 
   
    ++_nTau;
    ++_nL;
  }
		
  
  /*
   * refitting vertices displaced *********************************************************** 
   */
  // (Why double??) 
  unsigned iMu_plus=0;
  unsigned iMu_minus_mu=0;
  unsigned iMu_minus_e=0;
  unsigned iE_plus=_nMu;
  unsigned iE_minus_mu=_nMu;
  unsigned iE_minus_e=_nMu;
  cleanDileptonVertexArrays(_nVFit);

  for(const pat::Muon& mu_1 : *muons){ // for muons
    // Check if muon passes minimum criteria
    if(passMuonPreselection(mu_1)==false) continue;
    iMu_plus++;
    //+++++++++++++++    mu+
    if (mu_1.charge() < 0) continue;
    const reco::Track&  tk_1 = (!mu_1.innerTrack().isNull()) ? *mu_1.innerTrack () :  *mu_1.outerTrack () ;

    // ------------------  loop mu-
    iMu_minus_mu=0;
    for(const pat::Muon& mu_2 : *muons){ 
      // Check if muon passes minimum criteria
      if(passMuonPreselection(mu_2)==false) continue;
      iMu_minus_mu++;
      if (mu_2.charge() > 0) continue;  // only opposite charge
      const reco::Track&  tk_2 = (!mu_2.innerTrack().isNull()) ? *mu_2.innerTrack () :  *mu_2.outerTrack () ;
      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(dilvtx.isValid()) { 
	fillDileptonVertexArrays(_nVFit, iMu_plus, iMu_minus_mu, dilvtx, tk_1, tk_2);
	++_nVFit;   
      }
      else {   
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
	/*std::cout << "              1: "
		  << (mu_1.innerTrack().isNull() ? "OUT" : "TRK")
		  << " (" << tk_1.eta() << ", " << tk_1.phi() << ", " << tk_1.pt() << "); 2: "
		  << (mu_2.innerTrack().isNull() ? "OUT" : "TRK")
		  << " (" << tk_2.eta() << ", " << tk_2.phi() << ", " << tk_2.pt() << ")"
		  << std::endl;*/
      }
    }// end loop mu-
          
    // ------------------  loop e-
    iE_minus_mu=_nMu;
    for(auto ele_2 = electrons->begin(); ele_2 != electrons->end(); ++ele_2){
      // Check if electron passes minimum criteria
      if(passElectronPreselection(*ele_2)==false) continue;
      iE_minus_mu++; // it is already _nMu
      if(ele_2->charge() > 0) continue; // only opposite charge
      const reco::Track&  tk_2 =  *ele_2->gsfTrack() ;
      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(dilvtx.isValid()) { 
	fillDileptonVertexArrays(_nVFit, iMu_plus, iE_minus_mu, dilvtx, tk_1, tk_2);
	++_nVFit;   
      } 
      else {
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
	/*std::cout << "              1: "
		  << (mu_1.innerTrack().isNull() ? "OUT" : "TRK")
		  << " (" << tk_1.eta() << ", " << tk_1.phi() << ", " << tk_1.pt() << "); 2: "
		  << "GSF"
		  << " (" << tk_2.eta() << ", " << tk_2.phi() << ", " << tk_2.pt() << ")"
		  << std::endl;*/
      }
    }// end loop e-
  }//end loop mu

  iMu_minus_e=0;
  iE_plus=_nMu;
  iE_minus_e=_nMu;

  for(auto ele_1 = electrons->begin(); ele_1 != electrons->end(); ++ele_1){ // for electrons
    // Check if electron passes minimum criteria
    if(passElectronPreselection(*ele_1)==false) continue;
    iE_plus++;
    //+++++++++++++++++++++ e+
    if(ele_1->charge() < 0) continue;
    const reco::Track&  tk_1 =  *ele_1->gsfTrack() ;

    //------------------  loop mu+
    iMu_minus_e=0;
    for(const pat::Muon& mu_2 : *muons){ 
      // Check if muon passes minimum criteria
      if(passMuonPreselection(mu_2)==false) continue;
      iMu_minus_e++;
      if (mu_2.charge() > 0) continue;  // only opposite charge
      const reco::Track&  tk_2 = (!mu_2.innerTrack().isNull()) ? *mu_2.innerTrack () :  *mu_2.outerTrack () ;
      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(dilvtx.isValid()) { 
	fillDileptonVertexArrays(_nVFit, iE_plus, iMu_minus_e, dilvtx, tk_1, tk_2);
	++_nVFit;
      } 
      else {
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
	/*std::cout << "              1: "
		  << "GSF"
		  << " (" << tk_1.eta() << ", " << tk_1.phi() << ", " << tk_1.pt() << "); 2: "
		  << (mu_2.innerTrack().isNull() ? "OUT" : "TRK")
		  << " (" << tk_2.eta() << ", " << tk_2.phi() << ", " << tk_2.pt() << ")"
		  << std::endl;*/
      } 
    }// end loop mu-
    
    
    //------------------  loop e+
    iE_minus_e=_nMu;
    for(auto ele_2 = electrons->begin(); ele_2 != electrons->end(); ++ele_2){
      // Check if electron passes minimum criteria
      if(passElectronPreselection(*ele_2)==false) continue;
      iE_minus_e++;
      if(ele_2->charge() > 0) continue; // only opposite charge
      const reco::Track&  tk_2 =  *ele_2->gsfTrack();  
      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(dilvtx.isValid()) { 
	fillDileptonVertexArrays(_nVFit, iE_plus, iE_minus_e, dilvtx, tk_1, tk_2);
	++_nVFit;   
      } 
      else {
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
	/*std::cout << "              1: "
		  << "GSF"
		  << " (" << tk_1.eta() << ", " << tk_1.phi() << ", " << tk_1.pt() << "); 2: "
		  << "GSF"
		  << " (" << tk_2.eta() << ", " << tk_2.phi() << ", " << tk_2.pt() << ")"
		  << std::endl;*/
      }
    }// end loop e+
  }//end electrons

  //if(multilepAnalyzer->skim == "trilep"    and (_nLight < 3 || _nGoodLeading < 1 || _nGoodDisplaced < 2) ) return false;
  if(multilepAnalyzer->skim == "trilep"      and (_nLight < 3 || _nGoodLeading < 1                       ) ) return false;
  if(multilepAnalyzer->skim == "displtrilep" and (_nLight < 3 || _nGoodLeading < 1 || _nGoodDisplaced < 2) ) return false;

  if(multilepAnalyzer->skim == "dilep"       and _nLight <  2) return false;
  if(multilepAnalyzer->skim == "ttg"         and _nLight <  2) return false;
  if(multilepAnalyzer->skim == "singlelep"   and _nLight <  1) return false;
  if(multilepAnalyzer->skim == "FR"          and _nLight != 1) return false;

  return true;
}



/*
 * //--// Refit dilepton vertex:
 * Provide a transientvertex 
 */
TransientVertex LeptonAnalyzer::dileptonVertex(const reco::Track& tk1, const reco::Track& tk2) {
  MagneticField *bfield = new OAEParametrizedMagneticField("3_8T"); 
  std::vector<reco::TransientTrack> ttks;
  ttks.push_back(reco::TransientTrack(tk1, bfield)); 
  ttks.push_back(reco::TransientTrack(tk2, bfield)); 
  KalmanVertexFitter vtxFitter;
  return vtxFitter.vertex(ttks); 
} 


void LeptonAnalyzer::cleanDileptonVertexArrays(unsigned nVFit){
   for (int i =0; i < 50; i++){
	  for (int j =0; j < 24 ; j++){
		  if (j < 12) _vertices[i][j] = 0;
		  _lDisplaced[i][j] = 0;
	  }
   }	
}
// Fill the arrays of displaced vertices and leptons 
void LeptonAnalyzer::fillDileptonVertexArrays(unsigned nVFit, unsigned iL_plus, unsigned iL_minus,
					      const TransientVertex& dvtx,
					      const reco::Track& tk1, const reco::Track& tk2) {
  _vertices[nVFit][0]  = iL_plus*100 + iL_minus;
  _vertices[nVFit][1]  = dvtx.position().x();
  _vertices[nVFit][2]  = dvtx.position().y();
  _vertices[nVFit][3]  = dvtx.position().z();
  _vertices[nVFit][4]  = dvtx.positionError().cxx();
  _vertices[nVFit][5]  = dvtx.positionError().cyy();
  _vertices[nVFit][6]  = dvtx.positionError().czz();
  _vertices[nVFit][7]  = dvtx.positionError().cyx();
  _vertices[nVFit][8]  = dvtx.positionError().czy();
  _vertices[nVFit][9]  = dvtx.positionError().czx();
  _vertices[nVFit][10] = dvtx.degreesOfFreedom();
  _vertices[nVFit][11] = dvtx.totalChiSquared();

  GlobalPoint  vtxpos(_vertices[nVFit][1], _vertices[nVFit][2], _vertices[nVFit][3]);

  GlobalPoint  l1r(tk1.vx(), tk1.vy(), tk1.vz());
  GlobalVector l1p(tk1.px(), tk1.py(), tk1.pz());
  GlobalTrajectoryParameters l1gtp(l1r, l1p, tk1.charge(), _bField.product()); 
  CurvilinearTrajectoryError l1cov(tk1.covariance());
  FreeTrajectoryState l1fts(l1gtp, l1cov);
  FreeTrajectoryState l1newfts = _shProp->propagate(l1fts, vtxpos);
  if(!l1newfts.hasCurvilinearError()) { // instead of isValid()... 
    std::cout << "Propagation of L1 to dilepton vertex (" << _vertices[nVFit][0] << ") failed!" << std::endl;
    //return false;
    l1newfts = l1fts;
  }
  // Position
  _lDisplaced[nVFit][0] = l1newfts.position().x();
  _lDisplaced[nVFit][1] = l1newfts.position().y();
  _lDisplaced[nVFit][2] = l1newfts.position().z();
  // Momentum
  _lDisplaced[nVFit][3] = l1newfts.momentum().x();
  _lDisplaced[nVFit][4] = l1newfts.momentum().y();
  _lDisplaced[nVFit][5] = l1newfts.momentum().z();
  // Position error
  _lDisplaced[nVFit][6]  = (l1newfts.cartesianError().matrix())(0, 0);
  _lDisplaced[nVFit][7]  = (l1newfts.cartesianError().matrix())(1, 1);
  _lDisplaced[nVFit][8]  = (l1newfts.cartesianError().matrix())(2, 2);
  // Momentum error
  _lDisplaced[nVFit][9]  = (l1newfts.cartesianError().matrix())(3, 3);
  _lDisplaced[nVFit][10] = (l1newfts.cartesianError().matrix())(4, 4);
  _lDisplaced[nVFit][11] = (l1newfts.cartesianError().matrix())(5, 5);

  GlobalPoint  l2r(tk2.vx(), tk2.vy(), tk2.vz());
  GlobalVector l2p(tk2.px(), tk2.py(), tk2.pz());
  GlobalTrajectoryParameters l2gtp(l2r, l2p, tk2.charge(), _bField.product()); 
  CurvilinearTrajectoryError l2cov(tk2.covariance());
  FreeTrajectoryState l2fts(l2gtp, l2cov);
  FreeTrajectoryState l2newfts = _shProp->propagate(l2fts, vtxpos);
  if(!l2newfts.hasCurvilinearError()) { // instead of isValid()... 
    std::cout << "Propagation of L2 to dilepton vertex (" << _vertices[nVFit][0] << ") failed!" << std::endl;
    //return false;
    l2newfts = l2fts;
  }
  // Position
  _lDisplaced[nVFit][12] = l2newfts.position().x();
  _lDisplaced[nVFit][13] = l2newfts.position().y();
  _lDisplaced[nVFit][14] = l2newfts.position().z();
  // Momentum
  _lDisplaced[nVFit][15] = l2newfts.momentum().x();
  _lDisplaced[nVFit][16] = l2newfts.momentum().y();
  _lDisplaced[nVFit][17] = l2newfts.momentum().z();
  // Position error
  _lDisplaced[nVFit][18]  = (l2newfts.cartesianError().matrix())(0, 0);
  _lDisplaced[nVFit][19]  = (l2newfts.cartesianError().matrix())(1, 1);
  _lDisplaced[nVFit][20]  = (l2newfts.cartesianError().matrix())(2, 2);
  // Momentum error
  _lDisplaced[nVFit][21]  = (l2newfts.cartesianError().matrix())(3, 3);
  _lDisplaced[nVFit][22] = (l2newfts.cartesianError().matrix())(4, 4);
  _lDisplaced[nVFit][23] = (l2newfts.cartesianError().matrix())(5, 5);

}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton){
  _lPt[_nL]     = lepton.pt();
  _lEta[_nL]    = lepton.eta();
  _lPhi[_nL]    = lepton.phi();
  _lE[_nL]      = lepton.energy();
  _lCharge[_nL] = lepton.charge();
}

void LeptonAnalyzer::fillLeptonIsoVars(const pat::Muon& mu, const double rho){
  _puCorr[_nL] = rho*muonsEffectiveAreas.getEffectiveArea(mu.eta());
  double pucorr= rho*muonsEffectiveAreas.getEffectiveArea(mu.eta());
  _absIso03 [_nL] = mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - pucorr);
  _absIso04 [_nL] = mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - pucorr);
  _sumNeutralHadronEt04 [_nL] = mu.pfIsolationR04().sumNeutralHadronEt;
  _sumChargedHadronPt04 [_nL] = mu.pfIsolationR04().sumChargedHadronPt;
  _sumPhotonEt04[_nL]   = mu.pfIsolationR04().sumPhotonEt ;
  _sumNeutralHadronEt03 [_nL] = mu.pfIsolationR03().sumNeutralHadronEt;
  _sumChargedHadronPt03 [_nL] = mu.pfIsolationR03().sumChargedHadronPt;
  _sumPhotonEt03[_nL]    = mu.pfIsolationR03().sumPhotonEt ;
  _trackIso[_nL]         = mu.trackIso	();
  _ecalIso[_nL]          = mu.ecalIso ()  ;
  _hcalIso[_nL]          = mu.hcalIso()  ;
  _deltaBIso[_nL]        = mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumPhotonEt + mu.pfIsolationR03().sumNeutralHadronEt - 0.5*pucorr);
  _ecalPFClusterIso[_nL] =  -1.;
  _hcalPFClusterIso[_nL] =  -1.; 
		
}


void LeptonAnalyzer::fillLeptonIsoVars(const pat::Electron& ele, const double rho){
  _puCorr[_nL] = rho*electronsEffectiveAreas.getEffectiveArea(ele.superCluster()->eta());
  double pucorr= rho*electronsEffectiveAreas.getEffectiveArea(ele.superCluster()->eta());
  _absIso03 [_nL] = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - pucorr);
  _absIso04 [_nL] = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - pucorr);
  _sumNeutralHadronEt04 [_nL] = ele.pfIsolationVariables().sumNeutralHadronEt;
  _sumChargedHadronPt04 [_nL] = ele.pfIsolationVariables().sumChargedHadronPt;
  _sumPhotonEt04[_nL]   = ele.pfIsolationVariables().sumPhotonEt ;
  _sumNeutralHadronEt03 [_nL] = ele.pfIsolationVariables().sumNeutralHadronEt;
  _sumChargedHadronPt03 [_nL] = ele.pfIsolationVariables().sumChargedHadronPt;
  _sumPhotonEt03[_nL]   = ele.pfIsolationVariables().sumPhotonEt ;
  _trackIso[_nL]        = ele.trackIso	();
  _ecalIso[_nL]         = ele.ecalIso ()  ;
  _hcalIso[_nL]         = ele.hcalIso()  ;
  _deltaBIso[_nL]       = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumPhotonEt +  ele.pfIsolationVariables().sumNeutralHadronEt - 0.5*pucorr);
  _ecalPFClusterIso[_nL]= ele.ecalPFClusterIso();
  _hcalPFClusterIso[_nL]= ele.hcalPFClusterIso();
}
template <typename Lepton> void LeptonAnalyzer::fillLeptonGenVars(const Lepton& lepton, GenMatching* genMatcher){
    genMatcher->fillMatchingVars(lepton);
    _lGenIndex[_nL] = genMatcher->genLIndex();
    _lIsPrompt[_nL] = genMatcher->promptMatch();
    _lIsPromptFinalState[_nL] = genMatcher->promptFinalStateMatch();
    _lIsPromptDecayed[_nL] = genMatcher->promptDecayedMatch();

    _lMatchPdgId[_nL] = genMatcher->pdgIdMatch();
    _lProvenance[_nL] = genMatcher->getProvenance();
    _lProvenanceCompressed[_nL] = genMatcher->getProvenanceCompressed();
    _lProvenanceConversion[_nL] = genMatcher->getProvenanceConversion();
    _lMatchPt[_nL] = genMatcher->getMatchPt();
    _lMatchEta[_nL] = genMatcher->getMatchEta();
    _lMatchPhi[_nL] = genMatcher->getMatchPhi();
    _lMatchVertexX[_nL] = genMatcher->getMatchVertexX();
    _lMatchVertexY[_nL] = genMatcher->getMatchVertexY();
    _lMatchVertexZ[_nL] = genMatcher->getMatchVertexZ();
}



/*
 * Impact parameters:
 * Provide PV to dxy/dz otherwise you get dxy/dz to the beamspot instead of the primary vertex
 * For taus: dxy is pre-computed with PV it was constructed with
 */
void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Electron& ele, const reco::Vertex& vertex){
  _dxy[_nL]     = ele.gsfTrack()->dxy(vertex.position());
  _dz[_nL]      = ele.gsfTrack()->dz(vertex.position());
  _3dIP[_nL]    = ele.dB(pat::Electron::PV3D);
  _3dIPSig[_nL] = ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D);
  _2dIP[_nL]    = ele.dB();
  _2dIPSig[_nL] = ele.dB()/ele.edB();
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon, const reco::Vertex& vertex){
		
  _dxy[_nL]     = (!muon.innerTrack().isNull()) ? muon.innerTrack()->dxy(vertex.position()) : muon.outerTrack()->dxy(vertex.position());
  _dz[_nL]      = (!muon.innerTrack().isNull()) ? muon.innerTrack()->dz(vertex.position()) : muon.outerTrack()->dz(vertex.position());
  _3dIP[_nL]    = muon.dB(pat::Muon::PV3D);
  _3dIPSig[_nL] = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D);
  _2dIP[_nL]    = muon.dB();
  _2dIPSig[_nL] = muon.dB()/muon.edB();
}




//Check if electron overlaps with loose muon
bool LeptonAnalyzer::eleMuOverlap(const pat::Electron& ele, const bool* loose) const {
  TLorentzVector eleV;
  eleV.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy());
  for(unsigned m = 0; m < _nMu; ++m){
      
    if(_lPOGLoose[m]){                   // changed from HNL loose
     TLorentzVector muV;
            muV.SetPtEtaPhiE(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
            if(eleV.DeltaR(muV) < 0.05) return true;
        }
    }
    return false;
}


//Check if tau overlaps with light lepton
bool LeptonAnalyzer::tauLightOverlap(const pat::Tau& tau, const bool* loose){
  TLorentzVector tauV;
  tauV.SetPtEtaPhiE(tau.pt(), tau.eta(), tau.phi(), tau.energy());
  for(unsigned l = 0; l < _nLight; ++l){
    if(loose[l]){
      TLorentzVector lightV;
      lightV.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
      if(tauV.DeltaR(lightV) < 0.4) return true;
    }
  }
  return false;
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Tau& tau, const reco::Vertex& vertex){
  _dxy[_nL]     = (double) tau.dxy();                                      // warning: float while dxy of tracks are double; could also return -1000
  _dz[_nL]      = tau_dz(tau, vertex.position());
  _3dIP[_nL]    = tau.ip3d();
  _3dIPSig[_nL] = tau.ip3d_Sig();
}

//Function returning tau dz
double LeptonAnalyzer::tau_dz(const pat::Tau& tau, const reco::Vertex::Point& vertex){
  const reco::Candidate::Point& tauVtx = tau.leadChargedHadrCand()->vertex();
  return (tauVtx.Z() - vertex.z()) - ((tauVtx.X() - vertex.x())*tau.px()+(tauVtx.Y()-vertex.y())*tau.py())/tau.pt()*tau.pz()/tau.pt();
}

// To synchronize lepton selection
bool LeptonAnalyzer::passElectronPreselection(const pat::Electron& elec) const {
  //// Copied from the electron loop
  // if(!elec.hasTrackDetails())                                                            return false;
  // if(elec.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)>2) return false;
  // if(!elec.passConversionVeto())                                                         return false;  
  // if(std::abs(_dxy[_nL])>0.05)                                                               return false;
  // if(std::abs(_dz[_nL])>0.1)                                                                 return false;
  // if(_relIso[_nL]>1)                                                                     return false;
  if(elec.gsfTrack().isNull())     return false; 
  if(elec.pt()<5.)                 return false;
  if(std::abs(elec.eta())>2.5)     return false;
  if(!isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&elec)) return false;
  if(eleMuOverlap(elec, _lPFMuon)) return false; // overlap muon-electron deltaR<0.05

  return true;
}
//
bool LeptonAnalyzer::passMuonPreselection(const pat::Muon& muon) const {
  //// Copied from the muon loop
  // if(muon.innerTrack().isNull())                     return false;
  // if(!(muon.isTrackerMuon() || muon.isGlobalMuon())) return false;
  // if(!muon.isMediumMuon())                           return false;
  // if(std::abs(_dxy[_nL])>0.05)                       return false;
  // if(std::abs(_dz[_nL])>0.1)                         return false;
  // if(!muon.hasTrackDetails())                        return false;
  // if(_relIso[_nL]>1)                                 return false;
  // if(!_lPOGLoose[_nL])                               return false;
  if(!muon.isPFMuon())         return false;
  if(muon.pt()<3)              return false;
  if(std::abs(muon.eta())>2.4) return false;

  return true;
}
//
bool LeptonAnalyzer::passTauPreselection(const pat::Tau& tauh, unsigned taupos) const {
  //// Copied from the tau loop
  if(tauh.pt()<20.)                   return false; // Minimum pt for tau reconstruction
  if(std::abs(tauh.eta())>2.3)        return false;
  if(!tauh.tauID("decayModeFinding")) return false;
  if(_dz[taupos]<0.4)                 return false; // tau dz cut used in ewkino

  return true;
}


void LeptonAnalyzer::fillLeptonJetVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::Jet>>& jets, const reco::Vertex& vertex, const double rho){
    //Make skimmed "close jet" collection
    std::vector<pat::Jet> selectedJetsAll;
    for(auto jet = jets->cbegin(); jet != jets->cend(); ++jet){
        double jetPt = jet->pt()*multilepAnalyzer->jec->jetCorrection(jet->correctedP4("Uncorrected").Pt(), jet->correctedP4("Uncorrected").Eta(), rho, jet->jetArea(), jecLevel); 
        if( jetPt > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
        //if( jet->pt() > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
    }
    // Find closest selected jet
    unsigned closestIndex = 0;
    for(unsigned j = 1; j < selectedJetsAll.size(); ++j){
        if(reco::deltaR(selectedJetsAll[j], lepton) < reco::deltaR(selectedJetsAll[closestIndex], lepton)) closestIndex = j;
    }
    const pat::Jet& jet = selectedJetsAll[closestIndex];
    if(selectedJetsAll.size() == 0 || reco::deltaR(jet, lepton) > 0.4){ //Now includes safeguard for 0 jet events
        _ptRatio[_nL] = 1;
        _ptRel[_nL] = 0;
        _closestJetCsvV2[_nL] = 0;
        _closestJetDeepCsv_b[_nL] = 0;
        _closestJetDeepCsv_bb[_nL] = 0;
        _selectedTrackMult[_nL] = 0;
    } else {
        
        double totalJEC = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), rho, jet.jetArea(), jecLevel);
        double l1JEC = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), rho, jet.jetArea(), "L1FastJet");
        TLorentzVector l1Jet;
        l1Jet.SetPtEtaPhiE(jet.correctedP4("Uncorrected").Pt()*l1JEC, jet.correctedP4("Uncorrected").Eta(), jet.correctedP4("Uncorrected").Phi(), jet.correctedP4("Uncorrected").E()*l1JEC);
        TLorentzVector l(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());
        float JEC = totalJEC/l1JEC;
        TLorentzVector lepAwareJet = (l1Jet - l)*JEC + l;
        
        /*
        auto  l1Jet       = jet.correctedP4("L1FastJet");
        float JEC         = jet.p4().E()/l1Jet.E();
        auto  l           = lepton.p4();
        auto  lepAwareJet = (l1Jet - l)*JEC + l;
        */

        TLorentzVector lV(l.Px(), l.Py(), l.Pz(), l.E());
        TLorentzVector jV(lepAwareJet.Px(), lepAwareJet.Py(), lepAwareJet.Pz(), lepAwareJet.E());
        _ptRatio[_nL]       = l.Pt()/lepAwareJet.Pt();
        _ptRel[_nL]         = lV.Perp((jV - lV).Vect());
        _closestJetCsvV2[_nL] = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        _closestJetDeepCsv_b[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        //compute selected track multiplicity of closest jet
        _selectedTrackMult[_nL] = 0;
        for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
            const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
            if(daughter->hasTrackDetails()){
                const reco::Track& daughterTrack = daughter->pseudoTrack();
                TLorentzVector trackVec(daughterTrack.px(), daughterTrack.py(), daughterTrack.pz(), daughterTrack.p());
                double daughterDeltaR            = trackVec.DeltaR(jV);
                bool goodTrack                   = daughterTrack.pt() > 1 && daughterTrack.charge() != 0 && daughterTrack.hitPattern().numberOfValidHits() > 7
                    && daughterTrack.hitPattern().numberOfValidPixelHits() > 1 && daughterTrack.normalizedChi2() < 5 && std::abs(daughterTrack.dz(vertex.position())) < 17
                    && std::abs(daughterTrack.dxy(vertex.position())) < 17;
                if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack) ++_selectedTrackMult[_nL];
            }
        }
    }
}

bool LeptonAnalyzer::matchSingleTrigger(bool isele, double aeta, double aphi, 
					const edm::TriggerNames &names, 
					edm::Handle<pat::TriggerObjectStandAloneCollection> objs) {
  for(pat::TriggerObjectStandAlone iobj : *objs) { // NOTE: not const nor by reference, because we need to 'unpackPathNames'
    if(reco::deltaR(iobj.eta(), iobj.phi(), aeta, aphi)<0.15) {
      iobj.unpackPathNames(names);
      std::vector<std::string> &singletrigs = isele ? singleEleTrigs : singleMuoTrigs;
      for(std::string& itrig : singletrigs) {
	if(iobj.hasPathName(itrig.c_str(), true, true)) {
	  return true;
	}
      }
    }
  }

  return false;
}
