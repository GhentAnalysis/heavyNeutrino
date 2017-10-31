#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TLorentzVector.h"

LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  electronsEffectiveAreas((iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas")).fullPath()),
  muonsEffectiveAreas(    (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath()),
  multilepAnalyzer(multilepAnalyzer)
{
  leptonMvaComputer = std::make_shared<LeptonMvaHelper>(iConfig);
};


void LeptonAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nL",                           &_nL,                           "_nL/b");
  outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
  outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
  outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");
  outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/b");
  outputTree->Branch("_nVFit",                        &_nVFit,                        "_nVFit/b");
  outputTree->Branch("_nGoodLeading",                 &_nGoodLeading,                 "_nGoodLeading/b");
  outputTree->Branch("_lIndex",                       &_lIndex,                       "_lIndex[_nL]/D");
  outputTree->Branch("_vertices",                     &_vertices,                     "_vertices[12][_nVFit]/D");
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
  outputTree->Branch("_lElectronMva",                 &_lElectronMva,                 "_lElectronMva[_nLight]/F");
  outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
  outputTree->Branch("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", &_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, "_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL]/O");
  outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
  outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
  outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
  outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
  outputTree->Branch("_lpassConversionVeto",          &_lpassConversionVeto,          "_lpassConversionVeto[_nL]/O");
  outputTree->Branch("_muNumberInnerHits",            &_muNumberInnerHits,            "_muNumberInnerHits[_nL]/D");
  outputTree->Branch("_eleNumberInnerHitsMissing",    &_eleNumberInnerHitsMissing,    "_eleNumberInnerHitsMissing[_nL]/D");
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
  outputTree->Branch("_ecalIso",                      &_ecalIso ,                     " _ecalIso[_nL]/D");
  outputTree->Branch("_hcalIso",                      &_hcalIso ,                     " _hcalIso[_nL]/D");
  outputTree->Branch("_deltaBIso",                    &_deltaBIso,                    " _deltaBIso[_nL]/D");
  outputTree->Branch("_ecalPFClusterIso",             &_ecalPFClusterIso ,            " _ecalPFClusterIso[_nL]/D");
  outputTree->Branch("_hcalPFClusterIso",             &_hcalPFClusterIso ,            " _hcalPFClusterIso[_nL]/D");
  outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
  outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
  outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
  outputTree->Branch("_muonSegComp",                  &_muonSegComp,                  "_muonSegComp[_nMu]/D");
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
  outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
  outputTree->Branch("_muonSegComp",                  &_muonSegComp,                  "_muonSegComp[_nMu]/D");
  if(!multilepAnalyzer->isData){
    outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
    outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
  }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const reco::Vertex& primaryVertex){
  edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
  edm::Handle<edm::ValueMap<float>> electronsMva;                  iEvent.getByToken(multilepAnalyzer->eleMvaToken,                       electronsMva);
  //edm::Handle<edm::ValueMap<float>> electronsMvaHZZ;               iEvent.getByToken(multilepAnalyzer->eleMvaHZZToken,                    electronsMvaHZZ);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedVeto;          iEvent.getByToken(multilepAnalyzer->eleCutBasedVetoToken,              electronsCutBasedVeto);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedLoose;         iEvent.getByToken(multilepAnalyzer->eleCutBasedLooseToken,             electronsCutBasedLoose);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedMedium;        iEvent.getByToken(multilepAnalyzer->eleCutBasedMediumToken,            electronsCutBasedMedium);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedTight;         iEvent.getByToken(multilepAnalyzer->eleCutBasedTightToken,             electronsCutBasedTight);
  edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
  edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(multilepAnalyzer->tauToken,                          taus);
  edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
  edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);
  edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);

  _nL     = 0;
  _nLight = 0;
  _nMu    = 0;
  _nEle   = 0;
  _nTau   = 0;
  _nVFit  = 0;
  _nGoodLeading = 0;
  
  //bool good_leading=false; // to check 1 leading-well_isolated lepton
  int counter_index_leptons   = 0;
  int counter_number_vertices = 0;
  
  //loop over muons
  for(const pat::Muon& mu : *muons){
   // if(mu.innerTrack().isNull()) continue;
    if(!mu.isPFMuon()) continue;
   // if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue; // loose POG muon
    if(mu.pt() < 3)              continue;                   // from 5 to 3 GeV
    if(fabs(mu.eta()) > 2.4)     continue;
    //if(!mu.isMediumMuon()) continue;
    
    counter_index_leptons++  ;                               // unique index to identify the 2 tracks for each vertex
    _lIndex[_nL] = counter_index_leptons;
    
    _lPFMuon[_nL]=  mu.isPFMuon();
    // ===>  if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue; // loose POG muon
    fillLeptonImpactParameters(mu, primaryVertex);

    // if(fabs(_dxy[_nL]) > 0.05) continue;                   // no impact parameter cuts
    // if(fabs(_dz[_nL]) > 0.1) continue;                     // no impact parameter cuts
    fillLeptonKinVars(mu);
    fillLeptonGenVars(mu.genParticle());
    fillLeptonIsoVars(mu, *rho);	
    fillLeptonJetVariables(mu, jets, primaryVertex);

    _lFlavor[_nL]        = 1;
    _muonSegComp[_nL]    = mu.segmentCompatibility();
    _relIso[_nL]         = getRelIso03(mu, *rho);                                               // Isolation variables
  
     // TODO: this is a possible solution to the missing trackRef, but maybe not what you want 
    _muNumberInnerHits[_nL]= (!mu.globalTrack().isNull()) ?   mu.globalTrack()->hitPattern().numberOfValidMuonHits() : (!mu.outerTrack().isNull() ? mu.outerTrack()->hitPattern().numberOfValidMuonHits() : -1);
    _lPOGVeto[_nL]   = mu.isLooseMuon();
    _lPOGLoose[_nL]  = mu.isLooseMuon();
    if ( mu.isLooseMuon()) _lPOGMedium[_nL] = mu.isMediumMuon();
    if ( mu.isLooseMuon()) _lPOGTight[_nL]  = mu.isTightMuon(primaryVertex);

    _eleNumberInnerHitsMissing[_nL] =-1;
    _lpassConversionVeto[_nL] = false;
    _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL] = false;
    
    if (mu.pt() > 20 && fabs(_dxy[_nL]) < 0.05 && fabs(_dz[_nL])< 0.1 && getRelIso04(mu, *rho) < 0.3 && !mu.innerTrack().isNull() && (mu.isTrackerMuon() || mu.isGlobalMuon()) ) ++_nGoodLeading;

    ++_nMu;
    ++_nL;
    ++_nLight;
  }

  // Loop over electrons (note: using iterator we can easily get the ref too)
  for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
    auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele - electrons->begin()));
    if(ele->gsfTrack().isNull()) continue; 
    if(ele->pt() < 5)           continue; // from 10 to 6
    if(fabs(ele->eta()) > 2.5)   continue;
    // ---->  loose requirements about number of hits and VetoConversion
        //if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2) continue;
        //if(!ele->passConversionVeto()) continue;  
    if(!isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele)) continue; // check the loooong name
    if(eleMuOverlap(*ele, _lPFMuon))  continue;              // overlap muon-electron deltaR  ==  0.05
    
    counter_index_leptons++  ;                               // unique index to identify the 2 tracks for each vertex
    _lIndex[_nL] = counter_index_leptons;
    
    
    _eleNumberInnerHitsMissing[_nL]=ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    _muNumberInnerHits[_nL] =-1;
    fillLeptonImpactParameters(*ele, primaryVertex);
    fillLeptonIsoVars(*ele, *rho);
	  
    _lpassConversionVeto[_nL] = ele->passConversionVeto();


    //if(fabs(_dxy[_nL]) > 0.05) continue;                   // no impact parameter cuts
    // if(fabs(_dz[_nL]) > 0.1) continue;                   // no impact parameter cuts
    fillLeptonKinVars(*ele);
    fillLeptonGenVars(ele->genParticle());
    fillLeptonJetVariables(*ele, jets, primaryVertex);
    _lFlavor[_nL]      = 0;
    _lEtaSC[_nL]       = ele->superCluster()->eta();
    _relIso[_nL]       = getRelIso03(*ele, *rho);
    

    _lElectronMva[_nL] = (*electronsMva)[electronRef];
    _lElectronPassEmu[_nL] = passTriggerEmulationDoubleEG(&*ele);
   
   
    _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL] = isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele);
    _lPOGVeto[_nL]     = (*electronsCutBasedVeto)[electronRef];
    _lPOGLoose[_nL]    = (*electronsCutBasedLoose)[electronRef]; // cut-based cuts
    _lPOGMedium[_nL]   = (*electronsCutBasedMedium)[electronRef];
    _lPOGTight[_nL]    = (*electronsCutBasedTight)[electronRef];             // Actually in SUS-17-001 we applied addtionaly lostHists==0, probably not a big impact

    
   if (ele->pt() > 20 && fabs(_dxy[_nL]) < 0.05 && fabs(_dz[_nL])< 0.1 && _relIso[_nL] < 0.3 && !ele->gsfTrack().isNull() && _eleNumberInnerHitsMissing[_nL] <=2 && ele->passConversionVeto()) ++_nGoodLeading;

    
    ++_nEle;
    ++_nL;
    ++_nLight;
  }

   //loop over taus
  for(const pat::Tau& tau : *taus){
    if(tau.pt() < 20)         continue;          // Minimum pt for tau reconstruction
    if(fabs(tau.eta()) > 2.3) continue;
    if(!tau.tauID("decayModeFinding")) continue;
    fillLeptonKinVars(tau);
    fillLeptonGenVars(tau.genParticle());
    fillLeptonImpactParameters(tau, primaryVertex);
    if(_dz[_nL] < 0.4)        continue;         //tau dz cut used in ewkino

    _lFlavor[_nL]  = 2;
    _tauMuonVeto[_nL] = tau.tauID("againstMuonLoose3");                                        //Light lepton vetos
    _tauEleVeto[_nL] = tau.tauID("againstElectronLooseMVA6");

    _lPOGVeto[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");                        //old tau ID
    _lPOGLoose[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    _lPOGMedium[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    _lPOGTight[_nL] = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    _tauVTightMvaOld[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

    _decayModeFindingNew[_nL] = tau.tauID("decayModeFindingNewDMs");                           //new Tau ID 
    _tauVLooseMvaNew[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
    _tauLooseMvaNew[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
    _tauMediumMvaNew[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
    _tauTightMvaNew[_nL] = tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT");
    _tauVTightMvaNew[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT");

   
    ++_nTau;
    ++_nL;
  }
		
  
  /*
   * refitting vertices displaced *********************************************************** 
   */  
  
  float iMu_plus=0;
  float iMu_minus_mu=0;
  float iMu_minus_e=0;
  float iE_plus=_nMu;
  float iE_minus_mu=_nMu;
  float iE_minus_e=_nMu;
  
  for(const pat::Muon& mu_1 : *muons){ // for µ
    if(mu_1.pt() < 3 || fabs(mu_1.eta()) > 2.4 || !mu_1.isPFMuon())              continue;   
	  	  //    if(!mu_1.isLooseMuon()) continue;

    iMu_plus++;
    // +++++++++++++++    µ+
    if (mu_1.charge() < 0) continue;
    const reco::Track&  tk_1 = (!mu_1.innerTrack().isNull()) ? *mu_1.innerTrack () :  *mu_1.outerTrack () ;
    //std::cout<<"µ+ : "<<mu_1.pt() <<"  charge: "<<mu_1.charge()<<"  track pos: ("<<mu_1.outerTrack ()->outerX()<<","<<mu_1.outerTrack ()->outerY()<<","<<mu_1.outerTrack ()->outerZ()<<")"<<std::endl;
  //  std::cout<<"µ+ : "<<mu_1.pt() <<"  charge: "<<mu_1.charge()<<std::endl;

    //------------------  loop µ-
    for(const pat::Muon& mu_2 : *muons){ 
	    	//  	      if(!mu_2.isLooseMuon()) continue;

      if(mu_2.pt() < 3 || fabs(mu_2.eta()) > 2.4 || !mu_2.isPFMuon())              continue;   
      iMu_minus_mu++;
      if (mu_2.charge() > 0) continue;  // only opposite charge

      const reco::Track&  tk_2 = (!mu_2.innerTrack().isNull()) ? *mu_2.innerTrack () :  *mu_2.outerTrack () ;
      //std::cout<<"µ- : "<<mu_2.pt() <<"  charge: "<<mu_2.charge()<<"  track pos: ("<<mu_2.outerTrack ()->outerX()<<","<<mu_2.outerTrack ()->outerY()<<","<<mu_2.outerTrack ()->outerZ()<<")"<<std::endl;
 
     // std::cout<<"µ- : "<<mu_2.pt() <<"  charge: "<<mu_2.charge()<<std::endl;
    
      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(!dilvtx.isValid()) { 
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
      } 
      else {   
	 //std::cout<<"--> indices: "<<iMu_plus*100 + iMu_minus_mu<<" ("<<iMu_plus<<" - "<<iMu_minus_mu<< ")  pos:  "<<dilvtx.position().x()<<" , "<<dilvtx.position().y()<<" , "<<dilvtx.position().z()<<std::endl;

	_vertices[0][_nVFit] = iMu_plus*100 + iMu_minus_mu;                   
	_vertices[1][_nVFit] = dilvtx.position().x(); 
	_vertices[2][_nVFit] = dilvtx.position().y(); 
	_vertices[3][_nVFit] = dilvtx.position().z(); 
	_vertices[4][_nVFit] = dilvtx.positionError().cxx(); 
	_vertices[5][_nVFit] = dilvtx.positionError().cyy(); 
	_vertices[6][_nVFit] = dilvtx.positionError().czz(); 
	_vertices[7][_nVFit] = dilvtx.positionError().cyx(); 
	_vertices[8][_nVFit] = dilvtx.positionError().czy(); 
	_vertices[9][_nVFit] = dilvtx.positionError().czx(); 
	_vertices[10][_nVFit] = dilvtx.degreesOfFreedom(); 
	_vertices[11][_nVFit] = dilvtx.totalChiSquared(); 
	++_nVFit;   
      } 
    }// end loop µ-
          
    //------------------  loop e-
    for(auto ele_2 = electrons->begin(); ele_2 != electrons->end(); ++ele_2){
      auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele_2 - electrons->begin()));
      if(ele_2->gsfTrack().isNull()) continue;
      if(ele_2->pt() < 5 || fabs(ele_2->eta()) > 2.5 || !isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele_2) || eleMuOverlap(*ele_2, _lPFMuon) )           continue; // from 10 to 6
      iE_minus_mu++; // it is already _nMu
      if(ele_2->charge() > 0) continue; // only opposite charge

      const reco::Track&  tk_2 =  *ele_2->gsfTrack() ;
      //std::cout<<"e- : "<<ele_2->pt() <<"  charge: "<<ele_2->charge()<<"  track pos: ("<<ele_2->gsfTrack ()->outerX()<<","<<ele_2->gsfTrack ()->outerY()<<","<<ele_2->gsfTrack ()->outerZ()<<")"<<std::endl;
      //std::cout<<"e- : "<<ele_2->pt() <<"  charge: "<<ele_2->charge()<<std::endl;

      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(!dilvtx.isValid()) { 
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
      } 
      else {      
	 //std::cout<<"--> indices: "<<iMu_plus*100 + iE_minus_mu<<" ("<<iMu_plus<<" - "<<iE_minus_mu<< ")  pos:  "<<dilvtx.position().x()<<" , "<<dilvtx.position().y()<<" , "<<dilvtx.position().z()<<std::endl;

	_vertices[0][_nVFit] = iMu_plus*100 + iE_minus_mu;          
	_vertices[1][_nVFit] = dilvtx.position().x(); 
	_vertices[2][_nVFit] = dilvtx.position().y(); 
	_vertices[3][_nVFit] = dilvtx.position().z(); 
	_vertices[4][_nVFit] = dilvtx.positionError().cxx(); 
	_vertices[5][_nVFit] = dilvtx.positionError().cyy(); 
	_vertices[6][_nVFit] = dilvtx.positionError().czz(); 
	_vertices[7][_nVFit] = dilvtx.positionError().cyx(); 
	_vertices[8][_nVFit] = dilvtx.positionError().czy(); 
	_vertices[9][_nVFit] = dilvtx.positionError().czx(); 
	_vertices[10][_nVFit] = dilvtx.degreesOfFreedom(); 
	_vertices[11][_nVFit] = dilvtx.totalChiSquared(); 
	std::cout<<"---> "<<     _vertices[2][_nVFit]<<std::endl;
	      
	      
	      
	++_nVFit;   
      } 
    }// end loop e-
  }//end loop µ
  
  
  iMu_minus_e=0;
  iE_plus=_nMu;
  iE_minus_e=_nMu;
  
      
  for(auto ele_1 = electrons->begin(); ele_1 != electrons->end(); ++ele_1){ // for electrons
    auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele_1 - electrons->begin()));
    if(ele_1->gsfTrack().isNull()) continue;
    if(ele_1->pt() < 5 || fabs(ele_1->eta()) > 2.5 || !isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele_1) || eleMuOverlap(*ele_1, _lPFMuon) )           continue; // from 10 to 6
    iE_plus++;
    //+++++++++++++++++++++ e+
    if(ele_1->charge() < 0) continue; 
    
    const reco::Track&  tk_1 =  *ele_1->gsfTrack() ;
    //std::cout<<"****  e+ : "<<ele_1->pt() <<"  charge: "<<ele_1->charge()<<std::endl;

    //------------------  loop µ+
    for(const pat::Muon& mu_2 : *muons){ 
	    	  	//      if(!mu_2.isLooseMuon()) continue;

      if(mu_2.pt() < 3 || fabs(mu_2.eta()) > 2.4 || !mu_2.isPFMuon())              continue;   
      iMu_minus_e++;
      if (mu_2.charge() > 0) continue;  // only opposite charge
        
      const reco::Track&  tk_2 = (!mu_2.innerTrack().isNull()) ? *mu_2.innerTrack () :  *mu_2.outerTrack () ; 
	//std::cout<<"****  µ- : "<<mu_2.pt() <<"  charge: "<<mu_2.charge()<<std::endl;

      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(!dilvtx.isValid()) { 
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
      } 
      else {       
	 //std::cout<<"--> indices: "<<iE_plus*100 + iMu_minus_e<<" ("<<iE_plus<<" - "<<iMu_minus_e<< ")  pos:  "<<dilvtx.position().x()<<" , "<<dilvtx.position().y()<<" , "<<dilvtx.position().z()<<std::endl;

	_vertices[0][_nVFit] = iE_plus*100+iMu_minus_e;    
	_vertices[1][_nVFit] = dilvtx.position().x(); 
	_vertices[2][_nVFit] = dilvtx.position().y(); 
	_vertices[3][_nVFit] = dilvtx.position().z(); 
	_vertices[4][_nVFit] = dilvtx.positionError().cxx(); 
	_vertices[5][_nVFit] = dilvtx.positionError().cyy(); 
	_vertices[6][_nVFit] = dilvtx.positionError().czz(); 
	_vertices[7][_nVFit] = dilvtx.positionError().cyx(); 
	_vertices[8][_nVFit] = dilvtx.positionError().czy(); 
	_vertices[9][_nVFit] = dilvtx.positionError().czx(); 
	_vertices[10][_nVFit] = dilvtx.degreesOfFreedom(); 
	_vertices[11][_nVFit] = dilvtx.totalChiSquared(); 
	++_nVFit;   
      } 
    }// end loop µ-
    
    
    //------------------  loop e+
    for(auto ele_2 = electrons->begin(); ele_2 != electrons->end(); ++ele_2){
      auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele_2 - electrons->begin()));
      if(ele_2->gsfTrack().isNull()) continue;
      if(ele_2->pt() < 5 || fabs(ele_2->eta()) > 2.5 || !isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele_2) || eleMuOverlap(*ele_2, _lPFMuon) )           continue; // from 10 to 6
      iE_minus_e++;
        
      if(ele_2->charge() > 0) continue; // only opposite charge
        
      const reco::Track&  tk_2 =  *ele_2->gsfTrack();  
      //std::cout<<"****  e- : "<<ele_2->pt() <<"  charge: "<<ele_2->charge()<<std::endl;

      TransientVertex dilvtx = dileptonVertex(tk_1, tk_2);
      if(!dilvtx.isValid()) { 
	std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl; 
      } 
      else {    
	      //	 std::cout<<"--> indices: "<<iE_plus*100 + iE_minus_e<<" ("<<iE_plus<<" - "<<iE_minus_e<< ")  pos:  "<<dilvtx.position().x()<<" , "<<dilvtx.position().y()<<" , "<<dilvtx.position().z()<<std::endl;

	_vertices[0][_nVFit] = iE_plus*100 + iE_minus_e;          
	_vertices[1][_nVFit] = dilvtx.position().x(); 
	_vertices[2][_nVFit] = dilvtx.position().y(); 
	_vertices[3][_nVFit] = dilvtx.position().z(); 
	_vertices[4][_nVFit] = dilvtx.positionError().cxx(); 
	_vertices[5][_nVFit] = dilvtx.positionError().cyy(); 
	_vertices[6][_nVFit] = dilvtx.positionError().czz(); 
	_vertices[7][_nVFit] = dilvtx.positionError().cyx(); 
	_vertices[8][_nVFit] = dilvtx.positionError().czy(); 
	_vertices[9][_nVFit] = dilvtx.positionError().czx(); 
	_vertices[10][_nVFit] = dilvtx.degreesOfFreedom(); 
	_vertices[11][_nVFit] = dilvtx.totalChiSquared(); 
	++_nVFit;   
      }
    }// end loop e+
    
    
  }//end electrons
   for (int i =0; i < 12; i++){ 
	
	
	std::cout<<i<< ")  ***** "<<    _vertices[i][0]<<std::endl;
	   
   }
  if(multilepAnalyzer->skim == "trilep"    and (_nLight < 3  || _nGoodLeading < 1)  ) return false;

  //if(multilepAnalyzer->skim == "trilep"    and (_nLight     < 3   ||   !good_leading)) return false;
  if(multilepAnalyzer->skim == "dilep"     and _nLight < 2) return false;
  if(multilepAnalyzer->skim == "ttg"       and _nLight < 2) return false;
  if(multilepAnalyzer->skim == "singlelep" and _nLight < 1) return false;
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

void LeptonAnalyzer::fillLeptonGenVars(const reco::GenParticle* genParticle){
  if(genParticle != nullptr){
    _lIsPrompt[_nL] = (genParticle)->isPromptFinalState();
    _lMatchPdgId[_nL] = (genParticle)->pdgId();
  } else{
    _lIsPrompt[_nL] = false;
    _lMatchPdgId[_nL] = 0;
  }
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
bool LeptonAnalyzer::eleMuOverlap(const pat::Electron& ele, const bool* loose){
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
















void LeptonAnalyzer::fillLeptonJetVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::Jet>>& jets, const reco::Vertex& vertex){
    //Make skimmed "close jet" collection
    std::vector<pat::Jet> selectedJetsAll;
    for(auto jet = jets->cbegin(); jet != jets->cend(); ++jet){
        if( jet->pt() > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
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
        auto  l1Jet       = jet.correctedP4("L1FastJet");
        float JEC         = jet.p4().E()/l1Jet.E();
        auto  l           = lepton.p4();
        auto  lepAwareJet = (l1Jet - l)*JEC + l;
        TLorentzVector lV = TLorentzVector(l.px(), l.py(), l.pz(), l.E());
        TLorentzVector jV = TLorentzVector(lepAwareJet.px(), lepAwareJet.py(), lepAwareJet.pz(), lepAwareJet.E());
        _ptRatio[_nL]       = l.pt()/lepAwareJet.pt();
        _ptRel[_nL]         = lV.Perp((jV - lV).Vect());
        _closestJetCsvV2[_nL] = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        _closestJetDeepCsv_b[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        //compute selected track multiplicity of closest jet
        _selectedTrackMult[_nL] = 0;
        for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
            const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
            try {                                                                                                     // In principle, from CMSSW_9_X you need to use if(daughter->hasTrackDetails()){ here, bus that function does not exist in CMSSW_8_X
                const reco::Track& daughterTrack = daughter->pseudoTrack();                                             // Using try {} catch (...){} the code compiles in both versions
                TLorentzVector trackVec          = TLorentzVector(daughterTrack.px(), daughterTrack.py(), daughterTrack.pz(), daughterTrack.p());
                double daughterDeltaR            = trackVec.DeltaR(jV);
                bool goodTrack                   = daughterTrack.pt() > 1 && daughterTrack.charge() != 0 && daughterTrack.hitPattern().numberOfValidHits() > 7
                                                   && daughterTrack.hitPattern().numberOfValidPixelHits() > 1 && daughterTrack.normalizedChi2() < 5 && fabs(daughterTrack.dz(vertex.position())) < 17
                                                   && fabs(daughterTrack.dxy(vertex.position())) < 17;
                if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack) ++_selectedTrackMult[_nL];
            } catch (...){}
        }
    }
  }


