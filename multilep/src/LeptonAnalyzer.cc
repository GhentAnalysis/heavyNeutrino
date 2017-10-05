#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TLorentzVector.h"
/*
 * Calculating all lepton-related variables
 * I know, code is still messy here, but it will improve
 */


LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  electronsEffectiveAreas((iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas")).fullPath()),
  muonsEffectiveAreas(    (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath()),
  multilepAnalyzer(multilepAnalyzer)
{};


void LeptonAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nL",                           &_nL,                           "_nL/b");
  outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
  outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
  outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");
  outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/b");

  outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
  outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
  outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
  outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
  outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
  outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/I");
  outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
  outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
  outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
  outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
  outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
  outputTree->Branch("_lElectronMva",                 &_lElectronMva,                 "_lElectronMva[_nLight]/F");
  outputTree->Branch("_lElectronPassEmu",	          &_lElectronPassEmu,	          "_lElectronPassEmu[_nLight]/O");
  outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
  outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
  outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
  outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
  outputTree->Branch("_leptonMva",                    &_leptonMva,                    "_leptonMva[_nL]/O");

  outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
  outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
  outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
  outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
  outputTree->Branch("_closestJetCsv",                &_closestJetCsv,                "_closestJetCsv[_nLight]/D");
  outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/D");
  if(!multilepAnalyzer->isData){
      outputTree->Branch("_lIsPrompt",                &_lIsPrompt,                    "_lIsPrompt[_nL]/D");
      outputTree->Branch("_lMatchPdgId",              &_lMatchPdgId,                  "_lMatchPdgId[_nL]/D");
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
  edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoTokenAll,                       rho);
  edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);

  _nL     = 0;
  _nLight = 0;
  _nMu    = 0;
  _nEle   = 0;
  _nTau   = 0;

  //loop over muons
  for(const pat::Muon& mu : *muons){
    if(_nL == nL_max)            continue;
    if(mu.innerTrack().isNull()) continue;
    if(mu.pt() < 3)              continue;                   // from 5 to 3 GeV
    if(fabs(mu.eta()) > 2.4)     continue;
    if(!mu.isPFMuon()) continue;
    if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue; // loose POG muon
    fillLeptonImpactParameters(mu, primaryVertex);
   // if(fabs(_dxy[_nL]) > 0.05) continue;                   // no impact parameter cuts
   // if(fabs(_dz[_nL]) > 0.1) continue;                     // no impact parameter cuts
    fillLeptonKinVars(mu);
    fillLeptonGenVars(mu.genParticle());
    fillLeptonJetVariables(mu, jets, primaryVertex);
    _lFlavor[_nL] = 1;
    //Isolation variables
    _relIso[_nL]  = getRelIso03(mu, *rho);
    _miniIso[_nL] = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho);
    _miniIsoCharged[_nL] = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho, true);

    _lPOGVeto[_nL]   = mu.isLooseMuon();
    _lPOGLoose[_nL]  = mu.isLooseMuon();
    _lPOGMedium[_nL] = mu.isMediumMuon();
    _lPOGTight[_nL]  = mu.isTightMuon(primaryVertex);
    _leptonMva[_nL]  = leptonMvaVal(mu);

    ++_nMu;
    ++_nL;
    ++_nLight;
  }

  // Loop over electrons (note: using iterator we can easily get the ref too)
  for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
    auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele - electrons->begin()));
    if(_nL == nL_max)            break;
    if(ele->gsfTrack().isNull()) continue;
    if(ele->pt() < 5)           continue; // from 10 to 5
    if(fabs(ele->eta()) > 2.5)   continue;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2) continue;
    if(!ele->passConversionVeto()) continue;
    fillLeptonImpactParameters(*ele, primaryVertex);
    //if(fabs(_dxy[_nL]) > 0.05) continue;                   // no impact parameter cuts
    // if(fabs(_dz[_nL]) > 0.1) continue;                   // no impact parameter cuts
    fillLeptonKinVars(*ele);
    fillLeptonGenVars(ele->genParticle());
    fillLeptonJetVariables(*ele, jets, primaryVertex);
    _lFlavor[_nL]      = 0;
    _lEtaSC[_nL]       = ele->superCluster()->eta();
    _relIso[_nL]       = getRelIso03(*ele, *rho);
    _miniIso[_nL]      = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho);
    _miniIsoCharged[_nL] = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho, true);

    _lElectronMva[_nL] = (*electronsMva)[electronRef];
    _lElectronPassEmu[_nL] = passTriggerEmulationDoubleEG(&*ele);
   
   
    _lPOGVeto[_nL]     = (*electronsCutBasedVeto)[electronRef];
    _lPOGLoose[_nL]    = (*electronsCutBasedLoose)[electronRef]; // cut-based cuts
    _lPOGMedium[_nL]   = (*electronsCutBasedMedium)[electronRef];
    _lPOGTight[_nL]    = (*electronsCutBasedTight)[electronRef];             // Actually in SUS-17-001 we applied addtionaly lostHists==0, probably not a big impact
    _leptonMva[_nL]    = leptonMvaVal(*ele);

    ++_nEle;
    ++_nL;
    ++_nLight;
  }

  //loop over taus
  for(const pat::Tau& tau : *taus){
    if(_nL == nL_max)         continue;
    if(tau.pt() < 20)         continue;          // Minimum pt for tau reconstruction
    if(fabs(tau.eta()) > 2.3) continue;
    fillLeptonKinVars(tau);
    fillLeptonGenVars(tau.genParticle());
    fillLeptonImpactParameters(tau, primaryVertex);
    _lFlavor[_nL]  = 2;
    _lHNLoose[_nL] = false;                      // TO BE IMPLEMENTED
    _lHNFO[_nL]    = false;
    _lHNTight[_nL] = false;
    ++_nTau;
    ++_nL;
  }

  if(multilepAnalyzer->skim == "trilep"    and _nLight < 3) return false;
  if(multilepAnalyzer->skim == "dilep"     and _nLight < 2) return false;
  if(multilepAnalyzer->skim == "ttg"       and _nLight < 2) return false;
  if(multilepAnalyzer->skim == "singlelep" and _nLight < 1) return false;
  return true;
}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton){
  _lPt[_nL]     = lepton.pt();
  _lEta[_nL]    = lepton.eta();
  _lPhi[_nL]    = lepton.phi();
  _lE[_nL]      = lepton.energy();
  _lCharge[_nL] = lepton.charge();
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
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon, const reco::Vertex& vertex){
  _dz[_nL]      = muon.innerTrack()->dz(vertex.position());
  _3dIP[_nL]    = muon.dB(pat::Muon::PV3D);
  _3dIPSig[_nL] = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D);
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Tau& tau, const reco::Vertex& vertex){
  _dxy[_nL]     = (double) tau.dxy();                                      // warning: float while dxy of tracks are double; could also return -1000
  _dz[_nL]      =  tau_dz(tau, vertex.position());
  _3dIP[_nL]    = tau.ip3d();
  _3dIPSig[_nL] = tau.ip3d_Sig(); 
}

//Function returning tau dz
double LeptonAnalyzer::tau_dz(const pat::Tau& tau, const reco::Vertex::Point& vertex){
    const reco::Candidate::Point& tauVtx = tau.leadChargedHadrCand()->vertex();
    return (tauVtx.Z() - vertex.z()) - ((tauVtx.X() - vertex.x())*tau.px()+(tauVtx.Y()-vertex.y())*tau.py())/tau.pt()*tau.pz()/tau.pt();
}

//Check if electron overlaps with loose muon
bool LeptonAnalyzer::eleMuOverlap(const pat::Electron& ele){
  TLorentzVector eleV;
  eleV.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy());
  for(unsigned m = 0; m < _nMu; ++m){
    if(_lHNLoose[m]){
      TLorentzVector muV;
      muV.SetPtEtaPhiE(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
      if(eleV.DeltaR(muV) < 0.05) return true;
    }
  }
  return false;	
}


void LeptonAnalyzer::fillLeptonJetVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::Jet>>& jets, const reco::Vertex& vertex){
  //Make skimmed "close jet" collection
  std::vector<pat::Jet> selectedJetsAll;
  for(auto jet = jets->cbegin(); jet != jets->cend(); ++jet){
      if( jet->pt() > 5 && fabs( jet->eta() ) < 3){ 
          selectedJetsAll.push_back(*jet);
      }
  }
  // Find closest selected jet
  unsigned closestIndex = 0;
  for(unsigned j = 1; j < selectedJetsAll.size(); ++j){
      if(reco::deltaR(selectedJetsAll[j], lepton) < reco::deltaR(selectedJetsAll[closestIndex], lepton)){
          closestIndex = j;
       }
   }
   const pat::Jet& jet = selectedJetsAll[closestIndex];
   if(reco::deltaR(jet, lepton) > 0.4){
       _ptRatio[_nL] = 1;
       _ptRel[_nL]   = 0;
   } else {
       //WARNING, these jets currently remain uncorrected!!
       auto  l1Jet       = jet.correctedP4("L1FastJet");
       float JEC         = jet.p4().E()/l1Jet.E();
       auto  l           = lepton.p4();
       auto  lepAwareJet = (l1Jet - l)*JEC + l;

       TLorentzVector lV = TLorentzVector(l.px(), l.py(), l.pz(), l.E());
       TLorentzVector jV = TLorentzVector(lepAwareJet.px(), lepAwareJet.py(), lepAwareJet.pz(), lepAwareJet.E());

       _ptRatio[_nL] = l.pt()/lepAwareJet.pt();
       _ptRel[_nL]   = lV.Perp((jV - lV).Vect());
       _closestJetCsv[_nL] = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
       //compute selected track multiplicity of closest jet
       _selectedTrackMult[_nL] = 0;
       for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
           const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
           const reco::Track& daughterTrack = daughter->pseudoTrack();
           TLorentzVector trackVec = TLorentzVector(daughterTrack.px(), daughterTrack.py(), daughterTrack.pz(), daughterTrack.p());
           double daughterDeltaR = trackVec.DeltaR(jV);
           bool goodTrack = daughterTrack.pt() > 1 && daughterTrack.charge() != 0 && daughterTrack.hitPattern().numberOfValidHits() > 7 && daughterTrack.hitPattern().numberOfValidPixelHits() > 1 && daughterTrack.normalizedChi2() < 5 && fabs(daughterTrack.dz(vertex.position())) < 17 && fabs(daughterTrack.dxy(vertex.position())) < 17; 
           if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack){
               ++_selectedTrackMult[_nL];
           } 
       }
   }
}

