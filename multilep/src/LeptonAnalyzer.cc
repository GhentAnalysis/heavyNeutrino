#include "../interface/LeptonAnalyzer.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"

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
  outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
  outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
  outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/I");
  outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
  outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
  outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
  outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
  outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
  outputTree->Branch("_lElectronMva",                 &_lElectronMva,                 "_lElectronMva[_nLight]/F");
  outputTree->Branch("_lHNLoose",                     &_lHNLoose,                     "_lHNLoose[_nL]/O");
  outputTree->Branch("_lHNFO",                        &_lHNFO,                        "_lHNFO[_nL]/O");
  outputTree->Branch("_lHNTight",                     &_lHNTight,                     "_lHNTight[_nL]/O");

  outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
  outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
}

void LeptonAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
  edm::Handle<edm::ValueMap<float>> electronsMva;                  iEvent.getByToken(multilepAnalyzer->eleMvaToken,                       electronsMva);
  edm::Handle<edm::ValueMap<float>> electronsMvaHZZ;               iEvent.getByToken(multilepAnalyzer->eleMvaHZZToken,                    electronsMvaHZZ);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedTight;         iEvent.getByToken(multilepAnalyzer->eleCutBasedTightToken,             electronsCutBasedTight);
  edm::Handle<edm::ValueMap<bool>> electronsCutBasedMedium;        iEvent.getByToken(multilepAnalyzer->eleCutBasedMediumToken,            electronsCutBasedMedium);
  edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
  edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(multilepAnalyzer->tauToken,                          taus);
  edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
  edm::Handle<double> rhoJets;                                     iEvent.getByToken(multilepAnalyzer->rhoToken,                          rhoJets);

  _nL     = 0;
  _nLight = 0;
  _nMu    = 0;
  _nEle   = 0;
  _nTau   = 0;

  //loop over muons
  for(const pat::Muon& mu : *muons){
    if(_nL == nL_max)            continue;
    if(mu.innerTrack().isNull()) continue;
    if(mu.pt() < 5)              continue;
    if(fabs(mu.eta()) > 2.4)     continue;
    fillLeptonKinVars(mu);
    fillLeptonGenVars(mu.genParticle());
    _lFlavor[_nL] = 1;
    //Vertex variables // better move this too in fillLeptonIdVars
    _dxy[_nL]     = mu.innerTrack()->dxy(); // Still need to provide PV !!!!!!!!!!!!
    _dz[_nL]      = mu.innerTrack()->dz();
    _3dIP[_nL]    = mu.dB(pat::Muon::PV3D);
    _3dIPSig[_nL] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);
    //Isolation variables
    _relIso[_nL]  = getRelIso03(mu, *rhoJets);
    _miniIso[_nL] = getMiniIso(mu, *packedCands, 0.2, *rhoJets);

    _lHNLoose[_nL] = isHNLoose(mu);
    _lHNFO[_nL]    = isHNFO(mu);    // don't change order, they rely on above variables
    _lHNTight[_nL] = isHNTight(mu);

    ++_nMu;
    ++_nL;
    ++_nLight;
  }

  // Loop over electrons (note: using iterator we can easily get the ref too)
  for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
    auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele - electrons->begin()));
    if(_nL == nL_max)            continue;
    if(ele->gsfTrack().isNull()) continue;
    if(ele->pt() < 10)           continue;
    if(fabs(ele->eta()) > 2.5)   continue;
    fillLeptonKinVars(*ele);
    fillLeptonGenVars(ele->genParticle());
    _lFlavor[_nL]  = 0;
    //Vertex varitables
    _dxy[_nL]           = ele->gsfTrack()->dxy();
    _dz[_nL]            = ele->gsfTrack()->dz();
    _3dIP[_nL]          = ele->dB(pat::Electron::PV3D);
    _3dIPSig[_nL]       = ele->dB(pat::Electron::PV3D)/ele->edB(pat::Electron::PV3D);

    //isolation variables
    _relIso[_nL]        = getRelIso03(*ele, *rhoJets);
    _miniIso[_nL]       = getMiniIso(*ele, *packedCands, 0.2, *rhoJets);

  //float_mvaValue_HZZ = (*electronsMva)[electronRef];
    _lElectronMva[_nL] = (*electronsMvaHZZ)[electronRef];
    _lHNLoose[_nL]     = isHNLoose(*ele);
    _lHNFO[_nL]        = isHNFO(*ele);
    _lHNTight[_nL]     = isHNTight(*ele);

    ++_nEle;
    ++_nL;
    ++_nLight;
  }

  //loop over taus
  for(const pat::Tau& tau : *taus){
    if(_nL == nL_max)         continue;
    if(tau.pt() < 20)         continue; //investigate up to what Pt threshold taus can be properly reconstructed
    if(fabs(tau.eta()) > 2.3) continue;
    fillLeptonKinVars(tau);
    _lFlavor[_nL]  = 2;
    _dxy[_nL]      = tau.dxy();
  //_dz[_nL]       = tau.dz();
    _dz[_nL]       = 0;
    _3dIP[_nL]     = tau.ip3d();
    _3dIPSig[_nL]  = tau.ip3d_Sig();
    _lHNLoose[_nL] = false; // TO BE IMPLEMENTED
    _lHNFO[_nL]    = false;
    _lHNTight[_nL] = false;
    ++_nTau;
    ++_nL;
  }
}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton){
    //kinematics
    _lPt[_nL]     = lepton.pt();
    _lEta[_nL]    = lepton.eta();
    _lPhi[_nL]    = lepton.phi();
    _lE[_nL]      = lepton.energy();
    _lCharge[_nL] = lepton.charge();
}

void LeptonAnalyzer::fillLeptonGenVars(const reco::GenParticle* genParticle){
    if(genParticle != nullptr) _isPrompt[_nL] = (genParticle)->isPromptFinalState();
    else                       _isPrompt[_nL] = false;
}


//Check if electron overlaps with loose muon
//There's a better way to do this
bool LeptonAnalyzer::eleMuOverlap(const pat::Electron& ele){
  TLorentzVector eleV, muV;
  eleV.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy());
  for(unsigned m = 0; m < _nMu; ++m){
    if(_lHNLoose[m]){
      TLorentzVector muV(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
      if(eleV.DeltaR(muV) < 0.05) return true;
    }
  }
  return false;	
}
