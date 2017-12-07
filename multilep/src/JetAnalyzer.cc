#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  jecUnc((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath()),
  multilepAnalyzer(multilepAnalyzer)
{};

// Note that here only the uncertainty is saved, the Up and Down variations still need to be calculated later, probably easier to do on python level
// Also b-tagging up/down is easier on python level (see https://github.com/GhentAnalysis/StopsDilepton/blob/leptonSelectionUpdate_80X/tools/python/btagEfficiency.py)
// Storing here jet id variables, id itself also to be implemented at python level https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016, or maybe still here as the loose wp does not change often?
// WARNING: the _nJets is number of stored jets (i.e. including those where JECUp/JERUp passes the cut), do not use as selection
void JetAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/b");
  outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
  outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
  outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
  outputTree->Branch("_jetPt_JERUp",               &_jetPt_JERUp,              "_jetPt_JERUp[_nJets]/D");
  outputTree->Branch("_jetPt_JERDown",             &_jetPt_JERDown,            "_jetPt_JERDown[_nJets]/D");
  outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
  outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
  outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
  outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
  outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
  outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
  outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
  outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
//outputTree->Branch("_jetDeepCsv_cc",             &_jetDeepCsv_cc,            "_jetDeepCsv_cc[_nJets]/D");
  outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
  outputTree->Branch("_jetId",                     &_jetId,                    "_jetId[_nJets]/i");
  outputTree->Branch("_nJetswithDaughters",	   &_nJetswithDaughters,       "_nJetswithDaughters/b");
  outputTree->Branch("_nDaughters",		   &_nDaughters,	       "_nDaughters/I");
  outputTree->Branch("_jet_tag_for_daughters",	   &_jet_tag_for_daughters,    "_jet_tag_for_daughters[_nDaughters]/I");
  outputTree->Branch("_jet_daughter_pdgid",	   &_jet_daughter_pdgid,       "_jet_daughter_pdgid[_nDaughters]/I");
  outputTree->Branch("_jet_daughter_pt",	   &_jet_daughter_pt,          "_jet_daughter_pt[_nDaughters]/D");
  outputTree->Branch("_jet_daughter_eta",	   &_jet_daughter_eta,         "_jet_daughter_eta[_nDaughters]/D");
  outputTree->Branch("_jet_daughter_phi",	   &_jet_daughter_phi,         "_jet_daughter_phi[_nDaughters]/D");
  outputTree->Branch("_jet_daughter_energy",	   &_jet_daughter_energy,      "_jet_daughter_energy[_nDaughters]/D");
}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<pat::Jet>> jets;            iEvent.getByToken(multilepAnalyzer->jetToken,            jets);
  edm::Handle<std::vector<pat::Jet>> jetsSmeared;     iEvent.getByToken(multilepAnalyzer->jetSmearedToken,     jetsSmeared);
  edm::Handle<std::vector<pat::Jet>> jetsSmearedUp;   iEvent.getByToken(multilepAnalyzer->jetSmearedUpToken,   jetsSmearedUp);
  edm::Handle<std::vector<pat::Jet>> jetsSmearedDown; iEvent.getByToken(multilepAnalyzer->jetSmearedDownToken, jetsSmearedDown);

  _nJets = 0;
  _nJetswithDaughters = 0;
  _nDaughters = 0;

  auto jet            = jets->begin();
  auto jetSmeared     = jetsSmeared->begin();
  auto jetSmearedUp   = jetsSmearedUp->begin();
  auto jetSmearedDown = jetsSmearedDown->begin();
  for(; jet != jets->end(); ++jet, ++jetSmeared, ++jetSmearedUp, ++jetSmearedDown){
    if(_nJets == nJets_max) break;
    jecUnc.setJetEta(jetSmeared->eta());
    jecUnc.setJetPt(jetSmeared->pt());
    double unc = jecUnc.getUncertainty(true);

    if(std::max((1+unc)*jetSmeared->pt(), std::max(jetSmearedUp->pt(), jetSmearedDown->pt())) < 25) continue;

    _jetPt[_nJets]                    = jetSmeared->pt();
    _jetPt_JECDown[_nJets]            = jetSmeared->pt()*(1-unc);
    _jetPt_JECUp[_nJets]              = jetSmeared->pt()*(1+unc);
    _jetPt_JERDown[_nJets]            = jetSmearedDown->pt();
    _jetPt_JERUp[_nJets]              = jetSmearedUp->pt();
    _jetEta[_nJets]                   = jet->eta();
    _jetPhi[_nJets]                   = jet->phi();
    _jetE[_nJets]                     = jet->energy();
    //Old csvV2 b-tagger
    _jetCsvV2[_nJets]                 = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //new DeepFlavour tagger
    _jetDeepCsv_udsg[_nJets]          = jet->bDiscriminator("pfDeepCSVJetTags:probudsg");
    _jetDeepCsv_b[_nJets]             = jet->bDiscriminator("pfDeepCSVJetTags:probb");
    _jetDeepCsv_c[_nJets]             = jet->bDiscriminator("pfDeepCSVJetTags:probc");
    _jetDeepCsv_bb[_nJets]            = jet->bDiscriminator("pfDeepCSVJetTags:probbb");
//  _jetDeepCsv_cc[_nJets]            = jet->bDiscriminator("pfDeepCSVJetTags:probcc");
    _jetHadronFlavor[_nJets]         = jet->hadronFlavour();
    _jetId[_nJets]                    = jetId(*jet, false) + jetId(*jet, true); // 1: loose, 2: tight
    if(jet->numberOfDaughters() > 0) ++_nJetswithDaughters;
    
    for(unsigned d = 0; d < jet->numberOfDaughters(); ++d){
      if(_nDaughters > 100) std::cout << _nDaughters;
      const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet->daughter(d);
      _jet_tag_for_daughters[_nDaughters] = _nJets;
      _jet_daughter_pdgid[_nDaughters] 	  = daughter->pdgId();
      _jet_daughter_pt[_nDaughters] 	  = daughter->pt();
      _jet_daughter_eta[_nDaughters] 	  = daughter->eta();
      _jet_daughter_phi[_nDaughters] 	  = daughter->phi();
      _jet_daughter_energy[_nDaughters]   = daughter->energy();
      ++_nDaughters;
    }
    ++_nJets;
  }
  std::cout << std::endl;
  if(multilepAnalyzer->skim == "singlejet" and _nJets < 1) return false;
  if(multilepAnalyzer->skim == "FR" and _nJets < 1) return false;
  return true;
}

// Jet Id (https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016)
bool JetAnalyzer::jetId(const pat::Jet& j, bool tight){
  if(fabs(j.eta()) < 2.7){
    if(j.neutralHadronEnergyFraction() >= (tight? 0.90 : 0.99)) return false;
    if(j.neutralEmEnergyFraction() >= (tight? 0.90 : 0.99))     return false;
    if(j.chargedMultiplicity()+j.neutralMultiplicity() <= 1)    return false;
    if(fabs(j.eta()) < 2.4){
      if(j.chargedHadronEnergyFraction() <= 0)                  return false;
      if(j.chargedMultiplicity() <= 0)                          return false;
      if(j.chargedEmEnergyFraction() >= 0.99)                   return false;
    }
  } else if(fabs(j.eta()) < 3.0){
    if(j.neutralHadronEnergyFraction() >= 0.98)                 return false;
    if(j.neutralEmEnergyFraction() <= 0.01)                     return false;
    if(j.neutralMultiplicity() <= 2)                            return false;
  } else {
    if(j.neutralEmEnergyFraction() >= 0.90)                     return false;
    if(j.neutralMultiplicity() <= 10)                           return false;
  }
  return true;
}
