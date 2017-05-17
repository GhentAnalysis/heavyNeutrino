#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer)
{};


void JetAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nJets",                        &_nJets,                        "_nJets/b");
  outputTree->Branch("_jetPt",                        &_jetPt,                        "_jetPt[_nJets]/D");
  outputTree->Branch("_jetEta",                       &_jetEta,                       "_jetEta[_nJets]/D");
  outputTree->Branch("_jetPhi",                       &_jetPhi,                       "_jetPhi[_nJets]/D");
  outputTree->Branch("_jetE",                         &_jetE,                         "_jetE[_nJets]/D");
}

void JetAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);

  _nJets = 0;
  for(const pat::Jet& jet: *jets){
    if(jet.pt() < 25)         continue;
    if(fabs(jet.eta()) > 2.4) continue;
    _jetPt[_nJets]  = jet.pt();
    _jetEta[_nJets] = jet.eta();
    _jetPhi[_nJets] = jet.phi();
    _jetE[_nJets]   = jet.energy();
    ++_nJets;
  }
}
