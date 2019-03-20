//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//include ROOT classes
#include "TLorentzVector.h"

#include "heavyNeutrino/multilep/interface/ParticleLevelAnalyzer.h"

/*
 * Storing plerator particles
 */

ParticleLevelAnalyzer::ParticleLevelAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer){};

void ParticleLevelAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_pl_met",                   &_pl_met,                   "_pl_met/D");
    outputTree->Branch("_pl_metPhi",                &_pl_metPhi,                "_pl_metPhi/D");
    outputTree->Branch("_pl_nPh",                   &_pl_nPh,                   "_pl_nPh/i");
    outputTree->Branch("_pl_phPt",                  &_pl_phPt,                  "_pl_phPt[_pl_nPh]/D");
    outputTree->Branch("_pl_phEta",                 &_pl_phEta,                 "_pl_phEta[_pl_nPh]/D");
    outputTree->Branch("_pl_phPhi",                 &_pl_phPhi,                 "_pl_phPhi[_pl_nPh]/D");
    outputTree->Branch("_pl_phE",                   &_pl_phE,                   "_pl_phE[_pl_nPh]/D");
    outputTree->Branch("_pl_nL",                    &_pl_nL,                    "_pl_nL/i");
    outputTree->Branch("_pl_lPt",                   &_pl_lPt,                   "_pl_lPt[_pl_nL]/D");
    outputTree->Branch("_pl_lEta",                  &_pl_lEta,                  "_pl_lEta[_pl_nL]/D");
    outputTree->Branch("_pl_lPhi",                  &_pl_lPhi,                  "_pl_lPhi[_pl_nL]/D");
    outputTree->Branch("_pl_lE",                    &_pl_lE,                    "_pl_lE[_pl_nL]/D");
    outputTree->Branch("_pl_lFlavor",               &_pl_lFlavor,               "_pl_lFlavor[_pl_nL]/i");
    outputTree->Branch("_pl_lCharge",               &_pl_lCharge,               "_pl_lCharge[_pl_nL]/I");
    outputTree->Branch("_pl_nJet",                  &_pl_nJet,                  "_pl_nJet/i");
    outputTree->Branch("_pl_jetPt",                 &_pl_jetPt,                 "_pl_jetPt[_pl_nJet]/D");
    outputTree->Branch("_pl_jetEta",                &_pl_jetEta,                "_pl_jetEta[_pl_nJet]/D");
    outputTree->Branch("_pl_jetPhi",                &_pl_jetPhi,                "_pl_jetPhi[_pl_nJet]/D");
    outputTree->Branch("_pl_jetE",                  &_pl_jetE,                  "_pl_jetE[_pl_nJet]/D");
    outputTree->Branch("_pl_jetHadronFlavor",       &_pl_jetHadronFlavor,       "_pl_jetHadronFlavor[_pl_nJet]/i");

}

void ParticleLevelAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<reco::GenParticle>> photons; iEvent.getByToken(multilepAnalyzer->particleLevelPhotonsToken, photons);
    edm::Handle<std::vector<reco::GenJet>> leptons;      iEvent.getByToken(multilepAnalyzer->particleLevelLeptonsToken, leptons);
    edm::Handle<std::vector<reco::GenJet>> jets;         iEvent.getByToken(multilepAnalyzer->particleLevelJetsToken, jets);
    edm::Handle<std::vector<reco::MET>> mets;            iEvent.getByToken(multilepAnalyzer->particleLevelMetsToken, mets);

    _pl_met    = mets->at(0).pt();
    _pl_metPhi = mets->at(0).phi();

    _pl_nPh = 0;
    for(const reco::GenParticle& p : *photons){
      _pl_phPt[_pl_nPh]  = p.pt();
      _pl_phEta[_pl_nPh] = p.eta();
      _pl_phPhi[_pl_nPh] = p.phi();
      _pl_phE[_pl_nPh]   = p.energy();
      ++_pl_nPh;
    }

    _pl_nL = 0;
    for(const reco::GenJet& p : *leptons){
      _pl_lPt[_pl_nL]     = p.pt();
      _pl_lEta[_pl_nL]    = p.eta();
      _pl_lPhi[_pl_nL]    = p.phi();
      _pl_lE[_pl_nL]      = p.energy();
      _pl_lCharge[_pl_nL] = p.charge();

      if(abs(p.pdgId()) == 11)      _pl_lFlavor[_pl_nL] = 0;
      else if(abs(p.pdgId()) == 13) _pl_lFlavor[_pl_nL] = 1;
      else                          _pl_lFlavor[_pl_nL] = 2;
      ++_pl_nL;
    }

    _pl_nJet = 0;
    for(const reco::GenJet& p : *jets){
      _pl_jetPt[_pl_nJet]           = p.pt();
      _pl_jetEta[_pl_nJet]          = p.eta();
      _pl_jetPhi[_pl_nJet]          = p.phi();
      _pl_jetE[_pl_nJet]            = p.energy();
      _pl_jetHadronFlavor[_pl_nJet] = p.pdgId();
      ++_pl_nJet;
    }
}
