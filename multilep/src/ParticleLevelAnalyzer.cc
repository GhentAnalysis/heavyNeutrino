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
    outputTree->Branch("_pl_nJets",                 &_pl_nJets,                 "_pl_nJets/i");
    outputTree->Branch("_pl_jetPt",                 &_pl_jetPt,                 "_pl_jetPt[_pl_nJet]/D");
    outputTree->Branch("_pl_jetEta",                &_pl_jetEta,                "_pl_jetEta[_pl_nJet]/D");
    outputTree->Branch("_pl_jetPhi",                &_pl_jetPhi,                "_pl_jetPhi[_pl_nJet]/D");
    outputTree->Branch("_pl_jetE",                  &_pl_jetE,                  "_pl_jetE[_pl_nJet]/D");
    outputTree->Branch("_pl_jetHadronFlavor",       &_pl_jetHadronFlavor,       "_pl_jetHadronFlavor[_pl_nJet]/i");

}

bool ParticleLevelAnalyzer::analyze(const edm::Event& iEvent){
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

    _pl_nJets = 0;
    for(const reco::GenJet& p : *jets){
      _pl_jetPt[_pl_nJets]           = p.pt();
      _pl_jetEta[_pl_nJets]          = p.eta();
      _pl_jetPhi[_pl_nJets]          = p.phi();
      _pl_jetE[_pl_nJets]            = p.energy();
      _pl_jetHadronFlavor[_pl_nJets] = p.pdgId();
      ++_pl_nJets;
    }

    if(multilepAnalyzer->skim == "trilep"       and _pl_nL < 3)                    return false;
    if(multilepAnalyzer->skim == "dilep"        and _pl_nL < 3)                    return false;
    if(multilepAnalyzer->skim == "singlelep"    and _pl_nL < 1)                    return false;
    if(multilepAnalyzer->skim == "FR"           and (_pl_nL < 1 or _pl_nJets < 1)) return false;
    if(multilepAnalyzer->skim == "singlephoton" and _pl_nPh < 1)                   return false;
    if(multilepAnalyzer->skim == "diphoton"     and _pl_nPh < 2)                   return false;
    if(multilepAnalyzer->skim == "singlejet"    and _pl_nJets < 1)                 return false;
    return true;
}
