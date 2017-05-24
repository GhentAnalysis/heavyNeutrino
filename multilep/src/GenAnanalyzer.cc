#include "heavyNeutrino/multilep/interface/GenAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
/*
 * Storing generator particles
 * Currently storing everything such that we have it available for studies on tuple level
 * Might consider trimming it down to save space
 */


GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer)
{};


void GenAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nGen",                      &_nGen,                      "_nGen/I");  // do not use /b here, could save you a lot of debugging time
  outputTree->Branch("_genPt",                     &_genPt,                     "_genPt[_nGen]/D");
  outputTree->Branch("_genEta",                    &_genEta,                    "_genEta[_nGen]/D");
  outputTree->Branch("_genPhi",                    &_genPhi,                    "_genPhi[_nGen]/D");
  outputTree->Branch("_genMass",                   &_genMass,                   "_genMass[_nGen]/D");
  outputTree->Branch("_genCharge",                 &_genCharge,                 "_genCharge[_nGen]/I");
  outputTree->Branch("_genPdgId",                  &_genPdgId,                  "_genPdgId[_nGen]/I");
  outputTree->Branch("_genStatus",                 &_genStatus,                 "_genStatus[_nGen]/I");
  outputTree->Branch("_genFromHardProcess",        &_genFromHardProcess,        "_genFromHardProcess[_nGen]/O");
  outputTree->Branch("_genIsPrompt",               &_genIsPrompt,               "_genIsPrompt[_nGen]/O");
  outputTree->Branch("_genMotherIndex1",           &_genMotherIndex1,           "_genMotherIndex1[_nGen]/I");
  outputTree->Branch("_genMotherIndex2",           &_genMotherIndex2,           "_genMotherIndex2[_nGen]/I");
  outputTree->Branch("_genDaughterIndex1",         &_genDaughterIndex1,         "_genDaughterIndex1[_nGen]/I");
  outputTree->Branch("_genDaughterIndex2",         &_genDaughterIndex2,         "_genDaughterIndex2[_nGen]/I");
}

void GenAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<reco::GenParticle>> genParticles; iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);

  if(!genParticles.isValid()) return; 

  _nGen = 0;
  for(auto p = genParticles->begin(); p != genParticles->end() and _nGen < nGen_max; ++p){
    _genPt[_nGen]                 = p->pt();
    _genEta[_nGen]                = p->eta();
    _genPhi[_nGen]                = p->phi();
    _genMass[_nGen]               = p->mass();
    _genCharge[_nGen]             = p->charge();
    _genPdgId[_nGen]              = p->pdgId();
    _genStatus[_nGen]             = p->status();
    _genFromHardProcess[_nGen]    = p->fromHardProcessFinalState();
    _genIsPrompt[_nGen]           = p->isPromptFinalState() or p->isDirectPromptTauDecayProductFinalState() or p->isHardProcess();
    _genMotherIndex1[_nGen]       = p->numberOfMothers() > 0 ?   (p->motherRef(0).key()) : -1;
    _genMotherIndex2[_nGen]       = p->numberOfMothers() > 1 ?   (p->motherRef(1).key()) : -1;
    _genDaughterIndex1[_nGen]     = p->numberOfDaughters() > 0 ? (p->daughterRef(0).key()) : -1;
    _genDaughterIndex2[_nGen]     = p->numberOfDaughters() > 1 ? (p->daughterRef(1).key()) : -1;

    ++_nGen;
  }
}
