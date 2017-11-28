#ifndef GenMatching_h
#define GenMatching_h

//include other parts of framework
#include "heavyNeutrino/multilep/plugins/multilep.h"

class multilep;

class GenMatching{
  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    
    void setGenParticles(const edm::Event&);    
    reco::GenParticle const* findGenMatch(const reco::Candidate&, const bool differentId = false);
    bool toConsider(const reco::Candidate&, const reco::GenParticle&);
  public:
    GenMatching(const edm::ParameterSet& iConfig, multilep*);
    ~GenMatching(){};
    //check if given reco candidate is prompt
    bool isPrompt(const reco::Candidate&);
};
#endif
