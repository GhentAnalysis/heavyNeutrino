#ifndef GenMatching_h
#define GenMatching_h

//include other parts of framework
#include "heavyNeutrino/multilep/plugins/multilep.h"

class multilep;

class GenMatching{
  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    
    reco::GenParticle const* findGenMatch(const reco::Candidate&, const bool differentId = false) const;
    bool toConsider(const reco::Candidate&, const reco::GenParticle&, const bool differentId = false) const;
    int matchPdgId;
    bool matchIsPrompt;
    unsigned provenance;
    unsigned origin;
    unsigned originReduced;
  public:
    GenMatching(const edm::ParameterSet& iConfig, multilep*);
    ~GenMatching(){};
    //check if given reco candidate is prompt
    bool isPrompt(const reco::Candidate&, const reco::GenParticle&) const;
    void setGenParticles(const edm::Event&);    
    //fill match variables
    void fillMatchingVars(const reco::Candidate&);
    int pdgIdMatch() const {return matchPdgId; }
    bool promptMatch() const {return matchIsPrompt;}
    unsigned getProvenance() const {return provenance;}
};
#endif
