#ifndef GenMatching_h
#define GenMatching_h

//include other parts of framework
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"


class multilep;

class GenMatching{
  public:
    GenMatching(const edm::ParameterSet& iConfig, multilep*);
    ~GenMatching(){};

    //check if given reco candidate is prompt
    bool isPrompt(const reco::Candidate&, const reco::GenParticle&) const;
    void setGenParticles(const edm::Event&);    

    template<typename Lepton> const reco::GenParticle* findGenMatch(const Lepton& lepton) const{
        const reco::GenParticle* match = lepton.genParticle();
        if(match and match->pdgId() == lepton.pdgId()) return match;                  // Take the match found by genParticle() if it exists and its pdgId is the same
        else                                           return geometricMatch(lepton); // Else do geometric match
    }


  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    const reco::GenParticle* geometricMatch(const reco::Candidate&, const bool differentId = false) const;
    bool considerForMatching(const reco::Candidate&, const reco::GenParticle&, const bool differentId = false) const;
};
#endif
