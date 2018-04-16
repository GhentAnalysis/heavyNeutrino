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

    //fill match variables
    template <typename Lepton> void fillMatchingVars(const Lepton&);

    //return values
    unsigned genLIndex() {return genLindex;}
    // unsigned genPhIndex() {return genPhindex;}
    int pdgIdMatch() const {return matchPdgId; }
    bool promptMatch() const {return matchIsPrompt;}
    unsigned getProvenance() const {return provenance;}
    unsigned getProvenanceCompressed() const{ return provenanceCompressed; }
    unsigned getProvenanceConversion() const{ return provenanceConversion; }

  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    const reco::GenParticle* geometricMatch(const reco::Candidate&, const bool differentId = false) const;

    template<typename Lepton> const reco::GenParticle* findGenMatch(const Lepton& lepton) const{
        const reco::GenParticle* match = lepton.genParticle();
        //short circuit assumed here!
        if(match == nullptr || match->pdgId() != lepton.pdgId()){
            return geometricMatch(lepton);
        }
        return match;
    }

    bool considerForMatching(const reco::Candidate&, const reco::GenParticle&, const bool differentId = false) const;
    bool sameParticle(const reco::Candidate&, const reco::GenParticle&) const;
    unsigned genLindex;
    // unsigned genPhindex;
    int matchPdgId;
    bool matchIsPrompt;
    unsigned provenance;
    unsigned provenanceCompressed;
    unsigned provenanceConversion;
};
#endif
