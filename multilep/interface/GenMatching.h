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
    bool promptFinalStateMatch() const {return matchIsPromptFinalState;}
    bool promptDecayedMatch() const {return matchIsPromptDecayed;}

    unsigned getProvenance() const {return provenance;}
    unsigned getProvenanceCompressed() const{ return provenanceCompressed; }
    unsigned getProvenanceConversion() const{ return provenanceConversion; }
    double getMatchPt() const{ return matchPt; }
    double getMatchEta() const{ return matchEta; }
    double getMatchPhi() const{ return matchPhi; }
    double getMatchVertexX() const{ return matchXvtx; }
    double getMatchVertexY() const{ return matchYvtx; }
    double getMatchVertexZ() const{ return matchZvtx; }

  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    const reco::GenParticle* geometricMatch(const reco::Candidate&, const bool differentId = false) const;

    template<typename Lepton> const reco::GenParticle* findGenMatch(const Lepton& lepton, const bool allowallids = false) const{
        const reco::GenParticle* match = lepton.genParticle();
        //short circuit assumed here!
        if(match == nullptr || (allowallids==false && match->pdgId()!=lepton.pdgId())){
	    return geometricMatch(lepton, allowallids);
        }
        return match;
    }

    bool considerForMatching(const reco::Candidate&, const reco::GenParticle&, const bool differentId = false) const;
    bool sameParticle(const reco::Candidate&, const reco::GenParticle&) const;
    unsigned genLindex;
    // unsigned genPhindex;
    int matchPdgId;
    int momPdgId;
    bool matchIsPrompt;
    unsigned provenance;
    unsigned provenanceCompressed;
    unsigned provenanceConversion;
    double matchPt;
    double matchEta;
    double matchPhi;
    double matchXvtx;
    double matchYvtx;
    double matchZvtx;
    bool allowMatchToAllIds;
};
#endif
