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
    template <typename Lepton> void fillMatchingVars(const Lepton& reco){
        const reco::GenParticle* match = findGenMatch(reco);
        if(match != nullptr){
            matchIsPrompt = isPrompt(reco, *match);
            matchPdgId = match->pdgId();
            provenance = GenTools::provenance(*match, *genParticles);
            provenanceCompressed = (matchIsPrompt ? 0 : GenTools::provenanceCompressed(*match, *genParticles) );
            partonPt = GenTools::getParton(*match, *genParticles) == 0 ? 0 : GenTools::getParton(*match, *genParticles)->pt();
            provenanceConversion = GenTools::provenanceConversion(*match, *genParticles);
            momPdgId = ( GenTools::getMother(*match, *genParticles) )->pdgId();
        } else{
            matchIsPrompt = false;
            matchPdgId = 0;
            provenanceCompressed = 4;
            provenance = 18;
            partonPt = 0;
            provenanceConversion = 99;
            momPdgId = 0;
        }
    }

    //return values
    int pdgIdMatch() const { return matchPdgId; }
    int pdgIdMom() const { return momPdgId; }
    bool promptMatch() const { return matchIsPrompt; }
    unsigned getProvenance() const { return provenance; }
    unsigned getProvenanceCompressed() const{ return provenanceCompressed; }
    double getPartonPt() const {return partonPt;}
    unsigned getProvenanceConversion() const{ return provenanceConversion; }

  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    const reco::GenParticle* geometricMatch(const reco::Candidate&, const bool differentId = false) const;

    template<typename Lepton> const reco::GenParticle* findGenMatch(const Lepton& lepton) const{
        const reco::GenParticle* match = lepton.genParticle();
        //short circuit assumed here!
        if(match == nullptr || match->pdgId() != lepton.pdgId()){
            //std::cout << "gen matching is running for lepton with pt: " << lepton.pt() << std::endl;
            return geometricMatch(lepton);
        }
        //std::cout << "cmssw matching is running for lepton with pt " << lepton.pt() << std::endl;
        return match;
    }

    bool considerForMatching(const reco::Candidate&, const reco::GenParticle&, const bool differentId = false) const;
    bool sameParticle(const reco::Candidate&, const reco::GenParticle&) const;
    int matchPdgId;
    int momPdgId;
    bool matchIsPrompt;
    unsigned provenance;
    unsigned provenanceCompressed;
    double partonPt;
    unsigned provenanceConversion;
};
#endif
