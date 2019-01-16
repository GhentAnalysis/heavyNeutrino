#ifndef GenMatching_h
#define GenMatching_h

//include other parts of framework
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"


class multilep;

class GenMatching{
  public:
    typedef std::pair<const reco::GenParticle*, unsigned>   GenTypeMatch;
    typedef std::pair<const reco::GenParticle*, float>      GenDrMatch;
    typedef std::vector<GenDrMatch>                         GenDrMatches;
    typedef std::pair<const reco::Candidate*, GenDrMatches> LepToGenDrMatches;
    typedef std::vector<LepToGenDrMatches>                  LepToGenDrMatchesVector;
    typedef std::pair<const reco::Candidate*, GenTypeMatch> LepToGenTypeMatch;
    typedef std::vector<LepToGenTypeMatch>                  LepToGenTypeMatchVector;

    GenMatching(const edm::ParameterSet& iConfig, multilep*);
    ~GenMatching(){};

    //check if given reco candidate is prompt
    bool isPrompt(const reco::Candidate&, const reco::GenParticle&) const;

    void resetGenMatchingVector() {recogenmatchlist.clear();}
    void setGenParticles(const edm::Event&);
    void setPatParticles(std::vector<const pat::Electron*> &eleVec, std::vector<const pat::Muon*> &muoVec, std::vector<const pat::Tau*> &tauVec){     
      patElectrons = eleVec;
      patMuons     = muoVec;
      patTaus      = tauVec;
   }

    void matchGenToReco();
    template <typename Lepton> void individualGenToRecoMatch(const Lepton*, LepToGenDrMatchesVector&);
    const reco::GenParticle* returnGenMatch(const reco::Candidate*, unsigned&) const;

    //fill match variables
    template <typename Lepton> void fillMatchingVars(const Lepton&, const reco::GenParticle*, unsigned, const unsigned lepcnt = 0);
    template <typename Lepton> void fillMatchingVars(const Lepton&);

    //return values
    unsigned genLIndex() {return genLindex;}
    // unsigned genPhIndex() {return genPhindex;}
    unsigned typeMatch() {return matchType;}
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
    std::vector<const pat::Electron*> patElectrons;
    std::vector<const pat::Muon*    > patMuons;
    std::vector<const pat::Tau*     > patTaus;

    LepToGenTypeMatchVector recogenmatchlist;

    template<typename Lepton> const reco::GenParticle* findGenMatch(const Lepton& lepton, const bool allowallids = false) const{
        const reco::GenParticle* match = lepton.genParticle();
        //short circuit assumed here!
        if(match == nullptr || (allowallids==false && match->pdgId()!=lepton.pdgId())){
          return GenTools::geometricMatch(lepton, *genParticles, allowallids);
        }
        return match;
    }

    unsigned genLindex;
    // unsigned genPhindex;
    unsigned matchType;
    int matchPdgId;
    int momPdgId;
    bool matchIsPrompt;
    bool matchIsPromptFinalState;
    bool matchIsPromptDecayed;
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
