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

  private:
    multilep* multilepAnalyzer;
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    std::vector<const pat::Electron*> patElectrons;
    std::vector<const pat::Muon*    > patMuons;
    std::vector<const pat::Tau*     > patTaus;

    LepToGenTypeMatchVector recogenmatchlist;

    bool allowMatchToAllIds;
};
#endif
