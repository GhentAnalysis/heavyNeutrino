#ifndef GenMatching_h
#define GenMatching_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include <utility>

class GenMatching{
  public:
    typedef std::pair<const reco::GenParticle*, unsigned>   GenTypeMatch;
    typedef std::pair<const reco::GenParticle*, float>      GenDrMatch;
    typedef std::vector<GenDrMatch>                         GenDrMatches;
    typedef std::pair<const reco::Candidate*, GenDrMatches> LepToGenDrMatches;
    typedef std::vector<LepToGenDrMatches>                  LepToGenDrMatchesVector;
    typedef std::pair<const reco::Candidate*, GenTypeMatch> LepToGenTypeMatch;
    typedef std::vector<LepToGenTypeMatch>                  LepToGenTypeMatchVector;

    GenMatching(const edm::ParameterSet& iConfig);
    ~GenMatching(){};

    void matchGenToReco(const std::vector<reco::GenParticle>&, std::vector<const pat::Electron*>&, std::vector<const pat::Muon*>&, std::vector<const pat::Tau*>&);
    template <typename Lepton> void individualGenToRecoMatch(const std::vector<reco::GenParticle>& genParticles, const Lepton*, LepToGenDrMatchesVector&);
    const reco::GenParticle* returnGenMatch(const reco::Candidate&, unsigned&) const;

  private:
    LepToGenTypeMatchVector recogenmatchlist;
    bool                    allowMatchToAllIds;
};
#endif
