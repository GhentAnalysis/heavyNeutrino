#ifndef TauTools_H
#define TauTools_H

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

namespace TauTools{

    void getNextDaughters(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list);

    //Check whether a tau decayed hadronically
    bool decayedHadronically(const reco::GenParticle&, const std::vector<reco::GenParticle>&);

    const reco::GenParticle* findMatch(const pat::Tau& reco, const std::vector<reco::GenParticle>& genParticles);

    const bool considerForMatching(const pat::Tau&, const reco::GenParticle&, const std::vector<reco::GenParticle>& genParticles);

    const unsigned int tauGenStatus(const reco::GenParticle*);
}

#endif
