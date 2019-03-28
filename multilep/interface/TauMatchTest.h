#ifndef TAU_MATCH_TEST_H
#define TAU_MATCH_TEST_H

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

namespace TauMatchTest{

    void getNextDaughters(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list);

    //Check whether a tau decayed hadronically
    bool decayedHadronically(const reco::GenParticle&, const std::vector<reco::GenParticle>&);

    const reco::GenParticle* findJetMatch(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles);

    const reco::GenParticle* findLeptonMatch(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles);
    
    const reco::GenParticle* findMatch(const pat::Tau& reco, const std::vector<reco::GenParticle>& genParticles);

    const reco::GenParticle* matchGenJetToGenPart(const reco::GenJet& jet, const std::vector<reco::GenParticle>& genParticles);

    const unsigned int tauGenStatus(const reco::GenParticle*);
}

#endif
