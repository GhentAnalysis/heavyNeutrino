#ifndef GenTools_H
#define GenTools_H

#include <set>

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

namespace GenTools{
    const reco::GenParticle* getFirstMother(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    const int getFirstMotherIndex(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    const reco::GenParticle* getMother(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    int getMotherPdgId(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    //return decay chain for a particle;
    void setDecayChain(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list);
    bool particleInChain(const reco::GenParticle&, const std::vector<reco::GenParticle>&, const reco::GenParticle&);
    bool hasOnlyIncomingGluonsInChain(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles);

    //scan decay chain for certain types of particles
    bool bosonInChain(const std::set<int>&);
    bool bBaryonInChain(const std::set<int>&);
    bool bMesonInChain(const std::set<int>&);
    bool cBaryonInChain(const std::set<int>&);
    bool cMesonInChain(const std::set<int>&);
    bool sBaryonInChain(const std::set<int>&);
    bool sMesonInChain(const std::set<int>&);
    bool lightBaryonInChain(const std::set<int>&);
    bool lightMesonInChain(const std::set<int>&);
    bool pi0InChain(const std::set<int>&);
    bool photonInChain(const std::set<int>&);
    bool udsInChain(const std::set<int>&);
    bool tauInChain(const std::set<int>&);
    //find the provenance of a particle using the contents of its decayChain
    unsigned provenance(const reco::GenParticle*, const std::vector<reco::GenParticle>&);
    unsigned provenanceCompressed(const reco::GenParticle*, const std::vector<reco::GenParticle>&, bool isPrompt);

    //check whether photon comes from ME in conversion
    unsigned provenanceConversion(const reco::GenParticle*, const std::vector<reco::GenParticle>&);

    bool isPrompt(const reco::GenParticle&, const std::vector<reco::GenParticle>&); //function to check if particle is prompt TO BE USED INSTEAD OF CMSSW BUILTIN
    bool passParentage(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles);
    bool noMesonsInChain(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles);
    bool phoAndPiNear(const pat::Photon& photon, const std::vector<reco::GenParticle>& genParticles);
    double getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles, float ptCut=5);
    double getMinDeltaRTTG(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles, float ptCut=5);

    const reco::GenParticle* geometricMatch(const reco::Candidate& reco, const std::vector<reco::GenParticle>& genParticles, const bool differentId=false);
    bool considerForMatching(const reco::Candidate& reco, const reco::GenParticle& gen, const bool differentId);


}
#endif
