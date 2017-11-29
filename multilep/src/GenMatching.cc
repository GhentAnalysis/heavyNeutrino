#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "TLorentzVector.h" 

GenMatching::GenMatching(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer): multilepAnalyzer(multilepAnalyzer){};

void GenMatching::setGenParticles(const edm::Event& iEvent){
    iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);
}

reco::GenParticle const* GenMatching::findGenMatch(const reco::Candidate& reco, const bool differentId){
    reco::GenParticle const* match = nullptr;
    TLorentzVector recoV(reco.px(), reco.py(), reco.pz(), reco.energy());
    double minDeltaR = 99999.;
    for(std::vector<reco::GenParticle>::const_iterator genIt = genParticles->cbegin(); genIt != genParticles->cend(); ++genIt){
        if(toConsider(reco, *genIt) && (differentId || abs(genIt->pdgId()) == abs(reco.pdgId()) ) ){
            TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
            double deltaR = recoV.DeltaR(genV);
            if(deltaR < minDeltaR){
                minDeltaR = deltaR;
                match = &*genIt;
            }
        }
    } 
    if(minDeltaR > 0.2) match = findGenMatch(reco, true);
    return match;
};

bool GenMatching::toConsider(const reco::Candidate& reco, const reco::GenParticle& gen){
    if(abs(reco.pdgId()) == 15 && abs(gen.pdgId()) == 15) return gen.status() == 2 && gen.isLastCopy();
    return gen.status() == 1;
}

bool GenMatching::isPrompt(const reco::Candidate& reco){
    const reco::GenParticle& match = *findGenMatch(reco);
    if(fabs(reco.pdgId()) == fabs(match.pdgId()) || match.pdgId() == 22) match.isPromptFinalState();
    return false;
}
