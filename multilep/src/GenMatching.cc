#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "heavyNeutrino/multilep/interface/GenParticleManager.h"
#include "TLorentzVector.h" 

GenMatching::GenMatching(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer): multilepAnalyzer(multilepAnalyzer){};

void GenMatching::setGenParticles(const edm::Event& iEvent){
    iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);
}

reco::GenParticle const* GenMatching::findGenMatch(const reco::Candidate& reco, const bool differentId) const{
    reco::GenParticle const* match = nullptr;
    TLorentzVector recoV(reco.px(), reco.py(), reco.pz(), reco.energy());
    double minDeltaR = 99999.;
    for(std::vector<reco::GenParticle>::const_iterator genIt = genParticles->cbegin(); genIt != genParticles->cend(); ++genIt){
        if(toConsider(reco, *genIt, differentId) ){
            TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
            double deltaR = recoV.DeltaR(genV);
            if(deltaR < minDeltaR){
                minDeltaR = deltaR;
                match = &*genIt;
            }
        }
    } 
    if(minDeltaR > 0.2 && !differentId) match = findGenMatch(reco, true);
    return match;
}

bool GenMatching::toConsider(const reco::Candidate& reco, const reco::GenParticle& gen, const bool differentId) const{
    if(!differentId && (abs(reco.pdgId()) != abs(gen.pdgId())) ) return false;
    if(abs(reco.pdgId()) == 15 && abs(gen.pdgId()) == 15) return gen.status() == 2 && gen.isLastCopy();
    return gen.status() == 1;
}

bool GenMatching::isPrompt(const reco::Candidate& reco, const reco::GenParticle& match) const{
    if(abs(reco.pdgId()) == abs(match.pdgId()) || match.pdgId() == 22) return match.isPromptFinalState();
    return false;
}

void GenMatching::fillMatchingVars(const reco::Candidate& reco){
    const reco::GenParticle* match = findGenMatch(reco);
    if(match != nullptr){
        matchIsPrompt = isPrompt(reco, *match);
        matchPdgId = match->pdgId();
        provenance = GenTools::provenance(*match, *genParticles);
        origin = GenParticleManager::origin(*&match);
        originReduced = GenParticleManager::originReduced(GenParticleManager::origin(match));
    } else{
        matchIsPrompt = false;
        matchPdgId = 0;
        provenance = 4.;
        origin = -1;
        originReduced = -1;

    }
}




