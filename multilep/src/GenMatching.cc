#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "TLorentzVector.h" 

GenMatching::GenMatching(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer): multilepAnalyzer(multilepAnalyzer){};

void GenMatching::setGenParticles(const edm::Event& iEvent){
    iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);
}

const reco::GenParticle* GenMatching::geometricMatch(const reco::Candidate& reco, const bool differentId) const{
    reco::GenParticle const* match = nullptr;
    TLorentzVector recoV(reco.px(), reco.py(), reco.pz(), reco.energy());
    double minDeltaR = 99999.;
    for(std::vector<reco::GenParticle>::const_iterator genIt = genParticles->cbegin(); genIt != genParticles->cend(); ++genIt){
        if(considerForMatching(reco, *genIt, differentId) ){
            TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
            double deltaR = recoV.DeltaR(genV);
            if(deltaR < minDeltaR){
                minDeltaR = deltaR;
                match = &*genIt;
            }
        }
    } 
    if(minDeltaR > 0.2){
        if(!differentId){
            match = geometricMatch(reco, true);
        } else{
            //no decent match found!
            return nullptr;
        }
    }
    return match;
}

bool GenMatching::considerForMatching(const reco::Candidate& reco, const reco::GenParticle& gen, const bool differentId) const{
    //if gen particle is not of same as reco particle
    if(!sameParticle(reco, gen)){
        if(!differentId){
            return false;
        } else {
            //allow matching to photons
            if (abs(gen.pdgId()) != 22) return false;
        }
    }
    if(abs(reco.pdgId()) == 15 && abs(gen.pdgId()) == 15) return gen.status() == 2 && gen.isLastCopy();
    return gen.status() == 1;
}

bool GenMatching::sameParticle(const reco::Candidate& reco, const reco::GenParticle& gen) const{
    return ( abs(reco.pdgId()) == abs(gen.pdgId()) );
}

bool GenMatching::isPrompt(const reco::Candidate& reco, const reco::GenParticle& match) const{
    if(abs(reco.pdgId()) == abs(match.pdgId()) || match.pdgId() == 22) return GenTools::isPrompt(match, *genParticles);
    return false;
}
