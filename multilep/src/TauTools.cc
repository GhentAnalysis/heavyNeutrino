#include "heavyNeutrino/multilep/interface/TauTools.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//include ROOT classes
#include "TLorentzVector.h"


bool TauTools::decayedHadronically(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    std::set<int> decayProducts;
    if(abs(gen.pdgId()) == 11 or abs(gen.pdgId()) == 13) return false;
    getNextDaughters(gen, genParticles, decayProducts);
    
    if(decayProducts.find(11) != decayProducts.cend() || decayProducts.find(-11) != decayProducts.cend() || decayProducts.find(13) != decayProducts.cend() || decayProducts.find(-13) != decayProducts.cend()){
        return false;
    }   
    return true;
    
}

void TauTools::getNextDaughters(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list){
    for(unsigned d = 0; d < gen.numberOfDaughters(); ++d){
        if(list.empty() or list.find(gen.pdgId())==list.end()) list.insert(genParticles[gen.daughterRef(d).key()].pdgId());
    }      
}

//If match to genJet, look for the corresponding gen tau, else try to see if match with electron or muon 
const bool TauTools::considerForMatching(const pat::Tau& tau, const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    if(tau.genJet()){
        if(abs(gen.pdgId()) != 15 or !TauTools::decayedHadronically(gen, genParticles)) return false;
        else return true;
    }
    else{
        if(abs(gen.pdgId()) != 11 and abs(gen.pdgId()) != 13) return false;
        if(gen.pt() < 8)     return false;
        else return true;
    }
    return false; 
}

const reco::GenParticle* TauTools::findMatch(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles){
    reco::GenParticle const*match{nullptr};
    TLorentzVector recoV(tau.px(), tau.py(), tau.pz(), tau.energy());
    
    double minDeltaR = 0.2;
    for(auto genIt = genParticles.cbegin(); genIt != genParticles.cend(); ++genIt){
        if(!TauTools::considerForMatching(tau, *genIt, genParticles)) continue;
    
        TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
        
        double deltaR = recoV.DeltaR(genV);
        if(deltaR < minDeltaR){
            minDeltaR = deltaR;
            match = &*genIt;
        }
    }
    return match;
}

const reco::GenJet* TauTools::findMatchedGenJet(const reco::GenParticle& genTau, const std::vector<reco::GenJet>& genJets){
    reco::GenJet const*match{nullptr};
    TLorentzVector recoV(genTau.px(), genTau.py(), genTau.pz(), genTau.energy());
    
    double minDeltaR = 0.2;
    for(auto genIt = genJets.cbegin(); genIt != genJets.cend(); ++genIt){    
        TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
        
        double deltaR = recoV.DeltaR(genV);
        if(deltaR < minDeltaR){
            minDeltaR = deltaR;
            match = &*genIt;
        }
    }
    return match;
}


//As in https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#MC_Matching
const unsigned int TauTools::tauGenStatus(const reco::GenParticle* match){

    if(!match) return 6;
    
    if(abs(match->pdgId()) == 15){
        return 5;
    }

    else if(abs(match->pdgId()) == 11){
        if(match->statusFlags().isDirectPromptTauDecayProduct()) return 3;
        else if(match->statusFlags().isPrompt()) return 1;
    }

    else if(abs(match->pdgId()) == 13){
        if(match->statusFlags().isDirectPromptTauDecayProduct()) return 4;
        else if(match->statusFlags().isPrompt()) return 2;

    }
    return 6;

}
