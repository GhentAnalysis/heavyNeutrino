#include "heavyNeutrino/multilep/interface/TauMatchTest.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//include ROOT classes
#include "TLorentzVector.h"


bool TauMatchTest::decayedHadronically(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    std::set<int> decayProducts;
    if(abs(gen.pdgId()) == 11 or abs(gen.pdgId()) == 13) return false;
    getNextDaughters(gen, genParticles, decayProducts);
    
    if(decayProducts.find(11) != decayProducts.cend() || decayProducts.find(-11) != decayProducts.cend() || decayProducts.find(13) != decayProducts.cend() || decayProducts.find(-13) != decayProducts.cend()){
        return false;
    }   
    return true;
    
}

void TauMatchTest::getNextDaughters(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list){
    for(unsigned d = 0; d < gen.numberOfDaughters(); ++d){
        if(list.empty() or list.find(gen.pdgId())==list.end()) list.insert(genParticles[gen.daughterRef(d).key()].pdgId());
    }      
}


const reco::GenParticle* TauMatchTest::matchGenJetToGenPart(const reco::GenJet& jet, const std::vector<reco::GenParticle>& genParticles){
    reco::GenParticle const* match = nullptr;
    //std::cout << "Clone " << jet.eta() << " " << jet.phi() << std::endl;
    TLorentzVector recoV(jet.px(), jet.py(), jet.pz(), jet.energy());
    double minDeltaR = 99999.;
    for(auto genIt = genParticles.cbegin(); genIt != genParticles.cend(); ++genIt){
        if(abs(genIt->pdgId()) != 15 or !TauMatchTest::decayedHadronically(*genIt, genParticles))   continue;
        TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
        //std::cout << "genPart " << genIt->eta() << " " << genIt->phi() << std::endl;
        double deltaR = recoV.DeltaR(genV);
        //std::cout << deltaR << std::endl;
        if(deltaR < minDeltaR){
            minDeltaR = deltaR;
            match = &*genIt;
        }
    }
    if(match) std::cout << "GenJET MATCH" << match->eta() << " " << match->phi() << " " << std::endl;
      
    return match;          

}

const reco::GenParticle* TauMatchTest::findJetMatch(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles){
    reco::GenJet const* match = tau.genJet();
    if(match) return TauMatchTest::matchGenJetToGenPart(*match, genParticles);
    return nullptr;
    
}

const reco::GenParticle* TauMatchTest::findLeptonMatch(const pat::Tau& reco, const std::vector<reco::GenParticle>& genParticles){
    reco::GenParticle const* match = nullptr;
    TLorentzVector recoV(reco.px(), reco.py(), reco.pz(), reco.energy());
    double minDeltaR = 99999.;
    for(auto genIt = genParticles.cbegin(); genIt != genParticles.cend(); ++genIt){
        if(abs(genIt->pdgId()) != 11 and abs(genIt->pdgId()) != 13) continue;
        if(genIt->pt() < 8)     continue;
        TLorentzVector genV(genIt->px(), genIt->py(), genIt->pz(), genIt->energy());
        double deltaR = recoV.DeltaR(genV);
        if(deltaR < minDeltaR){
            minDeltaR = deltaR;
            match = &*genIt;
        }
    }

    if(minDeltaR > 0.2) return nullptr;
    
    if(match) std::cout << "Lepton MATCH" << match->eta() << " " << match->phi() << " " << std::endl;
    return match;    
}


const reco::GenParticle* TauMatchTest::findMatch(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles){
    reco::GenParticle const* match = TauMatchTest::findJetMatch(tau, genParticles);
    if(!match){
        match = TauMatchTest::findLeptonMatch(tau, genParticles);          
    }
    return match;
}

const unsigned int TauMatchTest::tauGenStatus(const reco::GenParticle* match){

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
