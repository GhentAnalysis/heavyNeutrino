#include <iostream>
#include "heavyNeutrino/multilep/interface/GenTools.h"

const reco::GenParticle* GenTools::getFirstMother(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    return (gen.numberOfMothers() == 0) ? nullptr : &genParticles[gen.motherRef(0).key()];
}

const reco::GenParticle* GenTools::getMother(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* mom = getFirstMother(gen, genParticles);
    if(!mom) return 0;
    else if(mom->pdgId() == gen.pdgId()) return getMother(*mom, genParticles);
    else return mom;
}

void GenTools::setDecayChain(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list){
   //  for(auto& i : list) std::cout << "list---> "<<i << std::endl;
    if((list.empty() or list.find(gen.pdgId())==list.end()) and gen.pdgId() != 2212) list.insert(gen.pdgId());
    //std::cout<<"in setdecaychain: gen.pdgId()"<<gen.pdgId()<<std::endl;
    if(gen.numberOfMothers() > 1) setDecayChain(genParticles[gen.motherRef(1).key()], genParticles, list);
    if(gen.numberOfMothers() > 0) setDecayChain(genParticles[gen.motherRef(0).key()], genParticles, list);
  //  for(auto& i : list) std::cout << "list  2   ---> "<<i << std::endl;

}

bool GenTools::bosonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ return (abs(entry) > 22 && abs(entry) < 26) || (abs(entry) == 9900012); }) != chain.cend();
}
 
bool GenTools::bBaryonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ return (abs(entry)/1000)%10 == 5; }) != chain.cend();
}

bool GenTools::bMesonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ unsigned mod = abs(entry)%1000; return mod >= 500 && mod < 600; }) != chain.cend();
}

bool GenTools::cBaryonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){return (abs(entry)/1000)%10 == 4; }) != chain.cend();
}

bool GenTools::cMesonInChain(const std::set<int>& chain){
    return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ unsigned mod = abs(entry)%1000; return mod >= 400 && mod < 500; }) != chain.cend();
}

bool GenTools::sBaryonInChain(const std::set<int>& chain){
    return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ return (abs(entry)/1000)%10 == 3; }) != chain.cend(); 
}

bool GenTools::lightMesonInChain(const std::set<int>& chain){
    return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ unsigned mod = abs(entry)%1000; return (mod >= 100 && mod < 400) || entry == 21; }) != chain.cend();
}

bool GenTools::lightBaryonInChain(const std::set<int>& chain){
    return std::find_if(chain.cbegin(), chain.cend(), 
            [](const int entry){ 
                if(abs(entry) == 2212) return false;
                unsigned red = (abs(entry)/1000)%10; 
                return (red == 1 || red == 2); 
            }) != chain.cend();
}

bool GenTools::pi0InChain(const std::set<int>& chain){
    return chain.find(111) != chain.cend();  
}

bool GenTools::photonInChain(const std::set<int>& chain){
    return chain.find(22) != chain.cend();
}

bool GenTools::tauInChain(const std::set<int>& chain){
    return chain.find(15) != chain.cend() || chain.find(-15) != chain.cend();
}

bool GenTools::udsInChain(const std::set<int>& chain){
    if(sBaryonInChain(chain))       return true;
    if(lightMesonInChain(chain))    return true;
    if(lightBaryonInChain(chain))   return true;
    /*
     to consider protons here too? 
    */
    return false;
}

//enumerated type to specify decay
enum decayType {
    W_L,
    W_T_L,
    W_B_L,
    W_B_C_L,
    W_B_C_T_L,
    W_B_T_L,
    W_C_L,
    W_C_T_L,
    B_L,
    B_C_L,
    B_C_T_L,
    B_T_L,
    C_L,
    C_T_L,
    B_Baryon,
    C_Baryon,
    pi_0,
    photon_,
    F_L
};

unsigned GenTools::provenance(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    std::set<int> decayChain;
    setDecayChain(gen, genParticles, decayChain);
    //first consider decays involving a boson
    if(bosonInChain(decayChain)){
        if(bMesonInChain(decayChain)){
            if(cMesonInChain(decayChain)){
                if(tauInChain(decayChain)){
                    return W_B_C_T_L;
                }
                return W_B_C_L;
            } else if(tauInChain(decayChain)){
                return W_B_T_L;
            }
            return W_B_L;
        } else if(cMesonInChain(decayChain)){
            if(tauInChain(decayChain)){
                return W_C_T_L;
            }
            return W_C_L;
        } else if(udsInChain(decayChain)){
            return pi_0;
        } else if(tauInChain(decayChain)){
            return W_T_L;
        }
        return W_L;
    } else if(bMesonInChain(decayChain)){
        if(cMesonInChain(decayChain)){
            if(tauInChain(decayChain)){
                return B_C_T_L;
            }
            return B_C_L;
        } else if(tauInChain(decayChain)){
            return B_T_L;
        }
        return B_L;
   } else if(cMesonInChain(decayChain)){
        if(tauInChain(decayChain)){
            return C_T_L;
        }
        return C_L;
   } else if(bBaryonInChain(decayChain)){
       return B_Baryon;
   } else if(cBaryonInChain(decayChain)){
       return C_Baryon;
   } else if(udsInChain(decayChain)){
       return pi_0;
   } else if(photonInChain(decayChain)){
       return photon_;
   }
   return F_L;
}

unsigned GenTools::provenanceCompressed(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    std::set<int> decayChain;
    setDecayChain(gen, genParticles, decayChain);
    if(bMesonInChain(decayChain) || bBaryonInChain(decayChain) ) return 1;          //lepton from heavy flavor decay
    if(cMesonInChain(decayChain) || cBaryonInChain(decayChain) ) return 2;          //lepton from c flavor decay
    if(bosonInChain(decayChain) ) return 0;                                         //lepton from boson
    if(!decayChain.empty()) return 3;                                               //light flavor fake
    return 4;                                                                       //unkown origin
}

bool GenTools::isPrompt(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* mom = getMother(gen, genParticles);
    if(abs(mom->pdgId()) == 15 && mom->isPromptDecayed()) return true;
    return (gen.isPromptFinalState() || gen.isPromptDecayed());
}

double GenTools::getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles){
    double minDeltaR = 10;
    for(auto& q : genParticles){
        if(q.pt() < 5)                                            continue;
        if(p.pt()-q.pt() < 0.0001)                                continue; // same particle
        if(q.status() != 1)                                       continue;
        if(q.pdgId() == 12 or q.pdgId() == 14 or q.pdgId() == 16) continue;
        minDeltaR = std::min(minDeltaR, deltaR(p.eta(), p.phi(), q.eta(), q.phi()));
    }
    return minDeltaR;
}
