#include "heavyNeutrino/multilep/interface/GenTools.h"

//include ROOT classes
#include "TLorentzVector.h"

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
    if((list.empty() or list.find(gen.pdgId())==list.end()) and gen.pdgId() != 2212) list.insert(gen.pdgId());
    if(gen.numberOfMothers() > 1) setDecayChain(genParticles[gen.motherRef(1).key()], genParticles, list);
    if(gen.numberOfMothers() > 0) setDecayChain(genParticles[gen.motherRef(0).key()], genParticles, list);
}

void GenTools::setDecayChainVector(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::vector<int>& list){
    if((list.empty() or gen.pdgId()!=list.back()) and gen.pdgId() != 2212) list.push_back(gen.pdgId());
    if(gen.numberOfMothers() > 1) setDecayChainVector(genParticles[gen.motherRef(1).key()], genParticles, list);
    if(gen.numberOfMothers() > 0) setDecayChainVector(genParticles[gen.motherRef(0).key()], genParticles, list);
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

unsigned GenTools::provenanceConversion(const reco::GenParticle& photon, const std::vector<reco::GenParticle>& genParticles) {
    //https://hypernews.cern.ch/HyperNews/CMS/get/susy-interpretations/192.html
    //99: not a photon
    //0: direct prompt photon (prompt and delta R with ME parton > 0.05)
    //1: fragmentation photon (prompt and delta R with ME parton < 0.05)
    //2: non-prompt photon
    if (photon.pdgId() != 22) return 99;
    if (!photon.isPromptFinalState()) return 2;
    if (photon.pt() < 10) return 1;

    TLorentzVector photonVec(photon.px(), photon.py(), photon.pz(), photon.energy() );
    for(auto& parton : genParticles){

        //only compare photon to ME partons
        if(parton.status() != 23) continue;

        //make sure parton is a parton
        unsigned partonId = abs( parton.pdgId() );
        if( ! ( (partonId == 21) || (partonId > 0 && partonId < 7) ) ) continue;
        
        //check separation of photon to parton
        TLorentzVector partonVec(parton.px(), parton.py(), parton.pz(), parton.energy() );
        if( photonVec.DeltaR(partonVec) < 0.05){
            return 1;
        }
    }
    return 0;
}

bool GenTools::isPrompt(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* mom = getMother(gen, genParticles);
    if(abs(mom->pdgId()) == 15 && mom->isPromptDecayed()) return true;
    return (gen.isPromptFinalState() || gen.isPromptDecayed());
}


/*
 * Returns true if an incoming gluon is found in the parental chain
 * (i.e. ignore gluons which have top, W-boson or Z-boson as parent)
 * Works for top-process [might need a more general implementation]
 */
bool GenTools::parentGluonIsIncoming(std::vector<int>& list){
  bool gluonEncountered = false;
  for(auto d : list){
    if(d==21) gluonEncountered = true;
    if(gluonEncountered and (abs(d)==6 or d==23 or d==24)) return false;
  }
  return true;
}

/*
 * Minimum deltaR between a gen particle and other gen particles with pt > ptCut
 * This could be used in trying to select the madgraph phase space on pythia level
 * The madgraph run card often contains deltaR cuts, such that events with getMinDeltaR(ptCut=5)<0.2
 * are typically out of the phase space of the generated sample
 * [but this is based on tuning and agreement with other groups, so maybe room for more studies/tuning]
 */
double GenTools::getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles, float ptCut){
    double minDeltaR = 10;
    for(auto& q : genParticles){
        if(q.pt() < ptCut)                                                       continue;
        if(q.status() != 1)                                                      continue;
        if(fabs(p.pt()-q.pt()) < 0.0001)                                         continue; // same particle
        if(abs(q.pdgId()) == 12 or abs(q.pdgId()) == 14 or abs(q.pdgId()) == 16) continue;
        minDeltaR = std::min(minDeltaR, deltaR(p.eta(), p.phi(), q.eta(), q.phi()));
    }
    return minDeltaR;
}
