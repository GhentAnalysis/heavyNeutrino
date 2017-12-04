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
    if(list.empty() && gen.pdgId() != 2212 && list.find(gen.pdgId()) == list.end()) list.insert(gen.pdgId());
    if(gen.numberOfMothers() > 1) setDecayChain(genParticles[gen.motherRef(1).key()], genParticles, list);
    if(gen.numberOfMothers() > 0) setDecayChain(genParticles[gen.motherRef(0).key()], genParticles, list);
}

bool GenTools::bosonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ return abs(entry) > 22 && abs(entry) < 26; }) != chain.cend();
}
 
bool GenTools::bBaryonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ return (entry/1000)%10 == 5; }) != chain.cend();
}

bool GenTools::bMesonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ unsigned mod = entry%1000; return mod >= 500 && mod <= 600; }) != chain.cend();
}

bool GenTools::cBaryonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){return (entry/1000)%10 == 4; }) != chain.cend();
}

bool GenTools::cMesonInChain(const std::set<int>& chain){
   return std::find_if(chain.cbegin(), chain.cend(), [](const int entry){ unsigned mod = entry%1000; return mod >= 400 && mod <= 500; }) != chain.cend();
}

unsigned GenTools::provenance(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles){
    std::set<int> decayChain;
    setDecayChain(gen, genParticles, decayChain);
    if(bMesonInChain(decayChain) || bBaryonInChain(decayChain) ) return 1;          //lepton from heavy flavor decay
    if(cMesonInChain(decayChain) || cBaryonInChain(decayChain) ) return 2;          //lepton from c flavor decay
    if(bosonInChain(decayChain) ) return 0;                                         //lepton from boson
    if(!decayChain.empty()) return 3;                                               //light flavor fake
    return 4;                                                                       //unkown origin
}
