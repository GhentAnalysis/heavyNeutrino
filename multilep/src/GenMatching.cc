#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "heavyNeutrino/multilep/interface/GenParticleManager.h"
#include "TLorentzVector.h" 

GenMatching::GenMatching(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer): multilepAnalyzer(multilepAnalyzer){
    allowMatchToAllIds = iConfig.existsAs<bool>("allowMatchingToAllIds") ? iConfig.getParameter<bool>("allowMatchingToAllIds") : false;
};

void GenMatching::setGenParticles(const edm::Event& iEvent){
    iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);
}

//fill match variables
template <typename Lepton> void GenMatching::fillMatchingVars(const Lepton& reco) {
  const reco::GenParticle* match = findGenMatch(reco, allowMatchToAllIds);
  if(match != nullptr){
    //genLindex = genPhindex = 0;
    genLindex = 0;
    for(; genLindex<multilepAnalyzer->genAnalyzer->gen_nL_max; ++genLindex)
      if(multilepAnalyzer->genAnalyzer->_gen_lRefs[genLindex]==match) break;
    // for(; genPhindex<multilepAnalyzer->genAnalyzer->gen_nPh_max; ++genPhindex)
    //     if(multilepAnalyzer->genAnalyzer->_gen_phRefs[genPhindex]==match) break;
    matchIsPrompt = isPrompt(reco, *match);
    matchIsPromptFinalState = isPromptFinalState(reco, *match);
    matchIsPromptDecayed = isPromptDecayed(reco, *match);

    matchPdgId = match->pdgId();
    provenance = GenTools::provenance(*match, *genParticles);
    provenanceCompressed = (matchIsPrompt ? 0 : GenTools::provenanceCompressed(*match, *genParticles) );
    matchPt = match->pt();
    matchEta = match->eta();
    matchPhi = match->phi();
    matchXvtx = match->vertex().x();
    matchYvtx = match->vertex().y();
    matchZvtx = match->vertex().z();
  } else{
    genLindex = multilepAnalyzer->genAnalyzer->gen_nL_max; // out of range
    // genPhindex = multilepAnalyzer->genAnalyzer->gen_nPh_max; // out of range
    matchIsPrompt = false;
    matchIsPromptFinalState = false;
    matchIsPromptDecayed = false;

    matchPdgId = 0;
    provenanceCompressed = 4;
    provenance = 18;
    matchPt = 0.;
    matchEta = 0.;
    matchPhi = 0.;
    matchXvtx = 0.;
    matchYvtx = 0.;
    matchZvtx = 0.;
  }
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

/********************************************************
   Declare explicitly all the specific instances of
   template <typename Lepton> void GenMatching::fillMatchingVars(const Lepton& reco)
   This is to avoid the following linking error in compilation:
--------------------
tmp/slc6_amd64_gcc630/src/heavyNeutrino/multilep/src/heavyNeutrinomultilep/LeptonAnalyzer.o: In function `LeptonAnalyzer::analyze(edm::Event const&, edm::EventSetup const&, reco::Vertex const&)':
LeptonAnalyzer.cc:(.text+0x696a): undefined reference to `void GenMatching::fillMatchingVars<pat::Muon>(pat::Muon const&)'
LeptonAnalyzer.cc:(.text+0x752f): undefined reference to `void GenMatching::fillMatchingVars<pat::Electron>(pat::Electron const&)'
LeptonAnalyzer.cc:(.text+0x7bd6): undefined reference to `void GenMatching::fillMatchingVars<pat::Tau>(pat::Tau const&)'
collect2: error: ld returned 1 exit status
--------------------
   This can also be avoided by moving the implementation of fillMatchingVars(...) 
   into the header file GenMatching.h. But this conflicts with using multilepAnalyzer 
   inside method fillMatchingVars(...) (error: invalid use of incomplete type 'class multilep')
********************************************************/
template void GenMatching::fillMatchingVars<pat::Muon>(pat::Muon const&);
template void GenMatching::fillMatchingVars<pat::Electron>(pat::Electron const&);
template void GenMatching::fillMatchingVars<pat::Tau>(pat::Tau const&);
