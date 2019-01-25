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

void GenMatching::matchGenToReco() {
  LepToGenDrMatchesVector recogenmatches;

  for(auto iele : patElectrons) { individualGenToRecoMatch(iele, recogenmatches); }
  for(auto imuo : patMuons    ) { individualGenToRecoMatch(imuo, recogenmatches); }
  for(auto itau : patTaus     ) { individualGenToRecoMatch(itau, recogenmatches); }

  // Loop on recogenmatches (i.e. on pat::Leptons)
  bool ambig = true;
  size_t niter = 0; // safety escape (to avoid infinite loops)
  while(ambig && niter<50) { // keep "trimming" matches until there are no double-matched GenParticles left
    ++niter;
    ambig = false;
    for(size_t imatch=0; imatch<recogenmatches.size(); ++imatch) {
      size_t imtchsize = recogenmatches[imatch].second.size();
      if(imtchsize==0) continue;
      for(size_t jmatch=imatch+1; jmatch<recogenmatches.size(); ++jmatch) {
        if(imtchsize==0 or recogenmatches[jmatch].second.size()==0) continue; // need to check imtchsize again (can change in between loops)
        // Pick the best gen-matches for leptons i and j (i.e. the last in the vector)
        if(recogenmatches[imatch].second.back().first == recogenmatches[jmatch].second.back().first) {
          ambig = true; // there are still ambiguous matches => do another iteration
          if(recogenmatches[imatch].second.back().second > recogenmatches[jmatch].second.back().second) {
            recogenmatches[imatch].second.pop_back();
            --imtchsize;
          }
          else {
            recogenmatches[jmatch].second.pop_back();
          }
        }
      } // end for(size_t jmatch=imatch+1; jmatch<recogenmatches.size(); ++jmatch)
    } // end for(size_t imatch=0; imatch<recogenmatches.size(); ++imatch)
  } // end while

  if(niter==50){
    std::cout << " *** WARNING[GenMatching]: reached the limit of 50 iterations, with ambig = " << (ambig ? "true" : "false") << " ***" << std::endl;
  }

  for(auto& imatch : recogenmatches) {
    if(imatch.second.size()==0) continue;
    unsigned group;
    if(imatch.second.back().second<0.5)      group = 1; // Group 1.A
    else if(imatch.second.back().second<1.5) group = 2; // Group 1.B
    else if(imatch.second.back().second<2.5) group = 3; // Group 2
    else if(imatch.second.back().second<3.5) group = 4; // Group 3.A
    else                                     group = 5; // Group 3.B
    recogenmatchlist.push_back(std::make_pair(imatch.first, std::make_pair(imatch.second.back().first, group)));
  }

  return;
}

template <typename Lepton> void GenMatching::individualGenToRecoMatch(const Lepton* lep, LepToGenDrMatchesVector& recgenmatches){
  GenDrMatches tmpgenmatches;

  if(lep->genParticle()!=nullptr) {  // match by reference
    recogenmatchlist.push_back(std::make_pair(lep, std::make_pair(lep->genParticle(), 0)));
  }
  else { // match by DR
    TLorentzVector recp4(lep->px(), lep->py(), lep->pz(), lep->energy());
    auto recid = lep->pdgId();
    for(auto&& gp : *genParticles) {
      TLorentzVector genp4(gp.px(), gp.py(), gp.pz(), gp.energy());
      double rgdr = recp4.DeltaR(genp4);
      if(rgdr<0.2) {
        // Skip if already matched by reference
        // -- Match option 1 --
        bool matchedbyref = false;
        for(auto&& imatch : recogenmatchlist) {if(imatch.second.first == &gp) {matchedbyref = true; break;}}
        if(matchedbyref) continue;
        // -- Match option 2 --
        // LGTvec_ci mit = std::find_if(recogenmatchlist.begin(), recogenmatchlist.end(),
        //                  [&](const LepToGenTypeMatch& ltgtm) {return (ltgtm.second.first == &gp);});
        // if(mit==recogenmatchlist.end()) continue;
        //
        if(gp.pdgId()==recid && gp.status()==1){                         // * Group 1: status 1, same PDG ID
          if(std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.2)      rgdr += 0.;  //    * Group 1.A: pT within 20%
          else if(std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.5) rgdr += 1.;  //    * Group 1.B: pT within 50%  (increase DR by 1.0, to give it less priority than group 1.A)
        }
        else if(std::abs(recid)==11 && gp.pdgId()==22 && (gp.isPromptFinalState() || gp.isPromptDecayed())) rgdr += 2.; // * Group 2: photon conversions to electrons (increase DR by 2.0, to give it less priority than group 1)
        else if(gp.pdgId()!=recid && gp.status()==1){                    // * Group 3: status 1, different PDG ID
          if(std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.2)      rgdr += 3.;  //    * Group 3.A: pT within 20% (increase DR by 3.0, to give it lower priority than groups 1 and 2)
          else if(std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.5) rgdr += 4.;  //    * Group 3.B: pT within 50% (increase DR by 4.0, to give it lower priority than groups 1, 2, and 3.A)
        }
        tmpgenmatches.push_back(std::make_pair(&gp, rgdr));
      }
    } // end for(auto&& gp : *genParticles)
  } // end match by DR

  // Now order all the matches by DR -- note that it is always group-1 < group-2 < group-3
  for(size_t imtch=0; imtch<tmpgenmatches.size(); ++imtch) {
    for(size_t jmtch=imtch+1; jmtch<tmpgenmatches.size(); ++jmtch) {
      if(tmpgenmatches[jmtch].second>tmpgenmatches[imtch].second) { // NB: DECREASING DR order!!! (i.e. the best match is the last!!!)
        GenDrMatch auxmatch = tmpgenmatches[imtch];
        tmpgenmatches[imtch] = tmpgenmatches[jmtch];
        tmpgenmatches[jmtch] = auxmatch;
      }
    }
  }

  // Finally, fill recgenmatches
  recgenmatches.push_back(std::make_pair(lep, tmpgenmatches));
  return;
}


const reco::GenParticle* GenMatching::returnGenMatch(const reco::Candidate* reclep, unsigned& mchtype) const {
  for(auto&& imatch : recogenmatchlist) {
    if(imatch.first == reclep) {
      mchtype = imatch.second.second;
      return imatch.second.first;
    }
  }
  mchtype = 6;
  return nullptr;
}


// Fill match variables
template <typename Lepton> void GenMatching::fillMatchingVars(const Lepton& reco, const reco::GenParticle* match, unsigned mtchtype){
  if(match != nullptr) {
    genLindex = 0;
    for(; genLindex<multilepAnalyzer->genAnalyzer->gen_nL_max; ++genLindex) {
      if(multilepAnalyzer->genAnalyzer->_gen_lRefs[genLindex]==match) break;
    }
    matchType               = mtchtype;
  } else {
    genLindex               = multilepAnalyzer->genAnalyzer->gen_nL_max; // out of range
    matchType               = mtchtype==5 ? 7 : 6;
  }
  matchIsPrompt           = match ? isPrompt(reco, *match) : false;
  matchIsPromptFinalState = match ? match->isPromptFinalState(): false;
  matchIsPromptDecayed    = match ? match->isPromptDecayed() : false;

  matchPdgId              = match ? match->pdgId() : 0;
  provenance              = match ? GenTools::provenance(match, *genParticles) : 18;
  provenanceCompressed    = match ? GenTools::provenanceCompressed(match, *genParticles, matchIsPrompt) : 4;
  matchPt                 = match ? match->pt() : 0;
  matchEta                = match ? match->eta() : 0;
  matchPhi                = match ? match->phi() : 0;
  matchXvtx               = match ? match->vertex().x() : 0;
  matchYvtx               = match ? match->vertex().y() : 0;
  matchZvtx               = match ? match->vertex().z() : 0;
}

bool GenMatching::isPrompt(const reco::Candidate& reco, const reco::GenParticle& match) const{
    if(abs(reco.pdgId()) == abs(match.pdgId()) || match.pdgId() == 22) return GenTools::isPrompt(match, *genParticles);
    return false;
}

/********************************************************
   Declare explicitly all the specific instances of
   template <typename Lepton> void GenMatching::fillMatchingVars<Lepton>(const Lepton&, reco::GenParticle const*, unsigned int)
   This is to avoid the following linking error in compilation:
--------------------
LeptonAnalyzer.cc:(.text+0x959d): undefined reference to `void GenMatching::fillMatchingVars<pat::Muon>(pat::Muon const&, reco::GenParticle const*, unsigned int)'
LeptonAnalyzer.cc:(.text+0xa544): undefined reference to `void GenMatching::fillMatchingVars<pat::Electron>(pat::Electron const&, reco::GenParticle const*, unsigned int)'
LeptonAnalyzer.cc:(.text+0xb03c): undefined reference to `void GenMatching::fillMatchingVars<pat::Tau>(pat::Tau const&, reco::GenParticle const*, unsigned int)'
collect2: error: ld returned 1 exit status
--------------------
   This can also be avoided by moving the implementation of fillMatchingVars(...)
   into the header file GenMatching.h. But this conflicts with using multilepAnalyzer
   inside method fillMatchingVars(...) (error: invalid use of incomplete type 'class multilep')
********************************************************/
template void GenMatching::fillMatchingVars<pat::Muon>(pat::Muon const&, reco::GenParticle const*, unsigned);
template void GenMatching::fillMatchingVars<pat::Electron>(pat::Electron const&, reco::GenParticle const*, unsigned);
template void GenMatching::fillMatchingVars<pat::Tau>(pat::Tau const&, reco::GenParticle const*, unsigned);
