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

  if(niter==50) {
    std::cout << " *** WARNING[GenMatching]: reached the limit of 50 iterations, with ambig = "
	      << (ambig ? "true" : "false") << " ***" << std::endl;
  }

  for(auto& imatch : recogenmatches) {
    if(imatch.second.size()==0) continue;
    //
    // Group 1.A
    if(imatch.second.back().second<0.5)
      recogenmatchlist.push_back(std::make_pair(imatch.first, std::make_pair(imatch.second.back().first, 1)));
    // Group 1.B
    else if(imatch.second.back().second<1.5)
      recogenmatchlist.push_back(std::make_pair(imatch.first, std::make_pair(imatch.second.back().first, 2)));
    // Group 2
    else if(imatch.second.back().second<2.5)
      recogenmatchlist.push_back(std::make_pair(imatch.first, std::make_pair(imatch.second.back().first, 3)));
    // Group 3.A
    else if(imatch.second.back().second<3.5)
      recogenmatchlist.push_back(std::make_pair(imatch.first, std::make_pair(imatch.second.back().first, 4)));
    // Group 3.B
    else
      recogenmatchlist.push_back(std::make_pair(imatch.first, std::make_pair(imatch.second.back().first, 5)));
  }

  // // OLD!
  // // Loop on recogenmatches (i.e. on pat::Leptons)
  // bool ambig = true;
  // size_t niter = 0;
  // while(ambig && niter<50) {
  //   ++niter;
  //   ambig = false;
  //   for(size_t imatch=0; imatch<recogenmatches.size(); ++imatch) {
  //     size_t imtchsize = recogenmatches[imatch].second.size();
  //     if(imtchsize==0) continue;
  //     for(size_t jmatch=imatch+1; jmatch<recogenmatches.size(); ++jmatch) {
  // 	size_t jmtchsize = recogenmatches[jmatch].second.size();
  // 	if(jmtchsize==0) continue;
  // 	size_t ii(0), ji(0);
  // 	// Pick the best gen-match for lepton i
  // 	for(; ii<imtchsize; ++ii) {if(recogenmatches[imatch].second[ii].first!=nullptr) break;}
  // 	if(ii==imtchsize) continue;
  // 	// Pick the best gen-match for lepton j
  // 	for(; ji<jmtchsize; ++ji) {if(recogenmatches[jmatch].second[ji].first!=nullptr) break;}
  // 	if(ji==jmtchsize) continue;
  // 	if(recogenmatches[imatch].second[ii].first == recogenmatches[jmatch].second[ji].first) {
  // 	  iambig = true;
  // 	  f(recogenmatches[imatch].second[ii].second > recogenmatches[jmatch].second[ji].second) {
  // 	    recogenmatches[imatch].second[ii].first = nullptr;
  // 	    recogenmatches[imatch].second[ii].second = 9999.;
  // 	  }
  // 	  else {
  // 	    recogenmatches[jmatch].second[ji].first = nullptr;
  // 	    recogenmatches[jmatch].second[ji].second = 9999.;
  // 	  }
  // 	}
  //     } // end for(size_t jmatch=imatch+1; jmatch<recogenmatches.size(); ++jmatch)
  //   } // end for(size_t imatch=0; imatch<recogenmatches.size(); ++imatch)
  // } // end while


  // std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Matchings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  // for(auto&& imatch : recogenmatchlist) {
  //   printf("  REC: (%+5d, %+5.2f, %+5.2f, %7.2f);  GEN(%+5d, %+5.2f, %+5.2f, %7.2f;  %2d)\n",
  // 	   imatch.first       ->pdgId(), imatch.first       ->eta(), imatch.first       ->phi(), imatch.first       ->pt(), 
  // 	   imatch.second.first->pdgId(), imatch.second.first->eta(), imatch.second.first->phi(), imatch.second.first->pt(), imatch.second.second);
  // }
  // std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  return;
}

template <typename Lepton>
void GenMatching::individualGenToRecoMatch(const Lepton* lep, LepToGenDrMatchesVector& recgenmatches) {
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
	// 			     [&](const LepToGenTypeMatch& ltgtm) {return (ltgtm.second.first == &gp);});
	// if(mit==recogenmatchlist.end()) continue;
	//
	if( gp.pdgId()==recid && gp.status()==1 ) { 	           // * Group 1: status 1, same PDG ID
	  if( std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.2 ) {  	   //    * Group 1.A: status 1, same PDG ID, pT within 20%
	    tmpgenmatches.push_back(std::make_pair(&gp, rgdr));
	  }
	  else if( std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.5) {	   //    * Group 1.B: status 1, same PDG ID, pT within 50%
	    tmpgenmatches.push_back(std::make_pair(&gp, rgdr+1.)); //      (increase DR by 1.0, to give it less priority than group 1.A)
	  }
	}
	else if( std::abs(recid)==11 && gp.pdgId()==22 && 
	    (gp.isPromptFinalState() || gp.isPromptDecayed()) ) {  // * Group 2: photon conversions to electrons 
	  tmpgenmatches.push_back(std::make_pair(&gp, rgdr+2.));   //   (increase DR by 2.0, to give it less priority than group 1)
	}
	else if( gp.pdgId()!=recid && gp.status()==1 ) {	   // * Group 3: status 1, different PDG ID
	  if( std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.2 ) {	   //    * Group 3.A: status 1, different PDG ID, pT within 20%
	    tmpgenmatches.push_back(std::make_pair(&gp, rgdr+3.)); //      (increase DR by 3.0, to give it lower priority than groups 1 and 2)
	  }
	  else if( std::abs(1.-(genp4.Pt()/recp4.Pt()))<0.5 ) {	   //    * Group 3.B: status 1, different PDG ID, pT within 20%
	    tmpgenmatches.push_back(std::make_pair(&gp, rgdr+4.)); //      (increase DR by 4.0, to give it lower priority than groups 1, 2, and 3.A)
	  }
	}
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
//
// (1) to be used with matchGenToReco()
template <typename Lepton> void GenMatching::fillMatchingVars(const Lepton& reco, const reco::GenParticle* match, unsigned mtchtype, const unsigned lepcnt) {
  if(match != nullptr) {
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) " << 0 << " --- " 
    //	      << &reco << " (" << reco.genParticle() << ")" << " - " << match << " - " << mtchtype << std::endl;
    //genLindex = genPhindex = 0;
    genLindex = 0;
    for(; genLindex<multilepAnalyzer->genAnalyzer->gen_nL_max; ++genLindex) {
      //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- genLindex (before): " << genLindex << std::endl;
      if(multilepAnalyzer->genAnalyzer->_gen_lRefs[genLindex]==match) break;
      //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- genLindex (after ): " << genLindex << std::endl;
    }
    // for(; genPhindex<multilepAnalyzer->genAnalyzer->gen_nPh_max; ++genPhindex)
    //     if(multilepAnalyzer->genAnalyzer->_gen_phRefs[genPhindex]==match) break;
    matchType = mtchtype;
    matchIsPrompt = isPrompt(reco, *match);
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after isPrompt: " << matchIsPrompt << std::endl;
    matchIsPromptFinalState = isPromptFinalState(reco, *match);
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after isPromptFinalState: " << matchIsPromptFinalState << std::endl;
    matchIsPromptDecayed = isPromptDecayed(reco, *match);
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after isPromptDecayed: " << matchIsPromptDecayed << std::endl;

    matchPdgId = match->pdgId();
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after pdgId: " << matchPdgId << std::endl;
    provenance = GenTools::provenance(match, *genParticles);
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after provenance: " << provenance << std::endl;
    provenanceCompressed = GenTools::provenanceCompressed(match, *genParticles, matchIsPrompt);
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after provenanceCompressed: " << provenanceCompressed << std::endl;
    matchPt = match->pt();
    matchEta = match->eta();
    matchPhi = match->phi();
    matchXvtx = match->vertex().x();
    matchYvtx = match->vertex().y();
    matchZvtx = match->vertex().z();
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - FULL match) -- after ALL: " << genLindex << std::endl;
  } else {
    //std::cout << " ~~~ DEBUG: ELE (gen - fill - NULL match) " << 0 << " --- " 
    //	      << &reco << " - " << match << " - " << mtchtype << std::endl;
    genLindex = multilepAnalyzer->genAnalyzer->gen_nL_max; // out of range
    // genPhindex = multilepAnalyzer->genAnalyzer->gen_nPh_max; // out of range
    matchType = 6;
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

  // printf("  >>>> %2d REC: (%+5d, %+5.2f, %+5.2f, %7.2f);  GEN(%+5d, %+5.2f, %+5.2f, %7.2f;  %2d;  %d, %d, %d;  %d, %2d) <<<<\n",
  // 	 lepcnt, reco.pdgId(), reco.eta(), reco.phi(), reco.pt(),
  // 	 matchPdgId, matchEta, matchPhi, matchPt, matchType,
  // 	 matchIsPrompt, matchIsPromptFinalState, matchIsPromptDecayed,
  // 	 provenanceCompressed, provenance);
}
//
// (2) to be used with findGenMatch()
template <typename Lepton> void GenMatching::fillMatchingVars(const Lepton& reco) {
  const reco::GenParticle* match = findGenMatch(reco, allowMatchToAllIds);
  if(match != nullptr){
    //genLindex = genPhindex = 0;
    genLindex = 0;
    for(; genLindex<multilepAnalyzer->genAnalyzer->gen_nL_max; ++genLindex)
      if(multilepAnalyzer->genAnalyzer->_gen_lRefs[genLindex]==match) break;
    // for(; genPhindex<multilepAnalyzer->genAnalyzer->gen_nPh_max; ++genPhindex)
    //     if(multilepAnalyzer->genAnalyzer->_gen_phRefs[genPhindex]==match) break;
    matchType = 5;
    matchIsPrompt = isPrompt(reco, *match);
    matchIsPromptFinalState = isPromptFinalState(reco, *match);
    matchIsPromptDecayed = isPromptDecayed(reco, *match);

    matchPdgId = match->pdgId();
    provenance = GenTools::provenance(match, *genParticles);
    provenanceCompressed = GenTools::provenanceCompressed(match, *genParticles, matchIsPrompt);
    matchPt = match->pt();
    matchEta = match->eta();
    matchPhi = match->phi();
    matchXvtx = match->vertex().x();
    matchYvtx = match->vertex().y();
    matchZvtx = match->vertex().z();
  } else{
    genLindex = multilepAnalyzer->genAnalyzer->gen_nL_max; // out of range
    // genPhindex = multilepAnalyzer->genAnalyzer->gen_nPh_max; // out of range
    matchType = 7;
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
    if(minDeltaR > 0.3){
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
bool GenMatching::isPromptFinalState(const reco::Candidate& reco, const reco::GenParticle& match) const{
    return GenTools::isPromptFinalState(match, *genParticles);
}
bool GenMatching::isPromptDecayed(const reco::Candidate& reco, const reco::GenParticle& match) const{
    return GenTools::isPromptDecayed(match, *genParticles);
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
template void GenMatching::fillMatchingVars<pat::Muon>(pat::Muon const&, reco::GenParticle const*, unsigned, unsigned const);
template void GenMatching::fillMatchingVars<pat::Electron>(pat::Electron const&, reco::GenParticle const*, unsigned, unsigned const);
template void GenMatching::fillMatchingVars<pat::Tau>(pat::Tau const&, reco::GenParticle const*, unsigned, unsigned const);

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
