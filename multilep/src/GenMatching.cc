#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "TLorentzVector.h"
#include <algorithm>

GenMatching::GenMatching(const edm::ParameterSet& iConfig){
    allowMatchToAllIds = iConfig.existsAs<bool>("allowMatchingToAllIds") ? iConfig.getParameter<bool>("allowMatchingToAllIds") : false;
};

void GenMatching::matchGenToReco(const std::vector<reco::GenParticle>& genParticles, std::vector<const pat::Electron*>& patElectrons, std::vector<const pat::Muon*>& patMuons, std::vector<const pat::Tau*>& patTaus){
  recogenmatchlist.clear();

  LepToGenDrMatchesVector recogenmatches;
  for(auto iele : patElectrons) { individualGenToRecoMatch(genParticles, iele, recogenmatches); }
  for(auto imuo : patMuons    ) { individualGenToRecoMatch(genParticles, imuo, recogenmatches); }
  for(auto itau : patTaus     ) { individualGenToRecoMatch(genParticles, itau, recogenmatches); }

  // Loop on recogenmatches (i.e. on pat::Leptons)
  bool ambig = true;
  size_t niter = 0;
  while(ambig) { // keep "trimming" matches until there are no double-matched GenParticles left
    ++niter;
    ambig = false;
    for(size_t imatch=0; imatch<recogenmatches.size(); ++imatch){                                        // loop over leptons
      for(size_t jmatch=imatch+1; jmatch<recogenmatches.size(); ++jmatch){                               // loop over all next leptons
        auto& genMatchesI = recogenmatches[imatch].second;
        auto& genMatchesJ = recogenmatches[jmatch].second;
        if(!genMatchesI.size() or !genMatchesJ.size()) continue;                                         // skip when no matches for lepton i or j (i.e. nothing to compare/resolve here)
        if(genMatchesI.back().first == genMatchesJ.back().first){                                        // If the current best gen-matches for leptons i and j are the same (i.e. the last in the vector)
          ambig = true;                                                                                  // there are still ambiguous matches to the same genparticle => do another iteration
          if(genMatchesI.back().second > genMatchesJ.back().second) genMatchesI.pop_back();              // Keep the one with best deltaR/group measure, delete the other
          else                                                      genMatchesJ.pop_back();
        }
      }
    }
    if(niter==50 and ambig){ // safety escape (to avoid infinite loops --> are infinite loops even possible with the above code?)
      std::cout << " *** WARNING[GenMatching]: reached the limit of 50 iterations" << std::endl;
      break;
    }
  }

  for(auto& imatch : recogenmatches){
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

template <typename Lepton> void GenMatching::individualGenToRecoMatch(const std::vector<reco::GenParticle>& genParticles, const Lepton* lep, LepToGenDrMatchesVector& recgenmatches){
  // If a reference to a genparticle is stored for the lepton use this
  if(lep->genParticle()){
    recogenmatchlist.push_back(std::make_pair(lep, std::make_pair(lep->genParticle(), 0)));
    return;
  }

  // Else use deltaR matching/groups
  GenDrMatches tmpgenmatches;
  TLorentzVector recp4(lep->px(), lep->py(), lep->pz(), lep->energy());
  auto recid = lep->pdgId();
  for(auto& gp : genParticles){ // TODO: maybe you want to exclude some genparticles of matching?? currently we match even to incoming protons
    // Skip if the genparticle is already matched by reference
    if(std::any_of(recogenmatchlist.begin(), recogenmatchlist.end(), [&gp](auto& match){return match.second.first == &gp;})) continue;

    TLorentzVector genp4(gp.px(), gp.py(), gp.pz(), gp.energy());
    double rgdr = recp4.DeltaR(genp4);
    if(rgdr<0.2){
      if(gp.pdgId()==recid && gp.status()==1){                         // * Group 1: status 1, same PDG ID
        if(fabs(1.-(genp4.Pt()/recp4.Pt()))<0.2)      rgdr += 0.;      //    * Group 1.A: pT within 20%
        else if(fabs(1.-(genp4.Pt()/recp4.Pt()))<0.5) rgdr += 1.;      //    * Group 1.B: pT within 50%  (increase DR by 1.0, to give it less priority than group 1.A)
        else continue;
      }
      else if(abs(recid)==11 && gp.pdgId()==22 && (gp.isPromptFinalState() || gp.isPromptDecayed())) rgdr += 2.; // * Group 2: photon conversions to electrons (increase DR by 2.0, to give it less priority than group 1)
      else if(gp.pdgId()!=recid && gp.status()==1){                    // * Group 3: status 1, different PDG ID
        if(fabs(1.-(genp4.Pt()/recp4.Pt()))<0.2)      rgdr += 3.;      //    * Group 3.A: pT within 20% (increase DR by 3.0, to give it lower priority than groups 1 and 2)
        else if(fabs(1.-(genp4.Pt()/recp4.Pt()))<0.5) rgdr += 4.;      //    * Group 3.B: pT within 50% (increase DR by 4.0, to give it lower priority than groups 1, 2, and 3.A)
        else continue;
      }
      else continue;
      tmpgenmatches.push_back(std::make_pair(&gp, rgdr));
    }
  }

  // Now order all the matches by DR -- note that it is always group-1 < group-2 < group-3
  std::sort(tmpgenmatches.begin(), tmpgenmatches.end(), [](auto &left, auto &right){return left.second > right.second;}); // decreasing deltaR order (i.e. the best match is the last)

  // Finally, fill recgenmatches
  recgenmatches.push_back(std::make_pair(lep, tmpgenmatches));
}


const reco::GenParticle* GenMatching::returnGenMatch(const reco::Candidate& reclep, unsigned& mchtype) const {
  for(auto&& imatch : recogenmatchlist) {
    if(&*(imatch.first) == &reclep) {
      mchtype = imatch.second.second;
      return imatch.second.first;
    }
  }
  mchtype = 6;
  return nullptr;
}
