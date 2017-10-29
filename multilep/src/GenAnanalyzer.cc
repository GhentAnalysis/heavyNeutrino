#include "heavyNeutrino/multilep/interface/GenAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//TlorentzVector, useful for calculating generator MET
#include "TLorentzVector.h"
/*
 * Storing generator particles
 * Currently storing everything such that we have it available for studies on tuple level
 * Might consider trimming it down to save space
 */


GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer)
{};


void GenAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_ttgEventType",              &_ttgEventType,              "_ttgEventType/b");
  outputTree->Branch("_zgEventType",               &_zgEventType,               "_zgEventType/b");
  outputTree->Branch("_gen_met",                   &_gen_met,                   "_gen_met/D");
  outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
  outputTree->Branch("_gen_nPh",                   &_gen_nPh,                   "_gen_nPh/b");
  outputTree->Branch("_gen_phPt",                  &_gen_phPt,                  "_gen_phPt[_gen_nPh]/D");
  outputTree->Branch("_gen_phEta",                 &_gen_phEta,                 "_gen_phEta[_gen_nPh]/D");
  outputTree->Branch("_gen_phPhi",                 &_gen_phPhi,                 "_gen_phPhi[_gen_nPh]/D");
  outputTree->Branch("_gen_phE",                   &_gen_phE,                   "_gen_phE[_gen_nPh]/D");
  outputTree->Branch("_gen_phMomPdg",              &_gen_phMomPdg,              "_gen_phMomPdg[_gen_nPh]/I");
  outputTree->Branch("_gen_phIsPrompt",            &_gen_phIsPrompt,            "_gen_phIsPrompt[_gen_nPh]/O");
  outputTree->Branch("_gen_phMinDeltaR",           &_gen_phMinDeltaR,           "_gen_phMinDeltaR[_gen_nPh]/D");
  outputTree->Branch("_gen_phPassParentage",       &_gen_phPassParentage,       "_gen_phPassParentage[_gen_nPh]/O");
  outputTree->Branch("_gen_nL",                    &_gen_nL,                    "_gen_nL/b");
  outputTree->Branch("_gen_lPt",                   &_gen_lPt,                   "_gen_lPt[_gen_nL]/D");
  outputTree->Branch("_gen_lEta",                  &_gen_lEta,                  "_gen_lEta[_gen_nL]/D");
  outputTree->Branch("_gen_lPhi",                  &_gen_lPhi,                  "_gen_lPhi[_gen_nL]/D");
  outputTree->Branch("_gen_lE",                    &_gen_lE,                    "_gen_lE[_gen_nL]/D");
  outputTree->Branch("_gen_lFlavor",               &_gen_lFlavor,               "_gen_lFlavor[_gen_nL]/i");
  outputTree->Branch("_gen_lCharge",               &_gen_lCharge,               "_gen_lCharge[_gen_nL]/I");
  outputTree->Branch("_gen_lMomPdg",               &_gen_lMomPdg,               "_gen_lMomPdg[_gen_nL]/I");
  outputTree->Branch("_gen_lIsPrompt",             &_gen_lIsPrompt,             "_gen_lIsPrompt[_gen_nL]/O");

  //jet stuff
  outputTree->Branch("_gen_nW",			           &_gen_nW,			        "_gen_nW/b");
  outputTree->Branch("_gen_WMomPdg",		       &_gen_WMomPdg,		        "_gen_WMomPdg[_gen_nW]/i");
  outputTree->Branch("_gen_nWfromN",		       &_gen_nWfromN,		        "_gen_nWfromN/b");
  outputTree->Branch("_gen_nq",			           &_gen_nq,			        "_gen_nq[6]/i");
  outputTree->Branch("_gen_nN",			           &_gen_nN,			        "_gen_nN/b");
  outputTree->Branch("_gen_nNdaughters",	       &_gen_nNdaughters,		    "_gen_nNdaughters/b");
  outputTree->Branch("_gen_Ndaughters_pdg",   	   &_gen_Ndaughters_pdg,	    "_gen_Ndaughters_pdg[_gen_nNdaughters]/i");
  outputTree->Branch("_gen_nstatus23",		       &_gen_nstatus23,		        "_gen_nstatus23/b");
  outputTree->Branch("_gen_nstatus23_fromNorW",	   &_gen_nstatus23_fromNorW,	"_gen_nstatus23_fromNorW/b");
  outputTree->Branch("_gen_nstatus23_fromNorW",	   &_gen_nstatus23_fromNorW,	"_gen_nstatus23_fromNorW/b");
  outputTree->Branch("_gen_nstatus23_fromN",	   &_gen_nstatus23_fromN,	    "_gen_nstatus23_fromN/b");
  outputTree->Branch("_gen_nstatus23_fromW",	   &_gen_nstatus23_fromW,	    "_gen_nstatus23_fromW/b");
  outputTree->Branch("_gen_status23_pdg",	       &_gen_status23_pdg,		        "_gen_status23_pdg[_gen_nstatus23]/I");
  outputTree->Branch("_gen_status23_fromNorW_mompdg", &_gen_status23_fromNorW_mompdg, "_gen_status23_fromNorW_mompdg[_gen_nstatus23_fromNorW]/i");
  outputTree->Branch("_gen_status23_fromNorW_pdg", &_gen_status23_fromNorW_pdg, "_gen_status23_fromNorW_pdg[_gen_nstatus23_fromNorW]/i");
  outputTree->Branch("_gen_status23_fromN_pdg",    &_gen_status23_fromN_pdg, 	"_gen_status23_fromN_pdg[_gen_nstatus23_fromN]/i");
  outputTree->Branch("_gen_status23_fromW_pdg",    &_gen_status23_fromW_pdg, 	"_gen_status23_fromW_pdg[_gen_nstatus23_fromW]/i");
  outputTree->Branch("_gen_nq23",		   &_gen_nq23,			"_gen_nq23/b");
  outputTree->Branch("_gen_qPt",		   &_gen_qPt,			"_gen_qPt[_gen_nq23]/D");
  outputTree->Branch("_gen_qEta",		   &_gen_qEta,			"_gen_qPt[_gen_nq23]/D");
  outputTree->Branch("_gen_qPhi",		   &_gen_qPhi,			"_gen_qPt[_gen_nq23]/D");
  outputTree->Branch("_gen_qE",		   	   &_gen_qE,			"_gen_qPt[_gen_nq23]/D");

}

void GenAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<reco::GenParticle>> genParticles; iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);

  if(!genParticles.isValid()) return;

  _ttgEventType = ttgEventType(*genParticles, 13., 3.0);
  _zgEventType  = ttgEventType(*genParticles, 10., 2.6);

  _gen_nW = 0;
  _gen_nWfromN = 0;
  _gen_nN = 0;
  _gen_nNdaughters = 0;
  _gen_nstatus23 = 0;
  _gen_nstatus23_fromNorW = 0;
  _gen_nstatus23_fromN = 0;
  _gen_nstatus23_fromW = 0;
  _gen_nq23 = 0;
  int pdgid;
  for(unsigned int i = 0; i < 6; i++) _gen_nq[i] = 0;
  
  _gen_nL = 0;
  _gen_nPh = 0;
  TLorentzVector genMetVector(0,0,0,0);
  for(const reco::GenParticle& p : *genParticles){
    //Calculate generator level MET
    if(p.status() == 1){
      if(abs(p.pdgId()) == 12 || abs(p.pdgId()) == 14 || abs(p.pdgId()) == 16 || (multilepAnalyzer->isSUSY &&  abs(p.pdgId()) == 1000022) ){
        TLorentzVector nuVect;
        nuVect.SetPtEtaPhiE(p.pt(), p.eta(), p.phi(), p.energy());
        genMetVector += nuVect;
      }
    }

    //store generator level lepton info
    if((p.status() == 1 && (abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)) || (p.status() == 2 && p.isLastCopy() && abs(p.pdgId()) == 15)){
      if(_gen_nL == gen_nL_max) break;
      _gen_lPt[_gen_nL]       = p.pt();
      _gen_lEta[_gen_nL]      = p.eta();
      _gen_lPhi[_gen_nL]      = p.phi();
      _gen_lE[_gen_nL]        = p.energy();
      _gen_lCharge[_gen_nL]   = p.charge();
      _gen_lIsPrompt[_gen_nL] = (p.isPromptDecayed() || p.isPromptFinalState());
      _gen_lMomPdg[_gen_nL]   = getMotherPdgId(p, *genParticles);
      if(abs(p.pdgId()) == 11)      _gen_lFlavor[_gen_nL] = 0;
      else if(abs(p.pdgId()) == 13) _gen_lFlavor[_gen_nL] = 1;
      else                          _gen_lFlavor[_gen_nL] = 2;
      ++_gen_nL;
    }

    //store generator level photon info
    if(p.status() == 1 && abs(p.pdgId()) == 22){
      if(_gen_nPh == gen_nPh_max) break;
      std::vector<int> motherList = {};
      getMotherList(p, *genParticles, motherList);
      _gen_phPt[_gen_nPh]            = p.pt();
      _gen_phEta[_gen_nPh]           = p.eta();
      _gen_phPhi[_gen_nPh]           = p.phi();
      _gen_phE[_gen_nPh]             = p.energy();
      _gen_phIsPrompt[_gen_nPh]      = p.isPromptFinalState();
      _gen_phMomPdg[_gen_nPh]        = getMotherPdgId(p, *genParticles);
      _gen_phMinDeltaR[_gen_nPh]     = getMinDeltaR(p, *genParticles);
      _gen_phPassParentage[_gen_nPh] = !(*(std::max_element(std::begin(motherList), std::end(motherList))) > 37 or *(std::min_element(std::begin(motherList), std::end(motherList))) < -37);
      ++_gen_nPh;
    }

    //attempt to store generator level jet/quark info
    if(abs(p.pdgId()) == 24 && p.isLastCopy()){ //Multiple iterations of same W will otherwise be stored!
      _gen_WMomPdg[_gen_nW] = getMotherPdgId(p, *genParticles);
      if(abs(_gen_WMomPdg[_gen_nW]) == 9900012) ++_gen_nWfromN;
      ++_gen_nW;
    }
    for(unsigned i = 0; i < 6; ++i){
        if(abs(p.pdgId()) == (i + 1) ) ++_gen_nq[i];
    } 
    if(abs(p.pdgId()) == 9900012){
      ++_gen_nN;
    }
    if(abs(getMotherPdgId(p, *genParticles)) == 9900012) {_gen_Ndaughters_pdg[_gen_nNdaughters] = abs(p.pdgId()); ++_gen_nNdaughters;}
    //only hard scatter:
    if(p.status() == 23){
      _gen_status23_pdg[_gen_nstatus23] = abs(p.pdgId()); 
      ++_gen_nstatus23;
      pdgid = getMotherPdgId(p, *genParticles);
      if((abs(pdgid) == 9900012 || abs(pdgid) == 24)){
        if(abs(pdgid) == 9900012) _gen_status23_fromNorW_mompdg[_gen_nstatus23_fromNorW] = 30; //What is 30???
        else _gen_status23_fromNorW_mompdg[_gen_nstatus23_fromNorW] = abs(pdgid);
        _gen_status23_fromNorW_pdg[_gen_nstatus23_fromNorW] = abs(p.pdgId());
        ++_gen_nstatus23_fromNorW;
      }
      if(abs(pdgid) == 9900012){
        _gen_status23_fromN_pdg[_gen_nstatus23_fromN] = abs(p.pdgId());
        ++_gen_nstatus23_fromN;
        if(abs(p.pdgId()) >= 1 && abs(p.pdgId()) <= 6){
          _gen_qPt[_gen_nq23]  = p.pt();
          _gen_qEta[_gen_nq23] = p.eta();
          _gen_qPhi[_gen_nq23] = p.phi();
          _gen_qE[_gen_nq23]   = p.energy();
          ++_gen_nq23;
        }
      }
      if(abs(pdgid) == 24){
        _gen_status23_fromW_pdg[_gen_nstatus23_fromW] = abs(p.pdgId());
        ++_gen_nstatus23_fromW;
      }
    }
  }
  _gen_met    = genMetVector.Pt();
  _gen_metPhi = genMetVector.Phi();
}

double GenAnalyzer::getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles){
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


const reco::GenParticle* GenAnalyzer::getMother(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles){
    return (p.numberOfMothers() == 0) ? nullptr : &genParticles[p.motherRef(0).key()];
}

const int GenAnalyzer::getMotherPdgId(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* mom = getMother(p, genParticles);
    if(!mom)                         return 0;
    else if(mom->pdgId()==p.pdgId()) return getMotherPdgId(*mom, genParticles);
    else                             return mom->pdgId();
}

// Make a (recursive) list of ancestors of a particle, taking out copies and protons
void GenAnalyzer::getMotherList(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles, std::vector<int>& list){
  if((list.empty() or p.pdgId() != list.back()) and p.pdgId() != 2212) list.push_back(p.pdgId());
  if(p.numberOfMothers() > 1) getMotherList(genParticles[p.motherRef(1).key()], genParticles, list);
  if(p.numberOfMothers() > 0) getMotherList(genParticles[p.motherRef(0).key()], genParticles, list);
}


bool GenAnalyzer::inMotherList(std::vector<int>& list, int i){
  return (std::find(list.begin(), list.end(), i) != list.end());
}


/*
 * Some event categorization in order to understand/debug/apply overlap removal for TTG <--> TTJets
 */
unsigned GenAnalyzer::ttgEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut){
  int type = 0;
  for(auto p = genParticles.begin(); p != genParticles.end(); ++p){
    if(p->status()<0)         continue;
    if(p->pdgId()!=22)        continue;
    type = std::max(type, 1);                                                            // Type 1: final state photon found in genparticles with generator level cuts
    if(p->pt()<ptCut)         continue;
    if(fabs(p->eta())>etaCut) continue;
    type = std::max(type, 2);                                                            // Type 2: photon from pion or other meson

    std::vector<int> motherList = {};
    getMotherList(*p, genParticles, motherList);

    if(*(std::max_element(std::begin(motherList), std::end(motherList))) > 37)  continue;
    if(*(std::min_element(std::begin(motherList), std::end(motherList))) < -37) continue;

    // Everything below is *signal*
    if(inMotherList(motherList, 24) or inMotherList(motherList, -24)){                   // If a W-boson in ancestry
      if(abs(getMotherPdgId(*p, genParticles)) == 24)     type =std::max(type, 6);       // Type 6: photon directly from W or decay products which are part of ME
      else if(abs(getMotherPdgId(*p, genParticles)) <= 6) type =std::max(type, 4);       // Type 4: photon from quark from W (photon from pythia, rarely)
      else                                                type =std::max(type, 5);       // Type 5: photon from lepton from W (photon from pythia)
    } else {
      if(abs(getMotherPdgId(*p, genParticles)) == 6)      type = std::max(type, 7);      // Type 7: photon from top
      else if(abs(getMotherPdgId(*p, genParticles)) == 5) type = std::max(type, 3);      // Type 3: photon from b
      else                                                type = std::max(type, 8);      // Type 8: photon from ME
    }
  }
  return type;
}
