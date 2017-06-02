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
  outputTree->Branch("_nTrueInt",                   &_nTrueInt,                 "_nTrueInt/F");

  outputTree->Branch("_ttgEventType",              &_ttgEventType,              "_ttgEventType/I");
  outputTree->Branch("_gen_met",                   &_gen_met,                   "_gen_met/D");
  outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
  outputTree->Branch("_gen_nPh",                   &_gen_nPh,                   "_gen_nPh/b");
  outputTree->Branch("_gen_phPt",                  &_gen_phPt,                  "_gen_phPt[_gen_nPh]/D");
  outputTree->Branch("_gen_phEta",                 &_gen_phEta,                 "_gen_phEta[_gen_nPh]/D");
  outputTree->Branch("_gen_phPhi",                 &_gen_phPhi,                 "_gen_phPhi[_gen_nPh]/D");
  outputTree->Branch("_gen_phE",                   &_gen_phE,                   "_gen_phE[_gen_nPh]/D");
  outputTree->Branch("_gen_phMomPdg",              &_gen_phMomPdg,              "_gen_phMomPdg[_gen_nPh]/I");
  outputTree->Branch("_gen_phIsPrompt",            &_gen_phIsPrompt,            "_gen_phIsPrompt[_gen_nPh]/O");
  outputTree->Branch("_gen_nL",                    &_gen_nL,                    "_gen_nL/b");
  outputTree->Branch("_gen_lPt",                   &_gen_lPt,                   "_gen_lPt[_gen_nL]/D");
  outputTree->Branch("_gen_lEta",                  &_gen_lEta,                  "_gen_lEta[_gen_nL]/D");
  outputTree->Branch("_gen_lPhi",                  &_gen_lPhi,                  "_gen_lPhi[_gen_nL]/D");
  outputTree->Branch("_gen_lE",                    &_gen_lE,                    "_gen_lE[_gen_nL]/D");
  outputTree->Branch("_gen_lFlavor",               &_gen_lFlavor,               "_gen_lFlavor[_gen_nL]/b");
  outputTree->Branch("_gen_lCharge",               &_gen_lCharge,               "_gen_lCharge[_gen_nL]/I");
  outputTree->Branch("_gen_lMomPdg",               &_gen_lMomPdg,               "_gen_lMomPdg[_gen_nL]/I");
  outputTree->Branch("_gen_lIsPrompt",             &_gen_lIsPrompt,             "_gen_lIsPrompt[_gen_nL]/O");
}

void GenAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<reco::GenParticle>> genParticles; iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);
  edm::Handle<std::vector<PileupSummaryInfo>>  pileUpInfo;  iEvent.getByToken(multilepAnalyzer->pileUpToken,      pileUpInfo);

  _nTrueInt = pileUpInfo->begin()->getTrueNumInteractions(); // getTrueNumInteractions should be the same for all bunch crosssings

  if(!genParticles.isValid()) return;

  _ttgEventType = ttgEventType(*genParticles);

  _gen_nL = 0;
  _gen_nPh = 0;
  TLorentzVector genMetVector;
  for(const reco::GenParticle& p : *genParticles){
      //Calculate generator level MET
      if(p.status() == 1){
          if(abs(p.pdgId()) == 12 || abs(p.pdgId()) == 14 || abs(p.pdgId()) == 16){
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
      else if( p.status() == 1 && abs(p.pdgId()) == 22){
          if(_gen_nPh == gen_nPh_max) break;
          _gen_phPt[_gen_nPh]       = p.pt();
          _gen_phEta[_gen_nPh]      = p.eta();
          _gen_phPhi[_gen_nPh]      = p.phi();
          _gen_phE[_gen_nPh]        = p.energy();
          _gen_phIsPrompt[_gen_nPh] = p.isPromptFinalState();
          _gen_phMomPdg[_gen_nPh]   = getMotherPdgId(p, *genParticles);
          ++_gen_nPh;
      }
  }
  _gen_met    = genMetVector.Pt();
  _gen_metPhi = genMetVector.Phi();
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
int GenAnalyzer::ttgEventType(const std::vector<reco::GenParticle>& genParticles){
  int type = 0;
  for(auto p = genParticles.begin(); p != genParticles.end(); ++p){
    if(p->status()<0)      continue;
    if(p->pdgId()!=22)     continue;
    type = std::max(type, 1);                                                            // Type 1: final state photon found in genparticles
    if(p->pt()<10)         continue;
    if(fabs(p->eta())>2.6) continue;
    type = std::max(type, 2);                                                            // Type 2: photon with generator level cuts

    std::vector<int> motherList = {};
    getMotherList(*p, genParticles, motherList);

    if(*(std::max_element(std::begin(motherList), std::end(motherList))) > 37)  continue;
    if(*(std::min_element(std::begin(motherList), std::end(motherList))) < -37) continue;
    type = std::max(type, 3);                                                            // Type 3: photon probably from pion (or other meson)

    if(inMotherList(motherList, 24) or inMotherList(motherList, -24)){                   // If a W-boson in ancestry
      if(abs(getMotherPdgId(*p, genParticles)) == 24)     type =std::max(type, 6);       // Type 6: photon directly from W
      else if(abs(getMotherPdgId(*p, genParticles)) <= 6) type =std::max(type, 4);       // Type 4: photon from quark from W
      else                                                type =std::max(type, 5);       // Type 5: photon from lepton from W
     
      continue;
    }

    if(abs(getMotherPdgId(*p, genParticles)) == 6) type = std::max(type, 7);             // Type 7: photon from top
    else                                           type = std::max(type, 8);             // Type 8: photon from other quark or gluon or top not recorded (like in the TTGJets sample)
  }
  return type;
}
