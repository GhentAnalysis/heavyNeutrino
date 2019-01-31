//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//include ROOT classes
#include "TLorentzVector.h"

//include other parts of code 
#include "heavyNeutrino/multilep/interface/GenAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"

/*
 * Storing generator particles
 */

GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer){};

void GenAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_ttgEventType",              &_ttgEventType,              "_ttgEventType/i");
    outputTree->Branch("_zgEventType",               &_zgEventType,               "_zgEventType/i");
    outputTree->Branch("_gen_met",                   &_gen_met,                   "_gen_met/D");
    outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
    outputTree->Branch("_gen_nPh",                   &_gen_nPh,                   "_gen_nPh/i");
    outputTree->Branch("_gen_phStatus",              &_gen_phStatus,              "_gen_phStatus[_gen_nPh]/i");
    outputTree->Branch("_gen_phPt",                  &_gen_phPt,                  "_gen_phPt[_gen_nPh]/D");
    outputTree->Branch("_gen_phEta",                 &_gen_phEta,                 "_gen_phEta[_gen_nPh]/D");
    outputTree->Branch("_gen_phPhi",                 &_gen_phPhi,                 "_gen_phPhi[_gen_nPh]/D");
    outputTree->Branch("_gen_phE",                   &_gen_phE,                   "_gen_phE[_gen_nPh]/D");
    outputTree->Branch("_gen_phMomPdg",              &_gen_phMomPdg,              "_gen_phMomPdg[_gen_nPh]/I");
    outputTree->Branch("_gen_phIsPrompt",            &_gen_phIsPrompt,            "_gen_phIsPrompt[_gen_nPh]/O");
    outputTree->Branch("_gen_phMinDeltaR",           &_gen_phMinDeltaR,           "_gen_phMinDeltaR[_gen_nPh]/D");
    outputTree->Branch("_gen_phPassParentage",       &_gen_phPassParentage,       "_gen_phPassParentage[_gen_nPh]/O");
    outputTree->Branch("_gen_nL",                    &_gen_nL,                    "_gen_nL/i");
    outputTree->Branch("_gen_pdgID",                 &_gen_pdgID,                 "_gen_pdgID[_gen_nL]/D");
    outputTree->Branch("_gen_lPt",                   &_gen_lPt,                   "_gen_lPt[_gen_nL]/D");
    outputTree->Branch("_gen_lEta",                  &_gen_lEta,                  "_gen_lEta[_gen_nL]/D");
    outputTree->Branch("_gen_lPhi",                  &_gen_lPhi,                  "_gen_lPhi[_gen_nL]/D");
    outputTree->Branch("_gen_lE",                    &_gen_lE,                    "_gen_lE[_gen_nL]/D");
    outputTree->Branch("_gen_lFlavor",               &_gen_lFlavor,               "_gen_lFlavor[_gen_nL]/i");
    outputTree->Branch("_gen_lCharge",               &_gen_lCharge,               "_gen_lCharge[_gen_nL]/I");
    outputTree->Branch("_gen_lMomPdg",               &_gen_lMomPdg,               "_gen_lMomPdg[_gen_nL]/I");
    outputTree->Branch("_gen_vertex_x",              &_gen_vertex_x,              "_gen_vertex_x[_gen_nL]/D");
    outputTree->Branch("_gen_vertex_y",              &_gen_vertex_y,              "_gen_vertex_y[_gen_nL]/D");
    outputTree->Branch("_gen_vertex_z",              &_gen_vertex_z,              "_gen_vertex_z[_gen_nL]/D");
    outputTree->Branch("_gen_lIsPrompt",             &_gen_lIsPrompt,             "_gen_lIsPrompt[_gen_nL]/O");
    outputTree->Branch("_gen_lMinDeltaR",            &_gen_lMinDeltaR,            "_gen_lMinDeltaR[_gen_nL]/D");
    outputTree->Branch("_gen_lPassParentage",        &_gen_lPassParentage,        "_gen_lPassParentage[_gen_nL]/O");
}

void GenAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<reco::GenParticle>> genParticles; iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);

    if(!genParticles.isValid()) return;

    // TODO: when applying overlap for new photon samples: check the pt and eta cuts of the photon
    _ttgEventType = overlapEventType(*genParticles, 13., 3.0); // for TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8
    _zgEventType  = overlapEventType(*genParticles, 15., 2.6); // for ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8

    _gen_nL = 0;
    _gen_nPh = 0;
    for(unsigned ig=0; ig<gen_nL_max ; ++ig) _gen_lRefs[ig]  = nullptr;
    TLorentzVector genMetVector(0,0,0,0);
    for(const reco::GenParticle& p : *genParticles){
        int absId = abs(p.pdgId());

        //Calculate generator level MET
        if(p.status() == 1){
            if(absId == 12 or absId == 14 or absId == 16 or absId == 1000022){
                TLorentzVector nuVect;
                nuVect.SetPtEtaPhiE(p.pt(), p.eta(), p.phi(), p.energy());
                genMetVector += nuVect;
            }
        }

        //store generator level lepton info
        if((p.status() == 1 and (absId == 11 or absId == 13)) and (p.status() == 2 and p.isLastCopy() and absId == 15)){
    //    if(p.status() == 1 || (p.status() == 2 && p.isLastCopy() && std::abs(p.pdgId()) == 15)){  // TODO: confusing, strange code in displaced branch; if there is a good reason to use it instead of the above, make a good comment
    //        if(p.pdgId()== 2212) continue;
            if(_gen_nL != gen_nL_max){
                _gen_lRefs[_gen_nL]          = &p;
                _gen_pdgID[_gen_nL]          = p.pdgId();
                _gen_lPt[_gen_nL]            = p.pt();
                _gen_lEta[_gen_nL]           = p.eta();
                _gen_lPhi[_gen_nL]           = p.phi();
                _gen_lE[_gen_nL]             = p.energy();
                _gen_lCharge[_gen_nL]        = p.charge();
                _gen_lIsPrompt[_gen_nL]      = GenTools::isPrompt(p, *genParticles);
                _gen_lMomPdg[_gen_nL]        = GenTools::getMother(p, *genParticles)->pdgId();
                _gen_lMinDeltaR[_gen_nL]     = GenTools::getMinDeltaR(p, *genParticles);
                _gen_lPassParentage[_gen_nL] = GenTools::passParentage(p, *genParticles);

                _gen_vertex_x[_gen_nL]       = p.vertex().x();
                _gen_vertex_y[_gen_nL]       = p.vertex().y();
                _gen_vertex_z[_gen_nL]       = p.vertex().z();

                if(absId == 11)      _gen_lFlavor[_gen_nL] = 0;
                else if(absId == 13) _gen_lFlavor[_gen_nL] = 1;
                else                 _gen_lFlavor[_gen_nL] = 2;
                ++_gen_nL;
            }
        }

        //store generator level photon info
        if((p.status() == 1 or p.status() == 71) and absId == 22){
            if(_gen_nPh != gen_nPh_max){
                _gen_phStatus[_gen_nPh]        = p.status();
                _gen_phPt[_gen_nPh]            = p.pt();
                _gen_phEta[_gen_nPh]           = p.eta();
                _gen_phPhi[_gen_nPh]           = p.phi();
                _gen_phE[_gen_nPh]             = p.energy();
                _gen_phIsPrompt[_gen_nPh]      = p.isPromptFinalState();
                _gen_phMomPdg[_gen_nPh]        = GenTools::getMother(p, *genParticles)->pdgId();
                _gen_phMinDeltaR[_gen_nPh]     = GenTools::getMinDeltaR(p, *genParticles);
                _gen_phPassParentage[_gen_nPh] = GenTools::passParentage(p, *genParticles);
                ++_gen_nPh;
            } 
        }
    }
    _gen_met    = genMetVector.Pt();
    _gen_metPhi = genMetVector.Phi();
}



/*
 * Some event categorization in order to understand/debug/apply overlap removal for TTGamma <--> TTJets and similar photon samples
 */
unsigned GenAnalyzer::overlapEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut) const{
    int type = 0;
    for(auto p = genParticles.cbegin(); p != genParticles.cend(); ++p){
        if(p->status()<0)         continue;
        if(p->pdgId()!=22)        continue;
        type = std::max(type, 1);                                                            // Type 1: final state photon found in genparticles with generator level cuts
        if(p->pt()<ptCut)         continue;
        if(fabs(p->eta())>etaCut) continue;
        type = std::max(type, 2);                                                            // Type 2: photon from pion or other meson

        if(GenTools::getMinDeltaR(*p, genParticles) < 0.2) continue;
        if(not GenTools::passParentage(*p, genParticles))  continue;

        // Everything below is *signal*
        std::set<int> decayChain;
        GenTools::setDecayChain(*p, genParticles, decayChain);
        const reco::GenParticle* mom = GenTools::getMother(*p, genParticles);
        if(std::any_of(decayChain.cbegin(), decayChain.cend(), [](const int entry){ return abs(entry) == 24;})){
            if(abs(mom->pdgId()) == 24)     type = std::max(type, 6);      // Type 6: photon directly from W or decay products which are part of ME
            else if(abs(mom->pdgId()) <= 6) type = std::max(type, 4);      // Type 4: photon from quark from W (photon from pythia, rarely)
            else                            type = std::max(type, 5);      // Type 5: photon from lepton from W (photon from pythia)
        } else {
            if(abs(mom->pdgId()) == 6)      type = std::max(type, 7);      // Type 7: photon from top
            else if(abs(mom->pdgId()) == 5) type = std::max(type, 3);      // Type 3: photon from b
            else                            type = std::max(type, 8);      // Type 8: photon from ME
        }
    }

    return type;
}

unsigned GenAnalyzer::getGenLeptonIndex(const reco::GenParticle* match){
  if(match == nullptr) return gen_nL_max; // out of range
  unsigned genLindex = 0;
  for(; genLindex<gen_nL_max; ++genLindex){
    if(_gen_lRefs[genLindex]==match) break;
  }
  return genLindex;
}
