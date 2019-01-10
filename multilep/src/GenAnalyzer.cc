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
    outputTree->Branch("_ttgEventType",              &_ttgEventType,              "_ttgEventType/b");
    outputTree->Branch("_zgEventType",               &_zgEventType,               "_zgEventType/b");
    outputTree->Branch("_gen_met",                   &_gen_met,                   "_gen_met/D");
    outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
    outputTree->Branch("_gen_nPh",                   &_gen_nPh,                   "_gen_nPh/b");
    outputTree->Branch("_gen_phStatus",              &_gen_phStatus,              "_gen_phStatus[_gen_nPh]/i");
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
    outputTree->Branch("_gen_lMinDeltaR",            &_gen_lMinDeltaR,            "_gen_lMinDeltaR[_gen_nL]/D");
    outputTree->Branch("_gen_lPassParentage",        &_gen_lPassParentage,        "_gen_lPassParentage[_gen_nL]/O");
    outputTree->Branch("_gen_HT",                    &_gen_HT,                    "_gen_HT/D");
}

void GenAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<reco::GenParticle>> genParticles; iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);

    if(!genParticles.isValid()) return;

    // TODO: when applying overlap for new photon samples: check the pt and eta cuts of the photon
    _ttgEventType = overlapEventType(*genParticles, 13., 3.0); // for TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8
    _zgEventType  = overlapEventType(*genParticles, 15., 2.6); // for ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8

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
            if(_gen_nL != gen_nL_max){
                _gen_lPt[_gen_nL]       = p.pt();
                _gen_lEta[_gen_nL]      = p.eta();
                _gen_lPhi[_gen_nL]      = p.phi();
                _gen_lE[_gen_nL]        = p.energy();
                _gen_lCharge[_gen_nL]   = p.charge();
                _gen_lIsPrompt[_gen_nL] = GenTools::isPrompt(p, *genParticles); //(p.isPromptDecayed() || p.isPromptFinalState());
                _gen_lMomPdg[_gen_nL]   = GenTools::getMother(p, *genParticles)->pdgId();

                std::set<int> decayChain;
                GenTools::setDecayChain(p, *genParticles, decayChain);
                _gen_lMinDeltaR[_gen_nL]     = GenTools::getMinDeltaR(p, *genParticles);
                _gen_lPassParentage[_gen_nL] = !(*(std::max_element(std::begin(decayChain), std::end(decayChain))) > 37 or *(std::min_element(std::begin(decayChain), std::end(decayChain))) < -37);

                if(abs(p.pdgId()) == 11)      _gen_lFlavor[_gen_nL] = 0;
                else if(abs(p.pdgId()) == 13) _gen_lFlavor[_gen_nL] = 1;
                else                          _gen_lFlavor[_gen_nL] = 2;
                ++_gen_nL;
            }
        }

        //store generator level photon info
        if((p.status() == 1 || p.status() == 71) && abs(p.pdgId()) == 22){
            if(_gen_nPh != gen_nPh_max){
                _gen_phStatus[_gen_nPh]        = p.status();
                _gen_phPt[_gen_nPh]            = p.pt();
                _gen_phEta[_gen_nPh]           = p.eta();
                _gen_phPhi[_gen_nPh]           = p.phi();
                _gen_phE[_gen_nPh]             = p.energy();
                _gen_phIsPrompt[_gen_nPh]      = p.isPromptFinalState();
                _gen_phMomPdg[_gen_nPh]        = GenTools::getMother(p, *genParticles)->pdgId();
                _gen_phMinDeltaR[_gen_nPh]     = GenTools::getMinDeltaR(p, *genParticles);
                std::set<int> decayChain;
                GenTools::setDecayChain(p, *genParticles, decayChain);
                _gen_phPassParentage[_gen_nPh] = !(*(std::max_element(std::begin(decayChain), std::end(decayChain))) > 37 or *(std::min_element(std::begin(decayChain), std::end(decayChain))) < -37);
                ++_gen_nPh;
            } 
        }
    }
    _gen_met    = genMetVector.Pt();
    _gen_metPhi = genMetVector.Phi();

    //compute gen HT as the sum of all status 23 partons
    _gen_HT = 0;
    for(const reco::GenParticle& p: *genParticles){
        if(p.status() == 23){
            if(abs(p.pdgId()) > 0 && abs(p.pdgId()) < 7){
                _gen_HT += p.pt();
            }
        }
    }
}



/*
 * Some event categorization in order to understand/debug/apply overlap removal for TTGamma <--> TTJets and similar
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
/*
        std::set<int> decayChain;
        GenTools::setDecayChain(*p, genParticles, decayChain);
        if(*(std::max_element(std::begin(decayChain), std::end(decayChain))) > 37)  continue;
        if(*(std::min_element(std::begin(decayChain), std::end(decayChain))) < -37) continue;
*/
        std::vector<int> decayChain;
        GenTools::setDecayChainVector(*p, genParticles, decayChain);
        if(*(std::max_element(std::begin(decayChain), std::end(decayChain))) > 37)      continue;
        if(*(std::min_element(std::begin(decayChain), std::end(decayChain))) < -37)     continue;

        bool fromDecayGluon = false;
        for(auto d : decayChain){
          if(d==21 and not GenTools::parentGluonIsIncoming(decayChain)) fromDecayGluon = true;
        }
        if(fromDecayGluon){
          for(auto d : decayChain) std::cout << d << " ";
          std::cout << "\t" << GenTools::parentGluonIsIncoming(decayChain) << std::endl;
          continue;
        }

        if(GenTools::getMinDeltaR(*p, genParticles) < 0.2)                          continue;



        // Everything below is *signal*
        const reco::GenParticle* mom = GenTools::getMother(*p, genParticles);
        if(std::find_if(decayChain.cbegin(), decayChain.cend(), [](const int entry) { return abs(entry) == 24; }) != decayChain.cend() ){
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
