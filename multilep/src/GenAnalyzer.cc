//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

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
    outputTree->Branch("_gen_HT",                    &_gen_HT,                    "_gen_HT/D");

  //jet stuff
  outputTree->Branch("_gen_nN",		               &_gen_nN,			        "_gen_nN/i");
  outputTree->Branch("_gen_NPt",		           &_gen_NPt,			        "_gen_NPt/D");
  outputTree->Branch("_gen_NEta",		           &_gen_NEta,			        "_gen_NEta/D");
  outputTree->Branch("_gen_NPhi",		           &_gen_NPhi,			        "_gen_NPhi/D");
  outputTree->Branch("_gen_NE",		   	           &_gen_NE,			        "_gen_NE/D");
  outputTree->Branch("_gen_Nvertex_x",		   	   &_gen_Nvertex_x,			    "_gen_Nvertex_x/D");
  outputTree->Branch("_gen_Nvertex_y",		   	   &_gen_Nvertex_y,			    "_gen_Nvertex_y/D");
  outputTree->Branch("_gen_Nvertex_z",		   	   &_gen_Nvertex_z,			    "_gen_Nvertex_z/D");

  outputTree->Branch("_gen_nNPackedDtrs",		   &_gen_nNPackedDtrs,			"_gen_nNPackedDtrs/i");
  outputTree->Branch("_gen_NPackedDtrsPt",		   &_gen_NPackedDtrsPt,			"_gen_NPackedDtrsPt[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrsEta",		   &_gen_NPackedDtrsEta,		"_gen_NPackedDtrsEta[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrsPhi",		   &_gen_NPackedDtrsPhi,		"_gen_NPackedDtrsPhi[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrsE",		   &_gen_NPackedDtrsE,			"_gen_NPackedDtrsE[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrsPdgId",      &_gen_NPackedDtrsPdgId,	    "_gen_NPackedDtrsPdgId[_gen_nNPackedDtrs]/I");
  outputTree->Branch("_gen_NPackedDtrsCharge",     &_gen_NPackedDtrsCharge,	    "_gen_NPackedDtrsCharge[_gen_nNPackedDtrs]/I");
  outputTree->Branch("matches",                    &matches,	                "matches[_gen_nNPackedDtrs]/I");
  outputTree->Branch("_gen_NPackedDtrsmineta",     &_gen_NPackedDtrsmineta,	    "_gen_NPackedDtrsmineta[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrsminphi",     &_gen_NPackedDtrsminphi,	    "_gen_NPackedDtrsminphi[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrsminpt",      &_gen_NPackedDtrsminpt,	    "_gen_NPackedDtrsminpt[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchPt",   &_gen_NPackedDtrs_matchPt,	"_gen_NPackedDtrs_matchPt[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchEta",  &_gen_NPackedDtrs_matchEta,	"_gen_NPackedDtrs_matchEta[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchPhi",  &_gen_NPackedDtrs_matchPhi,	"_gen_NPackedDtrs_matchPhi[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchE",    &_gen_NPackedDtrs_matchE,	"_gen_NPackedDtrs_matchE[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchdxy",  &_gen_NPackedDtrs_matchdxy,	"_gen_NPackedDtrs_matchdxy[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchdz",   &_gen_NPackedDtrs_matchdz,	"_gen_NPackedDtrs_matchdz[_gen_nNPackedDtrs]/D");
  outputTree->Branch("_gen_NPackedDtrs_matchcharge",&_gen_NPackedDtrs_matchcharge,"_gen_NPackedDtrs_matchcharge[_gen_nNPackedDtrs]/I");
  
  outputTree->Branch("_gen_nNdaughters",	       &_gen_nNdaughters,		    "_gen_nNdaughters/i");
  outputTree->Branch("_gen_Ndaughters_pdg",   	   &_gen_Ndaughters_pdg,	    "_gen_Ndaughters_pdg[_gen_nNdaughters]/i");
  outputTree->Branch("_gen_nstatus23",		       &_gen_nstatus23,		        "_gen_nstatus23/i");
  outputTree->Branch("_gen_nstatus23_fromN",	   &_gen_nstatus23_fromN,	    "_gen_nstatus23_fromN/i");
  outputTree->Branch("_gen_nstatus23_fromW",	   &_gen_nstatus23_fromW,	    "_gen_nstatus23_fromW/i");
  outputTree->Branch("_gen_status23_pdg",	       &_gen_status23_pdg,		    "_gen_status23_pdg[_gen_nstatus23]/I");
  outputTree->Branch("_gen_status23_fromN_pdg",    &_gen_status23_fromN_pdg, 	"_gen_status23_fromN_pdg[_gen_nstatus23_fromN]/i");
  outputTree->Branch("_gen_status23_fromW_pdg",    &_gen_status23_fromW_pdg, 	"_gen_status23_fromW_pdg[_gen_nstatus23_fromW]/i");
  
  outputTree->Branch("_gen_nq",		               &_gen_nq,			        "_gen_nq/i");
  outputTree->Branch("_gen_qPt",		           &_gen_qPt,			        "_gen_qPt[_gen_nq]/D");
  outputTree->Branch("_gen_qEta",		           &_gen_qEta,			        "_gen_qEta[_gen_nq]/D");
  outputTree->Branch("_gen_qPhi",		           &_gen_qPhi,			        "_gen_qPhi[_gen_nq]/D");
  outputTree->Branch("_gen_qE",		   	           &_gen_qE,			        "_gen_qE[_gen_nq]/D");
 
  outputTree->Branch("_gen_nq1dtr",		           &_gen_nq1dtr,		        "_gen_nq1dtr/i");
  outputTree->Branch("_gen_q1dtr_status",	       &_gen_q1dtr_status,		    "_gen_q1dtr_status[_gen_nq1dtr]/I");
  outputTree->Branch("_gen_q1dtr_pdgid",	       &_gen_q1dtr_pdgid,		    "_gen_q1dtr_pdgid[_gen_nq1dtr]/I");
  outputTree->Branch("_gen_q1dtr_Pt",		       &_gen_q1dtr_Pt,		        "_gen_q1dtr_Pt[_gen_nq1dtr]/D");
  outputTree->Branch("_gen_q1dtr_Eta",		       &_gen_q1dtr_Eta,		        "_gen_q1dtr_Eta[_gen_nq1dtr]/D");
  outputTree->Branch("_gen_q1dtr_Phi",		       &_gen_q1dtr_Phi,		        "_gen_q1dtr_Phi[_gen_nq1dtr]/D");
  outputTree->Branch("_gen_q1dtr_E",		       &_gen_q1dtr_E,		        "_gen_q1dtr_E[_gen_nq1dtr]/D");
  outputTree->Branch("_gen_nq2dtr",		           &_gen_nq2dtr,		        "_gen_nq2dtr/i");
  outputTree->Branch("_gen_q2dtr_status",	       &_gen_q2dtr_status,		    "_gen_q2dtr_status[_gen_nq2dtr]/I");
  outputTree->Branch("_gen_q2dtr_pdgid",	       &_gen_q2dtr_pdgid,		    "_gen_q2dtr_pdgid[_gen_nq2dtr]/I");
  outputTree->Branch("_gen_q2dtr_Pt",		       &_gen_q2dtr_Pt,		        "_gen_q2dtr_Pt[_gen_nq2dtr]/D");
  outputTree->Branch("_gen_q2dtr_Eta",		       &_gen_q2dtr_Eta,		        "_gen_q2dtr_Eta[_gen_nq2dtr]/D");
  outputTree->Branch("_gen_q2dtr_Phi",		       &_gen_q2dtr_Phi,		        "_gen_q2dtr_Phi[_gen_nq2dtr]/D");
  outputTree->Branch("_gen_q2dtr_E",		       &_gen_q2dtr_E,		        "_gen_q2dtr_E[_gen_nq2dtr]/D");
}


void GenAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<reco::GenParticle>> genParticles; iEvent.getByToken(multilepAnalyzer->genParticleToken, genParticles);
    edm::Handle<std::vector<pat::PackedGenParticle>> packedGenParticles; iEvent.getByToken(multilepAnalyzer->packedGenParticleToken, packedGenParticles);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
    //std::cout << "begin genanalyzer" << std::endl;
    if(!genParticles.isValid()) return;
    if(!packedGenParticles.isValid()) return;

    // TODO: when applying overlap for new photon samples: check the pt and eta cuts of the photon
    _ttgEventType = overlapEventType(*genParticles, 13., 3.0); // for TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8
    _zgEventType  = overlapEventType(*genParticles, 15., 2.6); // for ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8

    _gen_nN = 0;
    _gen_nNPackedDtrs = 0;
    _gen_nNdaughters = 0;
    _gen_nstatus23 = 0;
    _gen_nstatus23_fromN = 0;
    _gen_nstatus23_fromW = 0;
    _gen_nq = 0;
    _gen_nq1dtr = 0;
    _gen_nq2dtr = 0;

    _gen_nL = 0;
    _gen_nPh = 0;
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
        if((p.status() == 1 and (absId == 11 or absId == 13)) or (p.status() == 2 and p.isLastCopy() and absId == 15)){
            if(_gen_nL != gen_nL_max){
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
	
        // HNL
        if(abs(p.pdgId()) == 9900012 && p.isLastCopy()){
            if(_gen_nN != 1){
	            _gen_NPt  = p.pt();
	            _gen_NEta = p.eta();
	            _gen_NPhi = p.phi();
	            _gen_NE   = p.energy();
                _gen_Nvertex_x = p.vertex().x();
                _gen_Nvertex_y = p.vertex().y();
                _gen_Nvertex_z = p.vertex().z();
                ++_gen_nN;
            }

            // get info on packed genparticle daughters from HNL ( = all status 1 genparticles)
            for(const pat::PackedGenParticle& packed : *packedGenParticles){
                const reco::Candidate * motherInPrunedCollection = packed.mother(0);
                if(packed.pt() > 0.4 && motherInPrunedCollection != nullptr && isAncestor( &p , motherInPrunedCollection)){
                    _gen_NPackedDtrsPt[_gen_nNPackedDtrs]       = packed.pt();
                    _gen_NPackedDtrsEta[_gen_nNPackedDtrs]      = packed.eta();
                    _gen_NPackedDtrsPhi[_gen_nNPackedDtrs]      = packed.phi();
                    _gen_NPackedDtrsE[_gen_nNPackedDtrs]        = packed.energy();
                    _gen_NPackedDtrsPdgId[_gen_nNPackedDtrs]    = packed.pdgId();
                    _gen_NPackedDtrsCharge[_gen_nNPackedDtrs]   = packed.charge();
                    //find track in PF jet matching to this genparticle, to get its information
                    getMatchingPackedPFCandidateInfo(packed, packedCands);                 
                    _gen_nNPackedDtrs++;
                }
            }
        }
	    // daughters of HNL
        if(abs(GenTools::getMotherPdgId(p, *genParticles)) == 9900012){
	        _gen_Ndaughters_pdg[_gen_nNdaughters] = abs(p.pdgId()); 
	        ++_gen_nNdaughters;
	    }
	    // hard process (status 23)
        int mompdgid = GenTools::getMotherPdgId(p, *genParticles);
        // only hard scatter:
        if(p.status() == 23){
            _gen_status23_pdg[_gen_nstatus23] = abs(p.pdgId()); 
            ++_gen_nstatus23;
            if(abs(mompdgid) == 9900012){
                _gen_status23_fromN_pdg[_gen_nstatus23_fromN] = abs(p.pdgId());
                ++_gen_nstatus23_fromN;
            }
            if(abs(mompdgid) == 24){
                _gen_status23_fromW_pdg[_gen_nstatus23_fromW] = abs(p.pdgId());
                ++_gen_nstatus23_fromW;
            }
            // quarks
	        if(abs(mompdgid) == 9900012 || abs(mompdgid) == 24){
                if(abs(p.pdgId()) >= 1 && abs(p.pdgId()) <= 6){
                    _gen_qPt[_gen_nq]  = p.pt();
                    _gen_qEta[_gen_nq] = p.eta();
                    _gen_qPhi[_gen_nq] = p.phi();
                    _gen_qE[_gen_nq]   = p.energy();
                    ++_gen_nq;
                }
	        }
        }


        // Daughters of quarks      DELETE THESE IN NEXT ITERATION IF IM SURE THAT I DONT NEED THEM
        std::vector<reco::GenParticle> daughterList1 = {};
        std::vector<reco::GenParticle> daughterList2 = {};
        std::vector<int> chain_ends1 = {};
        std::vector<int> chain_ends2 = {};
        if(p.status() == 23 && abs(p.pdgId()) >= 1 && abs(p.pdgId()) <= 6 && abs(mompdgid) == 9900012){ // find 2 quarks from HNL decay
          if(_gen_nq == 1){
            getDaughterList(p, *genParticles, daughterList1, chain_ends1);
            removeDoubleCountedDaughters(daughterList1);
            for(auto daughters : daughterList1){
              //if(daughters.status() != 1) continue;
              _gen_q1dtr_status[_gen_nq1dtr] = daughters.status();
              _gen_q1dtr_pdgid[_gen_nq1dtr]  = daughters.pdgId();
              _gen_q1dtr_Pt[_gen_nq1dtr] 	 = daughters.pt();
              _gen_q1dtr_Eta[_gen_nq1dtr] 	 = daughters.eta();
              _gen_q1dtr_Phi[_gen_nq1dtr] 	 = daughters.phi();
              _gen_q1dtr_E[_gen_nq1dtr] 	 = daughters.energy();
              _gen_nq1dtr++;
            }
          }else if(_gen_nq == 2){
            getDaughterList(p, *genParticles, daughterList2, chain_ends2);
            removeDoubleCountedDaughters(daughterList2);
            for(auto daughters : daughterList2){
              //if(daughters.status() != 1) continue;
              _gen_q2dtr_status[_gen_nq2dtr] = daughters.status();
              _gen_q2dtr_pdgid[_gen_nq2dtr]  = daughters.pdgId();
              _gen_q2dtr_Pt[_gen_nq2dtr] 	 = daughters.pt();
              _gen_q2dtr_Eta[_gen_nq2dtr] 	 = daughters.eta();
              _gen_q2dtr_Phi[_gen_nq2dtr] 	 = daughters.phi();
              _gen_q2dtr_E[_gen_nq2dtr] 	 = daughters.energy();
              _gen_nq2dtr++;
            }
          }
          //print statements to analyze daughters
          /*std::cout << "Pdg Id " << p.pdgId() << " pt: " << p.pt() << std::endl;
          for(unsigned int i = 0; i < p.numberOfDaughters(); i++){
            std::cout << "Daughter " << (*genParticles)[p.daughterRef(i).key()].pdgId() << " pt: " << (*genParticles)[p.daughterRef(i).key()].pt() << std::endl;
          }
          if(_gen_nq == 1){
            std::cout << "number of daughters: " << daughterList1.size() << std::endl;
            for(unsigned int i = 0; i < daughterList1.size(); i++){
              if(chain_ends1[i] == 0) std::cout << daughterList1[i].pdgId() << " -> ";
              if(chain_ends1[i] == 1) std::cout << daughterList1[i].pdgId() << " | ";
            } std::cout << std::endl;
            for(unsigned int i = 0; i < daughterList1.size(); i++){
              if(chain_ends1[i] == 0) std::cout << daughterList1[i].status() << " -> ";
              if(chain_ends1[i] == 1) std::cout << daughterList1[i].status() << " | ";
            } std::cout << std::endl;
            for(unsigned int i = 0; i < daughterList1.size(); i++){
              if(chain_ends1[i] == 0) std::cout << daughterList1[i].pt() << " -> ";
              if(chain_ends1[i] == 1) std::cout << daughterList1[i].pt() << " | ";
            } std::cout << std::endl;
          }else if(_gen_nq == 2){
            std::cout << "number of daughters: " << daughterList2.size() << std::endl;
            for(unsigned int i = 0; i < daughterList2.size(); i++){
              if(chain_ends2[i] == 0) std::cout << daughterList2[i].pdgId() << " -> ";
              if(chain_ends2[i] == 1) std::cout << daughterList2[i].pdgId() << " | ";
            } std::cout << std::endl;
            for(unsigned int i = 0; i < daughterList2.size(); i++){
              if(chain_ends2[i] == 0) std::cout << daughterList2[i].status() << " -> ";
              if(chain_ends2[i] == 1) std::cout << daughterList2[i].status() << " | ";
            } std::cout << std::endl;
            for(unsigned int i = 0; i < daughterList2.size(); i++){
              if(chain_ends2[i] == 0) std::cout << daughterList2[i].pt() << " -> ";
              if(chain_ends2[i] == 1) std::cout << daughterList2[i].pt() << " | ";
            } std::cout << std::endl << std::endl;
          }*/
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
    //std::cout << "end genanalyzer" << std::endl;
}


void GenAnalyzer::getDaughterList(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles, std::vector<reco::GenParticle>& list, std::vector<int>& chain_ends)
{
  if(list.empty() or p.pdgId() != list.back().pdgId()){
    list.push_back(p);
  }
  int n = p.numberOfDaughters();
  chain_ends.push_back(check_for_daughter(p, genParticles));
  for(int i = 0; i < n; i++){
    getDaughterList(genParticles[p.daughterRef(i).key()], genParticles, list, chain_ends);
  }
}

void GenAnalyzer::removeDoubleCountedDaughters(std::vector<reco::GenParticle>& list)
{
  for(unsigned int i = 1; i < list.size() - 1; i++){
    for(unsigned int j = i+1; j < list.size(); j++){
      if(((double)list[i].pt() - (double)list[j].pt()) < 0.0001){ list.erase(list.begin() + j); j--; }
    }
  }
}

int GenAnalyzer::check_for_daughter(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles){
  //returns 0 if particle is not the end of a chain of daughters, 1 if it is
  int n = p.numberOfDaughters();
  if(n > 1) return 0;
  else if(n == 1){
    if(p.pdgId() == genParticles[p.daughterRef(0).key()].pdgId()) return check_for_daughter(genParticles[p.daughterRef(0).key()],genParticles);
    else return 0;
  }
  else return 1; // if n == 0 this is the end of the chain
}


bool GenAnalyzer::inMotherList(std::vector<int>& list, int i){
    return (std::find(list.begin(), list.end(), i) != list.end());
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

bool GenAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
    //particle is already the ancestor
    if(ancestor == particle ) return true;

    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0;i< particle->numberOfMothers();i++){
        if(isAncestor(ancestor,particle->mother(i))) return true;
    }

    //if we did not return yet, then particle and ancestor are not relatives
    return false;
}


void GenAnalyzer::getMatchingPackedPFCandidateInfo(const pat::PackedGenParticle &packed, edm::Handle<std::vector<pat::PackedCandidate>>& packedCands){
    //std::vector<pat::PackedCandidate> packedCandidates;
    
    matches[_gen_nNPackedDtrs]                      = 0;
    _gen_NPackedDtrsmineta[_gen_nNPackedDtrs]       = 0.5;
    _gen_NPackedDtrsminphi[_gen_nNPackedDtrs]       = 0.5;
    _gen_NPackedDtrsminpt[_gen_nNPackedDtrs]        = 2;
    _gen_NPackedDtrs_matchPt[_gen_nNPackedDtrs]     = 0;
    _gen_NPackedDtrs_matchEta[_gen_nNPackedDtrs]    = 0;
    _gen_NPackedDtrs_matchPhi[_gen_nNPackedDtrs]    = 0;
    _gen_NPackedDtrs_matchE[_gen_nNPackedDtrs]      = 0;
    _gen_NPackedDtrs_matchdxy[_gen_nNPackedDtrs]    = 0;
    _gen_NPackedDtrs_matchdz[_gen_nNPackedDtrs]     = 0;
    _gen_NPackedDtrs_matchcharge[_gen_nNPackedDtrs] = 0;

    for(auto cand = packedCands->cbegin(); cand != packedCands->cend(); ++cand){
        if(cand->charge() != 0 and fabs(packed.eta() - cand->eta()) < 0.05 and fabs(packed.phi() - cand->phi()) < 0.05 and fabs(packed.pt() - cand->pt()) < 2){
            if(fabs(packed.eta() - cand->eta()) < _gen_NPackedDtrsmineta[_gen_nNPackedDtrs] and (fabs(packed.phi() - cand->phi()) < _gen_NPackedDtrsminphi[_gen_nNPackedDtrs] or fabs(packed.pt() - cand->pt()) < _gen_NPackedDtrsminpt[_gen_nNPackedDtrs])){
                _gen_NPackedDtrsminpt[_gen_nNPackedDtrs]        = fabs(packed.pt() - cand->pt());
                _gen_NPackedDtrsmineta[_gen_nNPackedDtrs]       = fabs(packed.eta() - cand->eta());
                _gen_NPackedDtrsminphi[_gen_nNPackedDtrs]       = fabs(packed.phi() - cand->phi());
                _gen_NPackedDtrs_matchPt[_gen_nNPackedDtrs]     = cand->pt();
                _gen_NPackedDtrs_matchEta[_gen_nNPackedDtrs]    = cand->eta();
                _gen_NPackedDtrs_matchPhi[_gen_nNPackedDtrs]    = cand->phi();
                _gen_NPackedDtrs_matchE[_gen_nNPackedDtrs]      = cand->energy();
                _gen_NPackedDtrs_matchdxy[_gen_nNPackedDtrs]    = cand->dxy();
                _gen_NPackedDtrs_matchdz[_gen_nNPackedDtrs]     = cand->dz();
                _gen_NPackedDtrs_matchcharge[_gen_nNPackedDtrs] = cand->charge();
            }
            matches[_gen_nNPackedDtrs]++;
        }
    }
} 
