#include "../interface/LeptonAnalyzer.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"


double LeptonAnalyzer::getRelIso03(const pat::Muon& mu, const double rho){ //Note: effective area correction is used instead of delta-beta correction
	double puCorr = rho*muonsEffectiveAreas.getEffectiveArea(mu.eta());
	double absIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - puCorr);
	return absIso/mu.pt();
}

double LeptonAnalyzer::getRelIso03(const pat::Electron& ele, const double rho){
	double puCorr = rho*electronsEffectiveAreas.getEffectiveArea(ele.superCluster()->eta());
	double absIso = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - puCorr);
	return absIso/ele.pt();
}

// This is used for the misolation below
double LeptonAnalyzer::getRelIso(const reco::RecoCandidate& ptcl, const std::vector<pat::PackedCandidate>& pfcands, const double coneSize, const double rho){
	if(ptcl.pt() < 5) return 99999.;
	if(!ptcl.isMuon() && !ptcl.isElectron()){
		std::cerr << "ERROR in getRelIso: function only defined for muons and electrons" << std::endl;
		return 99999.;
	}
	double puCorr;
  if(ptcl.isMuon()) puCorr = rho*muonsEffectiveAreas.getEffectiveArea(ptcl.eta());
  else              puCorr = rho*electronsEffectiveAreas.getEffectiveArea(ptcl.superCluster()->eta());
	double deadCone_nHadr = 0, deadCone_chHadr = 0, deadCone_ph = 0;
	//determine dead cone around particle
	if(ptcl.isElectron()){
		if(fabs(ptcl.eta()) > 1.479){
			deadCone_chHadr = 0.015;
			deadCone_ph = 0.08;
		}
	} else if(ptcl.isMuon()){
		deadCone_chHadr = 0.0001;
		deadCone_nHadr = 0.01;
		deadCone_ph = 0.01;
	} else{
		;
	}
	//pt threshold for isolation objects
	double ptThres = ptcl.isElectron() ? 0. : 0.5;
	//Calculate the different isolation components
	double phIso = 0, nHadrIso = 0, chHadrIso = 0;
	for(const pat::PackedCandidate& cand: pfcands){
		if(abs(cand.pdgId())<7) continue;
		double deltaR = reco::deltaR(ptcl, cand);
		if(deltaR < coneSize){
			if(cand.charge() == 0){
				if(cand.pt() > ptThres){
					if(cand.pdgId() == 130){
						if(deltaR > deadCone_nHadr){
							nHadrIso += cand.pt();
						}
					} else if(cand.pdgId() == 22){
						if(deltaR > deadCone_ph){
							phIso += cand.pt();
						}
					}
				}
			} else{
				if(cand.fromPV() > 1){//1 is recommended by miniAOD workbook
					if(fabs(cand.pdgId()) == 211){
						if(deltaR > deadCone_chHadr){
							chHadrIso += cand.pt();
						}
					}
				}
			}
		}
	}
	return (chHadrIso + std::max(0., nHadrIso + phIso - puCorr))/(ptcl.pt());
}
				
double LeptonAnalyzer::getMiniIso(const reco::RecoCandidate& ptcl, const std::vector<pat::PackedCandidate>& pfcands, const double maxCone, const double rho){
	double coneSize = std::max(0.05, std::min(maxCone, 10./ptcl.pt()));
	return LeptonAnalyzer::getRelIso(ptcl, pfcands, coneSize, rho);
}


