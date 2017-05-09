#include "testAnalyzer/multilep/interface/Tools.h"
#include <cmath>
#include <algorithm>

double Tools::getEffArea(const unsigned flavor, const double eta){
		double effA[2][7] = {{0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687}, {0.0735, 0.0619, 0.0465, 0.0433, 0.0577, 0, 0}};
		double etaBins[2][6] = { {1, 1.479, 2.0, 2.2, 2.3, 2.4}, {0.8, 1.3, 2.0, 2.2, 0, 0}};
		unsigned nEtaBins[2] = {7, 5};
		for(unsigned e = 0; e < nEtaBins[flavor]; ++e){
			if(fabs(eta) < etaBins[flavor][e]) return effA[flavor][e];
		}
		return effA[flavor][nEtaBins[flavor]];
}

double Tools::getRelIso03(const pat::Muon& mu, const double rho){ //Note: effective area correction is used instead of delta-beta correction
	double puCorr = rho*getEffArea(1, mu.eta());
	double absIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - puCorr);
	return absIso/mu.pt();
}

//double Tools::getIsoAlt(const pat::Muon& mu, double rho){
	//double puCorr = rho*getEffArea(1, mu.eta());
	//double absIso = mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - puCorr);
	//return absIso/mu.pt();
//}

double Tools::getRelIso03(const pat::Electron& ele, const double rho){
	double puCorr = rho*getEffArea(0, ele.eta());
	double absIso = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - puCorr);
	return absIso/ele.pt();
}

double Tools::getRelIso(const reco::RecoCandidate& ptcl, const std::vector<pat::PackedCandidate>& pfcands, const double coneSize, const double rho){
	if(ptcl.pt() < 5) return 99999.;
	if(!ptcl.isMuon() && !ptcl.isElectron()){
		std::cerr << "ERROR in getRelIso: function only defined for muons and electrons" << std::endl;
		return 99999.;
	}
	double puCorr = rho*getEffArea(ptcl.isMuon(), ptcl.eta());
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
				
double Tools::getMiniIso(const reco::RecoCandidate& ptcl, const std::vector<pat::PackedCandidate>& pfcands, const double maxCone, const double rho){
	double coneSize = std::max(0.05, std::min(maxCone, 10./ptcl.pt()));
	return Tools::getRelIso(ptcl, pfcands, coneSize, rho);
}


