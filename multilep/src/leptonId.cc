#include "testAnalyzer/multilep/plugins/multilep.h"

template <class leptonType> void multilep::fillLeptonIdVars(const leptonType& lepton, const std::vector<pat::PackedCandidate>& packedCands, const reco::ConversionCollection& conversions, const reco::BeamSpot& beamspot){
	static double eleMvaCuts[3][3]
	if(lepton.isMuon()){
		_isloose = 
		_isFO = 
		_istight = 
	} else if(lepton.isElectron())
		_isloose = 
		_isFO = 
		_istight = 
	} 
}

template <class leptonType> bool multilep::isLoose(const leptonType& lepton){
	if(fabs(_dxy[_nL]) < 0.05 || fabs(_dz[_nL]) < 0.1) return false;
	if(lepton.isMuon()){
		return (_relIso[_nL] < 0.6) && lepton.isPFMuon() && (lepton.isTrackerMuon() || lepton.isGlobalMuon()) && (lepton.pt > 5);
	} else if(lepton.isElectron()){
		if(iE.gsfTrack().hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) return false;
		if(ConversionTools::hasMatchedConversion(packedCands,_conversions, beamspot.position()) return false;
		return (_relIso[_nL] < 0.6) && (lepton.pt() > 10);
	} else if(lepton.isTau()){
		
	}


			

               // _islooseCut[leptonCounter] = fabs(_lEta[leptonCounter]) < 2.5 && _lPt[leptonCounter] > 10 && fabs(_ipPV[leptonCounter]) < 0.05 && fabs(_ipZPV[leptonCounter]) < 0.1 && _hitsNumber[leptonCounter] < 2 && !_vtxFitConversion[leptonCounter] && _clusterpass[leptonCounter] && HoverEpass && passedMVA[0] && _miniisolation[leptonCounter][0] < 0.4;
				//_islooseCut[leptonCounter] = fabs(_lEta[leptonCounter]) < 2.5 && _lPt[leptonCounter] > 10 && fabs(_ipPV[leptonCounter]) < 0.05 && fabs(_ipZPV[leptonCounter]) < 0.1 && _hitsNumber[leptonCounter] < 2 && !_vtxFitConversion[leptonCounter] && passedMVA[0] && _isolation[leptonCounter] < 0.6;
				//_islooseCut[leptonCounter] = fabs(_lEta[leptonCounter]) < 2.5 && _lPt[leptonCounter] > 10 && fabs(_ipPV[leptonCounter]) < 0.05 && fabs(_ipZPV[leptonCounter]) < 0.1 && _hitsNumber[leptonCounter] < 2 && !_vtxFitConversion[leptonCounter] && passedMVA_hnlFO && _isolation[leptonCounter] < 0.6;	
				_islooseCut[leptonCounter] = fabs(_lEta[leptonCounter]) < 2.5 && _lPt[leptonCounter] > 10 && fabs(_ipPV[leptonCounter]) < 0.05 && fabs(_ipZPV[leptonCounter]) < 0.1 && _hitsNumber[leptonCounter] < 2 && !_vtxFitConversion[leptonCounter] && (_isolation[leptonCounter] < 0.6); //< 0.6;







}

template <class leptonType> bool multilep::isFO(const leptonType& lepton){
	if(!_isloose[_nL]) return false;
	if(lepton.isMuon()){
		return lepton.isMediumMuon() && fabs(_3dIPSig[_nL]) < 4;
	} else if(lepton.isElectron()){
		if(iE.gsfTrack().hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) != 0) return false;
		if(!eleClusterPass(ele)) return false;
		return  fabs(_3dIPSig[_nL]) < 4;
	}

}

template <class leptonType> bool multilep::isTight(const lepttonType& lepton){
	if(!_isFO[_nL]) return false;

}

//Check if electron overlaps with loose muon
bool multilep::eleMuOverlap(const pat::Electron& ele){
	TLorentzVector eleV, muV;
	eleV.SetPtEtaPhiE.(ele.pt(), ele.eta(), ele.phi(), ele.energy());
	for(unsigned m = 0; m < _nMu; ++m){
		if(_isLoose[m]){
			TLorentzVector muV(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
			if(eleV.DeltaR(muV) < 0.05) return true;
		}
	}
	return false;	
}

//Check if electron passes trigger emulation requirements
bool multilep::eleClusterPass(const pat::Electron& ele){
	double etaSC =  electron.superCluster().eta();
	return fabs(electron.deltaEtaSuperClusterTrackAtVtx()) < (0.01 - 0.002*(fabs(etaSC) > 1.479)) && fabs(electron.deltaPhiSuperClusterTrackAtVtx()) < (0.04 + 0.03(fabs(etaSC) > 1.479)) && fabs(ele.full5x5_sigmaIetaIeta()) < (0.011 + 0.019*(fabs(etaSC) > 1.479)) && (1.0/ele.correctedEcalEnergy() - ele.eSuperClusterOverP()/ele.correctedEcalEnergy()) < (0.01 - 0.05*fabs(etaSC) > 1.479)) && (1.0/ele.correctedEcalEnergy() - ele.eSuperClusterOverP()/ele.correctedEcalEnergy()) < -0.05) && ele.hadronicOverEm() < 0.1 - fabs(etaSC > 1.479)*0.03;
}

