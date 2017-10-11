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

double LeptonAnalyzer::getRelIso04(const pat::Muon& mu, const double rho){ //Note: effective area correction is used instead of delta-beta correction
	double puCorr = rho*muonsEffectiveAreas.getEffectiveArea(mu.eta());
	double absIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - puCorr);
	return absIso/mu.pt();
}

double LeptonAnalyzer::getRelIso03(const pat::Electron& ele, const double rho){
	double puCorr = rho*electronsEffectiveAreas.getEffectiveArea(ele.superCluster()->eta());
	double absIso = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - puCorr);
	return absIso/ele.pt();
}

double LeptonAnalyzer::getMiniIsolation(const reco::RecoCandidate& ptcl, edm::Handle<pat::PackedCandidateCollection> pfcands,
                                        double r_iso_min, double r_iso_max, double kt_scale, double rho, const bool onlyCharged){
    bool deltaBeta   = false;
    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl.isElectron() and fabs(ptcl.superCluster()->eta() >1.479)){ deadcone_ch = 0.015;  deadcone_pu = 0.015; deadcone_ph = 0.08; deadcone_nh = 0;}
    else if(ptcl.isMuon())                                           { deadcone_ch = 0.0001; deadcone_pu = 0.01;  deadcone_ph = 0.01; deadcone_nh = 0.01;}

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh = ptcl.isElectron()? 0. : 0.5;

    double max_pt = kt_scale/r_iso_min;
    double min_pt = kt_scale/r_iso_max;
    double r_iso  = kt_scale/std::max(std::min(ptcl.pt(), max_pt), min_pt);

    for(const pat::PackedCandidate &pfc : *pfcands){
      if(fabs(pfc.pdgId()) < 7) continue;

      double dr = deltaR(pfc, ptcl);
      if(dr > r_iso) continue;

      if(pfc.charge()==0){		                                          						// Neutral
        if(pfc.pt()>ptThresh){
          if(fabs(pfc.pdgId())==22 and dr > deadcone_ph)        iso_ph += pfc.pt();	// Photons
          else if (fabs(pfc.pdgId())==130 and dr > deadcone_nh) iso_nh += pfc.pt();	// Neutral hadrons
        }
      } else if (pfc.fromPV()>1){
        if(fabs(pfc.pdgId())==211 and dr > deadcone_ch) iso_ch += pfc.pt();		      // Charged from PV
      } else if(pfc.pt()>ptThresh and dr > deadcone_pu) iso_pu += pfc.pt();		      // Charged from PU
    }

    double puCorr = 0;
    if(ptcl.isMuon()) puCorr = rho*muonsEffectiveAreas.getEffectiveArea(ptcl.eta());
    else              puCorr = rho*electronsEffectiveAreas.getEffectiveArea(ptcl.superCluster()->eta());

    double iso;
    if(onlyCharged)    iso = iso_ch;
    else if(deltaBeta) iso = iso_ch + std::max(0., iso_ph + iso_nh - 0.5*iso_pu);
    else               iso = iso_ch + std::max(0., iso_ph + iso_nh - puCorr*(r_iso*r_iso)/(0.3*0.3));

    return iso/ptcl.pt();
}
