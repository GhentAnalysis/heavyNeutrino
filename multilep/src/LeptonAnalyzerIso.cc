#include "../interface/LeptonAnalyzer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"


double LeptonAnalyzer::getRelIso04(const pat::Muon& mu) const{ //Note: effective area correction is used instead of delta-beta correction
    double puCorr = 0.5 * mu.pfIsolationR04().sumPUPt;
    double absIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - puCorr);
    return absIso/mu.pt();
}

double LeptonAnalyzer::getRelIso03(const pat::Muon& mu, const double rho) const{ //Note: effective area correction is used instead of delta-beta correction
    double puCorr = rho*muonsEffectiveAreas.getEffectiveArea(mu.eta());
    double absIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - puCorr);
    return absIso/mu.pt();
}

double LeptonAnalyzer::getRelIso03(const pat::Electron& ele, const double rho) const{
    double puCorr = rho*electronsEffectiveAreas.getEffectiveArea(ele.superCluster()->eta());
    double absIso = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - puCorr);
    return absIso/ele.pt();
}


