#include "../interface/LeptonAnalyzer.h"

double etaForEffectiveArea( const pat::Muon& muon ){
    return muon.eta();
}


double etaForEffectiveArea( const pat::Electron& electron ){
    return electron.superCluster()->eta();
}


double LeptonAnalyzer::getRelIso04(const pat::Muon& mu, const double rho, const EffectiveAreas& effectiveAreas, const bool DeltaBeta) const{
    double puCorr;
    if(!DeltaBeta) puCorr = rho*effectiveAreas.getEffectiveArea( etaForEffectiveArea( mu ) )*( 16. / 9. );
    else           puCorr = 0.5*mu.pfIsolationR04().sumPUPt;

    double absIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - puCorr);
    return absIso/mu.pt();
}


double LeptonAnalyzer::getRelIso03(const pat::Muon& mu, const double rho, const EffectiveAreas& effectiveAreas) const{ //Note: effective area correction is used instead of delta-beta correction
    double puCorr = rho*effectiveAreas.getEffectiveArea( etaForEffectiveArea( mu ) );
    double absIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - puCorr);
    return absIso/mu.pt();
}


double LeptonAnalyzer::getRelIso03(const pat::Electron& ele, const double rho, const EffectiveAreas& effectiveAreas) const{
    double puCorr = rho*effectiveAreas.getEffectiveArea( etaForEffectiveArea( ele ) );
    double absIso = ele.pfIsolationVariables().sumChargedHadronPt + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - puCorr);
    return absIso/ele.pt();
}


double LeptonAnalyzer::getRelIso04(const pat::Electron& ele, const double rho, const EffectiveAreas& effectiveAreas) const{

    //take into account that area is larger in 0.4 cone 
    double puCorr = rho*effectiveAreas.getEffectiveArea( etaForEffectiveArea( ele ) ) * ( 16./9. );
    double absIso = ele.chargedHadronIso() + std::max( 0., ele.neutralHadronIso() + ele.photonIso() - puCorr );
    return absIso/ele.pt();
}
