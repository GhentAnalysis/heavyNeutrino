#include "../interface/LeptonAnalyzer.h"
#include "TLorentzVector.h"

/*
 * Overlap of electrons with loose muons
 * Overlap of taus with loose electrons or muons
 */
bool LeptonAnalyzer::eleMuOverlap(const pat::Electron& ele, const bool* loose) const{
    TLorentzVector eleV(ele.px(), ele.py(), ele.pz(), ele.energy());
    for(unsigned m = 0; m < _nMu; ++m){
        if(loose[m]){
            TLorentzVector muV;
            muV.SetPtEtaPhiE(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
            if(eleV.DeltaR(muV) < 0.05) return true;
        }
    }
    return false;
}

bool LeptonAnalyzer::tauLightOverlap(const pat::Tau& tau, const bool* loose) const{
    TLorentzVector tauV(tau.px(), tau.py(), tau.pz(), tau.energy());
    for(unsigned l = 0; l < _nLight; ++l){
        if(loose[l]){
            TLorentzVector lightV;
            lightV.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
            if(tauV.DeltaR(lightV) < 0.4) return true;
        }
    }
    return false;
}

/*
 * Old cut based ID with some variables switched off
 * Recommendation to base yourself on the recommended IDs, or maybe tune one yourself, because what is the point of using official IDs if you remove half of the cuts
 */
float LeptonAnalyzer::dEtaInSeed(const pat::Electron* ele) const{
    if(ele->superCluster().isNonnull() and ele->superCluster()->seed().isNonnull()) return ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta();
    else                                                                            return std::numeric_limits<float>::max();
}

bool LeptonAnalyzer::isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(const pat::Electron* ele) const{
  if(not (ele->isEB() or ele->isEE())) return false;

  float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
  if(ele->full5x5_sigmaIetaIeta()                  >= (ele->isEB() ? 0.11    : 0.0314 ))       return false;
  if(fabs(dEtaInSeed(ele))                         >= (ele->isEB() ? 0.00477 : 0.00868))       return false;
  if(fabs(ele->deltaPhiSuperClusterTrackAtVtx())   >= (ele->isEB() ? 0.222   : 0.213  ))       return false;
  if(ele->hadronicOverEm()                         >= (ele->isEB() ? 0.298   : 0.101  ))       return false;
  if(eInvMinusPInv                                 >= (ele->isEB() ? 0.241   : 0.14   ))       return false;
  //if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1)    return false;
  //if(not ele->passConversionVeto())                                                            return false;
  return true;
}
  

/*
 * Trigger emulation for single electron triggers is available in VID
 * Trigger emulation for double electron triggers with CaloIdL_TrackIdL_IsoVL (is certainly outdated for everything past 2016)
 */
bool LeptonAnalyzer::passTriggerEmulationDoubleEG(const pat::Electron* ele, const bool hOverE) const{
    float eInvMinusPInv =  (1.0 - ele->eSuperClusterOverP())/ele->correctedEcalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                >= (ele->isEB() ? 0.011 : 0.030)) return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.04  : 0.07))  return false;
    if(fabs(ele->deltaEtaSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.01  : 0.008)) return false;
    if(eInvMinusPInv                               <= -0.05)                         return false;
    if(eInvMinusPInv                               >= (ele->isEB() ? 0.01  : 0.005)) return false;
    if(hOverE && (ele->hadronicOverEm()            >= (ele->isEB() ? 0.10  : 0.07))) return false;//Need to be able to check trigEmy without this cut for ewkino
    return true;
}

/*
 * SUSY POG MVA definitions [still here for dependencies in HN and EWK ID's, NEVER use them in new analyses]
 */
float LeptonAnalyzer::slidingCut(float pt, float low, float high) const{
    float slope = (high - low)/10.;
    return std::min(low, std::max(high, low + slope*(pt-15)));
}

bool LeptonAnalyzer::passingElectronMvaHZZ(const pat::Electron* ele, double mvaValueHZZ) const{
    if(fabs(ele->eta()) < 0.8)         return mvaValueHZZ > -0.3; 
    else if (fabs(ele->eta()) < 1.479) return mvaValueHZZ > -0.36;
    else                               return mvaValueHZZ > -0.63;
}

bool LeptonAnalyzer::passingElectronMvaLooseSusy(const pat::Electron* ele, double mvaValue, double mvaValueHZZ) const{
    if(ele->pt() < 10)                 return passingElectronMvaHZZ(ele, mvaValueHZZ);
    if(fabs(ele->eta()) < 0.8)         return mvaValue > slidingCut(ele->pt(), -0.86, -0.96);
    else if (fabs(ele->eta()) < 1.479) return mvaValue > slidingCut(ele->pt(), -0.85, -0.96);
    else                               return mvaValue > slidingCut(ele->pt(), -0.81, -0.95);
}

bool LeptonAnalyzer::passingElectronMvaTightSusy(const pat::Electron* ele, double mvaValue) const{
    if(ele->pt() < 10)                 return false; 
    if(fabs(ele->eta()) < 0.8)         return mvaValue > slidingCut(ele->pt(),  0.77,  0.52);
    else if (fabs(ele->eta()) < 1.479) return mvaValue > slidingCut(ele->pt(),  0.56,  0.11);
    else                               return mvaValue > slidingCut(ele->pt(),  0.48, -0.01);
}

double LeptonAnalyzer::leptonMvaVal(const pat::Muon& muon, LeptonMvaHelper* mvaHelper){
    return mvaHelper->leptonMvaMuon(_lPt[_nL],
            _lEta[_nL],
            _selectedTrackMult[_nL],
            _miniIsoCharged[_nL],
            _miniIso[_nL] - _miniIsoCharged[_nL],
            _ptRel[_nL],
            _ptRatio[_nL],
            _closestJetDeepCsv[_nL],
            _closestJetDeepFlavor[_nL],
            _3dIPSig[_nL],
            _dxy[_nL],
            _dz[_nL],
            _relIso[_nL],
            muon.segmentCompatibility()
            );
}

double LeptonAnalyzer::leptonMvaVal(const pat::Electron& electron, LeptonMvaHelper* mvaHelper){
    return mvaHelper->leptonMvaElectron(_lPt[_nL],
            _lEta[_nL],
            _selectedTrackMult[_nL],
            _miniIsoCharged[_nL],
            _miniIso[_nL] - _miniIsoCharged[_nL],
            _ptRel[_nL],
            _ptRatio[_nL],
            _closestJetDeepCsv[_nL],
            _closestJetDeepFlavor[_nL],
            _3dIPSig[_nL],
            _dxy[_nL],
            _dz[_nL],
            _relIso_ttH[_nL],
            _lElectronMvaSummer16GP[_nL],
            _lElectronMvaFall17v1NoIso[_nL],
            _lElectronMvaFall17NoIso[_nL]
            );
}
