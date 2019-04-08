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

/*
 * Own HeavyNeutrino FO tune [tuned on a very old electronMva, do NOT use them for new analyses]
 */
bool LeptonAnalyzer::passingElectronMvaHeavyNeutrinoFO(const pat::Electron* ele, double mvaValue) const{
    if(ele->pt() < 10)                 return false; 
    if(fabs(ele->eta()) < 0.8)         return mvaValue > -0.02;
    else                               return mvaValue > -0.52;
}

/*
 * Ewkino FO tune [tuned on a very old electronMva, do NOT use them for new analyses]
 */
bool LeptonAnalyzer::passElectronMvaEwkFO(const pat::Electron* ele, double mvaValue) const{
    if(ele->pt() < 10)               return false;
    if(fabs(ele->eta()) < 1.479)     return mvaValue > 0.0;
    else                             return mvaValue > 0.3;
}


double LeptonAnalyzer::leptonMvaVal(const pat::Muon& muon, LeptonMvaHelper* mvaHelper){
    return mvaHelper->leptonMvaMuon(_lPt[_nL],
            _lEta[_nL],
            _selectedTrackMult[_nL],
            _miniIsoCharged[_nL],
            _miniIso[_nL] - _miniIsoCharged[_nL],
            _ptRel[_nL],
            _ptRatio[_nL],
            _closestJetCsvV2[_nL],
            _closestJetDeepCsv_b[_nL] + _closestJetDeepCsv_bb[_nL],
            _3dIPSig[_nL],
            _dxy[_nL],
            _dz[_nL],
            _relIso[_nL],
            _relIso0p4[_nL],
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
            _closestJetCsvV2[_nL],
            _closestJetDeepCsv_b[_nL] + _closestJetDeepCsv_bb[_nL],
            _3dIPSig[_nL],
            _dxy[_nL],
            _dz[_nL],
            _relIso[_nL],
            _relIso0p4[_nL],
            _lElectronMvaSummer16GP[_nL],
            _lElectronMvaSummer16HZZ[_nL],
            _lElectronMvaFall17v1NoIso[_nL]
            );
}

bool LeptonAnalyzer::isEwkLoose(const pat::Muon& lep) const{
    if(fabs(_dxy[_nL]) >= 0.05 || fabs(_dz[_nL]) >= 0.1 || _3dIPSig[_nL] >= 8)      return false;
    if(_miniIso[_nL] >= 0.4)                                                        return false;
    if(_lPt[_nL] <= 5 || fabs(_lEta[_nL]) >= 2.4)                                   return false;
    return lep.isLooseMuon();
}

bool LeptonAnalyzer::isEwkLoose(const pat::Electron& lep) const{
    if(fabs(_dxy[_nL]) >= 0.05 || fabs(_dz[_nL]) >= 0.1 || _3dIPSig[_nL] >= 8)                      return false;
    if(_miniIso[_nL] >= 0.4)                                                                        return false;
    if(_lPt[_nL] <= 7 || fabs(_lEta[_nL]) >= 2.5)                                                   return false;
    if(lep.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1)     return false;
    if(eleMuOverlap(lep, _lEwkLoose))                                                               return false;
    return passingElectronMvaLooseSusy(&lep, _lElectronMvaSummer16GP[_nL], _lElectronMvaSummer16HZZ[_nL]); // TODO: consider replacing this by something better
}

bool LeptonAnalyzer::isEwkLoose(const pat::Tau& tau) const{
    if( _lPt[_nL] <= 20 || fabs(_lEta[_nL]) >= 2.3)     return false;
    if(!_lPOGVeto[_nL])                                 return false;
    if(!_tauEleVetoLoose[_nL])                          return false;
    return tauLightOverlap(tau, _lEwkLoose);
}

bool LeptonAnalyzer::isEwkFO(const pat::Muon& lep) const{
    if(!_lEwkLoose[_nL])        return false;
    if(_lPt[_nL] <= 10)         return false;
    if(!lep.isMediumMuon())     return false;
    return _leptonMvaSUSY[_nL] > -0.2 || (_ptRatio[_nL] > 0.3 && _closestJetCsvV2[_nL] < 0.3);
}

bool LeptonAnalyzer::isEwkFO(const pat::Electron& lep) const{
    if(!_lEwkLoose[_nL])                                                                            return false;
    if(_lPt[_nL] <= 10)                                                                             return false;
    if(!passTriggerEmulationDoubleEG(&lep, false))                                                  return false;
    if(lep.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) !=0)     return false;
    double ptCone = _lPt[_nL];
    if(_leptonMvaSUSY[_nL] <= 0.5){
        ptCone *= 0.85/_ptRatio[_nL];
    }
    if(ptCone >= 30 && lep.hadronicOverEm() >= (lep.isEB() ? 0.10  : 0.07) )                        return false;
    return _leptonMvaSUSY[_nL] > 0.5 || (passElectronMvaEwkFO(&lep, _lElectronMvaSummer16GP[_nL]) && _ptRatio[_nL] > 0.3 && _closestJetCsvV2[_nL] < 0.3);
}

bool LeptonAnalyzer::isEwkFO(const pat::Tau& tau) const{
    return _lEwkLoose[_nL];
}

bool LeptonAnalyzer::isEwkTight(const pat::Muon& lep) const{
    if(!_lEwkFO[_nL]) return false;
    return _leptonMvaSUSY[_nL] > -0.2;
}

bool LeptonAnalyzer::isEwkTight(const pat::Electron& lep) const{
    if(!_lEwkFO[_nL])                       return false;
    if(!passTriggerEmulationDoubleEG(&lep)) return false;
    if(!lep.passConversionVeto())           return false;
    return _leptonMvaSUSY[_nL] > 0.5;
}

bool LeptonAnalyzer::isEwkTight(const pat::Tau& tau) const{
    return _lEwkFO[_nL] && _lPOGTight[_nL];
}

