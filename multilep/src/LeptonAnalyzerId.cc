#include "../interface/LeptonAnalyzer.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"

/*
 * Manual electron Summer16 cut-based id without isolation
 * Based on https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedLeptonAnalyzerRun2#Offline_selection_criteria
 */
float LeptonAnalyzer::dEtaInSeed(const pat::Electron* ele){
    if(ele->superCluster().isNonnull() and ele->superCluster()->seed().isNonnull()) return ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta();
    else                                                                            return std::numeric_limits<float>::max();
}

bool LeptonAnalyzer::isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(const pat::Electron* ele){
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
  


bool LeptonAnalyzer::isLooseCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(!(ele->isEB() or ele->isEE()))                                                            return false;
    float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                  >= (ele->isEB() ? 0.11    : 0.0314 ))       return false;
    if(fabs(dEtaInSeed(ele))                         >= (ele->isEB() ? 0.00477 : 0.00868))       return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx())   >= (ele->isEB() ? 0.222   : 0.213  ))       return false;
    if(ele->hadronicOverEm()                         >= (ele->isEB() ? 0.298   : 0.101  ))       return false;
    if(eInvMinusPInv                                 >= (ele->isEB() ? 0.241   : 0.14   ))       return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1)    return false;
    if(!ele->passConversionVeto())                                                               return false;
    return true;
}

bool LeptonAnalyzer::isMediumCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(!(ele->isEB() or ele->isEE()))                                                            return false;
    float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                  >= (ele->isEB() ? 0.00998    : 0.0298 ))    return false;
    if(fabs(dEtaInSeed(ele))                         >= (ele->isEB() ? 0.00311    : 0.00609))    return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx())   >= (ele->isEB() ? 0.103      : 0.045  ))    return false;
    if(ele->hadronicOverEm()                         >= (ele->isEB() ? 0.253      : 0.0878  ))   return false;
    if(eInvMinusPInv                                 >= (ele->isEB() ? 0.134      : 0.13   ))    return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1)    return false;
    if(!ele->passConversionVeto())                                                               return false;
    return true;
}

bool LeptonAnalyzer::isTightCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(!(ele->isEB() or ele->isEE())) return false;
    float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                  >= (ele->isEB() ? 0.00998 : 0.0292))        return false;
    if(fabs(dEtaInSeed(ele))                         >= (ele->isEB() ? 0.00308 : 0.00605))       return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx())   >= (ele->isEB() ? 0.0816  : 0.0394))        return false;
    if(ele->hadronicOverEm()                         >= (ele->isEB() ? 0.0414  : 0.0641))        return false;
    if(eInvMinusPInv                                 >= (ele->isEB() ? 0.0129  : 0.0129))        return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1)    return false;
    if(!ele->passConversionVeto())                                                               return false;
    return true;
}

// Trigger emulation for single electron triggers is available in VID
// Trigger emulation for double electron triggers with CaloIdL_TrackIdL_IsoVL:
bool LeptonAnalyzer::passTriggerEmulationDoubleEG(const pat::Electron* ele, const bool hOverE){    
    float eInvMinusPInv =  (1.0 - ele->eSuperClusterOverP())/ele->correctedEcalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                >= (ele->isEB() ? 0.011 : 0.030)) return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.04  : 0.07))  return false;
    if(fabs(ele->deltaEtaSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.01  : 0.008)) return false;
    if(eInvMinusPInv                               <= -0.05)                         return false;
    if(eInvMinusPInv                               >= (ele->isEB() ? 0.01  : 0.005)) return false;
    if(hOverE && (ele->hadronicOverEm()            >= (ele->isEB() ? 0.10  : 0.07))) return false;//Need to be able to check trigEmy without this cut for ewkino
    return true;
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
            _3dIPSig[_nL],
            _dxy[_nL],
            _dz[_nL],
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
            _3dIPSig[_nL],
            _dxy[_nL],
            _dz[_nL],
            _lElectronMva[_nL],
            _lElectronMvaHZZ[_nL]
            );
}
<<<<<<< HEAD
=======

bool LeptonAnalyzer::isEwkLoose(const pat::Muon& lep){
    if(fabs(_dxy[_nL]) >= 0.05 || fabs(_dz[_nL]) >= 0.1 || _3dIPSig[_nL] >= 8) return false;
    if(_miniIso[_nL] >= 0.4) return false;
    if(_lPt[_nL] <= 5 || fabs(_lEta[_nL]) >= 2.4) return false;
    return lep.isLooseMuon();
}

bool LeptonAnalyzer::isEwkLoose(const pat::Electron& lep){
    if(fabs(_dxy[_nL]) >= 0.05 || fabs(_dz[_nL]) >= 0.1 || _3dIPSig[_nL] >= 8) return false;
    if(_miniIso[_nL] >= 0.4) return false;
    if(_lPt[_nL] <= 7 || fabs(_lEta[_nL]) >= 2.5) return false;
    if(lep.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) return false;
    if(eleMuOverlap(lep, _lEwkLoose)) return false;
    return passingElectronMvaLooseSusy(&lep, _lElectronMva[_nL], _lElectronMvaHZZ[_nL]);
}

bool LeptonAnalyzer::isEwkLoose(const pat::Tau& tau){
    if( _lPt[_nL] <= 20 || fabs(_lEta[_nL]) >= 2.3) return false;
    if(!_lPOGVeto[_nL]) return false;
    if(!_tauEleVeto[_nL]) return false;
    return tauLightOverlap(tau, _lEwkLoose);
}

bool LeptonAnalyzer::isEwkFO(const pat::Muon& lep){
    if(!_lEwkLoose[_nL]) return false;
    if(_lPt[_nL] <= 10) return false;
    if(!lep.isMediumMuon()) return false;
    return _leptonMvaSUSY[_nL] > -0.2 || (_ptRatio[_nL] > 0.3 && _closestJetCsvV2[_nL] < 0.3);
}

bool LeptonAnalyzer::isEwkFO(const pat::Electron& lep){
    if(!_lEwkLoose[_nL]) return false;
    if(_lPt[_nL] <= 10) return false;
    if(!passTriggerEmulationDoubleEG(&lep, false)) return false;
    if(lep.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) !=0) return false;
    double ptCone = _lPt[_nL];
    if(_leptonMvaSUSY[_nL] <= 0.5){
        ptCone *= 0.85/_ptRatio[_nL];
    }
    if(ptCone >= 30 && lep.hadronicOverEm() >= (lep.isEB() ? 0.10  : 0.07) ) return false;
    return _leptonMvaSUSY[_nL] > 0.5 || (passElectronMvaEwkFO(&lep, _lElectronMva[_nL]) && _ptRatio[_nL] > 0.3 && _closestJetCsvV2[_nL] < 0.3);
}

bool LeptonAnalyzer::isEwkFO(const pat::Tau& tau){
    return _lEwkLoose[_nL];
}

bool LeptonAnalyzer::isEwkTight(const pat::Muon& lep){
    if(!_lEwkFO[_nL]) return false;
    return _leptonMvaSUSY[_nL] > -0.2;
}

bool LeptonAnalyzer::isEwkTight(const pat::Electron& lep){
    if(!_lEwkFO[_nL]) return false;
    if(!passTriggerEmulationDoubleEG(&lep)) return false;
    if(!lep.passConversionVeto()) return false;
    return _leptonMvaSUSY[_nL] > 0.5;
}

bool LeptonAnalyzer::isEwkTight(const pat::Tau& tau){
    return _lEwkFO[_nL] && _lPOGTight[_nL];
}

>>>>>>> master
