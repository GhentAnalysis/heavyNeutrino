#include "DataFormats/PatCandidates/interface/Electron.h"

/*
 * Functions for electron identification
 */


namespace electronIdentification {
  /*
   * Manual electron Summer16 cut-based id without isolation
   * Based on https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
   */
  float dEtaInSeed(const pat::Electron* ele){
    if(ele->superCluster().isNonnull() and ele->superCluster()->seed().isNonnull()) return ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta();
    else                                                                            return std::numeric_limits<float>::max();
  }

  bool isLooseCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(not (ele->isEB() or ele->isEE())) return false;

    float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                  >= (ele->isEB() ? 0.11    : 0.0314 ))       return false;
    if(fabs(electronIdentification::dEtaInSeed(ele)) >= (ele->isEB() ? 0.00477 : 0.00868))       return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx())   >= (ele->isEB() ? 0.222   : 0.213  ))       return false;
    if(ele->hadronicOverEm()                         >= (ele->isEB() ? 0.298   : 0.101  ))       return false;
    if(eInvMinusPInv                                 >= (ele->isEB() ? 0.241   : 0.14   ))       return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1)    return false;
    if(not ele->passConversionVeto())                                                            return false;
    return true;
  }

  bool isTightCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(not (ele->isEB() or ele->isEE())) return false;

    float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                  >= (ele->isEB() ? 0.00998 : 0.0292))        return false;
    if(fabs(electronIdentification::dEtaInSeed(ele)) >= (ele->isEB() ? 0.00308 : 0.00605))       return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx())   >= (ele->isEB() ? 0.0816  : 0.0394))        return false;
    if(ele->hadronicOverEm()                         >= (ele->isEB() ? 0.0414  : 0.0641))        return false;
    if(eInvMinusPInv                                 >= (ele->isEB() ? 0.0129  : 0.0129))        return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1)    return false;
    if(not ele->passConversionVeto())                                                            return false;
    return true;
  }

  // Trigger emulation for single electron triggers is available in VID
  // Trigger emulation for double electron triggers with CaloIdL_TrackIdL_IsoVL:
  bool passTriggerEmulationDoubleEG(const pat::Electron* ele){    
    float eInvMinusPInv = fabs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()                >= (ele->isEB() ? 0.011 : 0.030)) return false;
    if(fabs(ele->deltaPhiSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.04  : 0.07))  return false;
    if(fabs(ele->deltaEtaSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.01  : 0.008)) return false;
    if(ele->hadronicOverEm()                       >= (ele->isEB() ? 0.10  : 0.07))  return false;
    if(eInvMinusPInv                               <= -0.05)                        return false;
    if(eInvMinusPInv                               >= (ele->isEB() ? 0.01  : 0.005)) return false;
    return true;
  }

  /*
   * SUSY POG MVA definitions
   */
  float slidingCut(float pt, float low, float high){
    float slope = (high - low)/10.;
    return std::min(low, std::max(high, low + slope*(pt-15)));
  }

  bool passed_MVA_HZZ(const pat::Electron* iE, double mvaValueHZZ){
    if(fabs(iE->eta()) < 0.8)         return mvaValueHZZ > -0.3; 
    else if (fabs(iE->eta()) < 1.479) return mvaValueHZZ > -0.36;
    else                              return mvaValueHZZ > -0.63;
  }

  bool passed_loose_MVA_FR_slidingCut(const pat::Electron* iE, double mvaValue, double mvaValueHZZ){
    if(iE->pt() < 10)                 return passed_MVA_HZZ(iE, mvaValueHZZ);
    if(fabs(iE->eta()) < 0.8)         return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.86, -0.96);
    else if (fabs(iE->eta()) < 1.479) return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.85, -0.96);
    else                              return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.81, -0.95);
  }

  bool passed_medium_MVA_FR_slidingCut(const pat::Electron* iE, double mvaValue){
    if(iE->pt() < 10)                 return false; 
    if(fabs(iE->eta()) < 0.8)         return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.86, -0.86);
    else if (fabs(iE->eta()) < 1.479) return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.85, -0.85);
    else                              return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.81, -0.81);
  }

  bool passed_tight_MVA_FR_slidingCut(const pat::Electron* iE, double mvaValue){
    if(iE->pt() < 10)                 return false; 
    if(fabs(iE->eta()) < 0.8)         return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.77, -0.52);
    else if (fabs(iE->eta()) < 1.479) return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.56, -0.11);
    else                              return mvaValue > electronIdentification::slidingCut(iE->pt(), -0.48, -0.01);
  }

}
