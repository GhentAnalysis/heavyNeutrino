#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "heavyNeutrino/multilep/interface/TauTools.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>


// TODO: we should maybe stop indentifying effective areas by year, as they are typically more connected to a specific ID than to a specific year
LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer),
    electronsEffectiveAreas( ( multilepAnalyzer->is2017() || multilepAnalyzer->is2018() ) ? iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreasFall17").fullPath() : iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreasRelIso2016").fullPath() ),
    electronsEffectiveAreasMiniIso( ( multilepAnalyzer->is2017() || multilepAnalyzer->is2018() ) ? iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreasFall17").fullPath() : iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreasMiniIso2016").fullPath() ),
    muonsEffectiveAreas    ((multilepAnalyzer->is2017() || multilepAnalyzer->is2018() )? (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath() )
{
    leptonMvaComputerTTH = new LeptonMvaHelper(iConfig, true, !multilepAnalyzer->is2016() );
    leptonMvaComputertZq = new LeptonMvaHelper(iConfig, false, !multilepAnalyzer->is2016() );
};

LeptonAnalyzer::~LeptonAnalyzer(){
    delete leptonMvaComputerTTH;
    delete leptonMvaComputertZq;
}

void LeptonAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nL",                           &_nL,                           "_nL/i");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/i");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/i");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/i");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/i");
    outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
    outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
    outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
    outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
    outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
    outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/i");
    outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
    outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
    outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
    outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
    outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
    outputTree->Branch("_lElectronSummer16MvaGP",       &_lElectronMvaSummer16GP,       "_lElectronMvaSummer16GP[_nLight]/F");
    outputTree->Branch("_lElectronSummer16MvaHZZ",      &_lElectronMvaSummer16HZZ,      "_lElectronMvaSummer16HZZ[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17v1NoIso",    &_lElectronMvaFall17v1NoIso,    "_lElectronMvaFall17v1NoIso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17Iso",        &_lElectronMvaFall17Iso,        "_lElectronMvaFall17Iso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17NoIso",      &_lElectronMvaFall17NoIso,      "_lElectronMvaFall17NoIso[_nLight]/F");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_lEleIsEB",                     &_lEleIsEB ,                    "_lEleIsEB[_nLight]/O");
    outputTree->Branch("_lEleIsEE",                     &_lEleIsEE ,                    "_lEleIsEE[_nLight]/O");
    outputTree->Branch("_lEleSuperClusterOverP",        &_lEleSuperClusterOverP ,       "_lEleSuperClusterOverP[_nLight]/D");
    outputTree->Branch("_lEleEcalEnergy",               &_lEleEcalEnergy ,              "_lEleEcalEnergy[_nLight]/D");
    outputTree->Branch("_lElefull5x5SigmaIetaIeta",     &_lElefull5x5SigmaIetaIeta ,    "_lElefull5x5SigmaIetaIeta[_nLight]/D");
    outputTree->Branch("_lEleDEtaInSeed",               &_lEleDEtaInSeed ,              "_lEleDEtaInSeed[_nLight]/D");
    outputTree->Branch("_lEleDeltaPhiSuperClusterTrackAtVtx", &_lEleDeltaPhiSuperClusterTrackAtVtx , "_lEleDeltaPhiSuperClusterTrackAtVtx[_nLight]/D");
    outputTree->Branch("_lElehadronicOverEm",           &_lElehadronicOverEm ,          "_lElehadronicOverEm[_nLight]/D");
    outputTree->Branch("_lEleInvMinusPInv",             &_lEleInvMinusPInv ,            "_lEleInvMinusPInv[_nLight]/D");
    outputTree->Branch("_eleNumberInnerHitsMissing",    &_eleNumberInnerHitsMissing,    "_eleNumberInnerHitsMissing[_nLight]/D");
    outputTree->Branch("_leptonMvaSUSY",                &_leptonMvaSUSY,                "_leptonMvaSUSY[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH",                 &_leptonMvaTTH,                 "_leptonMvaTTH[_nLight]/D");
    outputTree->Branch("_leptonMvatZq",                 &_leptonMvatZq,                 "_leptonMvatZq[_nLight]/D");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
    outputTree->Branch("_tauMuonVetoLoose",             &_tauMuonVetoLoose,             "_tauMuonVetoLoose[_nL]/O");
    outputTree->Branch("_tauEleVetoLoose",              &_tauEleVetoLoose,              "_tauEleVetoLoose[_nL]/O");
    outputTree->Branch("_tauDecayMode",                 &_tauDecayMode,                 "_tauDecayMode[_nL]/i");
    outputTree->Branch("_decayModeFinding",             &_decayModeFinding,             "_decayModeFinding[_nL]/O");
    outputTree->Branch("_tauPOGVTight2017v2",           &_tauPOGVTight2017v2,           "_tauPOGVTight2017v2[_nL]/O");
    outputTree->Branch("_tauAgainstElectronMVA6Raw",    &_tauAgainstElectronMVA6Raw,    "_tauAgainstElectronMVA6Raw[_nL]/D");
    outputTree->Branch("_tauCombinedIsoDBRaw3Hits",     &_tauCombinedIsoDBRaw3Hits,     "_tauCombinedIsoDBRaw3Hits[_nL]/D");
    outputTree->Branch("_tauIsoMVAPWdR03oldDMwLT",      &_tauIsoMVAPWdR03oldDMwLT,      "_tauIsoMVAPWdR03oldDMwLT[_nL]/D");
    outputTree->Branch("_tauIsoMVADBdR03oldDMwLT",      &_tauIsoMVADBdR03oldDMwLT,      "_tauIsoMVADBdR03oldDMwLT[_nL]/D");
    outputTree->Branch("_tauIsoMVADBdR03newDMwLT",      &_tauIsoMVADBdR03newDMwLT,      "_tauIsoMVADBdR03newDMwLT[_nL]/D");
    outputTree->Branch("_tauIsoMVAPWnewDMwLT",          &_tauIsoMVAPWnewDMwLT,          "_tauIsoMVAPWnewDMwLT[_nL]/D");
    outputTree->Branch("_tauIsoMVAPWoldDMwLT",          &_tauIsoMVAPWoldDMwLT,          "_tauIsoMVAPWoldDMwLT[_nL]/D");
    outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
    outputTree->Branch("_relIso0p4",                    &_relIso0p4,                    "_relIso0p4[_nLight]/D");
    outputTree->Branch("_relIso0p4MuDeltaBeta",         &_relIso0p4MuDeltaBeta,         "_relIso0p4MuDeltaBeta[_nMu]/D");
    outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
    outputTree->Branch("_miniIsoCharged",               &_miniIsoCharged,               "_miniIsoCharged[_nLight]/D");
    outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
    outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
    outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv",            &_closestJetDeepCsv,            "_closestJetDeepCsv[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_b",       &_closestJetDeepFlavor_b,       "_closestJetDeepFlavor_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_bb",      &_closestJetDeepFlavor_bb,      "_closestJetDeepFlavor_bb[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_lepb",    &_closestJetDeepFlavor_lepb,    "_closestJetDeepFlavor_lepb[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor",         &_closestJetDeepFlavor,         "_closestJetDeepFlavor[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    //outputTree->Branch("_lKVF_valid",			        &_lKVF_valid,			        "_lKVF_valid[_nLight]/O");
    //outputTree->Branch("_lKVF_x",			            &_lKVF_x,			            "_lKVF_x[_nLight]/D");
    //outputTree->Branch("_lKVF_y",			            &_lKVF_y,			            "_lKVF_y[_nLight]/D");
    //outputTree->Branch("_lKVF_z",			            &_lKVF_z,			            "_lKVF_z[_nLight]/D");
    //outputTree->Branch("_lKVF_cxx",			            &_lKVF_cxx,			            "_lKVF_cxx[_nLight]/D");
    //outputTree->Branch("_lKVF_cyy",			            &_lKVF_cyy,			            "_lKVF_cyy[_nLight]/D");
    //outputTree->Branch("_lKVF_czz",			            &_lKVF_czz,			            "_lKVF_czz[_nLight]/D");
    //outputTree->Branch("_lKVF_cyx",			            &_lKVF_cyx,			            "_lKVF_cyx[_nLight]/D");
    //outputTree->Branch("_lKVF_czy",			            &_lKVF_czy,			            "_lKVF_czy[_nLight]/D");
    //outputTree->Branch("_lKVF_czx",			            &_lKVF_czx,			            "_lKVF_czx[_nLight]/D");
    //outputTree->Branch("_lKVF_df",			            &_lKVF_df,			            "_lKVF_df[_nLight]/D");
    //outputTree->Branch("_lKVF_chi2",			        &_lKVF_chi2,			        "_lKVF_chi2[_nLight]/D");
    //outputTree->Branch("_lKVF_ntracks",		            &_lKVF_ntracks,		            "_lKVF_ntracks[_nLight]/i");
    //outputTree->Branch("_lKVF_dRcut",                   &_lKVF_dRcut,                   "_lKVF_dRcut[_nLight]/D");
    //outputTree->Branch("_lKVF_trackPt",                 &_lKVF_trackPt,                 "_lKVF_trackPt[_nLight][15]/D");
    //outputTree->Branch("_lKVF_trackEta",                &_lKVF_trackEta,                "_lKVF_trackEta[_nLight][15]/D");
    //outputTree->Branch("_lKVF_trackPhi",                &_lKVF_trackPhi,                "_lKVF_trackPhi[_nLight][15]/D");
    //outputTree->Branch("_lKVF_trackE",                  &_lKVF_trackE,                  "_lKVF_trackE[_nLight][15]/D");
    //outputTree->Branch("_lKVF_trackdR",                 &_lKVF_trackdR,                 "_lKVF_trackdR[_nLight][15]/D");
    //outputTree->Branch("_lKVF_trackdxy",                &_lKVF_trackdxy,                "_lKVF_trackdxy[_nLight][15]/D");
    //outputTree->Branch("_lKVF_trackdz",                 &_lKVF_trackdz,                 "_lKVF_trackdz[_nLight][15]/D");
    //outputTree->Branch("_IVF_nvertex",                  &_IVF_nvertex,                  "_IVF_nvertex/i");
    outputTree->Branch("_IVF_x",                        &_IVF_x,                        "_IVF_x[_nLight]/D");
    outputTree->Branch("_IVF_y",                        &_IVF_y,                        "_IVF_y[_nLight]/D");
    outputTree->Branch("_IVF_z",                        &_IVF_z,                        "_IVF_z[_nLight]/D");
    outputTree->Branch("_IVF_cx",                       &_IVF_cx,                       "_IVF_cx[_nLight]/D");
    outputTree->Branch("_IVF_cy",                       &_IVF_cy,                       "_IVF_cy[_nLight]/D");
    outputTree->Branch("_IVF_cz",                       &_IVF_cz,                       "_IVF_cz[_nLight]/D");
    outputTree->Branch("_IVF_df",                       &_IVF_df,                       "_IVF_df[_nLight]/D");
    outputTree->Branch("_IVF_chi2",                     &_IVF_chi2,                     "_IVF_chi2[_nLight]/D");
    outputTree->Branch("_IVF_pt",                       &_IVF_pt,                       "_IVF_pt[_nLight]/D");
    outputTree->Branch("_IVF_eta",                      &_IVF_eta,                      "_IVF_eta[_nLight]/D");
    outputTree->Branch("_IVF_phi",                      &_IVF_phi,                      "_IVF_phi[_nLight]/D");
    outputTree->Branch("_IVF_E",                        &_IVF_E,                        "_IVF_E[_nLight]/D");
    outputTree->Branch("_IVF_mass",                     &_IVF_mass,                     "_IVF_mass[_nLight]/D");
    outputTree->Branch("_IVF_ntracks",                  &_IVF_ntracks,                  "_IVF_ntracks[_nLight]/i");
    outputTree->Branch("_IVF_trackpt",                  &_IVF_trackpt,                  "_IVF_trackpt[_nLight][15]/D");
    outputTree->Branch("_IVF_tracketa",                 &_IVF_tracketa,                 "_IVF_tracketa[_nLight][15]/D");
    outputTree->Branch("_IVF_trackphi",                 &_IVF_trackphi,                 "_IVF_trackphi[_nLight][15]/D");
    outputTree->Branch("_IVF_trackE",                   &_IVF_trackE,                   "_IVF_trackE[_nLight][15]/D");
    outputTree->Branch("_IVF_trackdxy",                 &_IVF_trackdxy,                 "_IVF_trackdxy[_nLight][15]/D");
    outputTree->Branch("_IVF_trackdz",                  &_IVF_trackdz,                  "_IVF_trackdz[_nLight][15]/D");
    outputTree->Branch("_IVF_trackcharge",              &_IVF_trackcharge,              "_IVF_trackcharge[_nLight][15]/D");
    outputTree->Branch("_lIVF_match",                   &_lIVF_match,                   "_lIVF_match[_nLight]/O");
    outputTree->Branch("_lGlobalMuon",                  &_lGlobalMuon,                  "_lGlobalMuon[_nMu]/O");
    outputTree->Branch("_lTrackerMuon",                 &_lTrackerMuon,                 "_lTrackerMuon[_nMu]/O");
    outputTree->Branch("_lInnerTrackValidFraction",     &_lInnerTrackValidFraction,     "_lInnerTrackValidFraction[_nMu]/D");
    outputTree->Branch("_lGlobalTrackNormalizedChi2",   &_lGlobalTrackNormalizedChi2,   "_lGlobalTrackNormalizedChi2[_nMu]/D");
    outputTree->Branch("_lCQChi2Position",              &_lCQChi2Position,              "_lCQChi2Position[_nMu]/D");
    outputTree->Branch("_lCQTrackKink",                 &_lCQTrackKink,                 "_lCQTrackKink[_nMu]/D");
    outputTree->Branch("_lNumberOfMatchedStation",      &_lNumberOfMatchedStation,      "_lNumberOfMatchedStation[_nMu]/i");
    outputTree->Branch("_lNumberOfValidPixelHits",      &_lNumberOfValidPixelHits,      "_lNumberOfValidPixelHits[_nMu]/i");
    outputTree->Branch("_muNumberInnerHits",            &_muNumberInnerHits,            "_muNumberInnerHits[_nMu]/i");
    outputTree->Branch("_lTrackerLayersWithMeasurement",&_lTrackerLayersWithMeasurement,"_lTrackerLayersWithMeasurement[_nMu]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");
    if( multilepAnalyzer->isMC() ){
        outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
        outputTree->Branch("_lMatchCharge",               &_lMatchCharge,                 "_lMatchCharge[_nL]/I");
        outputTree->Branch("_tauGenStatus",               &_tauGenStatus,                 "_tauGenStatus[_nL]/i");
        outputTree->Branch("_lMomPdgId",                  &_lMomPdgId,                    "_lMomPdgId[_nL]/I");
        outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
        outputTree->Branch("_lProvenanceCompressed",      &_lProvenanceCompressed,        "_lProvenanceCompressed[_nL]/i");
        outputTree->Branch("_lProvenanceConversion",      &_lProvenanceConversion,        "_lProvenanceConversion[_nL]/i");
    }
    outputTree->Branch("_lPtCorr",                    &_lPtCorr,                      "_lPtCorr[_nLight]/D");
    outputTree->Branch("_lPtScaleUp",                 &_lPtScaleUp,                   "_lPtScaleUp[_nLight]/D");
    outputTree->Branch("_lPtScaleDown",               &_lPtScaleDown,                 "_lPtScaleDown[_nLight]/D");
    outputTree->Branch("_lPtResUp",                   &_lPtResUp,                     "_lPtResUp[_nLight]/D");
    outputTree->Branch("_lPtResDown",                 &_lPtResDown,                   "_lPtResDown[_nLight]/D");
    outputTree->Branch("_lECorr",                     &_lECorr,                       "_lECorr[_nLight]/D");
    outputTree->Branch("_lEScaleUp",                  &_lEScaleUp,                    "_lEScaleUp[_nLight]/D");
    outputTree->Branch("_lEScaleDown",                &_lEScaleDown,                  "_lEScaleDown[_nLight]/D");
    outputTree->Branch("_lEResUp",                    &_lEResUp,                      "_lEResUp[_nLight]/D");
    outputTree->Branch("_lEResDown",                  &_lEResDown,                    "_lEResDown[_nLight]/D");
    if(multilepAnalyzer->storeAllTauID){
      outputTree->Branch("_decayModeFindingNew",          &_decayModeFindingNew,          "_decayModeFindingNew[_nL]/O");
      outputTree->Branch("_tauPOGVLoose2015",             &_tauPOGVLoose2015,             "_tauPOGVLoose2015[_nL]/O");
      outputTree->Branch("_tauPOGLoose2015",              &_tauPOGLoose2015,              "_tauPOGLoose2015[_nL]/O");
      outputTree->Branch("_tauPOGMedium2015",             &_tauPOGMedium2015,             "_tauPOGMedium2015[_nL]/O");
      outputTree->Branch("_tauPOGTight2015",              &_tauPOGTight2015,              "_tauPOGTight2015[_nL]/O");
      outputTree->Branch("_tauPOGVTight2015",             &_tauPOGVTight2015,             "_tauPOGVTight2015[_nL]/O");
      outputTree->Branch("_tauVLooseMvaNew2015",          &_tauVLooseMvaNew2015,          "_tauVLooseMvaNew2015[_nL]/O");
      outputTree->Branch("_tauLooseMvaNew2015",           &_tauLooseMvaNew2015,           "_tauLooseMvaNew2015[_nL]/O");
      outputTree->Branch("_tauMediumMvaNew2015",          &_tauMediumMvaNew2015,          "_tauMediumMvaNew2015[_nL]/O");
      outputTree->Branch("_tauTightMvaNew2015",           &_tauTightMvaNew2015,           "_tauTightMvaNew2015[_nL]/O");
      outputTree->Branch("_tauVTightMvaNew2015",          &_tauVTightMvaNew2015,          "_tauVTightMvaNew2015[_nL]/O");
      outputTree->Branch("_tauPOGVVLoose2017v2",          &_tauPOGVVLoose2017v2,          "_tauPOGVVLoose2017v2[_nL]/O");
      outputTree->Branch("_tauPOGVVTight2017v2",          &_tauPOGVVTight2017v2,          "_tauPOGVVTight2017v2[_nL]/O");
      outputTree->Branch("_tauVLooseMvaNew2017v2",        &_tauVLooseMvaNew2017v2,        "_tauVLooseMvaNew2017v2[_nL]/O");
      outputTree->Branch("_tauLooseMvaNew2017v2",         &_tauLooseMvaNew2017v2,         "_tauLooseMvaNew2017v2[_nL]/O");
      outputTree->Branch("_tauMediumMvaNew2017v2",        &_tauMediumMvaNew2017v2,        "_tauMediumMvaNew2017v2[_nL]/O");
      outputTree->Branch("_tauTightMvaNew2017v2",         &_tauTightMvaNew2017v2,         "_tauTightMvaNew2017v2[_nL]/O");
      outputTree->Branch("_tauVTightMvaNew2017v2",        &_tauVTightMvaNew2017v2,        "_tauVTightMvaNew2017v2[_nL]/O");
      outputTree->Branch("_tauMuonVetoTight",             &_tauMuonVetoTight,             "_tauMuonVetoTight[_nL]/O");
      outputTree->Branch("_tauEleVetoVLoose",             &_tauEleVetoVLoose,             "_tauEleVetoVLoose[_nL]/O");
      outputTree->Branch("_tauEleVetoMedium",             &_tauEleVetoMedium,             "_tauEleVetoMedium[_nL]/O");
      outputTree->Branch("_tauEleVetoTight",              &_tauEleVetoTight,              "_tauEleVetoTight[_nL]/O");
      outputTree->Branch("_tauEleVetoVTight",             &_tauEleVetoVTight,             "_tauEleVetoVTight[_nL]/O");
     }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& primaryVertex){
    edm::Handle<std::vector<pat::Electron>> electrons          = getHandle(iEvent, multilepAnalyzer->eleToken);
    edm::Handle<std::vector<pat::Muon>> muons                  = getHandle(iEvent, multilepAnalyzer->muonToken);
    edm::Handle<std::vector<pat::Tau>> taus                    = getHandle(iEvent, multilepAnalyzer->tauToken);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands = getHandle(iEvent, multilepAnalyzer->packedCandidatesToken);
    edm::Handle<double> rho                                    = getHandle(iEvent, multilepAnalyzer->rhoToken);
    edm::Handle<std::vector<pat::Jet>> jets                    = getHandle(iEvent, multilepAnalyzer->jetToken);
  //edm::Handle<std::vector<pat::Jet>> jets                    = getHandle(iEvent, multilepAnalyzer->jetSmearedToken);  // Are we sure we do not want the smeared jets here???
    edm::Handle<std::vector<reco::GenParticle>> genParticles   = getHandle(iEvent, multilepAnalyzer->genParticleToken);
    edm::Handle<std::vector<reco::Vertex>> secVertices         = getHandle(iEvent, multilepAnalyzer->secondaryVerticesToken);

    iSetup.get<IdealMagneticFieldRecord>().get(_bField);
    iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", _shProp);

    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;

    //std::cout << std::endl << std::endl << std::endl << "=== EVENT ===" << std::endl;
    //std::cout << "sec vertices size: " << (*secVertices).size() << std::endl;
    //fillAllIVFVariables(*secVertices, primaryVertex);

    //loop over muons
    // muons need to be run first, because some ID's need to calculate a muon veto for electrons
    for(const pat::Muon& mu : *muons){
        if(_nL == nL_max)                              break;
        if(mu.innerTrack().isNull())                   continue;
        if(mu.pt() < 5)                                continue;
        if(fabs(mu.eta()) > 2.4)                       continue;
        if(!mu.isPFMuon())                             continue;
        if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue;
        fillLeptonImpactParameters(mu);
        //if(fabs(_dxy[_nL]) > 0.05)                     continue; this can be applied at root level, but not here to preserve displaced signal
        //if(fabs(_dz[_nL]) > 0.1)                       continue;
        fillLeptonKinVars(mu);
        if( multilepAnalyzer->isMC() ) fillLeptonGenVars(mu, *genParticles);
        fillLeptonJetVariables(mu, jets, primaryVertex, *rho);
        fillDisplacedIDVariables(mu);

	    //std::vector<reco::Track> vertex_tracks;
	    //vertex_tracks.push_back(*mu.bestTrack());
	    //fillLeptonKVFVariables(packedCands, vertex_tracks);
        //fillLeptonIVFVariables(*mu.bestTrack(), *secVertices);
        //std::cout << std::endl << "--- muon pt, eta, phi: " << (*mu.bestTrack()).pt() << " " << (*mu.bestTrack()).eta() << " " << (*mu.bestTrack()).phi() << " " << mu.numberOfSourceCandidatePtrs() << " "; 
        //if(mu.numberOfSourceCandidatePtrs() > 0) std::cout << mu.sourceCandidatePtr(0)->pt() << std::endl;
        fillMatchingIVFVariables(*secVertices, mu, primaryVertex);

        _lFlavor[_nL]        = 1;
        _lMuonSegComp[_nL]    = mu.segmentCompatibility();
        _lMuonTrackPt[_nL]    = mu.innerTrack()->pt();
        _lMuonTrackPtErr[_nL] = mu.innerTrack()->ptError();

        _relIso[_nL]         = getRelIso03(mu, *rho);                     // Isolation variables
        _relIso0p4[_nL]      = getRelIso04(mu, *rho);
        _relIso0p4MuDeltaBeta[_nL] = getRelIso04(mu, *rho, true);
        _miniIso[_nL]        = getMiniIsolation( mu, *rho, false );
        _miniIsoCharged[_nL] = getMiniIsolation( mu, *rho, true );

        _lPOGVeto[_nL]       = mu.passed(reco::Muon::CutBasedIdLoose); // no veto available, so we take loose here
        _lPOGLoose[_nL]      = mu.passed(reco::Muon::CutBasedIdLoose);
        _lPOGMedium[_nL]     = mu.passed(reco::Muon::CutBasedIdMedium);
        _lPOGTight[_nL]      = mu.passed(reco::Muon::CutBasedIdTight);

        //fillLeptonJetVariables MUST be called after setting isolation variables since _relIso0p4 is used to set ptRatio in the absence of a closest jet 
        //tZq lepton MVA uses old matching scheme, first set the lepton jet variables with old matching to compute this MVA
        fillLeptonJetVariables(mu, jets, primaryVertex, *rho, true);
        _leptonMvatZq[_nL]   = leptonMvaVal(mu, leptonMvaComputertZq);

        //the TTH MVA uses a newer matching scheme, so we recompute the lepton jet variables, THIS VERSION IS STORED IN THE NTUPLES
        fillLeptonJetVariables(mu, jets, primaryVertex, *rho, false);
        _leptonMvaTTH[_nL]   = leptonMvaVal(mu, leptonMvaComputerTTH);

        ++_nMu;
        ++_nL;
        ++_nLight;
    }

    // Loop over electrons (note: using iterator we can easily get the ref too)
    for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
        if(_nL == nL_max)                                                                               break;
        if(ele->gsfTrack().isNull())                                                                    continue;
        if(ele->pt() < 7)                                                                               continue;
        if(fabs(ele->eta()) > 2.5)                                                                      continue;
        //if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2) continue;
        fillLeptonImpactParameters(*ele);
        //if(fabs(_dxy[_nL]) > 0.05)                                                               continue;   //same argument as for muons, dont harm displaced samples and can be applied at root level
        //if(fabs(_dz[_nL]) > 0.1)                                                                 continue;
        fillLeptonKinVars(*ele);
        if( multilepAnalyzer->isMC() ) fillLeptonGenVars(*ele, *genParticles);
        fillLeptonJetVariables(*ele, jets, primaryVertex, *rho);
        fillDisplacedIDVariables(*ele);
	
	    //std::vector<reco::Track> vertex_tracks;
	    //vertex_tracks.push_back(*ele->gsfTrack());
	    //fillLeptonKVFVariables(packedCands, vertex_tracks); 
        //fillLeptonIVFVariables(*ele, secVertices);
        //std::cout << std::endl << "--- electron gsf pt, eta, phi, ptr: " << (*ele->gsfTrack()).pt() << " " << (*ele->gsfTrack()).eta() << " " << (*ele->gsfTrack()).phi() << " -" << ele->associatedPackedPFCandidates().size() << "- " << std::endl;
        //for(edm::Ref<pat::PackedCandidateCollection> cand : ele->associatedPackedPFCandidates()){
        //    std::cout << "--- sourcePFCand pt, eta, phi, charge: " << cand->pt() << " " << cand->eta() << " " << cand->phi() << " " << cand->charge() << std::endl; 
        //}
        fillMatchingIVFVariables(*secVertices, *ele, primaryVertex);

        _lFlavor[_nL]                   = 0;
        _lEtaSC[_nL]                    = ele->superCluster()->eta();

        _relIso[_nL]                    = getRelIso03(*ele, *rho);
        _relIso0p4[_nL]                 = getRelIso04(*ele, *rho);
        _miniIso[_nL]                   = getMiniIsolation( *ele, *rho, false );
        _miniIsoCharged[_nL]            = getMiniIsolation( *ele, *rho, true );
        _lElectronMvaSummer16GP[_nL]    = ele->userFloat("ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"); // OLD, do not use it
        _lElectronMvaSummer16HZZ[_nL]   = ele->userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values"); // OLD, do not use it
        _lElectronMvaFall17v1NoIso[_nL] = ele->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values"); // OLD, do not use it
        _lElectronMvaFall17Iso[_nL]     = ele->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
        _lElectronMvaFall17NoIso[_nL]   = ele->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values");
        _lElectronPassEmu[_nL]          = passTriggerEmulationDoubleEG(&*ele);                             // Keep in mind, this trigger emulation is for 2016 DoubleEG, the SingleEG trigger emulation is different
        _lElectronPassConvVeto[_nL]     = ele->passConversionVeto();
        _lElectronChargeConst[_nL]      = ele->isGsfCtfScPixChargeConsistent();
        _lElectronMissingHits[_nL]      = ele->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

        _lPOGVeto[_nL]                  = ele->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
        _lPOGLoose[_nL]                 = ele->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
        _lPOGMedium[_nL]                = ele->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
        _lPOGTight[_nL]                 = ele->electronID("cutBasedElectronID-Fall17-94X-V1-tight");

        //fillLeptonJetVariables MUST be called after setting isolation variables since _relIso0p4 is used to set ptRatio in the absence of a closest jet 
        //tZq lepton MVA uses old matching scheme, first set the lepton jet variables with old matching to compute this MVA
        fillLeptonJetVariables(*ele, jets, primaryVertex, *rho, true);
        _leptonMvatZq[_nL]              = leptonMvaVal(*ele, leptonMvaComputertZq);

        //the TTH MVA uses a newer matching scheme, so we recompute the lepton jet variables, THIS VERSION IS STORED IN THE NTUPLES
        fillLeptonJetVariables(*ele, jets, primaryVertex, *rho, false);
        _leptonMvaTTH[_nL]              = leptonMvaVal(*ele, leptonMvaComputerTTH);

        // Note: for the scale and smearing systematics we use the overall values, assuming we are not very sensitive to these systematics
        // In case these systematics turn out to be important, need to add their individual source to the tree (and propagate to their own templates):
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
        _lPtCorr[_nL]                   = ele->pt()*ele->userFloat("ecalTrkEnergyPostCorr")/ele->energy();
        _lPtScaleUp[_nL]                = ele->pt()*ele->userFloat("energyScaleUp")/ele->energy();
        _lPtScaleDown[_nL]              = ele->pt()*ele->userFloat("energyScaleDown")/ele->energy();
        _lPtResUp[_nL]                  = ele->pt()*ele->userFloat("energySigmaUp")/ele->energy();
        _lPtResDown[_nL]                = ele->pt()*ele->userFloat("energySigmaDown")/ele->energy();
        _lECorr[_nL]                    = ele->userFloat("ecalTrkEnergyPostCorr");
        _lEScaleUp[_nL]                 = ele->userFloat("energyScaleUp");
        _lEScaleDown[_nL]               = ele->userFloat("energyScaleDown");
        _lEResUp[_nL]                   = ele->userFloat("energySigmaUp");
        _lEResDown[_nL]                 = ele->userFloat("energySigmaDown");

        ++_nEle;
        ++_nL;
        ++_nLight;
    }

    //Initialize with default values for those electron-only arrays which weren't filled with muons [to allow correct comparison by the test script]
    for(auto array : {&_lEtaSC}) std::fill_n(*array, _nMu, 0.);
    for(auto array : {&_lElectronMvaSummer16GP, &_lElectronMvaSummer16HZZ, &_lElectronMvaFall17v1NoIso}) std::fill_n(*array, _nMu, 0.); // OLD, do not use them
    for(auto array : {&_lElectronMvaFall17Iso, &_lElectronMvaFall17NoIso}) std::fill_n(*array, _nMu, 0.);
    for(auto array : {&_lElectronPassEmu, &_lElectronPassConvVeto, &_lElectronChargeConst}) std::fill_n(*array, _nMu, false);
    for(auto array : {&_lElectronMissingHits}) std::fill_n(*array, _nMu, 0.);
    for(auto array : {&_lPtCorr, &_lPtScaleUp, &_lPtScaleDown, &_lPtResUp, &_lPtResDown}) std::fill_n(*array, _nMu, 0.);
    for(auto array : {&_lECorr, &_lEScaleUp, &_lEScaleDown, &_lEResUp, &_lEResDown}) std::fill_n(*array, _nMu, 0.);

    //loop over taus
    for(const pat::Tau& tau : *taus){
        if(_nL == nL_max)         break;
        if(tau.pt() < 20)         continue;          // Minimum pt for tau reconstruction
        if(fabs(tau.eta()) > 2.3) continue;
        fillLeptonKinVars(tau);
        
        if(multilepAnalyzer->isMC()) fillTauGenVars(tau, *genParticles);                    //Still needs to be tested
        fillLeptonImpactParameters(tau, primaryVertex);

        _lFlavor[_nL]  = 2;
        _tauDecayMode[_nL] = tau.decayMode();
        _tauMuonVetoLoose[_nL] = tau.tauID("againstMuonLoose3");                                        //Light lepton vetos
        _tauEleVetoLoose[_nL] = tau.tauID("againstElectronLooseMVA6");

        _decayModeFinding[_nL] = tau.tauID("decayModeFinding");                           //old tau ID

        _lPOGVeto[_nL] = tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
        _lPOGLoose[_nL] = tau.tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017");
        _lPOGMedium[_nL] = tau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017");
        _lPOGTight[_nL] = tau.tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017");
        _tauPOGVTight2017v2[_nL] = tau.tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017");
        
        _tauAgainstElectronMVA6Raw[_nL] = tau.tauID("againstElectronMVA6Raw");
        _tauCombinedIsoDBRaw3Hits[_nL]  = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        _tauIsoMVAPWdR03oldDMwLT[_nL]   = tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw");
        _tauIsoMVADBdR03oldDMwLT[_nL]   = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        _tauIsoMVADBdR03newDMwLT[_nL]   = tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        _tauIsoMVAPWnewDMwLT[_nL]       = tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        _tauIsoMVAPWoldDMwLT[_nL]       = tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");

        if(multilepAnalyzer->storeAllTauID){
            _tauMuonVetoTight[_nL] = tau.tauID("againstMuonTight3");                                        //Light lepton vetos
            _tauEleVetoVLoose[_nL] = tau.tauID("againstElectronVLooseMVA6");
            _tauEleVetoMedium[_nL] = tau.tauID("againstElectronMediumMVA6");
            _tauEleVetoTight[_nL] = tau.tauID("againstElectronTightMVA6");
            _tauEleVetoVTight[_nL] = tau.tauID("againstElectronVTightMVA6");
            
            _tauPOGVLoose2015[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");                        
            _tauPOGLoose2015[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
            _tauPOGMedium2015[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
            _tauPOGTight2015[_nL] = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
            _tauPOGVTight2015[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
            
            _tauPOGVVLoose2017v2[_nL] = tau.tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
            _tauPOGVVTight2017v2[_nL] = tau.tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017");
            
            _decayModeFindingNew[_nL]       = tau.tauID("decayModeFindingNewDMs");                   //new Tau ID
            _tauVLooseMvaNew2015[_nL]           = tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
            _tauLooseMvaNew2015[_nL]            = tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
            _tauMediumMvaNew2015[_nL]           = tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
            _tauTightMvaNew2015[_nL]            = tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT");
            _tauVTightMvaNew2015[_nL]           = tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT");

            _tauVLooseMvaNew2017v2[_nL]           = tau.tauID("byVLooseIsolationMVArun2017v2DBnewDMwLT2017");
            _tauLooseMvaNew2017v2[_nL]            = tau.tauID("byLooseIsolationMVArun2017v2DBnewDMwLT2017");
            _tauMediumMvaNew2017v2[_nL]           = tau.tauID("byMediumIsolationMVArun2017v2DBnewDMwLT2017");
            _tauTightMvaNew2017v2[_nL]            = tau.tauID("byTightIsolationMVArun2017v2DBnewDMwLT2017");
            _tauVTightMvaNew2017v2[_nL]           = tau.tauID("byVTightIsolationMVArun2017v2DBnewDMwLT2017");
    
        }

	    fillLeptonJetVariables(tau, jets, primaryVertex, *rho);

        ++_nTau;
        ++_nL;
    }
    
    
    //Preselection of one POGMedium lepton with most of the cuts (in order to reduce file size of the root tuples)
    //bool tightlepton = false;
    //for(unsigned i = 0; i < _nLight; i++){
    //    if(_lFlavor[i] == 0 and _lPt[i] > 24 and fabs(_lEta[i]) < 2.5 and _lPOGMedium[i] and _relIso[i] > 0.1 and fabs(_dxy[i]) < 0.05 and _3dIPSig[i] < 4){ tightlepton = true; break;}
    //    if(_lFlavor[i] == 1 and _lPt[i] > 24 and fabs(_lEta[i]) < 2.4 and _lPOGMedium[i] and _relIso[i] > 0.1 and fabs(_dxy[i]) < 0.05 and _3dIPSig[i] < 4){ tightlepton = true; break;}
    //}
    //
    //if(tightlepton) tightlepton = true;

    //Initialize with default values for those tau-only arrays which weren't filled with electrons and muons [to allow correct comparison by the test script]
    for(auto array : {&_tauMuonVetoLoose, &_tauEleVetoLoose, &_decayModeFinding, &_decayModeFindingNew}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauEleVetoVLoose, &_tauEleVetoMedium, &_tauEleVetoTight, &_tauEleVetoVTight, &_tauMuonVetoTight, &_decayModeFindingNew}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauVLooseMvaNew2015, &_tauLooseMvaNew2015, &_tauMediumMvaNew2015, &_tauTightMvaNew2015, &_tauVTightMvaNew2015}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauVLooseMvaNew2017v2, &_tauLooseMvaNew2017v2, &_tauMediumMvaNew2017v2, &_tauTightMvaNew2017v2, &_tauVTightMvaNew2017v2}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauPOGVLoose2015, &_tauPOGLoose2015, &_tauPOGMedium2015, &_tauPOGTight2015, &_tauPOGVTight2015}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauPOGVVLoose2017v2, &_tauPOGVTight2017v2, &_tauPOGVVTight2017v2}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauAgainstElectronMVA6Raw, &_tauCombinedIsoDBRaw3Hits, &_tauIsoMVAPWdR03oldDMwLT}) std::fill_n(*array, _nLight, 0.);
    for(auto array : {&_tauDecayMode}) std::fill_n(*array, _nLight, 0.);
    for(auto array : {&_tauIsoMVADBdR03oldDMwLT, &_tauIsoMVADBdR03newDMwLT, &_tauIsoMVAPWnewDMwLT, &_tauIsoMVAPWoldDMwLT}) std::fill_n(*array, _nLight, 0.);

    if(multilepAnalyzer->skim == "trilep"    &&  _nL     < 3) return false;
    if(multilepAnalyzer->skim == "dilep"     &&  _nLight < 2 /*|| !tightlepton)*/) return false;
    if(multilepAnalyzer->skim == "ttg"       &&  _nLight < 2) return false;
    if(multilepAnalyzer->skim == "singlelep" &&  _nL     < 1) return false;
    if(multilepAnalyzer->skim == "singletau" &&  _nTau   < 1) return false;
    if(multilepAnalyzer->skim == "FR"        &&  _nLight < 1) return false;
    
    return true;
}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton){
    _lPt[_nL]     = lepton.pt();
    _lEta[_nL]    = lepton.eta();
    _lPhi[_nL]    = lepton.phi();
    _lE[_nL]      = lepton.energy();
    _lCharge[_nL] = lepton.charge();
}

template <typename Lepton> void LeptonAnalyzer::fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* match = lepton.genParticle();
    if(!match || match->pdgId() != lepton.pdgId()) match = GenTools::geometricMatch(lepton, genParticles); // if no match or pdgId is different, try the geometric match

    _tauGenStatus[_nL]          = TauTools::tauGenStatus(match);        
    _lIsPrompt[_nL]             = match && (abs(lepton.pdgId()) == abs(match->pdgId()) || match->pdgId() == 22) && GenTools::isPrompt(*match, genParticles); // only when matched to its own flavor or a photon
    _lMatchPdgId[_nL]           = match != nullptr ? match->pdgId() : 0;
    _lMatchCharge[_nL]          = match != nullptr ? match->charge() : 0;
    _lProvenance[_nL]           = GenTools::provenance(match, genParticles);
    _lProvenanceCompressed[_nL] = GenTools::provenanceCompressed(match, genParticles, _lIsPrompt[_nL]);
    _lProvenanceConversion[_nL] = GenTools::provenanceConversion(match, genParticles);
    _lMomPdgId[_nL]             = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;
}

void LeptonAnalyzer::fillTauGenVars(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles){
    
    const reco::GenParticle* match = TauTools::findMatch(tau, genParticles);

    _tauGenStatus[_nL]          = TauTools::tauGenStatus(match);        
    _lIsPrompt[_nL]             = match && _tauGenStatus[_nL] != 6; 
    _lMatchPdgId[_nL]           = match ? match->pdgId() : 0;
    _lProvenance[_nL]           = GenTools::provenance(match, genParticles);
    _lProvenanceCompressed[_nL] = GenTools::provenanceCompressed(match, genParticles, _lIsPrompt[_nL]);
    _lProvenanceConversion[_nL] = GenTools::provenanceConversion(match, genParticles);
    _lMomPdgId[_nL]             = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;

}

/*
 * Impact parameters:
 * Note: dB function seems to be preferred and more accurate over track->dxy and dz functions 
 * as the latter ones have some simplified extrapolation used (leading to slightly different values or opposite sign)
 * For taus: dxy is pre-computed with PV it was constructed with
 */
void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Electron& ele){
    _dxy[_nL]     = ele.dB(pat::Electron::PV2D);
    _dz[_nL]      = ele.dB(pat::Electron::PVDZ);
    _3dIP[_nL]    = ele.dB(pat::Electron::PV3D);
    _3dIPSig[_nL] = fabs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D));
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon){
    _dxy[_nL]     = muon.dB(pat::Muon::PV2D);
    _dz[_nL]      = muon.dB(pat::Muon::PVDZ);
    _3dIP[_nL]    = muon.dB(pat::Muon::PV3D);
    _3dIPSig[_nL] = fabs(muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D));
}

void LeptonAnalyzer::fillDisplacedIDVariables(const pat::Electron& ele){
    _lEleIsEB [_nL] = ele.isEB();
    _lEleIsEE[_nL] = ele.isEE();
    _lEleSuperClusterOverP[_nL] = ele.eSuperClusterOverP();
    _lEleEcalEnergy[_nL]= ele.ecalEnergy();
    _lElefull5x5SigmaIetaIeta[_nL] = ele.full5x5_sigmaIetaIeta();
    _lEleDEtaInSeed[_nL] = std::abs(dEtaInSeed(&ele));
    _lEleDeltaPhiSuperClusterTrackAtVtx[_nL] = std::abs(ele.deltaPhiSuperClusterTrackAtVtx());
    _lElehadronicOverEm[_nL] = ele.hadronicOverEm();
    _lEleInvMinusPInv[_nL] = std::abs(1.0 - ele.eSuperClusterOverP())/ele.ecalEnergy();
    _eleNumberInnerHitsMissing[_nL]=ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
}   

void LeptonAnalyzer::fillDisplacedIDVariables(const pat::Muon& mu){
    _lGlobalMuon[_nL] = mu.isGlobalMuon();
    _lTrackerMuon[_nL]= mu.isTrackerMuon();
    _lInnerTrackValidFraction[_nL] = (!mu.innerTrack().isNull()) ?   mu.innerTrack()->validFraction()  : -1;
    _lGlobalTrackNormalizedChi2[_nL]= (!mu.globalTrack().isNull()) ?   mu.globalTrack()->normalizedChi2()  : -1;
    _lCQChi2Position[_nL] = mu.combinedQuality().chi2LocalPosition;
    _lCQTrackKink[_nL] = mu.combinedQuality().trkKink;
    _lNumberOfMatchedStation[_nL] = mu.numberOfMatchedStations();
    _lNumberOfValidPixelHits[_nL] = (!mu.innerTrack().isNull()) ?   mu.innerTrack()->hitPattern().numberOfValidPixelHits()  : 0; // cannot be -1 !!
    _lTrackerLayersWithMeasurement[_nL] = (!mu.innerTrack().isNull()) ?   mu.innerTrack()->hitPattern().trackerLayersWithMeasurement()  : 0; // cannot be -1 !! 
    _muNumberInnerHits[_nL]= (!mu.globalTrack().isNull()) ?   mu.globalTrack()->hitPattern().numberOfValidMuonHits() : (!mu.outerTrack().isNull() ? mu.outerTrack()->hitPattern().numberOfValidMuonHits() : 0); // cannot be -1 !!!
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Tau& tau, const reco::Vertex& vertex){
    _dxy[_nL]     = (double) tau.dxy();                                      // warning: float while dxy of tracks are double; could also return -1000
    _dz[_nL]      = tau_dz(tau, vertex.position());
    _3dIP[_nL]    = tau.ip3d();
    _3dIPSig[_nL] = tau.ip3d_Sig();
}

//Function returning tau dz
double LeptonAnalyzer::tau_dz(const pat::Tau& tau, const reco::Vertex::Point& vertex) const{
    const reco::Candidate::Point& tauVtx = tau.leadChargedHadrCand()->vertex();
    return (tauVtx.Z() - vertex.z()) - ((tauVtx.X() - vertex.x())*tau.px()+(tauVtx.Y()-vertex.y())*tau.py())/tau.pt()*tau.pz()/tau.pt();
}


//fucking disgusting lepton-jet matching based on obscure pointer defined somewhere in the abyss of CMSSW. (I feel ashamed for using this, but one has to be in sync with the POG - Willem )
template< typename T1, typename T2 > bool isSourceCandidatePtrMatch( const T1& lhs, const T2& rhs ){
    for( size_t lhsIndex = 0; lhsIndex < lhs.numberOfSourceCandidatePtrs(); ++lhsIndex ){
        auto lhsSourcePtr = lhs.sourceCandidatePtr( lhsIndex );
        for( size_t rhsIndex = 0; rhsIndex < rhs.numberOfSourceCandidatePtrs(); ++rhsIndex ){
            auto rhsSourcePtr = rhs.sourceCandidatePtr( rhsIndex );
            if( lhsSourcePtr == rhsSourcePtr ){
                return true;
            }
        }
    }
    return false;
}


const pat::Jet* findMatchedJet( const reco::Candidate& lepton, const edm::Handle< std::vector< pat::Jet > >& jets, const bool oldMatching ){
    
    //Look for jet that matches with lepton
    const pat::Jet* matchedJetPtr = nullptr;

    //old matching scheme looks for closest jet in terms of delta R, and required this to be within delta R 0.4 of the lepton
    if( oldMatching ){
        for( auto& jet : *jets ){
            if( jet.pt() <= 5 || fabs( jet.eta() ) >= 3 ) continue;
            if( ( matchedJetPtr == nullptr) || reco::deltaR( jet, lepton ) < reco::deltaR( *matchedJetPtr, lepton ) ){
                matchedJetPtr = &jet;
            }
        }
        if( matchedJetPtr != nullptr && reco::deltaR( lepton, *matchedJetPtr ) > 0.4 ){
            matchedJetPtr = nullptr;
        }
    } else {
        for( auto& jet : *jets ){
            if( isSourceCandidatePtrMatch( lepton, jet ) ){

                //immediately returning guarantees that the leading jet matched to the lepton is returned
                return &jet;
            }
        }
    }
    return matchedJetPtr;
}


//compute closest jet variables using new matching scheme 
void LeptonAnalyzer::fillLeptonJetVariables( const reco::Candidate& lepton, edm::Handle< std::vector< pat::Jet > >& jets, const reco::Vertex& vertex, const double rho, const bool oldMatching ){

    //find closest jet based on source candidate pointer matching
    const pat::Jet* matchedJetPtr = findMatchedJet( lepton, jets, oldMatching );

    if( matchedJetPtr == nullptr ){
        if( _lFlavor[_nL] == 1 ){
            _ptRatio[_nL] = ( oldMatching ? 1. : 1. / ( 1. + _relIso0p4MuDeltaBeta[_nL] ) );
        } else{
            _ptRatio[_nL] = ( oldMatching ? 1. : 1. / ( 1. + _relIso0p4[_nL] ) );
        }
        _ptRel[_nL] = 0;
        _selectedTrackMult[_nL] = 0;
        _closestJetDeepFlavor_b[_nL] = 0;
        _closestJetDeepFlavor_bb[_nL] = 0;
        _closestJetDeepFlavor_lepb[_nL] = 0;
        _closestJetDeepFlavor[_nL] = 0;
        _closestJetDeepCsv_b[_nL] = 0;
        _closestJetDeepCsv_bb[_nL] = 0;
        _closestJetDeepCsv[_nL] = 0;
        _closestJetCsvV2[_nL] = 0;
    } else {
        const pat::Jet& jet = *matchedJetPtr;

        auto rawJetP4 = jet.correctedP4("Uncorrected"); 
        auto leptonP4 = lepton.p4();

        bool leptonEqualsJet = ( ( rawJetP4 - leptonP4 ).P() < 1e-4 );

        //if lepton and jet vector are equal set _ptRatio, _ptRel and track multipliticy to defaults 
        if( leptonEqualsJet && !oldMatching ){
            _ptRatio[_nL] = 1;
            _ptRel[_nL] = 0;
            _selectedTrackMult[_nL] = 0;
        } else {

            //remove all corrections above L1 from the lepton
            auto L1JetP4 = jet.correctedP4("L1FastJet");
            double L2L3JEC = jet.pt()/L1JetP4.pt(); 
            auto lepAwareJetP4 = ( L1JetP4 - leptonP4 )*L2L3JEC + leptonP4;

            _ptRatio[_nL] = lepton.pt() / lepAwareJetP4.pt();

            //lepton momentum orthogonal to the jet axis
            //magnitude of cross-product between lepton and rest of jet 
            _ptRel[_nL ] = leptonP4.Vect().Cross( (lepAwareJetP4 - leptonP4 ).Vect().Unit() ).R();

            _selectedTrackMult[_nL] = 0;
            for( const auto daughterPtr : jet.daughterPtrVector() ){
                const pat::PackedCandidate& daughter = *( (const pat::PackedCandidate*) daughterPtr.get() );
            
                if( daughter.charge() == 0 ) continue;
                if( daughter.fromPV() < 2 ) continue;
                if( reco::deltaR( daughter, lepton ) > 0.4 ) continue;
                if( !daughter.hasTrackDetails() ) continue;

                auto daughterTrack = daughter.pseudoTrack();
                if( daughterTrack.pt() <= 1 ) continue;
                if( daughterTrack.hitPattern().numberOfValidHits() < 8 ) continue;
                if( daughterTrack.hitPattern().numberOfValidPixelHits() < 2 ) continue;
                if( daughterTrack.normalizedChi2() >= 5 ) continue;
                if( std::abs( daughterTrack.dz( vertex.position() ) ) >= 17 ) continue;
                if( std::abs( daughterTrack.dxy( vertex.position() ) ) >= 0.2 ) continue;
                ++_selectedTrackMult[_nL];
            }

        }

        //CSVv2 of closest jet
        _closestJetCsvV2[_nL]      = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        //DeepCSV of closest jet
        _closestJetDeepCsv_b[_nL]  = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _closestJetDeepCsv[_nL]    = _closestJetDeepCsv_b[_nL] + _closestJetDeepCsv_bb[_nL];
        if( std::isnan( _closestJetDeepCsv[_nL] ) ) _closestJetDeepCsv[_nL] = 0.;

        //DeepFlavor b-tag values of closest jet
        _closestJetDeepFlavor_b[_nL] = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
        _closestJetDeepFlavor_bb[_nL] = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
        _closestJetDeepFlavor_lepb[_nL] = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        _closestJetDeepFlavor[_nL] = _closestJetDeepFlavor_b[_nL] + _closestJetDeepFlavor_bb[_nL] + _closestJetDeepFlavor_lepb[_nL];
        if( std::isnan( _closestJetDeepFlavor[_nL] ) ) _closestJetDeepFlavor[_nL] = 0.;

        //compute selected track multiplicity of closest jet
        _selectedTrackMult[_nL] = 0;
        for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
            const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
            if(daughter->hasTrackDetails()){
                const reco::Track& daughterTrack = daughter->pseudoTrack();
                if(daughterTrack.pt() <= 1)                                 continue;
                if(daughterTrack.charge() == 0)                             continue;
                if(daughter->fromPV() < 2)                                  continue;
                if(daughterTrack.hitPattern().numberOfValidHits() < 8)      continue;
                if(daughterTrack.hitPattern().numberOfValidPixelHits() < 2) continue;
                if(daughterTrack.normalizedChi2() >= 5)                     continue;
                if(fabs(daughterTrack.dz(vertex.position())) >= 17)         continue;
                if(fabs(daughterTrack.dxy(vertex.position())) >= 0.2)       continue;


                //distance from jet core
                TLorentzVector trackVec(daughterTrack.px(), daughterTrack.py(), daughterTrack.pz(), daughterTrack.p());
                double daughterDeltaR = trackVec.DeltaR(jV);

                if(daughterDeltaR <= 0.4) ++_selectedTrackMult[_nL];
            } 
        }
    }
}


void LeptonAnalyzer::fillMatchingIVFVariables(const std::vector<reco::Vertex>& secVertices, const pat::Muon& muon, const reco::Vertex& primaryVertex){
    _lIVF_match[_nLight] = false;
    bool new_vtx = false;
    double ptdiff, normchi2;
    double minptdiff = 10, minnormchi2 = 10000;
    for(const reco::Vertex& vtx : secVertices){
        for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++){
            reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
            ptdiff   = fabs(muon.sourceCandidatePtr(0)->pt() - vtxTrack->pt());
            normchi2 = fabs(vtx.chi2()/vtx.ndof());
            if(ptdiff < 0.001 and (ptdiff < minptdiff or (ptdiff == minptdiff and normchi2 < minnormchi2))){
                new_vtx = true;
                _lIVF_match[_nLight] = true;
                minptdiff   = ptdiff;
                minnormchi2 = normchi2;
            }
        }
        if(new_vtx){
            _IVF_x[_nLight]    = vtx.x();
            _IVF_y[_nLight]    = vtx.y();
            _IVF_z[_nLight]    = vtx.z();
            _IVF_cx[_nLight]   = vtx.xError();
            _IVF_cy[_nLight]   = vtx.yError();
            _IVF_cz[_nLight]   = vtx.zError();
            _IVF_df[_nLight]   = vtx.ndof();
            _IVF_chi2[_nLight] = vtx.chi2();
            _IVF_pt[_nLight]   = vtx.p4().pt();
            _IVF_eta[_nLight]  = vtx.p4().eta();
            _IVF_phi[_nLight]  = vtx.p4().phi();
            _IVF_E[_nLight]    = vtx.p4().energy();
            _IVF_mass[_nLight] = vtx.p4().mass();
            
            _IVF_ntracks[_nLight] = 0;
            for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++){
                if(_IVF_ntracks[_nLight] == ntrack_max) break;
                reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
                _IVF_trackpt[_nLight][_IVF_ntracks[_nLight]]        = vtxTrack->pt();
                _IVF_tracketa[_nLight][_IVF_ntracks[_nLight]]       = vtxTrack->eta();
                _IVF_trackphi[_nLight][_IVF_ntracks[_nLight]]       = vtxTrack->phi();
                _IVF_trackE[_nLight][_IVF_ntracks[_nLight]]         = vtxTrack->p();
                _IVF_trackcharge[_nLight][_IVF_ntracks[_nLight]]    = vtxTrack->charge();
                _IVF_trackdxy[_nLight][_IVF_ntracks[_nLight]]       = std::abs(vtxTrack->dxy(primaryVertex.position()));
                _IVF_trackdz[_nLight][_IVF_ntracks[_nLight]]        = std::abs(vtxTrack->dz(primaryVertex.position()));
                _IVF_ntracks[_nLight]++;
            }
            new_vtx = false;
        }
    }
}

void LeptonAnalyzer::fillMatchingIVFVariables(const std::vector<reco::Vertex>& secVertices, const pat::Electron& electron, const reco::Vertex& primaryVertex){
    _lIVF_match[_nLight] = false;
    bool new_vtx = false;
    double dR, deta, normchi2;
    double mindR = 20, minnormchi2 = 10000;
    for(const reco::Vertex& vtx : secVertices){
        for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++){
            reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
            for(edm::Ref<pat::PackedCandidateCollection> cand : electron.associatedPackedPFCandidates()){
                dR       = reco::deltaR(cand->eta(), cand->phi(), vtxTrack->eta(), vtxTrack->phi());
                deta     = fabs(cand->eta() - vtxTrack->eta());
                normchi2 = fabs(vtx.chi2()/vtx.ndof());
                if((dR < 0.05 or (dR < 0.1 and deta < 0.03)) and (dR < mindR or (dR == mindR and normchi2 < minnormchi2))){
                    new_vtx = true;
                    _lIVF_match[_nLight] = true; 
                    mindR       = dR;
                    minnormchi2 = normchi2;
                }
            }
        }
        if(new_vtx){
            _IVF_x[_nLight]    = vtx.x();
            _IVF_y[_nLight]    = vtx.y();
            _IVF_z[_nLight]    = vtx.z();
            _IVF_cx[_nLight]   = vtx.xError();
            _IVF_cy[_nLight]   = vtx.yError();
            _IVF_cz[_nLight]   = vtx.zError();
            _IVF_df[_nLight]   = vtx.ndof();
            _IVF_chi2[_nLight] = vtx.chi2();
            _IVF_pt[_nLight]   = vtx.p4().pt();
            _IVF_eta[_nLight]  = vtx.p4().eta();
            _IVF_phi[_nLight]  = vtx.p4().phi();
            _IVF_E[_nLight]    = vtx.p4().energy();
            _IVF_mass[_nLight] = vtx.p4().mass();
            
            _IVF_ntracks[_nLight] = 0;
            for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++){
                if(_IVF_ntracks[_nLight] == ntrack_max) break;
                reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
                _IVF_trackpt[_nLight][_IVF_ntracks[_nLight]]        = vtxTrack->pt();
                _IVF_tracketa[_nLight][_IVF_ntracks[_nLight]]       = vtxTrack->eta();
                _IVF_trackphi[_nLight][_IVF_ntracks[_nLight]]       = vtxTrack->phi();
                _IVF_trackE[_nLight][_IVF_ntracks[_nLight]]         = vtxTrack->p();
                _IVF_trackcharge[_nLight][_IVF_ntracks[_nLight]]    = vtxTrack->charge();
                _IVF_trackdxy[_nLight][_IVF_ntracks[_nLight]]       = std::abs(vtxTrack->dxy(primaryVertex.position()));
                _IVF_trackdz[_nLight][_IVF_ntracks[_nLight]]        = std::abs(vtxTrack->dz(primaryVertex.position()));
                _IVF_ntracks[_nLight]++;
            }
            new_vtx = false;
        }
    }
}

//TransientVertex LeptonAnalyzer::constructKalmanVertex(std::vector<reco::Track>& tracks, MagneticField* bfield){
//    std::vector<reco::TransientTrack> ttks;
//    for(auto track : tracks){
//	ttks.push_back(reco::TransientTrack(track, bfield));
//    }
//    KalmanVertexFitter vtxFitter;
//    return vtxFitter.vertex(ttks);
//}
//
//void LeptonAnalyzer::fillLeptonKVFVariables(edm::Handle<std::vector<pat::PackedCandidate>>& packedCands, std::vector<reco::Track>& tracks){
//    MagneticField *bfield = new OAEParametrizedMagneticField("3_8T");
//
//    const GlobalPoint*  lepPoint  = new GlobalPoint(tracks[0].referencePoint().x(), tracks[0].referencePoint().y(), tracks[0].referencePoint().z()); // Calc. traj.param. for lepton, used to calculate PCA
//    const GlobalVector* lepVector = new GlobalVector(tracks[0].px(), tracks[0].py(), tracks[0].pz());
//    TrackCharge lepCharge         = tracks[0].charge();
//    GlobalTrajectoryParameters lepParam(*lepPoint, *lepVector, lepCharge, bfield);
//
//    double min_tracks_dR = 5;
//    int i_duplicated_lep = -1;
//    
//    std::vector<std::pair<double, reco::Track>> temp_tracks;
//    //for(; dRcut < 1.1 && tracks.size() == 1; dRcut += 0.1){
//        for(auto cand = packedCands->cbegin(); cand != packedCands->cend(); ++cand){
//            if(!cand->hasTrackDetails()) continue;
//            const reco::Track& candTrack = cand->pseudoTrack();
//
//            const GlobalPoint*  candPoint  = new GlobalPoint(candTrack.referencePoint().x(), candTrack.referencePoint().y(), candTrack.referencePoint().z());
//            const GlobalVector* candVector = new GlobalVector(candTrack.px(), candTrack.py(), candTrack.pz());
//            TrackCharge candCharge         = candTrack.charge();
//            GlobalTrajectoryParameters candParam(*candPoint, *candVector, candCharge, bfield);
//            TwoTrackMinimumDistance TTMinDist;
//                
//            double PCA_distance = 5;
//            double PCA_dphi     = 5;
//            if(TTMinDist.calculate(lepParam, candParam)){
//                PCA_distance = TTMinDist.distance();
//                PCA_dphi     = fabs(TTMinDist.firstAngle() - TTMinDist.secondAngle());
//                if(PCA_dphi > 3.14) PCA_dphi = 6.28 - PCA_dphi;
//            }
//            
//            double tracks_dz    = fabs(tracks[0].dz() - candTrack.dz());
//            double tracks_deta  = fabs(tracks[0].eta() - candTrack.eta());
//            double tracks_dR    = sqrt(tracks_deta*tracks_deta + PCA_dphi*PCA_dphi);
//
//            bool goodTrack      =   tracks_dR < 1 &&
//                                    PCA_distance < 1 && 
//	            			        candTrack.normalizedChi2() < 5 &&
//                                    candTrack.pt() > 0.95 &&
//                                    tracks_dz < 8;
//
//            if(goodTrack){
//                temp_tracks.push_back(std::pair<double, reco::Track>(tracks_dR, candTrack));
//            }
//
//            if(goodTrack && fabs(tracks[0].pt() - candTrack.pt()) < 0.5 && tracks_dR < 0.05 && tracks_dR < min_tracks_dR){ //find original lepton in PF candidates 
//                min_tracks_dR = tracks_dR;
//                i_duplicated_lep = tracks.size() - 1;
//            }
//        }
//    //}dRcut += -0.1;
//    if(i_duplicated_lep != -1) temp_tracks.erase(temp_tracks.begin() + i_duplicated_lep); //remove the duplicated original lepton
//
//    _lKVF_trackPt[_nL][0] = tracks[0].pt();
//    _lKVF_trackEta[_nL][0] = tracks[0].eta();
//    _lKVF_trackPhi[_nL][0] = tracks[0].phi();
//    _lKVF_trackE[_nL][0] = tracks[0].p();
//
//    _lKVF_trackdR[_nL][0] = 0;
//    _lKVF_trackdxy[_nL][0] = tracks[0].dxy();
//    _lKVF_trackdz[_nL][0] = tracks[0].dz();
//    double dRcut = 0.3;
//    while(tracks.size() < 2 && dRcut < 1){
//        dRcut += 0.1;
//        for(std::pair<double, reco::Track> track : temp_tracks){
//            if(track.first < dRcut){ 
//                tracks.push_back(track.second);
//                _lKVF_trackPt[_nL][tracks.size()-1] = track.second.pt();
//                _lKVF_trackEta[_nL][tracks.size()-1] = track.second.eta();
//                _lKVF_trackPhi[_nL][tracks.size()-1] = track.second.phi();
//                _lKVF_trackE[_nL][tracks.size()-1] = track.second.p();
//
//                _lKVF_trackdR[_nL][tracks.size()-1] = track.first;
//                _lKVF_trackdxy[_nL][tracks.size()-1] = track.second.dxy();
//                _lKVF_trackdz[_nL][tracks.size()-1] = track.second.dz();
//            }
//        }
//    }
//    if(tracks.size() < 2) dRcut = 0;
//
//    _lKVF_x[_nL]                 = 0;
//    _lKVF_y[_nL]                 = 0;
//    _lKVF_z[_nL]                 = 0;
//    _lKVF_cxx[_nL]               = 0;
//    _lKVF_cyy[_nL]               = 0;
//    _lKVF_czz[_nL]               = 0;
//    _lKVF_cyx[_nL]               = 0;
//    _lKVF_czy[_nL]               = 0;
//    _lKVF_czx[_nL]               = 0;
//    _lKVF_df[_nL]                = 0;
//    _lKVF_chi2[_nL]              = 0;
//    _lKVF_ntracks[_nL]           = tracks.size();
//    _lKVF_dRcut[_nL]             = dRcut;
//
//    //if(tracks.size() == 1) 
//    if(tracks.size() == 1){
//        _lKVF_valid[_nL] = false;
//        return;
//    } else{
//        TransientVertex vtx = constructKalmanVertex(tracks, bfield);
//        _lKVF_valid[_nL]                = vtx.isValid();
//        if(vtx.isValid()){
//            _lKVF_x[_nL]                 = vtx.position().x();
//            _lKVF_y[_nL]                 = vtx.position().y();
//            _lKVF_z[_nL]                 = vtx.position().z();
//            _lKVF_cxx[_nL]               = vtx.positionError().cxx();
//            _lKVF_cyy[_nL]               = vtx.positionError().cyy();
//            _lKVF_czz[_nL]               = vtx.positionError().czz();
//            _lKVF_cyx[_nL]               = vtx.positionError().cyx();
//            _lKVF_czy[_nL]               = vtx.positionError().czy();
//            _lKVF_czx[_nL]               = vtx.positionError().czx();
//            _lKVF_df[_nL]                = vtx.degreesOfFreedom();
//            _lKVF_chi2[_nL]              = vtx.totalChiSquared();
//        }
//    }
//}

//This function stores all vertices with their information: position, errors, chi2, dof, track kinematics
//void LeptonAnalyzer::fillAllIVFVariables(const std::vector<reco::Vertex>& secVertices, const reco::Vertex& primaryVertex){
//    _IVF_nvertex = 0;
//    for(const reco::Vertex& vtx : secVertices){
//        //store vertex info
//        if(_IVF_nvertex == nvtx_max) break;
//        _IVF_x[_IVF_nvertex]    = vtx.x();
//        _IVF_y[_IVF_nvertex]    = vtx.y();
//        _IVF_z[_IVF_nvertex]    = vtx.z();
//        _IVF_cx[_IVF_nvertex]   = vtx.xError();
//        _IVF_cy[_IVF_nvertex]   = vtx.yError();
//        _IVF_cz[_IVF_nvertex]   = vtx.zError();
//        _IVF_df[_IVF_nvertex]   = vtx.ndof();
//        _IVF_chi2[_IVF_nvertex] = vtx.chi2();
//        _IVF_pt[_IVF_nvertex]   = vtx.p4().pt();
//        _IVF_eta[_IVF_nvertex]  = vtx.p4().eta();
//        _IVF_phi[_IVF_nvertex]  = vtx.p4().phi();
//        _IVF_E[_IVF_nvertex]    = vtx.p4().energy();
//        _IVF_mass[_IVF_nvertex] = vtx.p4().mass();
//        
//        _IVF_ntracks[_IVF_nvertex] = 0;
//        for(reco::Vertex::trackRef_iterator vtxTrackref = vtx.tracks_begin(); vtxTrackref != vtx.tracks_end(); vtxTrackref++){
//            //store track info
//            if(_IVF_ntracks[_IVF_nvertex] == ntrack_max) break;
//            reco::TrackRef vtxTrack = vtxTrackref->castTo<reco::TrackRef>();
//            _IVF_trackpt[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]        = vtxTrack->pt();
//            _IVF_tracketa[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]       = vtxTrack->eta();
//            _IVF_trackphi[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]       = vtxTrack->phi();
//            _IVF_trackE[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]         = vtxTrack->p();
//            _IVF_trackcharge[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]    = vtxTrack->charge();
//            _IVF_trackdxy[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]       = std::abs(vtxTrack->dxy(primaryVertex.position()));
//            _IVF_trackdz[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]        = std::abs(vtxTrack->dz(primaryVertex.position()));
//            _IVF_ntracks[_IVF_nvertex]++;
//        }
//        _IVF_nvertex++;
//    }
//}
//void LeptonAnalyzer::OldfillMatchingIVFVariables(const pat::Muon& muon){
//    _lIVF_match[_nLight] = -1;
//    //const reco::Track& lepTrack = (*muon.bestTrack());
//    for(unsigned i_vtx = 0; i_vtx < _IVF_nvertex; i_vtx++){
//        for(unsigned i_track = 0; i_track < _IVF_ntracks[i_vtx]; i_track++){
//            if(fabs(muon.sourceCandidatePtr(0)->pt() - _IVF_trackpt[i_vtx][i_track]) < 0.001){
//                _lIVF_match[_nLight] = i_vtx;
//            }
//        }
//    }
//}
//
//void LeptonAnalyzer::OldfillMatchingIVFVariables(const pat::Electron& electron){
//    _lIVF_match[_nLight] = -1;
//    //const reco::Track lepTrack = (*electron.gsfTrack());
//    //std::cout << std::endl << "-------------------------" << std::endl;
//    //std::cout << "Electron pt, eta, phi, dxy, dz, charge: " << electron.pt() << " " << electron.eta() << " " << electron.phi() << " " << electron.charge() << std::endl;
//    //std::cout << "Cands pt, eta, phi, dxy, dz, charge: " << std::endl;
//    //for(edm::Ref<pat::PackedCandidateCollection> cand : electron.associatedPackedPFCandidates()){
//    //    std::cout << cand->pt() << " " << cand->eta() << " " << cand->phi() << " " << cand->dxy() << " " << cand->dz() << " " << cand->charge() << std::endl;
//    //}
//    //std::cout << "-- nvtx: " << _IVF_nvertex << std::endl;
//    double dR, deta;
//    double mindR = 20;
//    for(unsigned i_vtx = 0; i_vtx < _IVF_nvertex; i_vtx++){
//        //std::cout << "vtx " << _IVF_ntracks[i_vtx] << std::endl;
//        for(unsigned i_track = 0; i_track < _IVF_ntracks[i_vtx]; i_track++){
//            for(edm::Ref<pat::PackedCandidateCollection> cand : electron.associatedPackedPFCandidates()){
//                //if(fabs(cand->pt() - _IVF_trackpt[i_vtx][i_track]) < 0.1 and fabs(cand->eta() - _IVF_tracketa[i_vtx][i_track]) < 0.1 and cand->charge() != 0){
//                dR   = reco::deltaR(cand->eta(), cand->phi(), _IVF_tracketa[i_vtx][i_track], _IVF_trackphi[i_vtx][i_track]);
//                deta = fabs(cand->eta() - _IVF_tracketa[i_vtx][i_track]);
//                if(dR < 0.05 or (dR < 0.2 and deta < 0.01)){
//                    //std::cout << _IVF_trackpt[i_vtx][i_track] << " " << _IVF_tracketa[i_vtx][i_track] << " " << _IVF_trackphi[i_vtx][i_track] << " " << _IVF_trackdxy[i_vtx][i_track] << " " << _IVF_trackdz[i_vtx][i_track] << "      " << reco::deltaR(cand->eta(), cand->phi(), _IVF_tracketa[i_vtx][i_track], _IVF_trackphi[i_vtx][i_track]) << " " << fabs(cand->eta() - _IVF_tracketa[i_vtx][i_track]) << std::endl;
//                    if(dR < mindR){
//                        mindR = dR;
//                        _lIVF_match[_nLight] = i_vtx; 
//                    }
//                }
//            }
//        }
//    }
//}
//
