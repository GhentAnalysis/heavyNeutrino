#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"     // for displaced
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h" // for displaced
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"        // for displaced
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"      // for displaced
#include "TLorentzVector.h"
#include <algorithm>


// TODO: we should maybe stop indentifying effective areas by year, as they are typically more connected to a specific ID than to a specific year
LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer),
    electronsEffectiveAreas(iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas").fullPath()),
    muonsEffectiveAreas    ((multilepAnalyzer->is2017 || multilepAnalyzer->is2018)? (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath()),
    singleEleTrigs(multilepAnalyzer->is2017 ? iConfig.getParameter<std::vector<std::string> >("SingleEleTriggers2017") : (multilepAnalyzer->is2018 ? iConfig.getParameter<std::vector<std::string> >("SingleEleTriggers2018") : iConfig.getParameter<std::vector<std::string> >("SingleEleTriggers"))),
    singleMuoTrigs(multilepAnalyzer->is2017 ? iConfig.getParameter<std::vector<std::string> >("SingleMouTriggers2017") : (multilepAnalyzer->is2018 ? iConfig.getParameter<std::vector<std::string> >("SingleMuoTriggers2018") : iConfig.getParameter<std::vector<std::string> >("SingleMuoTriggers"))),
{
    leptonMvaComputerSUSY16   = new LeptonMvaHelper(iConfig, 0, false); // SUSY         // TODO: all of these really needed? could use some clean-up
    leptonMvaComputerTTH16    = new LeptonMvaHelper(iConfig, 1, false); // TTH
    leptonMvaComputerSUSY17   = new LeptonMvaHelper(iConfig, 0, true);  // SUSY
    leptonMvaComputerTTH17    = new LeptonMvaHelper(iConfig, 1, true);  // TTH
    leptonMvaComputertZqTTV16 = new LeptonMvaHelper(iConfig, 2, false); // tZq/TTV
    leptonMvaComputertZqTTV17 = new LeptonMvaHelper(iConfig, 2, true);  // tZq/TTV
    if(!multilepAnalyzer->isData) genMatcher = new GenMatching(iConfig); // displaced specific
};

LeptonAnalyzer::~LeptonAnalyzer(){
    delete leptonMvaComputerSUSY16;
    delete leptonMvaComputerTTH16;
    delete leptonMvaComputertZqTTV16;
    delete leptonMvaComputerSUSY17;
    delete leptonMvaComputerTTH17;
    delete leptonMvaComputertZqTTV17;
    if(!multilepAnalyzer->isData) delete genMatcher; // displaced specific
}

void LeptonAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nL",                           &_nL,                           "_nL/i");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/i");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/i");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/i");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/i");
    outputTree->Branch("_pvX",                          &_pvX,                          "_pvX/D");       // displaced specific [TODO: it seems these should be moved to multilep.cc or so]
    outputTree->Branch("_pvY",                          &_pvY,                          "_pvY/D");       // "
    outputTree->Branch("_pvZ",                          &_pvZ,                          "_pvZ/D");       // "
    outputTree->Branch("_pvXErr",                       &_pvXErr,                       "_pvXErr/D");    // "
    outputTree->Branch("_pvYErr",                       &_pvYErr,                       "_pvYErr/D");    // "
    outputTree->Branch("_pvZErr",                       &_pvZErr,                       "_pvZErr/D");    // "
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/b");
    outputTree->Branch("_nVFit",                        &_nVFit,                        "_nVFit/i");                   // displaced specific
    outputTree->Branch("_nGoodLeading",                 &_nGoodLeading,                 "_nGoodLeading/i");            // "
    outputTree->Branch("_nGoodDisplaced",               &_nGoodDisplaced,               "_nGoodDisplaced/i");          // "
//  outputTree->Branch("_lIndex",                       &_lIndex,                       "_lIndex[_nL]/i");             // DELETED (simply use your index+1 instead in case you have some crazy situation where you need your index to start from 1 instead of 0)
    outputTree->Branch("_vertices",                     &_vertices,                     "_vertices[_nVFit][12]/D");    // displaced specific
    outputTree->Branch("_lDisplaced",                   &_lDisplaced,                   "_lDisplaced[_nVFit][24]/D");  // "
    outputTree->Branch("_lHasTrigger",                  &_lHasTrigger,                  "_lHasTrigger[_nL]/i");        // "
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
    outputTree->Branch("_2dIP",                         &_2dIP,                         "_2dIP[_nL]/D");          // displaced specific
    outputTree->Branch("_2dIPSig",                      &_2dIPSig,                      "_2dIPSig[_nL]/D");       // "
    //outputTree->Branch("_lElectronSummer16MvaGP",       &_lElectronMvaSummer16GP,       "_lElectronMvaSummer16GP[_nLight]/F");
    //outputTree->Branch("_lElectronSummer16MvaHZZ",      &_lElectronMvaSummer16HZZ,      "_lElectronMvaSummer16HZZ[_nLight]/F");
    //outputTree->Branch("_lElectronMvaFall17v1NoIso",    &_lElectronMvaFall17v1NoIso,    "_lElectronMvaFall17v1NoIso[_nLight]/F");
    //outputTree->Branch("_lElectronMvaFall17Iso",        &_lElectronMvaFall17Iso,        "_lElectronMvaFall17Iso[_nLight]/F");
    //outputTree->Branch("_lElectronMvaFall17NoIso",      &_lElectronMvaFall17NoIso,      "_lElectronMvaFall17NoIso[_nLight]/F");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", &_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, "_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL]/O"); // displaced specific (still needed ???)
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    //outputTree->Branch("_leptonMvaSUSY16",              &_leptonMvaSUSY16,              "_leptonMvaSUSY16[_nLight]/D");
    //outputTree->Branch("_leptonMvaTTH16",               &_leptonMvaTTH16,               "_leptonMvaTTH16[_nLight]/D");
    //outputTree->Branch("_leptonMvaSUSY17",              &_leptonMvaSUSY17,              "_leptonMvaSUSY17[_nLight]/D");
    //outputTree->Branch("_leptonMvaTTH17",               &_leptonMvaTTH17,               "_leptonMvaTTH17[_nLight]/D");
    //outputTree->Branch("_leptonMvatZqTTV16",            &_leptonMvatZqTTV16,            "_leptonMvatZqTTV16[_nLight]/D");
    //outputTree->Branch("_leptonMvatZqTTV17",            &_leptonMvatZqTTV17,            "_leptonMvatZqTTV17[_nLight]/D");
    //outputTree->Branch("_lHNLoose",                     &_lHNLoose,                     "_lHNLoose[_nLight]/O");
    //outputTree->Branch("_lHNFO",                        &_lHNFO,                        "_lHNFO[_nLight]/O");
    //outputTree->Branch("_lHNTight",                     &_lHNTight,                     "_lHNTight[_nLight]/O");
    //outputTree->Branch("_lEwkLoose",                    &_lEwkLoose,                    "_lEwkLoose[_nL]/O");
    //outputTree->Branch("_lEwkFO",                       &_lEwkFO,                       "_lEwkFO[_nL]/O");
    //outputTree->Branch("_lEwkTight",                    &_lEwkTight,                    "_lEwkTight[_nL]/O");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
//  outputTree->Branch("_lpassConversionVeto",          &_lpassConversionVeto,          "_lpassConversionVeto[_nL]/O");             // DELETED (same as the original lElectronPassConvVeto)
//  outputTree->Branch("_eleNumberInnerHitsMissing",    &_eleNumberInnerHitsMissing,    "_eleNumberInnerHitsMissing[_nL]/D");       // DELETED (same as the original lElectronMissingHits)
    outputTree->Branch("_lGlobalMuon",                  &_lGlobalMuon,                  "_lGlobalMuon[_nMu]/O");                    // displaced specific
    outputTree->Branch("_lTrackerMuon",                 &_lTrackerMuon,                 "_lTrackerMuon[_nMu]/O");                   // "
    outputTree->Branch("_lInnerTrackValidFraction",     &_lInnerTrackValidFraction,     "_lInnerTrackValidFraction[_nMu]/D");       // "
    outputTree->Branch("_lGlobalTrackNormalizeChi2",    &_lGlobalTrackNormalizeChi2,    "_lGlobalTrackNormalizeChi2[_nMu]/D");      // "
    outputTree->Branch("_lCQChi2Position",              &_lCQChi2Position,              "_lCQChi2Position[_nMu]/D");                // "
    outputTree->Branch("_lCQTrackKink",                 &_lCQTrackKink,                 "_lCQTrackKink[_nMu]/D");                   // "
    outputTree->Branch("_lNumberOfMatchedStation",      &_lNumberOfMatchedStation,      "_lNumberOfMatchedStation[_nMu]/i");        // "
    outputTree->Branch("_lNumberOfValidPixelHits",      &_lNumberOfValidPixelHits,      "_lNumberOfValidPixelHits[_nMu]/i");        // "
    outputTree->Branch("_lTrackerLayersWithMeasurement",&_lTrackerLayersWithMeasurement,"_lTrackerLayersWithMeasurement[_nMu]/i");  // "
    outputTree->Branch("_lSimType",                     &_lSimType,                     "_lSimType[_nMu]/I");                       // "
    outputTree->Branch("_lSimExtType",                  &_lSimExtType,                  "_lSimExtType[_nMu]/I");                    // "
    outputTree->Branch("_lSimFlavour",                  &_lSimFlavour,                  "_lSimFlavour[_nMu]/I");                    // "
    outputTree->Branch("_muDTStationsWithValidHits",    &_muDTStationsWithValidHits,    "_muDTStationsWithValidHits[_nMu]/I");      // "
    outputTree->Branch("_muCSCStationsWithValidHits",   &_muCSCStationsWithValidHits,   "_muCSCStationsWithValidHits[_nMu]/I");     // "
    outputTree->Branch("_muRPCStationsWithValidHits",   &_muRPCStationsWithValidHits,   "_muRPCStationsWithValidHits[_nMu]/I");     // "
    outputTree->Branch("_muMuonStationsWithValidHits",  &_muMuonStationsWithValidHits,  "_muMuonStationsWithValidHits[_nMu]/I");    // "
    outputTree->Branch("_lMuRPCTimenDof",               &_lMuRPCTimenDof,               "_lMuRPCTimenDof[_nMu]/I");                 // "
    outputTree->Branch("_lMuTimenDof",                  &_lMuTimenDof,                  "_lMuTimenDof[_nMu]/I");                    // "
    outputTree->Branch("_lMuRPCTime",                   &_lMuRPCTime,                   "_lMuRPCTime[_nMu]/D");                     // "
    outputTree->Branch("_lMuRPCTimeErr",                &_lMuRPCTimeErr,                "_lMuRPCTimeErr[_nMu]/D");                  // "
    outputTree->Branch("_lMuTime",                      &_lMuTime,                      "_lMuTime[_nMu]/D");                        // "
    outputTree->Branch("_lMuTimeErr",                   &_lMuTimeErr,                   "_lMuTimeErr[_nMu]/D");                     // "
    outputTree->Branch("_muNumberInnerHits",            &_muNumberInnerHits,            "_muNumberInnerHits[_nMu]/i");              // "
    outputTree->Branch("_lEleIsEB",                     &_lEleIsEB ,                    "_lEleIsEB[_nLight]/O");                    // "
    outputTree->Branch("_lEleIsEE",                     &_lEleIsEE ,                    "_lEleIsEE[_nLight]/O");                    // "
    outputTree->Branch("_lEleSuperClusterOverP",        &_lEleSuperClusterOverP ,       "_lEleSuperClusterOverP[_nLight]/D");       // "
    outputTree->Branch("_lEleEcalEnergy",               &_lEleEcalEnergy ,              "_lEleEcalEnergy[_nLight]/D");              // "
    outputTree->Branch("_lElefull5x5SigmaIetaIeta",     &_lElefull5x5SigmaIetaIeta ,    "_lElefull5x5SigmaIetaIeta[_nLight]/D");    // "
    outputTree->Branch("_lEleDEtaInSeed",               &_lEleDEtaInSeed ,              "_lEleDEtaInSeed[_nLight]/D");              // "
    outputTree->Branch("_lEleDeltaPhiSuperClusterTrackAtVtx", &_lEleDeltaPhiSuperClusterTrackAtVtx , "_lEleDeltaPhiSuperClusterTrackAtVtx[_nLight]/D"); // "
    outputTree->Branch("_lElehadronicOverEm",           &_lElehadronicOverEm ,          "_lElehadronicOverEm[_nLight]/D");          // "
    outputTree->Branch("_lEleInvMinusPInv",             &_lEleInvMinusPInv ,            "_lEleInvMinusPInv[_nLight]/D");            // "
    outputTree->Branch("_puCorr",                       &_puCorr,                       "_puCorr[_nLight]/D");                      // "
    outputTree->Branch("_absIso03",                     &_absIso03,                     "_absIso03[_nL]/D");                        // "
    outputTree->Branch("_absIso04",                     &_absIso04,                     "_absIso04[_nMu]/D");                       // "
    outputTree->Branch("_sumNeutralHadronEt04",         &_sumNeutralHadronEt04,         "_sumNeutralHadronEt04[_nMu]/D");           // "
    outputTree->Branch("_sumChargedHadronPt04",         &_sumChargedHadronPt04,         "_sumChargedHadronPt04[_nMu]/D");           // "
    outputTree->Branch("_sumPhotonEt04",                &_sumPhotonEt04,                "_sumPhotonEt04[_nMu]/D");                  // "
    outputTree->Branch("_sumNeutralHadronEt03",         &_sumNeutralHadronEt03,         "_sumNeutralHadronEt03[_nLight]/D");        // "
    outputTree->Branch("_sumChargedHadronPt03",         &_sumChargedHadronPt03,         "_sumChargedHadronPt03[_nLight]/D");        // "
    outputTree->Branch("_sumPhotonEt03",                &_sumPhotonEt03,                "_sumPhotonEt03[_nLight]/D");               // "
    outputTree->Branch("_trackIso",                     &_trackIso ,                    "_trackIso[_nLight]/D");                    // "
    outputTree->Branch("_ecalIso",                      &_ecalIso ,                     "_ecalIso[_nLight]/D");                     // "
    outputTree->Branch("_hcalIso",                      &_hcalIso ,                     "_hcalIso[_nLight]/D");                     // "
//  outputTree->Branch("_deltaBIso",                    &_deltaBIso,                    "_deltaBIso[_nLight]/D");                   // DELETED (what was been filled was effective areas-based iso, not deltaBeta)
    outputTree->Branch("_ecalPFClusterIso",             &_ecalPFClusterIso ,            "_ecalPFClusterIso[_nLight]/D");            // displaced specific
    outputTree->Branch("_hcalPFClusterIso",             &_hcalPFClusterIso ,            "_hcalPFClusterIso[_nLight]/D");            // "
    outputTree->Branch("_tauMuonVeto",                  &_tauMuonVeto,                  "_tauMuonVeto[_nL]/O");
    outputTree->Branch("_tauEleVeto",                   &_tauEleVeto,                   "_tauEleVeto[_nL]/O");
    outputTree->Branch("_decayModeFindingNew",          &_decayModeFindingNew,          "_decayModeFindingNew[_nL]/O");
    outputTree->Branch("_tauVLooseMvaNew",              &_tauVLooseMvaNew,              "_tauVLooseMvaNew[_nL]/O");
    outputTree->Branch("_tauLooseMvaNew",               &_tauLooseMvaNew,               "_tauLooseMvaNew[_nL]/O");
    outputTree->Branch("_tauMediumMvaNew",              &_tauMediumMvaNew,              "_tauMediumMvaNew[_nL]/O");
    outputTree->Branch("_tauTightMvaNew",               &_tauTightMvaNew,               "_tauTightMvaNew[_nL]/O");
    outputTree->Branch("_tauVTightMvaNew",              &_tauVTightMvaNew,              "_tauVTightMvaNew[_nL]/O");
    outputTree->Branch("_tauVTightMvaOld",              &_tauVTightMvaOld,              "_tauVTightMvaOld[_nL]/O");
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
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");
    if(!multilepAnalyzer->isData){
      outputTree->Branch("_lGenIndex",                  &_lGenIndex,                    "_lGenIndex[_nL]/i");
      outputTree->Branch("_lMatchType",                 &_lMatchType,                   "_lMatchType[_nL]/i");
      outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
      outputTree->Branch("_lIsPromptFinalState",        &_lIsPromptFinalState,          "_lIsPromptFinalState[_nL]/O");
      outputTree->Branch("_lIsPromptDecayed",           &_lIsPromptDecayed,             "_lIsPromptDecayed[_nL]/O");
      outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
      outputTree->Branch("_lMomPdgId",                  &_lMomPdgId,                    "_lMomPdgId[_nL]/I");
      outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
      outputTree->Branch("_lProvenanceCompressed",      &_lProvenanceCompressed,        "_lProvenanceCompressed[_nL]/i");
      outputTree->Branch("_lProvenanceConversion",      &_lProvenanceConversion,        "_lProvenanceConversion[_nL]/i");
      outputTree->Branch("_lMatchPt",                   &_lMatchPt,                     "_lMatchPt[_nL]/D");
      outputTree->Branch("_lMatchEta",                  &_lMatchEta,                    "_lMatchEta[_nL]/D");
      outputTree->Branch("_lMatchPhi",                  &_lMatchPhi,                    "_lMatchPhi[_nL]/D");
      outputTree->Branch("_lMatchVertexX",              &_lMatchVertexX,                "_lMatchVertexX[_nL]/D");
      outputTree->Branch("_lMatchVertexY",              &_lMatchVertexY,                "_lMatchVertexY[_nL]/D");
      outputTree->Branch("_lMatchVertexZ",              &_lMatchVertexZ,                "_lMatchVertexZ[_nL]/D");
    }
    if(!multilepAnalyzer->is2018){
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
    }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& primaryVertex){
    edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
    edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
    edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(multilepAnalyzer->tauToken,                          taus);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
    edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);
    edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);
  //edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetSmearedToken,                   jets);  // Are we sure we do not want the smeared jets here???
    edm::Handle<std::vector<reco::GenParticle>> genParticles;        iEvent.getByToken(multilepAnalyzer->genParticleToken,                  genParticles);
    edm::Handle<edm::TriggerResults> trigBits;                       iEvent.getByToken(multilepAnalyzer->triggerToken,                      trigBits);
    edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjs;    iEvent.getByToken(multilepAnalyzer->trigObjToken,                      trigObjs);

    iSetup.get<IdealMagneticFieldRecord>().get(_bField);
    iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", _shProp);
    GsfPropagatorAdapter gsfPropagator(AnalyticalPropagator(&(*_bField), anyDirection));
    _gsfProp = new TransverseImpactPointExtrapolator(gsfPropagator);

    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;
    _nVFit  = 0;
    _nGoodLeading = 0;
    _nGoodDisplaced = 0;

    // Reco primary vertex coordinates and uncertainties // TODO: seems to not belong in LeptonAnalyzer
    _pvX    = primaryVertex.x();
    _pvXErr = primaryVertex.xError();
    _pvY    = primaryVertex.y();
    _pvYErr = primaryVertex.yError();
    _pvZ    = primaryVertex.z();
    _pvZErr = primaryVertex.zError();

    // Check which muons pass minimum criteria
    std::vector<const pat::Muon*> selmuons;
    for(const pat::Muon& mu : *muons){
      if(passMuonPreselection(mu, *rho)) {
        selmuons.push_back(&mu);
      }
    }

    // Check which electrons pass minimum criteria
    std::vector<const pat::Electron*> selelectrons;
    for(const pat::Electron& ele : *electrons){
      if(passElectronPreselection(ele, *rho)) {
        selelectrons.push_back(&ele);
      }
    }

    // Check which taus pass minimum criteria
    std::vector<const pat::Tau*> seltaus;
    for(const pat::Tau& tau : *taus){
      if(passTauPreselection(tau, primaryVertex.position())) {
        seltaus.push_back(&tau);
      }
    }

    // Run the RECO-GEN matching based on the genParticles and the selected pat::Lepton collections
    if(!multilepAnalyzer->isData){
      genMatcher->matchGenToReco(*genParticles, selelectrons, selmuons, seltaus);
    }


    // loop over muons
    // muons need to be run first, because some ID's need to calculate a muon veto for electrons
    for(const pat::Muon* muptr : selmuons) {
        if(_nL == nL_max) break;
        const pat::Muon& mu = (*muptr);

        _lPFMuon[_nL] = mu.isPFMuon();

        const reco::MuonTime cmb = mu.time();
        const reco::MuonTime rpc = mu.rpcTime();

        //csc + dt
        _lMuTimenDof[_nL] = cmb.nDof;
        _lMuTime[_nL]     = cmb.timeAtIpInOut;
        _lMuTimeErr[_nL]  = cmb.timeAtIpInOutErr;

        //RPC
        _lMuRPCTimenDof[_nL] = rpc.nDof;
        _lMuRPCTime[_nL]     = rpc.timeAtIpInOut;
        _lMuRPCTimeErr[_nL]  = rpc.timeAtIpInOutErr;

        fillLeptonImpactParameters(mu, primaryVertex);

        fillLeptonKinVars(mu);
        fillLeptonIsoVars(mu, *rho);
        if(!multilepAnalyzer->isData) fillLeptonGenVars(mu, *genParticles);

        fillLeptonJetVariables(mu, jets, primaryVertex, *rho);
        _lGlobalMuon[_nL]                   = mu.isGlobalMuon();
        _lTrackerMuon[_nL]                  = mu.isTrackerMuon();
        _lInnerTrackValidFraction[_nL]      = (!mu.innerTrack().isNull()) ? mu.innerTrack()->validFraction()  : -1;
        _lGlobalTrackNormalizeChi2[_nL]     = (!mu.globalTrack().isNull()) ? mu.globalTrack()->normalizedChi2()  : -1;
        _lCQChi2Position[_nL]               = mu.combinedQuality().chi2LocalPosition;
        _lCQTrackKink[_nL]                  = mu.combinedQuality().trkKink;
        _lNumberOfMatchedStation[_nL]       = mu.numberOfMatchedStations();
        _lNumberOfValidPixelHits[_nL]       = (!mu.innerTrack().isNull()) ? mu.innerTrack()->hitPattern().numberOfValidPixelHits() : 0; // cannot be -1 !!
        _lTrackerLayersWithMeasurement[_nL] = (!mu.innerTrack().isNull()) ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0; // cannot be -1 !!

        _lFlavor[_nL]        = 1;
        _lMuonSegComp[_nL]    = mu.segmentCompatibility();
        _lMuonTrackPt[_nL]    = mu.innerTrack()->pt();
        _lMuonTrackPtErr[_nL] = mu.innerTrack()->ptError();
        _lSimType[_nL]        = mu.simType();
        _lSimExtType[_nL]     = mu.simExtType();
        _lSimFlavour[_nL]     = mu.simFlavour();

        // TODO: this is a possible solution to the missing trackRef, but maybe not what you want
        _muNumberInnerHits[_nL]= (!mu.globalTrack().isNull()) ?   mu.globalTrack()->hitPattern().numberOfValidMuonHits() : (!mu.outerTrack().isNull() ? mu.outerTrack()->hitPattern().numberOfValidMuonHits() : 0); // cannot be -1 !!!

        _relIso[_nL]         = getRelIso03(mu, *rho);                     // Isolation variables
        _relIso0p4[_nL]      = getRelIso04(mu, *rho);
        _relIso0p4MuDeltaBeta[_nL] = getRelIso04(mu, *rho, true);
        _miniIso[_nL]        = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho, false); // TODO: check how this compares with the MiniIsoLoose,etc... booleans
        _miniIsoCharged[_nL] = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho, true);

        _lHNLoose[_nL]       = isHNLoose(mu);                                                       // ID variables
        _lHNFO[_nL]          = isHNFO(mu);                                                          // don't change order, they rely on above variables
        _lHNTight[_nL]       = isHNTight(mu);

        _lPOGVeto[_nL]       = mu.passed(reco::Muon::CutBasedIdLoose); // no veto available, so we take loose here
        _lPOGLoose[_nL]      = mu.passed(reco::Muon::CutBasedIdLoose);
        _lPOGMedium[_nL]     = mu.passed(reco::Muon::CutBasedIdMedium);
        _lPOGTight[_nL]      = mu.passed(reco::Muon::CutBasedIdTight);
        // TODO: consider to add muon MVA

        _leptonMvaSUSY16[_nL]  = leptonMvaVal(mu, leptonMvaComputerSUSY16);
        _leptonMvaTTH16[_nL]   = leptonMvaVal(mu, leptonMvaComputerTTH16);
        _leptonMvaSUSY17[_nL]  = leptonMvaVal(mu, leptonMvaComputerSUSY17);
        _leptonMvaTTH17[_nL]   = leptonMvaVal(mu, leptonMvaComputerTTH17);
        _leptonMvatZqTTV16[_nL] = leptonMvaVal(mu, leptonMvaComputertZqTTV16);
        _leptonMvatZqTTV17[_nL] = leptonMvaVal(mu, leptonMvaComputertZqTTV17);

        _lEwkLoose[_nL]      = isEwkLoose(mu);
        _lEwkFO[_nL]         = isEwkFO(mu);
        _lEwkTight[_nL]      = isEwkTight(mu);

        _muDTStationsWithValidHits[_nL]   = mu.bestTrack()->hitPattern().dtStationsWithValidHits();
        _muCSCStationsWithValidHits[_nL]  = mu.bestTrack()->hitPattern().cscStationsWithValidHits();
        _muRPCStationsWithValidHits[_nL]  = mu.bestTrack()->hitPattern().rpcStationsWithValidHits();
        _muMuonStationsWithValidHits[_nL] = mu.bestTrack()->hitPattern().muonStationsWithValidHits();

        // maybe put this simply in a function in LeptonAnalyzerId.cc???
        bool someIdWithoutName = (mu.pt() > 24 and _lPOGMedium[_nL]
                                  and fabs(_dxy[_nL]) < 0.05 and fabs(_dz[_nL])< 0.1
                                  and getRelIso03(mu, *rho) < 0.2
                                  and !mu.innerTrack().isNull()
                                  and (mu.isTrackerMuon() or mu.isGlobalMuon()));

        if(someIdWithoutName) ++_nGoodLeading;
        if((multilepAnalyzer->skim == "trilep" or multilepAnalyzer->skim == "displtrilep") && someIdWithoutName) _lHasTrigger[_nL] = matchSingleTrigger(false, _lEta[_nL], _lPhi[_nL], iEvent.triggerNames(*trigBits), trigObjs);
        else                                                    _lHasTrigger[_nL] = 0;

        ++_nMu;
        ++_nL;
        ++_nLight;
    }


    // Loop over electrons (note: using iterator we can easily get the ref too)
    for(size_t iele=0; iele<selelectrons.size(); ++iele){
        if(_nL == nL_max) break;
        const pat::Electron *ele = selelectrons[iele];

        fillLeptonKinVars(*ele);
        if(!multilepAnalyzer->isData) fillLeptonGenVars(*ele, *genParticles);
        fillLeptonJetVariables(*ele, jets, primaryVertex, *rho);

        fillLeptonImpactParameters(*ele, primaryVertex);
        fillLeptonIsoVars(*ele, *rho);


        _lFlavor[_nL]                   = 0;
        _lEtaSC[_nL]                    = ele->superCluster()->eta();

        _relIso[_nL]                    = getRelIso03(*ele, *rho);
        _relIso0p4[_nL]                 = getRelIso(*ele, packedCands, 0.4, *rho, false);
        _miniIso[_nL]                   = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho, false);
        _miniIsoCharged[_nL]            = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho, true);
        _lElectronMvaSummer16GP[_nL]    = ele->userFloat("ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"); // OLD, do not use it
        _lElectronMvaSummer16HZZ[_nL]   = ele->userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values"); // OLD, do not use it
        _lElectronMvaFall17v1NoIso[_nL] = ele->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values"); // OLD, do not use it
        _lElectronMvaFall17Iso[_nL]     = ele->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
        _lElectronMvaFall17NoIso[_nL]   = ele->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values");
        _lElectronPassEmu[_nL]          = passTriggerEmulationDoubleEG(&*ele);                             // Keep in mind, this trigger emulation is for 2016 DoubleEG, the SingleEG trigger emulation is different
        _lElectronPassConvVeto[_nL]     = ele->passConversionVeto();
        _lElectronChargeConst[_nL]      = ele->isGsfCtfScPixChargeConsistent();
        _lElectronMissingHits[_nL]      = ele->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

        // ID ele variables
        _lEleIsEB [_nL]                 = ele->isEB();
        _lEleIsEE[_nL]                  = ele->isEE();
        _lEleSuperClusterOverP[_nL]     = ele->eSuperClusterOverP();
        _lEleEcalEnergy[_nL]            = ele->ecalEnergy();
        _lElefull5x5SigmaIetaIeta[_nL]  = ele->full5x5_sigmaIetaIeta();
        _lEleDEtaInSeed[_nL]            = std::abs(dEtaInSeed(&*ele));
        _lEleDeltaPhiSuperClusterTrackAtVtx[_nL] = std::abs(ele->deltaPhiSuperClusterTrackAtVtx());
        _lElehadronicOverEm[_nL]        = ele->hadronicOverEm();
        _lEleInvMinusPInv[_nL]          = std::abs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();

        _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL] = isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele);

        _lHNLoose[_nL]                  = isHNLoose(*ele);
        _lHNFO[_nL]                     = isHNFO(*ele);
        _lHNTight[_nL]                  = isHNTight(*ele);

        _lPOGVeto[_nL]                  = ele->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
        _lPOGLoose[_nL]                 = ele->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
        _lPOGMedium[_nL]                = ele->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
        _lPOGTight[_nL]                 = ele->electronID("cutBasedElectronID-Fall17-94X-V1-tight");

        _leptonMvaSUSY16[_nL]           = leptonMvaVal(*ele, leptonMvaComputerSUSY16);
        _leptonMvaTTH16[_nL]            = leptonMvaVal(*ele, leptonMvaComputerTTH16);
        _leptonMvaSUSY17[_nL]           = leptonMvaVal(*ele, leptonMvaComputerSUSY17);
        _leptonMvaTTH17[_nL]            = leptonMvaVal(*ele, leptonMvaComputerTTH17);
        _leptonMvatZqTTV16[_nL]         = leptonMvaVal(*ele, leptonMvaComputertZqTTV16);
        _leptonMvatZqTTV17[_nL]         = leptonMvaVal(*ele, leptonMvaComputertZqTTV17);

        _lEwkLoose[_nL]                 = isEwkLoose(*ele);
        _lEwkFO[_nL]                    = isEwkFO(*ele);
        _lEwkTight[_nL]                 = isEwkTight(*ele);

        // Note: for the scale and smearing systematics we use the overall values, assuming we are not very sensitive to these systematics
        // In case these systematics turn out to be important, need to add their individual source to the tree (and propagate to their own templates):
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
        // Currently only available for 2016/2017
        if(!multilepAnalyzer->is2018){
          _lPtCorr[_nL]                 = ele->pt()*ele->userFloat("ecalTrkEnergyPostCorr")/ele->energy();
          _lPtScaleUp[_nL]              = ele->pt()*ele->userFloat("energyScaleUp")/ele->energy();
          _lPtScaleDown[_nL]            = ele->pt()*ele->userFloat("energyScaleDown")/ele->energy();
          _lPtResUp[_nL]                = ele->pt()*ele->userFloat("energySigmaUp")/ele->energy();
          _lPtResDown[_nL]              = ele->pt()*ele->userFloat("energySigmaDown")/ele->energy();
          _lECorr[_nL]                  = ele->userFloat("ecalTrkEnergyPostCorr");
          _lEScaleUp[_nL]               = ele->userFloat("energyScaleUp");
          _lEScaleDown[_nL]             = ele->userFloat("energyScaleDown");
          _lEResUp[_nL]                 = ele->userFloat("energySigmaUp");
          _lEResDown[_nL]               = ele->userFloat("energySigmaDown");
        }


        // maybe put this simply in a function in LeptonAnalyzerId.cc???
        bool someIdWithoutName = (ele->pt() > 27 and _lPOGMedium[_nL]
                                  and fabs(_dxy[_nL]) < 0.05 and fabs(_dz[_nL])< 0.1
                                  and _relIso[_nL] < 0.2
                                  and !ele->gsfTrack().isNull()
                                  and _lElectronMissingHits[_nL] <=2
                                  and ele->passConversionVeto());

        if(ele->pt() > 7 && fabs(_dxy[_nL]) > 0.02) ++_nGoodDisplaced;
        if(someIdWithoutName) ++_nGoodLeading;
        if((multilepAnalyzer->skim == "trilep" or multilepAnalyzer->skim == "displtrilep")  && someIdWithoutName) _lHasTrigger[_nL] = matchSingleTrigger(true, _lEta[_nL], _lPhi[_nL], iEvent.triggerNames(*trigBits), trigObjs);
        else                                                    _lHasTrigger[_nL] = 0;

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
    for(auto array : {&_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, &_lEleIsEB, &_lEleIsEE}) std::fill_n(*array, _nMu, false);       // displaced specific
    for(auto array : {&_lEleSuperClusterOverP, &_lEleEcalEnergy, &_lElefull5x5SigmaIetaIeta, &_lEleDEtaInSeed}) std::fill_n(*array, _nMu, -1.); // displaced speficic
    for(auto array : {&_lEleDeltaPhiSuperClusterTrackAtVtx, &_lElehadronicOverEm, &_lEleInvMinusPInv}) std::fill_n(*array, _nMu, -1.);          // displaced specific
    for(auto array : {&_ecalPFClusterIso, &_hcalPFClusterIso}) std::fill_n(*array, _nMu, -1.);                                                  // displaced specific

    //loop over taus
    for(const pat::Tau* tauptr : seltaus) {
        if(_nL == nL_max) break;
        const pat::Tau& tau = (*tauptr);

        fillLeptonKinVars(tau);
        if(!multilepAnalyzer->isData) fillLeptonGenVars(tau, *genParticles);
        fillLeptonImpactParameters(tau, primaryVertex);

        _lFlavor[_nL]                   = 2;
        _tauMuonVeto[_nL]               = tau.tauID("againstMuonLoose3");                        //Light lepton vetos
        _tauEleVeto[_nL]                = tau.tauID("againstElectronLooseMVA6");

        _lPOGVeto[_nL]                  = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");     //old tau ID
        _lPOGLoose[_nL]                 = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
        _lPOGMedium[_nL]                = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
        _lPOGTight[_nL]                 = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
        _tauVTightMvaOld[_nL]           = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

        _decayModeFindingNew[_nL]       = tau.tauID("decayModeFindingNewDMs");                   //new Tau ID
        _tauVLooseMvaNew[_nL]           = tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
        _tauLooseMvaNew[_nL]            = tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
        _tauMediumMvaNew[_nL]           = tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
        _tauTightMvaNew[_nL]            = tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT");
        _tauVTightMvaNew[_nL]           = tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT");

        _tauAgainstElectronMVA6Raw[_nL] = tau.tauID("againstElectronMVA6Raw");
        _tauCombinedIsoDBRaw3Hits[_nL]  = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        _tauIsoMVAPWdR03oldDMwLT[_nL]   = tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw");
        _tauIsoMVADBdR03oldDMwLT[_nL]   = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        _tauIsoMVADBdR03newDMwLT[_nL]   = tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        _tauIsoMVAPWnewDMwLT[_nL]       = tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        _tauIsoMVAPWoldDMwLT[_nL]       = tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");
        // TODO:  Should try also deepTau?

        _lEwkLoose[_nL] = isEwkLoose(tau);
        _lEwkFO[_nL]    = isEwkFO(tau);
        _lEwkTight[_nL] = isEwkTight(tau);
        ++_nTau;
        ++_nL;
    }


    /*
     * refitting vertices displaced
     * try all possible combinations of electrons and muons (further selections on opposite-charge is required in the fillDileptonVertexArray function)
     */
    cleanDileptonVertexArrays(_nVFit);
    for(unsigned i=0; i < selmuons.size(); ++i){
      for(unsigned j=0; j < selmuons.size(); ++j)     fillDileptonVertexArrays(i, j,      selmuons.at(i), selmuons.at(j));
      for(unsigned j=0; j < selelectrons.size(); ++j) fillDileptonVertexArrays(i, _nMu+j, selmuons.at(i), selelectrons.at(j));
    }
    for(unsigned i=0; i < selelectrons.size(); ++i){
      for(unsigned j=0; j < selmuons.size(); ++j)     fillDileptonVertexArrays(_nMu+i, j,      selelectrons.at(i), selmuons.at(j));
      for(unsigned j=0; j < selelectrons.size(); ++j) fillDileptonVertexArrays(_nMu+i, _nMu+j, selelectrons.at(i), selelectrons.at(j));
    }

    //Initialize with default values for those tau-only arrays which weren't filled with electrons and muons [to allow correct comparison by the test script]
    for(auto array : {&_tauMuonVeto, &_tauEleVeto, &_decayModeFindingNew, &_tauVLooseMvaNew, &_tauLooseMvaNew}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauMediumMvaNew, &_tauTightMvaNew, &_tauVTightMvaNew, &_tauVTightMvaOld}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauAgainstElectronMVA6Raw, &_tauCombinedIsoDBRaw3Hits, &_tauIsoMVAPWdR03oldDMwLT}) std::fill_n(*array, _nLight, 0.);
    for(auto array : {&_tauIsoMVADBdR03oldDMwLT, &_tauIsoMVADBdR03newDMwLT, &_tauIsoMVAPWnewDMwLT, &_tauIsoMVAPWoldDMwLT}) std::fill_n(*array, _nLight, 0.);


    if(multilepAnalyzer->skim == "trilep"      and (_nLight < 3 || _nGoodLeading < 1                       ) ) return false;
    if(multilepAnalyzer->skim == "displtrilep" and (_nLight < 3 || _nGoodLeading < 1 || _nGoodDisplaced < 2) ) return false;

    if(multilepAnalyzer->skim == "dilep"       and _nLight <  2) return false;
    if(multilepAnalyzer->skim == "ttg"         and _nLight <  2) return false;
    if(multilepAnalyzer->skim == "singlelep"   and _nLight <  1) return false;
    if(multilepAnalyzer->skim == "FR"          and _nLight != 1 ) return false; // TODO: this _lHasTrigger[0] seems to be an unsigned but the name sounds as if it is a boolean, in which case "not _lHasTrigger" would make more sense
    return true;
}



/*
 * //--// Refit dilepton vertex:
 * Provide a transientvertex
 */
TransientVertex LeptonAnalyzer::dileptonVertex(const reco::RecoCandidate* lep1, const reco::RecoCandidate* lep2){
  MagneticField *bfield = new OAEParametrizedMagneticField("3_8T");
  std::vector<reco::TransientTrack> ttks;
  ttks.push_back(reco::TransientTrack(getTrack(lep1), bfield));
  ttks.push_back(reco::TransientTrack(getTrack(lep2), bfield));
  KalmanVertexFitter vtxFitter;
  return vtxFitter.vertex(ttks);
}


void LeptonAnalyzer::cleanDileptonVertexArrays(unsigned nVFit){
  for (int i =0; i < 50; i++){
    for (int j =0; j < 24 ; j++){
      if (j < 12) _vertices[i][j] = 0;
      _lDisplaced[i][j] = 0;
    }
  }
}

const reco::Track& LeptonAnalyzer::getTrack(const reco::RecoCandidate* lep){
  if(lep->isMuon()){
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>(lep);
    if(!mu->innerTrack().isNull()) return *(mu->innerTrack());
    else                           return *(mu->outerTrack());
  } else {
    return *(dynamic_cast<const pat::Electron*>(lep)->gsfTrack());
  }
}

// Fill the arrays of displaced vertices and leptons
void LeptonAnalyzer::fillDileptonVertexArrays(unsigned iL_plus, unsigned iL_minus, const reco::RecoCandidate* lep1, const reco::RecoCandidate* lep2){
  if(lep1->charge() < 0) return; // ensure opposite charge
  if(lep2->charge() > 0) return;

  TransientVertex dvtx = dileptonVertex(lep1, lep2);
  if(!dvtx.isValid()){
    std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl;
    return;
  }

  _vertices[_nVFit][0]  = (iL_plus+1)*100 + (iL_minus+1); // indices i and j stored as i*100+j, but start counting at 1 instead of 0 (historical, backwards compatibility)
  _vertices[_nVFit][1]  = dvtx.position().x();
  _vertices[_nVFit][2]  = dvtx.position().y();
  _vertices[_nVFit][3]  = dvtx.position().z();
  _vertices[_nVFit][4]  = dvtx.positionError().cxx();
  _vertices[_nVFit][5]  = dvtx.positionError().cyy();
  _vertices[_nVFit][6]  = dvtx.positionError().czz();
  _vertices[_nVFit][7]  = dvtx.positionError().cyx();
  _vertices[_nVFit][8]  = dvtx.positionError().czy();
  _vertices[_nVFit][9]  = dvtx.positionError().czx();
  _vertices[_nVFit][10] = dvtx.degreesOfFreedom();
  _vertices[_nVFit][11] = dvtx.totalChiSquared();

  int i=0;
  for(auto lep : {lep1, lep2}){
    const reco::Track& track = getTrack(lep);
    GlobalPoint  globPoint(track.vx(), track.vy(), track.vz());
    GlobalVector globVector(track.px(), track.py(), track.pz());
    GlobalTrajectoryParameters globTrajParam(globPoint, globVector, track.charge(), _bField.product());
    FreeTrajectoryState tempFreeTrajectoryState(globTrajParam, track.covariance());
    FreeTrajectoryState freeTrajectoryState;

    if(lep->isMuon()) freeTrajectoryState = _shProp->propagate(tempFreeTrajectoryState, dvtx.position());
    else              freeTrajectoryState = *(_gsfProp->extrapolate(tempFreeTrajectoryState, dvtx.position()).freeState());

    if(!freeTrajectoryState.hasCurvilinearError()){
      std::cout << "Propagation of lepton to dilepton vertex (" << _vertices[_nVFit][0] << ") failed!" << std::endl;
      freeTrajectoryState = tempFreeTrajectoryState;
    }

    // fill _lDisplaced[nVFit] with 2x12=24 elements structured as 3 position + 3 momentum + 3 position error + 3 momentum error
    _lDisplaced[_nVFit][i++] = freeTrajectoryState.position().x();
    _lDisplaced[_nVFit][i++] = freeTrajectoryState.position().y();
    _lDisplaced[_nVFit][i++] = freeTrajectoryState.position().z();

    _lDisplaced[_nVFit][i++] = freeTrajectoryState.momentum().x();
    _lDisplaced[_nVFit][i++] = freeTrajectoryState.momentum().y();
    _lDisplaced[_nVFit][i++] = freeTrajectoryState.momentum().z();

    _lDisplaced[_nVFit][i++] = (freeTrajectoryState.cartesianError().matrix())(0, 0); //position error
    _lDisplaced[_nVFit][i++] = (freeTrajectoryState.cartesianError().matrix())(1, 1);
    _lDisplaced[_nVFit][i++] = (freeTrajectoryState.cartesianError().matrix())(2, 2);

    _lDisplaced[_nVFit][i++] = (freeTrajectoryState.cartesianError().matrix())(3, 3); // momentum error
    _lDisplaced[_nVFit][i++] = (freeTrajectoryState.cartesianError().matrix())(4, 4);
    _lDisplaced[_nVFit][i++] = (freeTrajectoryState.cartesianError().matrix())(5, 5);
  }

  _nVFit++;
}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton){
  _lPt[_nL]     = lepton.pt();
  _lEta[_nL]    = lepton.eta();
  _lPhi[_nL]    = lepton.phi();
  _lE[_nL]      = lepton.energy();
  _lCharge[_nL] = lepton.charge();
}

void LeptonAnalyzer::fillLeptonIsoVars(const pat::Muon& mu, const double rho){  // TODO: is all this stuff really still needed?
  _puCorr[_nL]               = rho*muonsEffectiveAreas.getEffectiveArea(mu.eta());
  _absIso03[_nL]             = getRelIso03(mu, rho)*mu.pt();
  _absIso04[_nL]             = getRelIso04(mu, rho, false)*mu.pt();
  _sumNeutralHadronEt04[_nL] = mu.pfIsolationR04().sumNeutralHadronEt;
  _sumChargedHadronPt04[_nL] = mu.pfIsolationR04().sumChargedHadronPt;
  _sumPhotonEt04[_nL]        = mu.pfIsolationR04().sumPhotonEt;
  _sumNeutralHadronEt03[_nL] = mu.pfIsolationR03().sumNeutralHadronEt;
  _sumChargedHadronPt03[_nL] = mu.pfIsolationR03().sumChargedHadronPt;
  _sumPhotonEt03[_nL]        = mu.pfIsolationR03().sumPhotonEt;
  _trackIso[_nL]             = mu.trackIso();
  _ecalIso[_nL]              = mu.ecalIso();
  _hcalIso[_nL]              = mu.hcalIso();
}

void LeptonAnalyzer::fillLeptonIsoVars(const pat::Electron& ele, const double rho){
  _puCorr[_nL]               = rho*electronsEffectiveAreas.getEffectiveArea(ele.superCluster()->eta());
  _absIso03[_nL]             = getRelIso03(ele, rho)*ele.pt();
  _sumNeutralHadronEt03[_nL] = ele.pfIsolationVariables().sumNeutralHadronEt;
  _sumChargedHadronPt03[_nL] = ele.pfIsolationVariables().sumChargedHadronPt;
  _sumPhotonEt03[_nL]        = ele.pfIsolationVariables().sumPhotonEt;
  _trackIso[_nL]             = ele.trackIso();
  _ecalIso[_nL]              = ele.ecalIso();
  _hcalIso[_nL]              = ele.hcalIso();
  _ecalPFClusterIso[_nL]     = ele.ecalPFClusterIso();
  _hcalPFClusterIso[_nL]     = ele.hcalPFClusterIso();
}

// Fill match variables
template <typename Lepton> void LeptonAnalyzer::fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles){
    auto match = genMatcher->returnGenMatch(lepton, _lMatchType[_nL]);

    _lGenIndex[_nL]             = multilepAnalyzer->genAnalyzer->getGenLeptonIndex(match);
    _lIsPrompt[_nL]             = match and (abs(lepton.pdgId()) == abs(match->pdgId()) || match->pdgId() == 22) and GenTools::isPrompt(*match, genParticles); // only when matched to its own flavor or a photon
    _lIsPromptFinalState[_nL]   = match ? match->isPromptFinalState(): false; 
    _lIsPromptDecayed[_nL]      = match ? match->isPromptDecayed() : false; 

    _lMatchPdgId[_nL]           = match ? match->pdgId() : 0; 
    _lProvenance[_nL]           = GenTools::provenance(match, genParticles); 
    _lProvenanceCompressed[_nL] = GenTools::provenanceCompressed(match, genParticles, _lIsPrompt[_nL]); 
    _lProvenanceConversion[_nL] = GenTools::provenanceConversion(match, genParticles);
    _lMomPdgId[_nL]             = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;

    _lMatchPt[_nL]              = match ? match->pt() : 0; 
    _lMatchEta[_nL]             = match ? match->eta() : 0;
    _lMatchPhi[_nL]             = match ? match->phi() : 0;
    _lMatchVertexX[_nL]         = match ? match->vertex().x() : 0;
    _lMatchVertexY[_nL]         = match ? match->vertex().y() : 0;
    _lMatchVertexZ[_nL]         = match ? match->vertex().z() : 0;
}

/*
 * Impact parameters:
 * Provide PV to dxy/dz otherwise you get dxy/dz to the beamspot instead of the primary vertex
 * For taus: dxy is pre-computed with PV it was constructed with
 */
void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Electron& ele, const reco::Vertex& vertex){
    _dxy[_nL]     = ele.gsfTrack()->dxy(vertex.position());
    _dz[_nL]      = ele.gsfTrack()->dz(vertex.position());
    _3dIP[_nL]    = ele.dB(pat::Electron::PV3D);
    _3dIPSig[_nL] = fabs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D));
    _2dIP[_nL]    = ele.dB();
    _2dIPSig[_nL] = ele.dB()/ele.edB();
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon, const reco::Vertex& vertex){
    _dxy[_nL]     = (!muon.innerTrack().isNull()) ? muon.innerTrack()->dxy(vertex.position()) : muon.outerTrack()->dxy(vertex.position());
    _dz[_nL]      = (!muon.innerTrack().isNull()) ? muon.innerTrack()->dz(vertex.position()) : muon.outerTrack()->dz(vertex.position());
    _3dIP[_nL]    = muon.dB(pat::Muon::PV3D);
    _3dIPSig[_nL] = fabs(muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D));
    _2dIP[_nL]    = muon.dB();
    _2dIPSig[_nL] = muon.dB()/muon.edB();
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

/*
 * The lepton selections [for tuple and genmatching]
 */
bool LeptonAnalyzer::passElectronPreselection(const pat::Electron& elec, const double rho) const {
  if(elec.gsfTrack().isNull())                                                                    return false;
  if(getRelIso03(elec, rho) >  1.5)                                                               return false;
  if(elec.pt()<5.)                                                                                return false;
  if(std::abs(elec.eta())>2.5)                                                                    return false;
  if(!isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&elec)) return false; // Note: should be reviewd, especially for 2017-2018
  if(eleMuOverlap(elec, _lPFMuon))                                                                return false; // overlap muon-electron deltaR<0.05, using PF muons
  return true;
}

bool LeptonAnalyzer::passMuonPreselection(const pat::Muon& muon, const double rho) const {
  if(!muon.isPFMuon())           return false;
  if(!muon.isLooseMuon())        return false;
  if(getRelIso03(muon, rho) > 1.5) return false;
  if(muon.pt() < 3.)             return false;
  if(std::abs(muon.eta()) > 2.4) return false;
  return true;
}

bool LeptonAnalyzer::passTauPreselection(const pat::Tau& tau, const reco::Vertex::Point& vertex) const {
  if(tau.pt() < 20.)            return false; // Minimum pt for tau reconstruction
  if(std::abs(tau.eta()) > 2.3) return false;
  if(tau_dz(tau, vertex) > 2.)  return false;
  return true;
}


void LeptonAnalyzer::fillLeptonJetVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::Jet>>& jets, const reco::Vertex& vertex, const double rho){
    //Make skimmed "close jet" collection
    std::vector<pat::Jet> selectedJetsAll;
    for(auto jet = jets->cbegin(); jet != jets->cend(); ++jet){
        if(jet->pt() > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
    }
    // Find closest selected jet
    unsigned closestIndex = 0;
    for(unsigned j = 1; j < selectedJetsAll.size(); ++j){
        if(reco::deltaR(selectedJetsAll[j], lepton) < reco::deltaR(selectedJetsAll[closestIndex], lepton)) closestIndex = j;
    }

    const pat::Jet& jet = selectedJetsAll[closestIndex];
    if(selectedJetsAll.size() == 0 || reco::deltaR(jet, lepton) > 0.4){ //Now includes safeguard for 0 jet events
        _ptRatio[_nL]              = 1;
        _ptRel[_nL]                = 0;
        _closestJetCsvV2[_nL]      = 0;
        _closestJetDeepCsv_b[_nL]  = 0;
        _closestJetDeepCsv_bb[_nL] = 0;
        _selectedTrackMult[_nL]    = 0;
        _closestJEC[_nL]           = 0;
        _closest_lepAwareJet[_nL]  =0;
        _closest_l1Jet[_nL]        = 0;
        _closest_lJetE [_nL]        = 0;
        _closest_lJetPx [_nL]        = 0;
        _closest_lJetPy [_nL]        = 0;
        _closest_lJetPz [_nL]        = 0;
    } else {
        auto  l1Jet       = jet.correctedP4("L1FastJet");
        float JEC         = jet.p4().E()/l1Jet.E();
        auto  l           = lepton.p4();
        auto  lepAwareJet = (l1Jet - l)*JEC + l;

        TLorentzVector lV(l.Px(), l.Py(), l.Pz(), l.E());
        TLorentzVector jV(lepAwareJet.Px(), lepAwareJet.Py(), lepAwareJet.Pz(), lepAwareJet.E());
        _ptRatio[_nL]              = l.Pt()/lepAwareJet.Pt();
        _ptRel[_nL]                = lV.Perp((jV - lV).Vect());
        _closestJetCsvV2[_nL]      = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        _closestJetDeepCsv_b[_nL]  = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _closestJEC[_nL]           = JEC;
        _closest_lepAwareJet[_nL]  =lepAwareJet;
        _closest_l1Jet[_nL]        = l1Jet;
        _closest_lJetE [_nL]        = jet.p4().E();
        _closest_lJetPx [_nL]        = jet.p4().Px();
        _closest_lJetPy [_nL]        = jet.p4().Py();
        _closest_lJetPz [_nL]        = jet.p4().Pz();


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

unsigned LeptonAnalyzer::matchSingleTrigger(bool isele, double aeta, double aphi, const edm::TriggerNames &names, edm::Handle<pat::TriggerObjectStandAloneCollection> objs){
  unsigned trigmask(0);
  for(pat::TriggerObjectStandAlone iobj : *objs) { // NOTE: not const nor by reference, because we need to 'unpackPathNames'
    if(reco::deltaR(iobj.eta(), iobj.phi(), aeta, aphi)<0.15) {
      iobj.unpackPathNames(names);
      std::vector<std::string> &singletrigs = isele ? singleEleTrigs : singleMuoTrigs;
      int ipath(-1);
      for(std::string& itrig : singletrigs) {
        ++ipath;
        if(iobj.hasPathName(itrig.c_str(), true, true)){
          trigmask |= (1<<ipath);
          //return true; // do not return, because we need the list of all paths that fired!
        }
      }
    }
  }

  return trigmask;
}
