#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"     // for displaced
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h" // for displaced
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"        // for displaced
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"      // for displaced
#include "heavyNeutrino/multilep/interface/TauTools.h"
#include "TLorentzVector.h"
#include <algorithm>


LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer),
    electronsEffectiveAreas(            iConfig.getParameter<edm::FileInPath>("electronsEffAreas").fullPath()),
    electronsEffectiveAreas_ttH_relIso( iConfig.getParameter<edm::FileInPath>("electronsEffAreas_ttH_relIso").fullPath()),
    electronsEffectiveAreas_ttH_miniIso(iConfig.getParameter<edm::FileInPath>("electronsEffAreas_ttH_miniIso").fullPath()),
    muonsEffectiveAreas(                iConfig.getParameter<edm::FileInPath>("muonsEffAreas").fullPath())
{
    leptonMvaComputerTTH = new LeptonMvaHelper(iConfig, true, !multilepAnalyzer->is2016() );
    leptonMvaComputertZq = new LeptonMvaHelper(iConfig, false, !multilepAnalyzer->is2016() );
    if(multilepAnalyzer->isMC()) genMatcher = new GenMatching(iConfig); // displaced specific
};

LeptonAnalyzer::~LeptonAnalyzer(){
    delete leptonMvaComputerTTH;
    delete leptonMvaComputertZq;
    if(multilepAnalyzer->isMC()) delete genMatcher; // displaced specific
}

void LeptonAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nL",                           &_nL,                           "_nL/i");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/i");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/i");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/i");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/i");
    outputTree->Branch("_rho",                          &_rho,                          "_rho/D");       // displaced specific: needed to manually apply the H/E id cut
    outputTree->Branch("_pvX",                          &_pvX,                          "_pvX/D");       // displaced specific [TODO: it seems these should be moved to multilep.cc or so]
    outputTree->Branch("_pvY",                          &_pvY,                          "_pvY/D");       // "
    outputTree->Branch("_pvZ",                          &_pvZ,                          "_pvZ/D");       // "
    outputTree->Branch("_pvXErr",                       &_pvXErr,                       "_pvXErr/D");    // "
    outputTree->Branch("_pvYErr",                       &_pvYErr,                       "_pvYErr/D");    // "
    outputTree->Branch("_pvZErr",                       &_pvZErr,                       "_pvZErr/D");    // "
    outputTree->Branch("_nVFit_os",                     &_nVFit_os,                     "_nVFit_os/i");                      // displaced specific
    outputTree->Branch("_nVFit",                        &_nVFit,                        "_nVFit/i");                         // displaced specific
    outputTree->Branch("_nGoodLeading",                 &_nGoodLeading,                 "_nGoodLeading/i");                  // "
    outputTree->Branch("_nGoodDisplaced",               &_nGoodDisplaced,               "_nGoodDisplaced/i");                // "
    outputTree->Branch("_vertices_os",                  &_vertices_os,                  "_vertices_os[_nVFit_os][12]/D");    // displaced specific
    outputTree->Branch("_lDisplaced_os",                &_lDisplaced_os,                "_lDisplaced_os[_nVFit_os][24]/D");  // "
    outputTree->Branch("_vertices",                     &_vertices,                     "_vertices[_nVFit][12]/D");          // displaced specific
    outputTree->Branch("_lDisplaced",                   &_lDisplaced,                   "_lDisplaced[_nVFit][24]/D");        // "
    outputTree->Branch("_lHasTrigger",                  &_lHasTrigger,                  "_lHasTrigger[_nL]/i");              // "
    outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
    outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
    outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
    outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
    outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
    outputTree->Branch("_lEnergySC",                    &_lEnergySC,                    "_lEnergySC[_nL]/D"); // displaced specific: supercluster energy
    outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/i");
    outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
    outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
    outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
    outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
    outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
    outputTree->Branch("_2dIP",                         &_2dIP,                         "_2dIP[_nL]/D");          // displaced specific
    outputTree->Branch("_2dIPSig",                      &_2dIPSig,                      "_2dIPSig[_nL]/D");       // "
    outputTree->Branch("_lElectronSummer16MvaGP",       &_lElectronMvaSummer16GP,       "_lElectronMvaSummer16GP[_nLight]/F");
    outputTree->Branch("_lElectronSummer16MvaHZZ",      &_lElectronMvaSummer16HZZ,      "_lElectronMvaSummer16HZZ[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17v1NoIso",    &_lElectronMvaFall17v1NoIso,    "_lElectronMvaFall17v1NoIso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17Iso",        &_lElectronMvaFall17Iso,        "_lElectronMvaFall17Iso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17NoIso",      &_lElectronMvaFall17NoIso,      "_lElectronMvaFall17NoIso[_nLight]/F");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", &_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, "_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nLight]/O"); // displaced specific (still needed ???)
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_lElectronPassMVAFall17NoIsoWP80",    &_lElectronPassMVAFall17NoIsoWP80,    "_lElectronPassMVAFall17NoIsoWP80[_nLight]/O");
    outputTree->Branch("_lElectronPassMVAFall17NoIsoWP90",    &_lElectronPassMVAFall17NoIsoWP90,    "_lElectronPassMVAFall17NoIsoWP90[_nLight]/O");
    outputTree->Branch("_lElectronPassMVAFall17NoIsoWPLoose", &_lElectronPassMVAFall17NoIsoWPLoose, "_lElectronPassMVAFall17NoIsoWPLoose[_nLight]/O");
    outputTree->Branch("_lElectronSigmaIetaIeta",               &_lElectronSigmaIetaIeta,               "_lElectronSigmaIetaIeta[_nLight]/D");
    outputTree->Branch("_lElectronDeltaPhiSuperClusterTrack",   &_lElectronDeltaPhiSuperClusterTrack,   "_lElectronDeltaPhiSuperClusterTrack[_nLight]/D");
    outputTree->Branch("_lElectronDeltaEtaSuperClusterTrack",   &_lElectronDeltaEtaSuperClusterTrack,   "_lElectronDeltaEtaSuperClusterTrack[_nLight]/D");
    outputTree->Branch("_lElectronEInvMinusPInv",               &_lElectronEInvMinusPInv,               "_lElectronEInvMinusPInv[_nLight]/D");
    outputTree->Branch("_lElectronHOverE",                      &_lElectronHOverE,                      "_lElectronHOverE[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH",                 &_leptonMvaTTH,                 "_leptonMvaTTH[_nLight]/D");
    outputTree->Branch("_leptonMvatZq",                 &_leptonMvatZq,                 "_leptonMvatZq[_nLight]/D");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
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
    outputTree->Branch("_ecalPFClusterIso",             &_ecalPFClusterIso ,            "_ecalPFClusterIso[_nLight]/D");            // displaced specific
    outputTree->Branch("_hcalPFClusterIso",             &_hcalPFClusterIso ,            "_hcalPFClusterIso[_nLight]/D");            // "
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
    outputTree->Branch("_closestJEC",          &_closestJEC,          "_closestJEC[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetE",          &_closest_lepAwareJetE,          "_closest_lepAwareJetE[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetPx",          &_closest_lepAwareJetPx,          "_closest_lepAwareJetPx[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetPy",          &_closest_lepAwareJetPy,          "_closest_lepAwareJetPy[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetPz",          &_closest_lepAwareJetPz,          "_closest_lepAwareJetPz[_nLight]/D");

    outputTree->Branch("_closest_l1JetE",          &_closest_l1JetE,          "_closest_l1JetE[_nLight]/D");
    outputTree->Branch("_closest_l1JetPx",          &_closest_l1JetPx,          "_closest_l1JetPx[_nLight]/D");
    outputTree->Branch("_closest_l1JetPy",          &_closest_l1JetPy,          "_closest_l1JetPy[_nLight]/D");
    outputTree->Branch("_closest_l1JetPz",          &_closest_l1JetPz,          "_closest_l1JetPz[_nLight]/D");

    outputTree->Branch("_closest_lJetE",          &_closest_lJetE,          "_closest_lJetE[_nLight]/D");
    outputTree->Branch("_closest_lJetPx",          &_closest_lJetPx,          "_closest_lJetPx[_nLight]/D");
    outputTree->Branch("_closest_lJetPy",          &_closest_lJetPy,          "_closest_lJetPy[_nLight]/D");
    outputTree->Branch("_closest_lJetPz",          &_closest_lJetPz,          "_closest_lJetPz[_nLight]/D");  
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv",            &_closestJetDeepCsv,            "_closestJetDeepCsv[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_b",       &_closestJetDeepFlavor_b,       "_closestJetDeepFlavor_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_bb",      &_closestJetDeepFlavor_bb,      "_closestJetDeepFlavor_bb[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_lepb",    &_closestJetDeepFlavor_lepb,    "_closestJetDeepFlavor_lepb[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor",         &_closestJetDeepFlavor,         "_closestJetDeepFlavor[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");
    if( multilepAnalyzer->isMC() ){
        outputTree->Branch("_lGenIndex",                  &_lGenIndex,                    "_lGenIndex[_nL]/i");
        outputTree->Branch("_lMatchType",                 &_lMatchType,                   "_lMatchType[_nL]/i");
        outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lIsPromptFinalState",        &_lIsPromptFinalState,          "_lIsPromptFinalState[_nL]/O");
        outputTree->Branch("_lIsPromptDecayed",           &_lIsPromptDecayed,             "_lIsPromptDecayed[_nL]/O");
        outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
        outputTree->Branch("_lMatchCharge",               &_lMatchCharge,                 "_lMatchCharge[_nLight]/I");
        outputTree->Branch("_tauGenStatus",               &_tauGenStatus,                 "_tauGenStatus[_nL]/i");
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
      outputTree->Branch("_decayModeFindingDeepTau",          &_decayModeFindingDeepTau,          "_decayModeFindingDeepTau[_nL]/O");
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
      outputTree->Branch("_tauDeepTauVsJetsRaw",          &_tauDeepTauVsJetsRaw,          "_tauDeepTauVsJetsRaw[_nL]/O");
      outputTree->Branch("_tauVVVLooseDeepTauVsJets",      &_tauVVVLooseDeepTauVsJets,      "_tauVVVLooseDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauVVLooseDeepTauVsJets",      &_tauVVLooseDeepTauVsJets,      "_tauVVLooseDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauVLooseDeepTauVsJets",       &_tauVLooseDeepTauVsJets,       "_tauVLooseDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauLooseDeepTauVsJets",        &_tauLooseDeepTauVsJets,        "_tauLooseDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauMediumDeepTauVsJets",       &_tauMediumDeepTauVsJets,       "_tauMediumDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauTightDeepTauVsJets",        &_tauTightDeepTauVsJets,        "_tauTightDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauVTightDeepTauVsJets",       &_tauVTightDeepTauVsJets,       "_tauVTightDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauVVTightDeepTauVsJets",      &_tauVVTightDeepTauVsJets,      "_tauVVTightDeepTauVsJets[_nL]/O");
      outputTree->Branch("_tauDeepTauVsEleRaw",           &_tauDeepTauVsEleRaw,           "_tauDeepTauVsEleRaw[_nL]/O");
      outputTree->Branch("_tauVVVLooseDeepTauVsEle",       &_tauVVVLooseDeepTauVsEle,       "_tauVVVLooseDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauVVLooseDeepTauVsEle",       &_tauVVLooseDeepTauVsEle,       "_tauVVLooseDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauVLooseDeepTauVsEle",        &_tauVLooseDeepTauVsEle,        "_tauVLooseDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauLooseDeepTauVsEle",         &_tauLooseDeepTauVsEle,         "_tauLooseDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauMediumDeepTauVsEle",        &_tauMediumDeepTauVsEle,        "_tauMediumDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauTightDeepTauVsEle",         &_tauTightDeepTauVsEle,         "_tauTightDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauVTightDeepTauVsEle",        &_tauVTightDeepTauVsEle,        "_tauVTightDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauVVTightDeepTauVsEle",       &_tauVVTightDeepTauVsEle,       "_tauVVTightDeepTauVsEle[_nL]/O");
      outputTree->Branch("_tauDeepTauMuRaw",              &_tauDeepTauVsMuRaw,            "_tauDeepTauVsEleMu[_nL]/O");
      outputTree->Branch("_tauVLooseDeepTauVsMu",         &_tauVLooseDeepTauVsMu,         "_tauVLooseDeepTauVsMu[_nL]/O");
      outputTree->Branch("_tauLooseDeepTauVsMu",          &_tauLooseDeepTauVsMu,          "_tauLooseDeepTauVsMu[_nL]/O");
      outputTree->Branch("_tauMediumDeepTauVsMu",         &_tauMediumDeepTauVsMu,         "_tauMediumDeepTauVsMu[_nL]/O");
      outputTree->Branch("_tauTightDeepTauVsMu",          &_tauTightDeepTauVsMu,          "_tauTightDeepTauVsMu[_nL]/O");
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

    iSetup.get<IdealMagneticFieldRecord>().get(_bField);
    iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", _shProp);
    GsfPropagatorAdapter gsfPropagator(AnalyticalPropagator(&(*_bField), anyDirection));
    _gsfProp = new TransverseImpactPointExtrapolator(gsfPropagator);

    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;
    _rho    = *rho;
    _nVFit  = 0;
    _nVFit_os = 0;
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
    if(multilepAnalyzer->isMC()){
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

        fillLeptonImpactParameters(mu);

        fillLeptonKinVars(mu);
        fillLeptonIsoVars(mu, *rho);
        if( multilepAnalyzer->isMC() ) fillLeptonGenVars(mu, *genParticles);

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

        _relIso[_nL]         = getRelIso03(mu, *rho, muonsEffectiveAreas);                     // Isolation variables
        _relIso0p4[_nL]      = getRelIso04(mu, *rho, muonsEffectiveAreas);
        _relIso0p4MuDeltaBeta[_nL] = getRelIso04(mu, *rho, muonsEffectiveAreas, true);
        _miniIso[_nL]        = getMiniIsolation(mu, *rho, muonsEffectiveAreas, false);
        _miniIsoCharged[_nL] = getMiniIsolation(mu, *rho, muonsEffectiveAreas, true);

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

        _muDTStationsWithValidHits[_nL]   = mu.bestTrack()->hitPattern().dtStationsWithValidHits();
        _muCSCStationsWithValidHits[_nL]  = mu.bestTrack()->hitPattern().cscStationsWithValidHits();
        _muRPCStationsWithValidHits[_nL]  = mu.bestTrack()->hitPattern().rpcStationsWithValidHits();
        _muMuonStationsWithValidHits[_nL] = mu.bestTrack()->hitPattern().muonStationsWithValidHits();

        // maybe put this simply in a function in LeptonAnalyzerId.cc???
        bool someIdWithoutName = (mu.pt() > 24 and _lPOGMedium[_nL]
                                  and fabs(_dxy[_nL]) < 0.05 and fabs(_dz[_nL])< 0.1
                                  and _relIso[_nL] < 0.2
                                  and !mu.innerTrack().isNull()
                                  and (mu.isTrackerMuon() or mu.isGlobalMuon()));

        if(someIdWithoutName) ++_nGoodLeading;
        if((multilepAnalyzer->skim == "trilep" or multilepAnalyzer->skim == "displtrilep") && someIdWithoutName) _lHasTrigger[_nL] = matchSingleTrigger(iEvent, false, _lEta[_nL], _lPhi[_nL]);
        else                                                                                                     _lHasTrigger[_nL] = 0;

        ++_nMu;
        ++_nL;
        ++_nLight;
    }


    // Loop over electrons (note: using iterator we can easily get the ref too)
    for(size_t iele=0; iele<selelectrons.size(); ++iele){
        if(_nL == nL_max) break;
        const pat::Electron *ele = selelectrons[iele];

        fillLeptonKinVars(*ele);
        if( multilepAnalyzer->isMC() ) fillLeptonGenVars(*ele, *genParticles);

        fillLeptonImpactParameters(*ele);
        fillLeptonIsoVars(*ele, *rho);


        _lFlavor[_nL]                   = 0;
        _lEtaSC[_nL]                    = ele->superCluster()->eta();

        _relIso[_nL]                    = getRelIso03(*ele, *rho, electronsEffectiveAreas);
        _relIso0p4[_nL]                 = getRelIso04(*ele, *rho, electronsEffectiveAreas);
        _relIso_ttH[_nL]                = getRelIso03(*ele, *rho, electronsEffectiveAreas_ttH_relIso);
        _relIso0p4_ttH[_nL]             = getRelIso04(*ele, *rho, electronsEffectiveAreas_ttH_relIso);
        _miniIso[_nL]                   = getMiniIsolation(*ele, *rho, electronsEffectiveAreas_ttH_miniIso, false);
        _miniIsoCharged[_nL]            = getMiniIsolation(*ele, *rho, electronsEffectiveAreas_ttH_miniIso, true);
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
        _lEnergySC[_nL]                 = ele->superCluster()->energy();

        _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL] = isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&*ele);

        _lPOGVeto[_nL]                  = ele->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
        _lPOGLoose[_nL]                 = ele->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
        _lPOGMedium[_nL]                = ele->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
        _lPOGTight[_nL]                 = ele->electronID("cutBasedElectronID-Fall17-94X-V1-tight");

        _lElectronPassMVAFall17NoIsoWPLoose[_nL] = ele->electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
        _lElectronPassMVAFall17NoIsoWP90[_nL] = ele->electronID("mvaEleID-Fall17-noIso-V2-wp90");
        _lElectronPassMVAFall17NoIsoWP80[_nL] = ele->electronID("mvaEleID-Fall17-noIso-V2-wp80");

        _lElectronSigmaIetaIeta[_nL] = ele->full5x5_sigmaIetaIeta();
        _lElectronDeltaPhiSuperClusterTrack[_nL] = fabs(ele->deltaPhiSuperClusterTrackAtVtx());
        _lElectronDeltaEtaSuperClusterTrack[_nL] = fabs(ele->deltaEtaSuperClusterTrackAtVtx());
        _lElectronEInvMinusPInv[_nL] = (1.0 - ele->eSuperClusterOverP())/ele->correctedEcalEnergy();
        _lElectronHOverE[_nL] = ele->hadronicOverEm();

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


        // maybe put this simply in a function in LeptonAnalyzerId.cc???
        bool someIdWithoutName = (ele->pt() > 27 
                                  and fabs(_dxy[_nL]) < 0.05 and fabs(_dz[_nL])< 0.1
                                  and _relIso[_nL] < 0.2
                                  and !ele->gsfTrack().isNull()
                                  and _lElectronMissingHits[_nL] <=2
                                  and ele->passConversionVeto());

        if(ele->pt() > 7 && fabs(_dxy[_nL]) > 0.02) ++_nGoodDisplaced;
        if(someIdWithoutName) ++_nGoodLeading;
        if((multilepAnalyzer->skim == "trilep" or multilepAnalyzer->skim == "displtrilep")  && someIdWithoutName) _lHasTrigger[_nL] = matchSingleTrigger(iEvent, true, _lEta[_nL], _lPhi[_nL]);
        else                                                                                                      _lHasTrigger[_nL] = 0;

        ++_nEle;
        ++_nL;
        ++_nLight;
    }

    //Initialize with default values for those electron-only arrays which weren't filled with muons [to allow correct comparison by the test script]
    for(auto array : {_lEtaSC, _lEnergySC}) std::fill_n(array, _nMu, 0.);
    for(auto array : {_lElectronMvaSummer16GP, _lElectronMvaSummer16HZZ, _lElectronMvaFall17v1NoIso}) std::fill_n(array, _nMu, 0.); // OLD, do not use them
    for(auto array : {_lElectronMvaFall17Iso, _lElectronMvaFall17NoIso}) std::fill_n(array, _nMu, 0.);
    for(auto array : {_lElectronPassMVAFall17NoIsoWPLoose, _lElectronPassMVAFall17NoIsoWP90, _lElectronPassMVAFall17NoIsoWP80}) std::fill_n(array, _nMu, false);
    for(auto array : {_lElectronPassEmu, _lElectronPassConvVeto, _lElectronChargeConst}) std::fill_n(array, _nMu, false);
    for(auto array : {_lElectronMissingHits}) std::fill_n(array, _nMu, 0.);
    for(auto array : {_lElectronSigmaIetaIeta, _lElectronDeltaPhiSuperClusterTrack, _lElectronDeltaEtaSuperClusterTrack, _lElectronEInvMinusPInv, _lElectronHOverE} ) std::fill_n( array, _nMu, 0. );
    for(auto array : {_lPtCorr, _lPtScaleUp, _lPtScaleDown, _lPtResUp, _lPtResDown}) std::fill_n(array, _nMu, 0.);
    for(auto array : {_lECorr, _lEScaleUp, _lEScaleDown, _lEResUp, _lEResDown}) std::fill_n(array, _nMu, 0.);
    for(auto array : {_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, _lEleIsEB, _lEleIsEE}) std::fill_n(array, _nMu, false);      // displaced specific
    for(auto array : {_lEleSuperClusterOverP, _lEleEcalEnergy, _lElefull5x5SigmaIetaIeta, _lEleDEtaInSeed}) std::fill_n(array, _nMu, -1.); // displaced speficic
    for(auto array : {_lEleDeltaPhiSuperClusterTrackAtVtx, _lElehadronicOverEm, _lEleInvMinusPInv}) std::fill_n(array, _nMu, -1.);         // displaced specific
    for(auto array : {_ecalPFClusterIso, _hcalPFClusterIso}) std::fill_n(array, _nMu, -1.);                                                // displaced specific

    //loop over taus
    for(const pat::Tau* tauptr : seltaus) {
        break; // displaced specific: do not consider tau's in the displaced branch [note: switching this back on without reviewing the code will put unitialized values in the tree, breaking the tests]
        if(_nL == nL_max) break;
        const pat::Tau& tau = (*tauptr);

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
            _decayModeFindingDeepTau[_nL]       = tau.tauID("decayModeFindingNewDMs") and _tauDecayMode[_nL] != 5 and _tauDecayMode[_nL] != 6;                   //new Tau ID
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
           
             
            _tauDeepTauVsJetsRaw[_nL]           = tau.tauID("byDeepTau2017v2p1VSjetraw");
            _tauVVVLooseDeepTauVsJets[_nL]           = tau.tauID("byVVVLooseDeepTau2017v2p1VSjet");
            _tauVVLooseDeepTauVsJets[_nL]           = tau.tauID("byVVLooseDeepTau2017v2p1VSjet");
            _tauVLooseDeepTauVsJets[_nL]           = tau.tauID("byVLooseDeepTau2017v2p1VSjet");
            _tauLooseDeepTauVsJets[_nL]           = tau.tauID("byLooseDeepTau2017v2p1VSjet");
            _tauMediumDeepTauVsJets[_nL]           = tau.tauID("byMediumDeepTau2017v2p1VSjet");
            _tauTightDeepTauVsJets[_nL]           = tau.tauID("byTightDeepTau2017v2p1VSjet");
            _tauVTightDeepTauVsJets[_nL]           = tau.tauID("byVTightDeepTau2017v2p1VSjet");
            _tauVVTightDeepTauVsJets[_nL]           = tau.tauID("byVVTightDeepTau2017v2p1VSjet");
            
            _tauDeepTauVsEleRaw[_nL]           = tau.tauID("byDeepTau2017v2p1VSeraw");
            _tauVVVLooseDeepTauVsEle[_nL]           = tau.tauID("byVVVLooseDeepTau2017v2p1VSe");
            _tauVVLooseDeepTauVsEle[_nL]           = tau.tauID("byVVLooseDeepTau2017v2p1VSe");
            _tauVLooseDeepTauVsEle[_nL]           = tau.tauID("byVLooseDeepTau2017v2p1VSe");
            _tauLooseDeepTauVsEle[_nL]           = tau.tauID("byLooseDeepTau2017v2p1VSe");
            _tauMediumDeepTauVsEle[_nL]           = tau.tauID("byMediumDeepTau2017v2p1VSe");
            _tauTightDeepTauVsEle[_nL]           = tau.tauID("byTightDeepTau2017v2p1VSe");
            _tauVTightDeepTauVsEle[_nL]           = tau.tauID("byVTightDeepTau2017v2p1VSe");
            _tauVVTightDeepTauVsEle[_nL]           = tau.tauID("byVVTightDeepTau2017v2p1VSe");
            
            _tauDeepTauVsMuRaw[_nL]           = tau.tauID("byDeepTau2017v2p1VSmuraw");
            _tauVLooseDeepTauVsMu[_nL]           = tau.tauID("byVLooseDeepTau2017v2p1VSmu");
            _tauLooseDeepTauVsMu[_nL]           = tau.tauID("byLooseDeepTau2017v2p1VSmu");
            _tauMediumDeepTauVsMu[_nL]           = tau.tauID("byMediumDeepTau2017v2p1VSmu");
            _tauTightDeepTauVsMu[_nL]           = tau.tauID("byTightDeepTau2017v2p1VSmu");
    
        }

        ++_nTau;
        ++_nL;
    }


    /*
     * refitting vertices displaced
     * try all possible combinations of electrons and muons (further selections on opposite-charge is required in the fillDileptonVertexArray function)
     */
    cleanDileptonVertexArrays(_nVFit);
    for(unsigned i=0; i < selmuons.size(); ++i){
      for(unsigned j=0; j < selmuons.size(); ++j)     fillDileptonVertexArrays(i, j,      selmuons.at(i), selmuons.at(j), true);
      for(unsigned j=0; j < selelectrons.size(); ++j) fillDileptonVertexArrays(i, _nMu+j, selmuons.at(i), selelectrons.at(j), true);
    }
    for(unsigned i=0; i < selelectrons.size(); ++i){
      for(unsigned j=0; j < selmuons.size(); ++j)     fillDileptonVertexArrays(_nMu+i, j,      selelectrons.at(i), selmuons.at(j), true);
      for(unsigned j=0; j < selelectrons.size(); ++j) fillDileptonVertexArrays(_nMu+i, _nMu+j, selelectrons.at(i), selelectrons.at(j), true);
    }
    
   for(unsigned i=0; i < selmuons.size(); ++i){
      for(unsigned j=0; j < selmuons.size(); ++j)     fillDileptonVertexArrays(i, j,      selmuons.at(i), selmuons.at(j));
      for(unsigned j=0; j < selelectrons.size(); ++j) fillDileptonVertexArrays(i, _nMu+j, selmuons.at(i), selelectrons.at(j));
    }
    for(unsigned i=0; i < selelectrons.size(); ++i){
      //for(unsigned j=0; j < selmuons.size(); ++j)     fillDileptonVertexArrays(_nMu+i, j,      selelectrons.at(i), selmuons.at(j));
      for(unsigned j=0; j < selelectrons.size(); ++j) fillDileptonVertexArrays(_nMu+i, _nMu+j, selelectrons.at(i), selelectrons.at(j));
    }

    //Initialize with default values for those tau-only arrays which weren't filled with electrons and muons [to allow correct comparison by the test script]
    for(auto array : {_tauMuonVetoLoose, _tauEleVetoLoose, _decayModeFinding, _decayModeFindingNew, _decayModeFindingDeepTau}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauEleVetoVLoose, _tauEleVetoMedium, _tauEleVetoTight, _tauEleVetoVTight, _tauMuonVetoTight, _decayModeFindingNew}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauVLooseMvaNew2015, _tauLooseMvaNew2015, _tauMediumMvaNew2015, _tauTightMvaNew2015, _tauVTightMvaNew2015}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauVLooseMvaNew2017v2, _tauLooseMvaNew2017v2, _tauMediumMvaNew2017v2, _tauTightMvaNew2017v2, _tauVTightMvaNew2017v2}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauVVVLooseDeepTauVsJets, _tauVVLooseDeepTauVsJets, _tauVLooseDeepTauVsJets, _tauMediumDeepTauVsJets}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauTightDeepTauVsJets, _tauVTightDeepTauVsJets, _tauVVTightDeepTauVsJets}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauVVVLooseDeepTauVsEle, _tauVVLooseDeepTauVsEle, _tauVLooseDeepTauVsEle, _tauMediumDeepTauVsEle}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauTightDeepTauVsEle, _tauVTightDeepTauVsEle, _tauVVTightDeepTauVsEle}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauVLooseDeepTauVsMu, _tauLooseDeepTauVsMu, _tauMediumDeepTauVsMu, _tauTightDeepTauVsMu}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauPOGVLoose2015, _tauPOGLoose2015, _tauPOGMedium2015, _tauPOGTight2015, _tauPOGVTight2015}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauPOGVVLoose2017v2, _tauPOGVTight2017v2, _tauPOGVVTight2017v2}) std::fill_n(array, _nLight, false);
    for(auto array : {_tauAgainstElectronMVA6Raw, _tauCombinedIsoDBRaw3Hits, _tauIsoMVAPWdR03oldDMwLT}) std::fill_n(array, _nLight, 0.);
    for(auto array : {_tauDecayMode}) std::fill_n(array, _nLight, 0);
    for(auto array : {_tauIsoMVADBdR03oldDMwLT, _tauIsoMVADBdR03newDMwLT, _tauIsoMVAPWnewDMwLT, _tauIsoMVAPWoldDMwLT}) std::fill_n(array, _nLight, 0.);

    if(multilepAnalyzer->skim == "trilep"      and (_nLight < 3 || _nGoodLeading < 1                       ) ) return false;
    if(multilepAnalyzer->skim == "displtrilep" and (_nLight < 3 || _nGoodLeading < 1 || _nGoodDisplaced < 2) ) return false;

    if(multilepAnalyzer->skim == "dilep"       and _nLight <  2) return false;
    if(multilepAnalyzer->skim == "ttg"         and _nLight <  2) return false;
    if(multilepAnalyzer->skim == "singlelep"   and _nLight <  1) return false;
    if(multilepAnalyzer->skim == "FR"          and _nLight != 1 ) return false;
    if(multilepAnalyzer->skim == "FRsM"        and _nLight != 1 ) return false;

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
      if (j < 12) {
          _vertices[i][j] = 0;
          _vertices_os[i][j] = 0;
      }
      _lDisplaced[i][j] = 0;
      _lDisplaced_os[i][j] = 0;
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
void LeptonAnalyzer::fillDileptonVertexArrays(unsigned iL_plus, unsigned iL_minus, const reco::RecoCandidate* lep1, const reco::RecoCandidate* lep2, const bool ensureOppositeSign){
  if(iL_plus==iL_minus) return;
  if(ensureOppositeSign and lep1->charge() < 0) return; // ensure opposite charge
  if(ensureOppositeSign and lep2->charge() > 0) return;
    
  TransientVertex dvtx = dileptonVertex(lep1, lep2);
  if(!dvtx.isValid()){
    std::cout << " *** WARNING: refitted dilepton vertex is not valid! " << std::endl;
    return;
  }

  auto& vertexArray = ensureOppositeSign ? _vertices_os : _vertices;
  auto& displArray  = ensureOppositeSign ? _lDisplaced_os : _lDisplaced;
  auto& vertexIndex = ensureOppositeSign ? _nVFit_os : _nVFit;

  vertexArray[vertexIndex][0]  = (iL_plus+1)*100 + (iL_minus+1); // indices i and j stored as i*100+j, but start counting at 1 instead of 0 (historical, backwards compatibility)
  vertexArray[vertexIndex][1]  = dvtx.position().x();
  vertexArray[vertexIndex][2]  = dvtx.position().y();
  vertexArray[vertexIndex][3]  = dvtx.position().z();
  vertexArray[vertexIndex][4]  = dvtx.positionError().cxx();
  vertexArray[vertexIndex][5]  = dvtx.positionError().cyy();
  vertexArray[vertexIndex][6]  = dvtx.positionError().czz();
  vertexArray[vertexIndex][7]  = dvtx.positionError().cyx();
  vertexArray[vertexIndex][8]  = dvtx.positionError().czy();
  vertexArray[vertexIndex][9]  = dvtx.positionError().czx();
  vertexArray[vertexIndex][10] = dvtx.degreesOfFreedom();
  vertexArray[vertexIndex][11] = dvtx.totalChiSquared();

  int i=0;
  for(auto lep : {lep1, lep2}){
    const reco::Track& track = getTrack(lep);
    GlobalPoint  globPoint(track.vx(), track.vy(), track.vz());
    GlobalVector globVector(track.px(), track.py(), track.pz());
    GlobalTrajectoryParameters globTrajParam(globPoint, globVector, track.charge(), _bField.product());
    FreeTrajectoryState tempFreeTrajectoryState(globTrajParam, track.covariance());
    FreeTrajectoryState freeTrajectoryState;

    if(lep->isMuon()){
      freeTrajectoryState = _shProp->propagate(tempFreeTrajectoryState, dvtx.position());
    } else {
      auto trajectoryOnSurface = _gsfProp->extrapolate(tempFreeTrajectoryState, dvtx.position());
      if(trajectoryOnSurface.isValid()) freeTrajectoryState = *(trajectoryOnSurface.freeState());
    }

    if(!freeTrajectoryState.hasCurvilinearError()){
      std::cout << "Propagation of lepton to dilepton vertex (" << _vertices[_nVFit][0] << ") failed!" << std::endl;
      freeTrajectoryState = tempFreeTrajectoryState;
    }

    // fill _lDisplaced[nVFit] with 2x12=24 elements structured as 3 position + 3 momentum + 3 position error + 3 momentum error
    displArray[vertexIndex][i++] = freeTrajectoryState.position().x();
    displArray[vertexIndex][i++] = freeTrajectoryState.position().y();
    displArray[vertexIndex][i++] = freeTrajectoryState.position().z();

    displArray[vertexIndex][i++] = freeTrajectoryState.momentum().x();
    displArray[vertexIndex][i++] = freeTrajectoryState.momentum().y();
    displArray[vertexIndex][i++] = freeTrajectoryState.momentum().z();

    displArray[vertexIndex][i++] = (freeTrajectoryState.cartesianError().matrix())(0, 0); //position error
    displArray[vertexIndex][i++] = (freeTrajectoryState.cartesianError().matrix())(1, 1);
    displArray[vertexIndex][i++] = (freeTrajectoryState.cartesianError().matrix())(2, 2);

    displArray[vertexIndex][i++] = (freeTrajectoryState.cartesianError().matrix())(3, 3); // momentum error
    displArray[vertexIndex][i++] = (freeTrajectoryState.cartesianError().matrix())(4, 4);
    displArray[vertexIndex][i++] = (freeTrajectoryState.cartesianError().matrix())(5, 5);
  }

  vertexIndex++;
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
  _absIso03[_nL]             = getRelIso03(mu, rho, muonsEffectiveAreas)*mu.pt();
  _absIso04[_nL]             = getRelIso04(mu, rho, muonsEffectiveAreas, false)*mu.pt();
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
  _absIso03[_nL]             = getRelIso03(ele, rho, electronsEffectiveAreas)*ele.pt();
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

    _lProvenance[_nL]           = GenTools::provenance(match, genParticles); 
    _lProvenanceCompressed[_nL] = GenTools::provenanceCompressed(match, genParticles, _lIsPrompt[_nL]); 
    _lProvenanceConversion[_nL] = GenTools::provenanceConversion(match, genParticles);
    _lMomPdgId[_nL]             = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;

    _tauGenStatus[_nL]          = TauTools::tauGenStatus(match);        

    _lMatchPdgId[_nL]           = match ? match->pdgId() : 0;
    _lMatchCharge[_nL]          = match ? match->charge() : 0;
    _lMatchPt[_nL]              = match ? match->pt() : 0; 
    _lMatchEta[_nL]             = match ? match->eta() : 0;
    _lMatchPhi[_nL]             = match ? match->phi() : 0;
    _lMatchVertexX[_nL]         = match ? match->vertex().x() : 0;
    _lMatchVertexY[_nL]         = match ? match->vertex().y() : 0;
    _lMatchVertexZ[_nL]         = match ? match->vertex().z() : 0;
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
    _2dIP[_nL]    = ele.dB();
    _2dIPSig[_nL] = ele.dB()/ele.edB();
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon){
    _dxy[_nL]     = muon.dB(pat::Muon::PV2D);
    _dz[_nL]      = muon.dB(pat::Muon::PVDZ);
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
  if(getRelIso03(elec, rho, electronsEffectiveAreas) >  2)                                        return false;
  if(elec.pt()<5.)                                                                                return false;
  if(std::abs(elec.eta())>2.5)                                                                    return false;
  //if(!isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(&elec)) return false; // Note: should be reviewd, especially for 2017-2018
  if(eleMuOverlap(elec, _lPFMuon))                                                                return false; // overlap muon-electron deltaR<0.05, using PF muons
  return true;
}

bool LeptonAnalyzer::passMuonPreselection(const pat::Muon& muon, const double rho) const {
  if(!muon.isPFMuon())                                return false;
  if(!muon.isLooseMuon())                             return false;
  if(getRelIso03(muon, rho, muonsEffectiveAreas) > 2) return false;
  if(muon.pt() < 3.)                                  return false;
  if(std::abs(muon.eta()) > 2.4)                      return false;
  return true;
}

bool LeptonAnalyzer::passTauPreselection(const pat::Tau& tau, const reco::Vertex::Point& vertex) const {
  if(tau.pt() < 20.)            return false; // Minimum pt for tau reconstruction
  if(std::abs(tau.eta()) > 2.3) return false;
  if(tau_dz(tau, vertex) > 2.)  return false;
  return true;
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
            _ptRatio[_nL] = ( oldMatching ? 1. : 1. / ( 1. + _relIso0p4_ttH[_nL] ) );
        }
        _ptRel[_nL] = 0;
        _selectedTrackMult[_nL] = 0;
        _closestJetDeepFlavor_b[_nL] = 0;
        _closestJetDeepFlavor_bb[_nL] = 0;
        _closestJetDeepFlavor_lepb[_nL] = 0;
        _closestJetDeepFlavor[_nL] = 0;
        _closestJetDeepCsv_b[_nL] = 0;
        _closestJetDeepCsv_bb[_nL] = 0;
        _selectedTrackMult[_nL]    = 0;
        _closestJEC[_nL]           = 0;
        _closest_lepAwareJetE[_nL]  =0;
        _closest_lepAwareJetPx[_nL]  =0;
        _closest_lepAwareJetPy[_nL]  =0;
        _closest_lepAwareJetPz[_nL]  =0;
        _closest_l1JetE[_nL]        = 0;
        _closest_l1JetPx[_nL]        = 0;
        _closest_l1JetPy[_nL]        = 0;
        _closest_l1JetPz[_nL]        = 0;
        _closest_lJetE [_nL]        = 0;
        _closest_lJetPx [_nL]        = 0;
        _closest_lJetPy [_nL]        = 0;
        _closest_lJetPz [_nL]        = 0;
        _closestJetDeepCsv[_nL] = 0;
        _closestJetCsvV2[_nL] = 0;
    } else {
        const pat::Jet& jet = *matchedJetPtr;

        auto rawJetP4 = jet.correctedP4("Uncorrected"); 
        auto leptonP4 = lepton.p4();

        bool leptonEqualsJet = ( ( rawJetP4 - leptonP4 ).P() < 1e-4 );

        //if lepton and jet vector are equal set _ptRatio, _ptRel and track multipliticy to defaults 
        if( leptonEqualsJet && !oldMatching){
            _ptRatio[_nL] = 1;
            _ptRel[_nL] = 0;
            _selectedTrackMult[_nL] = 0;
            //displaced specific (TODO: cleanup?)
            _closestJEC[_nL]            = 0;
            _closest_lepAwareJetE[_nL]  = 0;
            _closest_lepAwareJetPx[_nL] = 0;
            _closest_lepAwareJetPy[_nL] = 0;
            _closest_lepAwareJetPz[_nL] = 0;
            _closest_l1JetE[_nL]        = 0;
            _closest_l1JetPx[_nL]       = 0;
            _closest_l1JetPy[_nL]       = 0;
            _closest_l1JetPz[_nL]       = 0;
            _closest_lJetE[_nL]         = 0;
            _closest_lJetPx[_nL]        = 0;
            _closest_lJetPy[_nL]        = 0;
            _closest_lJetPz[_nL]        = 0;
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
            //displaced specific
            _closestJEC[_nL]            = L2L3JEC;
            _closest_lepAwareJetE[_nL]  = lepAwareJetP4.E();
            _closest_lepAwareJetPx[_nL] = lepAwareJetP4.Px();
            _closest_lepAwareJetPy[_nL] = lepAwareJetP4.Py();
            _closest_lepAwareJetPz[_nL] = lepAwareJetP4.Pz();
            _closest_l1JetE[_nL]        = L1JetP4.E();
            _closest_l1JetPx[_nL]       = L1JetP4.Px();
            _closest_l1JetPy[_nL]       = L1JetP4.Py();
            _closest_l1JetPz[_nL]       = L1JetP4.Pz();
            _closest_lJetE[_nL]         = jet.p4().E();
            _closest_lJetPx[_nL]        = jet.p4().Px();
            _closest_lJetPy[_nL]        = jet.p4().Py();
            _closest_lJetPz[_nL]        = jet.p4().Pz();
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
    }
}

unsigned LeptonAnalyzer::matchSingleTrigger(const edm::Event& iEvent, bool isele, double aeta, double aphi){
  auto triggerResults = getHandle(iEvent, multilepAnalyzer->triggerToken);
  auto triggerObjects = getHandle(iEvent, multilepAnalyzer->trigObjToken);

  std::vector<std::string> &singletrigs = isele ? singleEleTrigs : singleMuoTrigs;

  unsigned trigmask(0);
  for(pat::TriggerObjectStandAlone iobj : *triggerObjects) { // NOTE: not const nor by reference, because we need to 'unpackPathNames'
    if(reco::deltaR(iobj.eta(), iobj.phi(), aeta, aphi)<0.15) {
      iobj.unpackPathNames(iEvent.triggerNames(*triggerResults));
      iobj.unpackFilterLabels(iEvent, *triggerResults);
      int ipath(-1);
      for(std::string& itrig : singletrigs) {
        ++ipath;
        if(multilepAnalyzer->is2017() and itrig == "HLT_Ele32_WPTight_Gsf"){ // special case for the 2017 Ele32 trigger
          if(!iobj.hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter")) continue;
          if(!iobj.hasFilterLabel("hltEGL1SingleEGOrFilter")) continue;
        } else {
          if(!iobj.hasPathName(itrig.c_str(), true, true)) continue;
        }
        trigmask |= (1<<ipath);
      }
    }
  }
  return trigmask;
}
