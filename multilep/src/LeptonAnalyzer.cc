#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>


// TODO: we should maybe stop indentifying effective areas by year, as they are typically more connected to a specific ID than to a specific year
LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer),
    electronsEffectiveAreas(iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas").fullPath()),
    muonsEffectiveAreas    ((multilepAnalyzer->is2017 || multilepAnalyzer->is2018)? (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath() )
{
    leptonMvaComputerSUSY16   = new LeptonMvaHelper(iConfig, 0, false); // SUSY         // TODO: all of these really needed? could use some clean-up
    leptonMvaComputerTTH16    = new LeptonMvaHelper(iConfig, 1, false); // TTH
    leptonMvaComputerSUSY17   = new LeptonMvaHelper(iConfig, 0, true);  // SUSY
    leptonMvaComputerTTH17    = new LeptonMvaHelper(iConfig, 1, true);  // TTH
    leptonMvaComputertZqTTV16 = new LeptonMvaHelper(iConfig, 2, false); // tZq/TTV
    leptonMvaComputertZqTTV17 = new LeptonMvaHelper(iConfig, 2, true);  // tZq/TTV
};

LeptonAnalyzer::~LeptonAnalyzer(){
    delete leptonMvaComputerSUSY16;
    delete leptonMvaComputerTTH16;
    delete leptonMvaComputertZqTTV16;
    delete leptonMvaComputerSUSY17;
    delete leptonMvaComputerTTH17;
    delete leptonMvaComputertZqTTV17;
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
    outputTree->Branch("_leptonMvaSUSY16",              &_leptonMvaSUSY16,              "_leptonMvaSUSY16[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH16",               &_leptonMvaTTH16,               "_leptonMvaTTH16[_nLight]/D");
    outputTree->Branch("_leptonMvaSUSY17",              &_leptonMvaSUSY17,              "_leptonMvaSUSY17[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH17",               &_leptonMvaTTH17,               "_leptonMvaTTH17[_nLight]/D");
    outputTree->Branch("_leptonMvatZqTTV16",            &_leptonMvatZqTTV16,            "_leptonMvatZqTTV16[_nLight]/D");
    outputTree->Branch("_leptonMvatZqTTV17",            &_leptonMvatZqTTV17,            "_leptonMvatZqTTV17[_nLight]/D");
    outputTree->Branch("_lHNLoose",                     &_lHNLoose,                     "_lHNLoose[_nLight]/O");
    outputTree->Branch("_lHNFO",                        &_lHNFO,                        "_lHNFO[_nLight]/O");
    outputTree->Branch("_lHNTight",                     &_lHNTight,                     "_lHNTight[_nLight]/O");
    outputTree->Branch("_lEwkLoose",                    &_lEwkLoose,                    "_lEwkLoose[_nL]/O");
    outputTree->Branch("_lEwkFO",                       &_lEwkFO,                       "_lEwkFO[_nL]/O");
    outputTree->Branch("_lEwkTight",                    &_lEwkTight,                    "_lEwkTight[_nL]/O");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
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
    if(!multilepAnalyzer->isData){
      outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
      outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
      outputTree->Branch("_lMomPdgId",                  &_lMomPdgId,                    "_lMomPdgId[_nL]/I");
      outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
      outputTree->Branch("_lProvenanceCompressed",      &_lProvenanceCompressed,        "_lProvenanceCompressed[_nL]/i");
      outputTree->Branch("_lProvenanceConversion",      &_lProvenanceConversion,        "_lProvenanceConversion[_nL]/i");
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
    //std::cout << "begin leptonanalyzer" << std::endl;
    edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
    edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
    edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(multilepAnalyzer->tauToken,                          taus);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
    edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);
    edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);
    //edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetSmearedToken,                   jets);  // Are we sure we do not want the smeared jets here???
    edm::Handle<std::vector<reco::GenParticle>> genParticles;        iEvent.getByToken(multilepAnalyzer->genParticleToken,                  genParticles);
    edm::Handle<std::vector<reco::Vertex>> secVertices;              iEvent.getByToken(multilepAnalyzer->secondaryVerticesToken, secVertices);
    edm::Handle<std::vector<reco::Vertex>> primvertices;             iEvent.getByToken(multilepAnalyzer->vtxToken, primvertices);

    iSetup.get<IdealMagneticFieldRecord>().get(_bField);
    iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", _shProp);
    const reco::Vertex pv = *(primvertices->begin());


    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;

    //std::cout << std::endl << std::endl << std::endl << "=== EVENT ===" << std::endl;
    //std::cout << "sec vertices size: " << (*secVertices).size() << std::endl;
    //fillAllIVFVariables(*secVertices, pv);

    //loop over muons
    // muons need to be run first, because some ID's need to calculate a muon veto for electrons
    for(const pat::Muon& mu : *muons){
        if(_nL == nL_max)                              break;
        if(mu.innerTrack().isNull())                   continue;
        if(mu.pt() < 5)                                continue;
        if(fabs(mu.eta()) > 2.4)                       continue;
        if(!mu.isPFMuon())                             continue;
        if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue;
        fillLeptonImpactParameters(mu, primaryVertex);
        //if(fabs(_dxy[_nL]) > 0.05)                     continue; this can be applied at root level, but not here to preserve displaced signal
        //if(fabs(_dz[_nL]) > 0.1)                       continue;
        fillLeptonKinVars(mu);
        if(!multilepAnalyzer->isData) fillLeptonGenVars(mu, *genParticles);
        fillLeptonJetVariables(mu, jets, primaryVertex, *rho);
        fillDisplacedIDVariables(mu);

	    //std::vector<reco::Track> vertex_tracks;
	    //vertex_tracks.push_back(*mu.bestTrack());
	    //fillLeptonKVFVariables(packedCands, vertex_tracks);
        //fillLeptonIVFVariables(*mu.bestTrack(), *secVertices);
        //std::cout << std::endl << "--- muon pt, eta, phi: " << (*mu.bestTrack()).pt() << " " << (*mu.bestTrack()).eta() << " " << (*mu.bestTrack()).phi() << " " << mu.numberOfSourceCandidatePtrs() << " "; 
        //if(mu.numberOfSourceCandidatePtrs() > 0) std::cout << mu.sourceCandidatePtr(0)->pt() << std::endl;
        fillMatchingIVFVariables(*secVertices, mu, pv);

        _lFlavor[_nL]        = 1;
        _lMuonSegComp[_nL]    = mu.segmentCompatibility();
        _lMuonTrackPt[_nL]    = mu.innerTrack()->pt();
        _lMuonTrackPtErr[_nL] = mu.innerTrack()->ptError();

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
        fillLeptonImpactParameters(*ele, primaryVertex);
        //if(fabs(_dxy[_nL]) > 0.05)                                                               continue;   //same argument as for muons, dont harm displaced samples and can be applied at root level
        //if(fabs(_dz[_nL]) > 0.1)                                                                 continue;
        fillLeptonKinVars(*ele);
        if(!multilepAnalyzer->isData) fillLeptonGenVars(*ele, *genParticles);
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
        fillMatchingIVFVariables(*secVertices, *ele, pv);

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
        //if(!tau.tauID("decayModeFinding")) continue;
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
    for(auto array : {&_tauMuonVeto, &_tauEleVeto, &_decayModeFindingNew, &_tauVLooseMvaNew, &_tauLooseMvaNew}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauMediumMvaNew, &_tauTightMvaNew, &_tauVTightMvaNew, &_tauVTightMvaOld}) std::fill_n(*array, _nLight, false);
    for(auto array : {&_tauAgainstElectronMVA6Raw, &_tauCombinedIsoDBRaw3Hits, &_tauIsoMVAPWdR03oldDMwLT}) std::fill_n(*array, _nLight, 0.);
    for(auto array : {&_tauIsoMVADBdR03oldDMwLT, &_tauIsoMVADBdR03newDMwLT, &_tauIsoMVAPWnewDMwLT, &_tauIsoMVAPWoldDMwLT}) std::fill_n(*array, _nLight, 0.);


    if(multilepAnalyzer->skim == "trilep"    &&  _nL     < 3) return false;
    if(multilepAnalyzer->skim == "dilep"     &&  _nLight < 2 /*|| !tightlepton)*/) return false;
    if(multilepAnalyzer->skim == "ttg"       &&  _nLight < 2) return false;
    if(multilepAnalyzer->skim == "singlelep" &&  _nL     < 1) return false;
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
    if(!match or match->pdgId() != lepton.pdgId()) match = GenTools::geometricMatch(lepton, genParticles); // if no match or pdgId is different, try the geometric match

    _lIsPrompt[_nL]             = match and (abs(lepton.pdgId()) == abs(match->pdgId()) || match->pdgId() == 22) and GenTools::isPrompt(*match, genParticles); // only when matched to its own flavor or a photon
    _lMatchPdgId[_nL]           = match ? match->pdgId() : 0;
    _lProvenance[_nL]           = GenTools::provenance(match, genParticles);
    _lProvenanceCompressed[_nL] = GenTools::provenanceCompressed(match, genParticles, _lIsPrompt[_nL]);
    _lProvenanceConversion[_nL] = GenTools::provenanceConversion(match, genParticles);
    _lMomPdgId[_nL]             = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;
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
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon, const reco::Vertex& vertex){
    _dxy[_nL]     = muon.innerTrack()->dxy(vertex.position());
    _dz[_nL]      = muon.innerTrack()->dz(vertex.position());
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


void LeptonAnalyzer::fillMatchingIVFVariables(const std::vector<reco::Vertex>& secVertices, const pat::Muon& muon, const reco::Vertex& pv){
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
                _IVF_trackdxy[_nLight][_IVF_ntracks[_nLight]]       = std::abs(vtxTrack->dxy(pv.position()));
                _IVF_trackdz[_nLight][_IVF_ntracks[_nLight]]        = std::abs(vtxTrack->dz(pv.position()));
                _IVF_ntracks[_nLight]++;
            }
            new_vtx = false;
        }
    }
}

void LeptonAnalyzer::fillMatchingIVFVariables(const std::vector<reco::Vertex>& secVertices, const pat::Electron& electron, const reco::Vertex& pv){
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
                _IVF_trackdxy[_nLight][_IVF_ntracks[_nLight]]       = std::abs(vtxTrack->dxy(pv.position()));
                _IVF_trackdz[_nLight][_IVF_ntracks[_nLight]]        = std::abs(vtxTrack->dz(pv.position()));
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
//void LeptonAnalyzer::fillAllIVFVariables(const std::vector<reco::Vertex>& secVertices, const reco::Vertex& pv){
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
//            _IVF_trackdxy[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]       = std::abs(vtxTrack->dxy(pv.position()));
//            _IVF_trackdz[_IVF_nvertex][_IVF_ntracks[_IVF_nvertex]]        = std::abs(vtxTrack->dz(pv.position()));
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
