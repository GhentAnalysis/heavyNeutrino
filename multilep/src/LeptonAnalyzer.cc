#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenMatching.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TLorentzVector.h"

LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer), 
    electronsEffectiveAreas(multilepAnalyzer->is2017 ? (iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas")).fullPath() ),
    muonsEffectiveAreas    (multilepAnalyzer->is2017 ? (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreasFall17")).fullPath() : (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath() )
{
    leptonMvaComputerSUSY16 = new LeptonMvaHelper(iConfig, 0, false);     //SUSY
    leptonMvaComputerTTH16 = new LeptonMvaHelper(iConfig, 1, false);     //TTH
    leptonMvaComputerSUSY17 = new LeptonMvaHelper(iConfig, 0, true);     //SUSY
    leptonMvaComputerTTH17 = new LeptonMvaHelper(iConfig, 1, true);     //TTH
    leptonMvaComputertZqTTV16 = new LeptonMvaHelper(iConfig, 2, false);  //tZq/TTV
    if(!multilepAnalyzer->isData) genMatcher = new GenMatching(iConfig, multilepAnalyzer);
    if(multilepAnalyzer->isData){
        jecLevel = "L2L3Residual";
    } else {
        jecLevel = "L3Absolute";
    }
};

LeptonAnalyzer::~LeptonAnalyzer(){
    delete leptonMvaComputerSUSY16;
    delete leptonMvaComputerTTH16;
    delete leptonMvaComputertZqTTV16;
    delete leptonMvaComputerSUSY17;
    delete leptonMvaComputerTTH17;
    if(!multilepAnalyzer->isData) delete genMatcher;
}

void LeptonAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nL",                           &_nL,                           "_nL/b");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/b");
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
    outputTree->Branch("_lElectronMva",                 &_lElectronMva,                 "_lElectronMva[_nLight]/F");
    outputTree->Branch("_lElectronMvaHZZ",              &_lElectronMvaHZZ,              "_lElectronMvaHZZ[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17Iso",        &_lElectronMvaFall17Iso,        "_lElectronMvaFall17Iso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17NoIso",      &_lElectronMvaFall17NoIso,      "_lElectronMvaFall17NoIso[_nLight]/F");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_leptonMvaSUSY16",              &_leptonMvaSUSY16,              "_leptonMvaSUSY16[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH16",               &_leptonMvaTTH16,               "_leptonMvaTTH16[_nLight]/D");
    outputTree->Branch("_leptonMvaSUSY17",              &_leptonMvaSUSY17,              "_leptonMvaSUSY17[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH17",               &_leptonMvaTTH17,               "_leptonMvaTTH17[_nLight]/D");
    outputTree->Branch("_leptonMvatZqTTV16",            &_leptonMvatZqTTV16,            "_leptonMvatZqTTV16[_nLight]/D");
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
    outputTree->Branch("_lPOGLooseWOIso",               &_lPOGLooseWOIso,               "_lPOGLooseWOIso[_nLight]/O");
    outputTree->Branch("_lPOGMediumWOIso",              &_lPOGMediumWOIso,              "_lPOGMediumWOIso[_nLight]/O");
    outputTree->Branch("_lPOGTightWOIso",               &_lPOGTightWOIso,               "_lPOGTightWOIso[_nLight]/O");
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
        outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
        outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
        outputTree->Branch("_lProvenanceCompressed",      &_lProvenanceCompressed,        "_lProvenanceCompressed[_nL]/i");
        outputTree->Branch("_lProvenanceConversion",      &_lProvenanceConversion,        "_lProvenanceConversion[_nL]/i");
    }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const reco::Vertex& primaryVertex){
    edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
    edm::Handle<edm::ValueMap<float>> electronsMva;                  iEvent.getByToken(multilepAnalyzer->eleMvaToken,                       electronsMva);
    edm::Handle<edm::ValueMap<float>> electronsMvaHZZ;               iEvent.getByToken(multilepAnalyzer->eleMvaHZZToken,                    electronsMvaHZZ);
    edm::Handle<edm::ValueMap<float>> electronMvaFall17Iso;          iEvent.getByToken(multilepAnalyzer->eleMvaFall17IsoToken,              electronMvaFall17Iso);
    edm::Handle<edm::ValueMap<float>> electronMvaFall17NoIso;        iEvent.getByToken(multilepAnalyzer->eleMvaFall17NoIsoToken,            electronMvaFall17NoIso);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedVeto;          iEvent.getByToken(multilepAnalyzer->eleCutBasedVetoToken,              electronsCutBasedVeto);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedLoose;         iEvent.getByToken(multilepAnalyzer->eleCutBasedLooseToken,             electronsCutBasedLoose);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedMedium;        iEvent.getByToken(multilepAnalyzer->eleCutBasedMediumToken,            electronsCutBasedMedium);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedTight;         iEvent.getByToken(multilepAnalyzer->eleCutBasedTightToken,             electronsCutBasedTight);
    edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
    edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(multilepAnalyzer->tauToken,                          taus);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
    edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);
    edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);
    //edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetSmearedToken,                   jets);

    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;

    //set up generator matching
    if(!multilepAnalyzer->isData) genMatcher->setGenParticles(iEvent);

    //loop over muons
    for(const pat::Muon& mu : *muons){
        if(_nL == nL_max)                              break;
        if(mu.innerTrack().isNull())                   continue;
        if(mu.pt() < 5)                                continue;
        if(fabs(mu.eta()) > 2.4)                       continue;
        if(!mu.isPFMuon())                             continue;
        if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue;
        fillLeptonImpactParameters(mu, primaryVertex);
        if(fabs(_dxy[_nL]) > 0.05)                     continue;
        if(fabs(_dz[_nL]) > 0.1)                       continue;
        fillLeptonKinVars(mu);
        //fillLeptonGenVars(mu.genParticle());
        if(!multilepAnalyzer->isData) fillLeptonGenVars(mu, genMatcher);
        fillLeptonJetVariables(mu, jets, primaryVertex, *rho);

        _lFlavor[_nL]        = 1;
        _lMuonSegComp[_nL]    = mu.segmentCompatibility();
        _lMuonTrackPt[_nL]    = mu.innerTrack()->pt();
        _lMuonTrackPtErr[_nL] = mu.innerTrack()->ptError();

        _relIso[_nL]         = getRelIso03(mu, *rho);                     // Isolation variables
        _relIso0p4[_nL]      = getRelIso04(mu, *rho);                                                     
        _miniIso[_nL]        = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho, false);
        _miniIsoCharged[_nL] = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho, true);

        _lHNLoose[_nL]       = isHNLoose(mu);                                                       // ID variables
        _lHNFO[_nL]          = isHNFO(mu);                                                          // don't change order, they rely on above variables
        _lHNTight[_nL]       = isHNTight(mu);

        _lPOGVeto[_nL]       = mu.isLooseMuon();
        _lPOGLoose[_nL]      = mu.isLooseMuon();
        _lPOGMedium[_nL]     = mu.isMediumMuon();
        _lPOGTight[_nL]      = mu.isTightMuon(primaryVertex);

        _leptonMvaSUSY16[_nL]  = leptonMvaVal(mu, leptonMvaComputerSUSY16);
        _leptonMvaTTH16[_nL]   = leptonMvaVal(mu, leptonMvaComputerTTH16);
        _leptonMvaSUSY17[_nL]  = leptonMvaVal(mu, leptonMvaComputerSUSY17);
        _leptonMvaTTH17[_nL]   = leptonMvaVal(mu, leptonMvaComputerTTH17);
        _leptonMvatZqTTV16[_nL] = leptonMvaVal(mu, leptonMvaComputertZqTTV16);

        _lEwkLoose[_nL]      = isEwkLoose(mu);
        _lEwkFO[_nL]         = isEwkFO(mu);
        _lEwkTight[_nL]      = isEwkTight(mu);

        ++_nMu;
        ++_nL;
        ++_nLight;
    }

    // Loop over electrons (note: using iterator we can easily get the ref too)
    for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
        auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele - electrons->begin()));
        if(_nL == nL_max)                                                                               break;
        if(ele->gsfTrack().isNull())                                                                    continue;
        if(ele->pt() < 7)                                                                               continue;
        if(fabs(ele->eta()) > 2.5)                                                                      continue;
        if(ele->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 2)    continue;
        fillLeptonImpactParameters(*ele, primaryVertex);
        if(fabs(_dxy[_nL]) > 0.05)                                                                      continue;
        if(fabs(_dz[_nL]) > 0.1)                                                                        continue;
        fillLeptonKinVars(*ele);
        //fillLeptonGenVars(ele->genParticle());
        if(!multilepAnalyzer->isData) fillLeptonGenVars(*ele, genMatcher);
        fillLeptonJetVariables(*ele, jets, primaryVertex, *rho);

        _lFlavor[_nL]          = 0;
        _lEtaSC[_nL]           = ele->superCluster()->eta();

        _relIso[_nL]                    = getRelIso03(*ele, *rho);
        _relIso0p4[_nL]                 = getRelIso(*ele, packedCands, 0.4, *rho, false);
        _miniIso[_nL]                   = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho, false);
        _miniIsoCharged[_nL]            = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho, true);
        _lElectronMva[_nL]              = (*electronsMva)[electronRef];
        _lElectronMvaHZZ[_nL]           = (*electronsMvaHZZ)[electronRef];
        _lElectronMvaFall17Iso[_nL]     = (*electronMvaFall17Iso)[electronRef];
        _lElectronMvaFall17NoIso[_nL]   = (*electronMvaFall17NoIso)[electronRef];
        _lElectronPassEmu[_nL]          = passTriggerEmulationDoubleEG(&*ele);                             // Keep in mind, this trigger emulation is for 2016 DoubleEG, the SingleEG trigger emulation is different
        _lElectronPassConvVeto[_nL]     = ele->passConversionVeto();
        _lElectronChargeConst[_nL]      = ele->isGsfCtfScPixChargeConsistent();
        _lElectronMissingHits[_nL]      = ele->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

        _lHNLoose[_nL]                  = isHNLoose(*ele);
        _lHNFO[_nL]                     = isHNFO(*ele);
        _lHNTight[_nL]                  = isHNTight(*ele);

        _lPOGVeto[_nL]                  = (*electronsCutBasedVeto)[electronRef];
        _lPOGLoose[_nL]                 = (*electronsCutBasedLoose)[electronRef];
        _lPOGMedium[_nL]                = (*electronsCutBasedMedium)[electronRef];
        _lPOGTight[_nL]                 = (*electronsCutBasedTight)[electronRef];

        _lPOGLooseWOIso[_nL]            = isLooseCutBasedElectronWithoutIsolation(&*ele);
        _lPOGMediumWOIso[_nL]           = isMediumCutBasedElectronWithoutIsolation(&*ele);
        _lPOGTightWOIso[_nL]            = isTightCutBasedElectronWithoutIsolation(&*ele);

        _leptonMvaSUSY16[_nL]           = leptonMvaVal(*ele, leptonMvaComputerSUSY16);
        _leptonMvaTTH16[_nL]            = leptonMvaVal(*ele, leptonMvaComputerTTH16);
        _leptonMvaSUSY17[_nL]           = leptonMvaVal(*ele, leptonMvaComputerSUSY17);
        _leptonMvaTTH17[_nL]            = leptonMvaVal(*ele, leptonMvaComputerTTH17);
        _leptonMvatZqTTV16[_nL]           = leptonMvaVal(*ele, leptonMvaComputertZqTTV16);

        _lEwkLoose[_nL]                 = isEwkLoose(*ele);
        _lEwkFO[_nL]                    = isEwkFO(*ele);
        _lEwkTight[_nL]                 = isEwkTight(*ele);

        ++_nEle;
        ++_nL;
        ++_nLight;
    }

    //loop over taus
    for(const pat::Tau& tau : *taus){
        if(_nL == nL_max)         break;
        if(tau.pt() < 20)         continue;          // Minimum pt for tau reconstruction
        if(fabs(tau.eta()) > 2.3) continue;
        //if(!tau.tauID("decayModeFinding")) continue;
        fillLeptonKinVars(tau);
        //fillLeptonGenVars(tau.genParticle());
        //if(!multilepAnalyzer->isData) fillLeptonGenVars(tau, genMatcher);
        fillLeptonImpactParameters(tau, primaryVertex);
        if(_dz[_nL] < 0.4)        continue;         //tau dz cut used in ewkino

        _lFlavor[_nL]  = 2;
        _tauMuonVeto[_nL] = tau.tauID("againstMuonLoose3");                                        //Light lepton vetos
        _tauEleVeto[_nL] = tau.tauID("againstElectronLooseMVA6");

        _lPOGVeto[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");                        //old tau ID
        _lPOGLoose[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
        _lPOGMedium[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
        _lPOGTight[_nL] = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
        _tauVTightMvaOld[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

        _decayModeFindingNew[_nL] = tau.tauID("decayModeFindingNewDMs");                           //new Tau ID 
        _tauVLooseMvaNew[_nL] = tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
        _tauLooseMvaNew[_nL] = tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
        _tauMediumMvaNew[_nL] = tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
        _tauTightMvaNew[_nL] = tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT");
        _tauVTightMvaNew[_nL] = tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT");

        _tauAgainstElectronMVA6Raw[_nL] = tau.tauID("againstElectronMVA6Raw");
        _tauCombinedIsoDBRaw3Hits[_nL] = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        _tauIsoMVAPWdR03oldDMwLT[_nL] = tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw");
        _tauIsoMVADBdR03oldDMwLT[_nL] = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        _tauIsoMVADBdR03newDMwLT[_nL] = tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        _tauIsoMVAPWnewDMwLT[_nL] = tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        _tauIsoMVAPWoldDMwLT[_nL] = tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw"); 

        _lEwkLoose[_nL] = isEwkLoose(tau);
        _lEwkFO[_nL]    = isEwkFO(tau);
        _lEwkTight[_nL] = isEwkTight(tau);
        ++_nTau;
        ++_nL;
    }

    if(multilepAnalyzer->skim == "trilep"    &&  _nL     < 3) return false;
    if(multilepAnalyzer->skim == "dilep"     &&  _nL     < 2) return false;
    if(multilepAnalyzer->skim == "ttg"       &&  _nLight < 2) return false;
    if(multilepAnalyzer->skim == "singlelep" &&  _nL     < 1) return false;
    if(multilepAnalyzer->skim == "FR" &&  _nLight < 1) return false;
    return true;
}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton){
    _lPt[_nL]     = lepton.pt();
    _lEta[_nL]    = lepton.eta();
    _lPhi[_nL]    = lepton.phi();
    _lE[_nL]      = lepton.energy();
    _lCharge[_nL] = lepton.charge();
}

template <typename Lepton> void LeptonAnalyzer::fillLeptonGenVars(const Lepton& lepton, GenMatching* genMatcher){
    genMatcher->fillMatchingVars(lepton);
    _lIsPrompt[_nL] = genMatcher->promptMatch();
    _lMatchPdgId[_nL] = genMatcher->pdgIdMatch();
    _lProvenance[_nL] = genMatcher->getProvenance();
    _lProvenanceCompressed[_nL] = genMatcher->getProvenanceCompressed();
    _lProvenanceConversion[_nL] = genMatcher->getProvenanceConversion();
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

//Check if electron overlaps with loose muon
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

//Check if tau overlaps with light lepton
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


void LeptonAnalyzer::fillLeptonJetVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::Jet>>& jets, const reco::Vertex& vertex, const double rho){
    //Make skimmed "close jet" collection
    std::vector<pat::Jet> selectedJetsAll;
    for(auto jet = jets->cbegin(); jet != jets->cend(); ++jet){
        //double jetPt = jet->pt()*multilepAnalyzer->jec->jetCorrection(jet->correctedP4("Uncorrected").Pt(), jet->correctedP4("Uncorrected").Eta(), rho, jet->jetArea(), jecLevel); 
        //if( jetPt > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
        if( jet->pt() > 5 && fabs( jet->eta() ) < 3) selectedJetsAll.push_back(*jet);
    }
    // Find closest selected jet
    unsigned closestIndex = 0;
    for(unsigned j = 1; j < selectedJetsAll.size(); ++j){
        if(reco::deltaR(selectedJetsAll[j], lepton) < reco::deltaR(selectedJetsAll[closestIndex], lepton)) closestIndex = j;
    }
    const pat::Jet& jet = selectedJetsAll[closestIndex];
    if(selectedJetsAll.size() == 0 || reco::deltaR(jet, lepton) > 0.4){ //Now includes safeguard for 0 jet events
        _ptRatio[_nL] = 1;
        _ptRel[_nL] = 0;
        _closestJetCsvV2[_nL] = 0;
        _closestJetDeepCsv_b[_nL] = 0;
        _closestJetDeepCsv_bb[_nL] = 0;
        _selectedTrackMult[_nL] = 0;
    } else {
        /*
        double totalJEC = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), rho, jet.jetArea(), jecLevel);
        double l1JEC = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), rho, jet.jetArea(), "L1FastJet");
        TLorentzVector l1Jet;
        l1Jet.SetPtEtaPhiE(jet.correctedP4("Uncorrected").Pt()*l1JEC, jet.correctedP4("Uncorrected").Eta(), jet.correctedP4("Uncorrected").Phi(), jet.correctedP4("Uncorrected").E()*l1JEC);
        TLorentzVector l(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());
        TLorentzVector lepAwareJet = (l1Jet - l)*JEC + l;
        float JEC = totalJEC/l1JEC;
        */
        auto  l1Jet       = jet.correctedP4("L1FastJet");
        float JEC         = jet.p4().E()/l1Jet.E();
        auto  l           = lepton.p4();
        auto  lepAwareJet = (l1Jet - l)*JEC + l;


        TLorentzVector lV(l.Px(), l.Py(), l.Pz(), l.E());
        TLorentzVector jV(lepAwareJet.Px(), lepAwareJet.Py(), lepAwareJet.Pz(), lepAwareJet.E());
        _ptRatio[_nL]       = l.Pt()/lepAwareJet.Pt();
        _ptRel[_nL]         = lV.Perp((jV - lV).Vect());
        _closestJetCsvV2[_nL] = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        _closestJetDeepCsv_b[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        //compute selected track multiplicity of closest jet
        _selectedTrackMult[_nL] = 0;
        for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
            const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
            if(daughter->hasTrackDetails()){
                const reco::Track& daughterTrack = daughter->pseudoTrack();
                TLorentzVector trackVec(daughterTrack.px(), daughterTrack.py(), daughterTrack.pz(), daughterTrack.p());
                double daughterDeltaR            = trackVec.DeltaR(jV);
                bool goodTrack                   = daughterTrack.pt() > 1 && daughterTrack.charge() != 0 && daughterTrack.hitPattern().numberOfValidHits() > 7
                    && daughterTrack.hitPattern().numberOfValidPixelHits() > 1 && daughterTrack.normalizedChi2() < 5 && fabs(daughterTrack.dz(vertex.position())) < 17
                    && fabs(daughterTrack.dxy(vertex.position())) < 17;
                if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack) ++_selectedTrackMult[_nL];
            }
        }
    }
}
