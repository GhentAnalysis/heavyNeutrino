#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
//#include "heavyNeutrino/multilep/interface/GenMatching.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TLorentzVector.h"

LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    electronsEffectiveAreas((iConfig.getParameter<edm::FileInPath>("electronsEffectiveAreas")).fullPath()),
    muonsEffectiveAreas(    (iConfig.getParameter<edm::FileInPath>("muonsEffectiveAreas")).fullPath()),
    multilepAnalyzer(multilepAnalyzer)
{
    leptonMvaComputerSUSY = new LeptonMvaHelper(iConfig);           //SUSY
    leptonMvaComputerTTH = new LeptonMvaHelper(iConfig, false);     //TTH
    //if(!multilepAnalyzer->isData) genMatcher = new GenMatching(iConfig, multilepAnalyzer);
};

LeptonAnalyzer::~LeptonAnalyzer(){
    delete leptonMvaComputerSUSY;
    delete leptonMvaComputerTTH;
    //if(!multilepAnalyzer->isData) delete genMatcher;
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
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_leptonMvaSUSY",                &_leptonMvaSUSY,                "_leptonMvaSUSY[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH",                 &_leptonMvaTTH,                 "_leptonMvaTTH[_nLight]/D");
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
    outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
    outputTree->Branch("_relIso0p4Mu",                  &_relIso0p4Mu,                  "_relIso0p4Mu[_nMu]/D");
    outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
    outputTree->Branch("_miniIsoCharged",               &_miniIsoCharged,               "_miniIsoCharged[_nLight]/D");
    outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
    outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
    outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_TrackMult_pt0",        	&_TrackMult_pt0,        	"_TrackMult_pt0[_nLight]/i");
    outputTree->Branch("_TrackMult_pt1",        	&_TrackMult_pt1,        	"_TrackMult_pt1[_nL]/i");
    outputTree->Branch("_TrackMult_pt2",        	&_TrackMult_pt2,        	"_TrackMult_pt2[_nLight]/i");
    outputTree->Branch("_TrackMult_pt3",        	&_TrackMult_pt3,        	"_TrackMult_pt3[_nLight]/i");
    outputTree->Branch("_TrackMult_pt4",        	&_TrackMult_pt4,        	"_TrackMult_pt4[_nLight]/i");
    outputTree->Branch("_TrackMult_pt5",        	&_TrackMult_pt5,        	"_TrackMult_pt5[_nL]/i");
    outputTree->Branch("_TrackMult_noIP_pt0",        	&_TrackMult_noIP_pt0,        	"_TrackMult_noIP_pt0[_nLight]/i");
    outputTree->Branch("_TrackMult_noIP_pt1",        	&_TrackMult_noIP_pt1,        	"_TrackMult_noIP_pt1[_nLight]/i");
    outputTree->Branch("_TrackMult_noIP_pt2",        	&_TrackMult_noIP_pt2,        	"_TrackMult_noIP_pt2[_nLight]/i");
    outputTree->Branch("_TrackMult_noIP_pt3",        	&_TrackMult_noIP_pt3,        	"_TrackMult_noIP_pt3[_nLight]/i");
    outputTree->Branch("_TrackMult_noIP_pt4",        	&_TrackMult_noIP_pt4,        	"_TrackMult_noIP_pt4[_nLight]/i");
    outputTree->Branch("_TrackMult_noIP_pt5",        	&_TrackMult_noIP_pt5,        	"_TrackMult_noIP_pt5[_nLight]/i");
    outputTree->Branch("_Nutau_TrackMult_pt1",      	&_Nutau_TrackMult_pt1,		"_Nutau_TrackMult_pt1[_nL]/i");
    outputTree->Branch("_Nutau_TrackMult_pt5",      	&_Nutau_TrackMult_pt5,		"_Nutau_TrackMult_pt5[_nL]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");  
    if(!multilepAnalyzer->isData){
        outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
        //outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
    }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const reco::Vertex& primaryVertex){
    edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
    edm::Handle<edm::ValueMap<float>> electronsMva;                  iEvent.getByToken(multilepAnalyzer->eleMvaToken,                       electronsMva);
    edm::Handle<edm::ValueMap<float>> electronsMvaHZZ;               iEvent.getByToken(multilepAnalyzer->eleMvaHZZToken,                    electronsMvaHZZ);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedVeto;          iEvent.getByToken(multilepAnalyzer->eleCutBasedVetoToken,              electronsCutBasedVeto);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedLoose;         iEvent.getByToken(multilepAnalyzer->eleCutBasedLooseToken,             electronsCutBasedLoose);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedMedium;        iEvent.getByToken(multilepAnalyzer->eleCutBasedMediumToken,            electronsCutBasedMedium);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedTight;         iEvent.getByToken(multilepAnalyzer->eleCutBasedTightToken,             electronsCutBasedTight);
    edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
    edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(multilepAnalyzer->tauToken,                          taus);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
    edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);
    edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);
    edm::Handle<std::vector<reco::GenParticle>> genParticles; 	     iEvent.getByToken(multilepAnalyzer->genParticleToken, 		    genParticles);

    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;

    //set up generator matching
    //if(!multilepAnalyzer->isData) genMatcher->setGenParticles(iEvent);

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
        fillLeptonGenVars(mu.genParticle());
        //if(!multilepAnalyzer->isData) fillLeptonGenVars(mu, genMatcher);
	/*std::cout << std::endl << "Muon pt, phi and eta: " << mu.pt() << ", " << mu.phi() << ", " << mu.eta() << std::endl;
	if(mu.innerTrack().isNonnull()){
	    const reco::Track* innertrack = mu.innerTrack().get();
	    std::cout << "Muon innertrack pt, phi and eta: " << innertrack->pt() << ", " << innertrack->phi() << ", " << innertrack->eta() << std::endl;
	}
	if(mu.outerTrack().isNonnull()){
	    const reco::Track* outertrack = mu.outerTrack().get();
	    std::cout << "Muon outertrack pt, phi and eta: " << outertrack->pt() << ", " << outertrack->phi() << ", " << outertrack->eta() << std::endl;
	}
	if(mu.globalTrack().isNonnull()){
	    const reco::Track* globaltrack = mu.globalTrack().get();
	    std::cout << "Muon globaltrack pt, phi and eta: " << globaltrack->pt() << ", " << globaltrack->phi() << ", " << globaltrack->eta() << std::endl;
	}*/
        fillLeptonJetVariables(mu, jets, primaryVertex);
	fillLeptonTrackVariables(mu, packedCands, primaryVertex);

        _lFlavor[_nL]        = 1;
        _lMuonSegComp[_nL]    = mu.segmentCompatibility();
        _lMuonTrackPt[_nL]    = mu.innerTrack()->pt();
        _lMuonTrackPtErr[_nL] = mu.innerTrack()->ptError();

        _relIso[_nL]         = getRelIso03(mu, *rho);                                               // Isolation variables
        _relIso0p4Mu[_nL]    = getRelIso04(mu);                                                     
        _miniIso[_nL]        = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho);
        _miniIsoCharged[_nL] = getMiniIsolation(mu, packedCands, 0.05, 0.2, 10, *rho, true);

        _lHNLoose[_nL]       = isHNLoose(mu);                                                       // ID variables
        _lHNFO[_nL]          = isHNFO(mu);                                                          // don't change order, they rely on above variables
        _lHNTight[_nL]       = isHNTight(mu);

        _lPOGVeto[_nL]       = mu.isLooseMuon();
        _lPOGLoose[_nL]      = mu.isLooseMuon();
        _lPOGMedium[_nL]     = mu.isMediumMuon();
        _lPOGTight[_nL]      = mu.isTightMuon(primaryVertex);

        _leptonMvaSUSY[_nL]  = leptonMvaVal(mu, leptonMvaComputerSUSY);
        _leptonMvaTTH[_nL]   = leptonMvaVal(mu, leptonMvaComputerTTH);

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
        if(_nL == nL_max)                                                                        break;
        if(ele->gsfTrack().isNull())                                                             continue;
        if(ele->pt() < 7)                                                                        continue;
        if(fabs(ele->eta()) > 2.5)                                                               continue;
        if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2) continue;
        fillLeptonImpactParameters(*ele, primaryVertex);
        if(fabs(_dxy[_nL]) > 0.05)                                                               continue;
        if(fabs(_dz[_nL]) > 0.1)                                                                 continue;
        fillLeptonKinVars(*ele);
        fillLeptonGenVars(ele->genParticle());
        //if(!multilepAnalyzer->isData) fillLeptonGenVars(*ele, genMatcher);
	//std::cout << std::endl << "electron pt, phi and eta: " << ele->pt() << ", " << ele->phi() << ", " << ele->eta() << std::endl;
        fillLeptonJetVariables(*ele, jets, primaryVertex);
	fillLeptonTrackVariables(*ele, packedCands, primaryVertex);

        _lFlavor[_nL]          = 0;
        _lEtaSC[_nL]           = ele->superCluster()->eta();

        _relIso[_nL]                = getRelIso03(*ele, *rho);
        _miniIso[_nL]               = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho);
        _miniIsoCharged[_nL]        = getMiniIsolation(*ele, packedCands, 0.05, 0.2, 10, *rho, true);
        _lElectronMva[_nL]          = (*electronsMva)[electronRef];
        _lElectronMvaHZZ[_nL]       = (*electronsMvaHZZ)[electronRef];
        _lElectronPassEmu[_nL]      = passTriggerEmulationDoubleEG(&*ele);                             // Keep in mind, this trigger emulation is for 2016 DoubleEG, the SingleEG trigger emulation is different
        _lElectronPassConvVeto[_nL] = ele->passConversionVeto();
        _lElectronChargeConst[_nL]  = ele->isGsfCtfScPixChargeConsistent();
        _lElectronMissingHits[_nL]  = ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

        _lHNLoose[_nL]              = isHNLoose(*ele);
        _lHNFO[_nL]                 = isHNFO(*ele);
        _lHNTight[_nL]              = isHNTight(*ele);

        _lPOGVeto[_nL]              = (*electronsCutBasedVeto)[electronRef];
        _lPOGLoose[_nL]             = (*electronsCutBasedLoose)[electronRef];
        _lPOGMedium[_nL]            = (*electronsCutBasedMedium)[electronRef];
        _lPOGTight[_nL]             = (*electronsCutBasedTight)[electronRef];

        _lPOGLooseWOIso[_nL]        = isLooseCutBasedElectronWithoutIsolation(&*ele);
        _lPOGMediumWOIso[_nL]       = isMediumCutBasedElectronWithoutIsolation(&*ele);
        _lPOGTightWOIso[_nL]        = isTightCutBasedElectronWithoutIsolation(&*ele);

        _leptonMvaSUSY[_nL]         = leptonMvaVal(*ele, leptonMvaComputerSUSY);
        _leptonMvaTTH[_nL]          = leptonMvaVal(*ele, leptonMvaComputerTTH);

        _lEwkLoose[_nL]             = isEwkLoose(*ele);
        _lEwkFO[_nL]                = isEwkFO(*ele);
        _lEwkTight[_nL]             = isEwkTight(*ele);

        ++_nEle;
        ++_nL;
        ++_nLight;
    }

    //loop over taus
    for(const pat::Tau& tau : *taus){
        if(_nL == nL_max)         break;
        if(tau.pt() < 20)         continue;          // Minimum pt for tau reconstruction
        if(fabs(tau.eta()) > 2.3) continue;
        if(!tau.tauID("decayModeFinding")) continue;
        fillLeptonKinVars(tau);
        fillLeptonGenVars(tau.genParticle());
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

        _lEwkLoose[_nL] = isEwkLoose(tau);
        _lEwkFO[_nL]    = isEwkFO(tau);
        _lEwkTight[_nL] = isEwkTight(tau);

	// todo: find nutau and look around that instead of the tau!!!
	fillLeptonJetVariables(tau, jets, primaryVertex);
	fillLeptonTrackVariables(tau, packedCands, primaryVertex);

	//find nutau and get trackmult around its direction
	//match gen tau with tau and then look at the neutrino daughter
	//fillNutauTrackVariables(tau.genParticle(), packedCands, primaryVertex);
	/*std::cout << std::endl << "Tau pt, eta, phi: " << tau.pt() << ", " << tau.eta() << ", " << tau.phi() << std::endl;
	auto lep = tau.p4();
        TLorentzVector lepVec(lep.px(), lep.py(), lep.pz(), lep.E());
	const reco::Candidate* leadchargedcand = tau.leadChargedHadrCand().get();
	std::cout << "lead charged Hadron cand pt, eta, phi and dR: " << leadchargedcand->pt() << ", " << leadchargedcand->eta() << ", " << leadchargedcand->phi() << std::endl;
	const reco::Track& candTrack = *(leadchargedcand->bestTrack());
	TLorentzVector trackVec(candTrack.px(), candTrack.py(), candTrack.pz(), candTrack.p());
	double leptondeltaR = trackVec.DeltaR(lepVec);
	bool goodTrack		= leptondeltaR < 0.4 &&
				  candTrack.normalizedChi2() < 5 &&
				  candTrack.charge() != 0;
	bool goodIP		= candTrack.hitPattern().numberOfValidHits() > 7 &&
				  candTrack.hitPattern().numberOfValidPixelHits() > 1 &&
				  fabs(candTrack.dz(primaryVertex.position())) < 17 &&
				  fabs(candTrack.dxy(primaryVertex.position())) < 17;
	bool is_lepton		= fabs(lepVec.Pt() - trackVec.Pt()) < 0.1 &&
				  leptondeltaR < 0.05;
	std::cout << "goodtrack, goodIP, is_lepton: " << goodTrack << " " << goodIP << " " << is_lepton << std::endl;
	*/ //std::cout << "# of genparticles: " << genParticles->size() << std::endl;
	for(const reco::GenParticle& p : *genParticles){
	    //std::cout << "gen id: " << p.pdgId() << std::endl;
	    //std::cout << "gen pt: " << p.pt() << std::endl;
	    if(abs(p.pdgId()) == 15){
		//std::cout << "gen tau pt, phi, eta: " << p.pt() << ", " << p.phi() << ", " << p.eta() << std::endl;
	    	for(unsigned i = 0; i < p.numberOfDaughters(); ++i){
		    //std::cout << "gen Tau daughter pdgid: " << (*genParticles)[p.daughterRef(i).key()].pdgId() << std::endl;
		    if(abs((*genParticles)[p.daughterRef(i).key()].pdgId()) == 16){
			const reco::GenParticle& Nutau = (*genParticles)[p.daughterRef(i).key()];
			fillNutauTrackVariables(Nutau, packedCands, primaryVertex);
			//std::cout << "Nutau found: " << Nutau.pdgId() << std::endl;
			//std::cout << "Nutau pt, phi, eta: " << Nutau.pt() << ", " << Nutau.phi() << ", " << Nutau.eta() << std::endl;
		    }
		}
	    }
	}

	/*const reco::GenParticle* gen_tau = tau.genParticle();
	std::cout << "pointer: " << gen_tau << std::endl;
	if(gen_tau) std::cout << "gen tau id: " << gen_tau->pdgId() << std::endl;*/

        ++_nTau;
        ++_nL;
    }

    if(multilepAnalyzer->skim == "trilep"    &&  _nL     < 3) return false;
    if(multilepAnalyzer->skim == "dilep"     &&  _nLight < 2) return false;
    if(multilepAnalyzer->skim == "ttg"       &&  _nLight < 2) return false;
    if(multilepAnalyzer->skim == "singlelep" &&  _nLight < 1) return false;
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


 void LeptonAnalyzer::fillLeptonGenVars(const reco::GenParticle* genParticle){
   if(genParticle != nullptr){
   _lIsPrompt[_nL] = (genParticle)->isPromptFinalState();
   _lMatchPdgId[_nL] = (genParticle)->pdgId();
   } else{
   _lIsPrompt[_nL] = false;
   _lMatchPdgId[_nL] = 0;
   }
 }
   

/*void LeptonAnalyzer::fillLeptonGenVars(const reco::Candidate& lepton, GenMatching* genMatcher){
    genMatcher->fillMatchingVars(lepton);
    _lIsPrompt[_nL] = genMatcher->promptMatch();
    _lMatchPdgId[_nL] = genMatcher->pdgIdMatch();
    _lProvenance[_nL] = genMatcher->getProvenance();
}*/


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
double LeptonAnalyzer::tau_dz(const pat::Tau& tau, const reco::Vertex::Point& vertex){
    const reco::Candidate::Point& tauVtx = tau.leadChargedHadrCand()->vertex();
    return (tauVtx.Z() - vertex.z()) - ((tauVtx.X() - vertex.x())*tau.px()+(tauVtx.Y()-vertex.y())*tau.py())/tau.pt()*tau.pz()/tau.pt();
}

//Check if electron overlaps with loose muon
bool LeptonAnalyzer::eleMuOverlap(const pat::Electron& ele, const bool* loose){
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
bool LeptonAnalyzer::tauLightOverlap(const pat::Tau& tau, const bool* loose){
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


void LeptonAnalyzer::fillLeptonJetVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::Jet>>& jets, const reco::Vertex& vertex){
    //Make skimmed "close jet" collection
    std::vector<pat::Jet> selectedJetsAll;
    for(auto jet = jets->cbegin(); jet != jets->cend(); ++jet){
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
        /*_TrackMult_pt0[_nL] = 0;
        _TrackMult_pt1[_nL] = 0;
        _TrackMult_pt2[_nL] = 0;
        _TrackMult_pt3[_nL] = 0;
        _TrackMult_pt4[_nL] = 0;
        _TrackMult_pt5[_nL] = 0;
        _TrackMult_noIP_pt0[_nL] = 0;
        _TrackMult_noIP_pt1[_nL] = 0;
        _TrackMult_noIP_pt2[_nL] = 0;
        _TrackMult_noIP_pt3[_nL] = 0;
        _TrackMult_noIP_pt4[_nL] = 0;
        _TrackMult_noIP_pt5[_nL] = 0;*/
    } else {
        auto  l1Jet       = jet.correctedP4("L1FastJet");
        float JEC         = jet.p4().E()/l1Jet.E();
        auto  l           = lepton.p4();
        auto  lepAwareJet = (l1Jet - l)*JEC + l;
        TLorentzVector lV(l.px(), l.py(), l.pz(), l.E());
        TLorentzVector jV(lepAwareJet.px(), lepAwareJet.py(), lepAwareJet.pz(), lepAwareJet.E());
        _ptRatio[_nL]       = l.pt()/lepAwareJet.pt();
        _ptRel[_nL]         = lV.Perp((jV - lV).Vect());
        _closestJetCsvV2[_nL] = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        _closestJetDeepCsv_b[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL] = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        //compute selected track multiplicity of closest jet
        _selectedTrackMult[_nL] = 0;
        /*_TrackMult_pt0[_nL] = 0;
        _TrackMult_pt1[_nL] = 0;
        _TrackMult_pt2[_nL] = 0;
        _TrackMult_pt3[_nL] = 0;
        _TrackMult_pt4[_nL] = 0;
        _TrackMult_pt5[_nL] = 0;
        _TrackMult_noIP_pt0[_nL] = 0;
        _TrackMult_noIP_pt1[_nL] = 0;
        _TrackMult_noIP_pt2[_nL] = 0;
        _TrackMult_noIP_pt3[_nL] = 0;
        _TrackMult_noIP_pt4[_nL] = 0;
        _TrackMult_noIP_pt5[_nL] = 0;*/
	//std::cout << "jet daughters:" << jet.numberOfDaughters() << std::endl;
        for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
            const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
            try {                                                                                                     // In principle, from CMSSW_9_X you need to use if(daughter->hasTrackDetails()){ here, bus that function does not exist in CMSSW_8_X
                const reco::Track& daughterTrack = daughter->pseudoTrack();                                             // Using try {} catch (...){} the code compiles in both versions
                TLorentzVector trackVec(daughterTrack.px(), daughterTrack.py(), daughterTrack.pz(), daughterTrack.p());
                double daughterDeltaR            = trackVec.DeltaR(jV);
                bool goodTrack                   = daughterTrack.pt() > 1 && daughterTrack.charge() != 0 && daughterTrack.hitPattern().numberOfValidHits() > 7
                    && daughterTrack.hitPattern().numberOfValidPixelHits() > 1 && daughterTrack.normalizedChi2() < 5 && fabs(daughterTrack.dz(vertex.position())) < 17
                    && fabs(daughterTrack.dxy(vertex.position())) < 17;
                if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack) ++_selectedTrackMult[_nL];
		// track mult as a function of pt:
		/*bool goodTrack_nopt              = daughterTrack.charge() != 0 && daughterTrack.hitPattern().numberOfValidHits() > 7
                                                   && daughterTrack.hitPattern().numberOfValidPixelHits() > 1 && daughterTrack.normalizedChi2() < 5 && fabs(daughterTrack.dz(vertex.position())) < 17
                                                   && fabs(daughterTrack.dxy(vertex.position())) < 17;
		if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack_nopt){
		    //track multiplicity for several pt thresholds:
		    ++_TrackMult_pt0[_nL];
		    if(daughterTrack.pt() > 1) ++_TrackMult_pt1[_nL];
		    if(daughterTrack.pt() > 2) ++_TrackMult_pt2[_nL];
		    if(daughterTrack.pt() > 3) ++_TrackMult_pt3[_nL];
		    if(daughterTrack.pt() > 4) ++_TrackMult_pt4[_nL];
		    if(daughterTrack.pt() > 5) ++_TrackMult_pt5[_nL];
		}
		// track mult without ip requirements:
		bool goodTrack_nopt_noIP         = daughterTrack.charge() != 0 && daughterTrack.normalizedChi2() < 5;
		if(daughterDeltaR < 0.4 && daughter->fromPV() > 1 && goodTrack_nopt_noIP){
		    //track multiplicity for several pt thresholds:
		    ++_TrackMult_noIP_pt0[_nL];
		    if(daughterTrack.pt() > 1) ++_TrackMult_noIP_pt1[_nL];
		    if(daughterTrack.pt() > 2) ++_TrackMult_noIP_pt2[_nL];
		    if(daughterTrack.pt() > 3) ++_TrackMult_noIP_pt3[_nL];
		    if(daughterTrack.pt() > 4) ++_TrackMult_noIP_pt4[_nL];
		    if(daughterTrack.pt() > 5) ++_TrackMult_noIP_pt5[_nL];
		}*/
            } catch (...){}
        }
    }
}

void LeptonAnalyzer::fillLeptonTrackVariables(const reco::Candidate& lepton, edm::Handle<std::vector<pat::PackedCandidate>>& packedCands, const reco::Vertex& vertex){
    //packedCandidate collection
    std::vector<pat::PackedCandidate> packedCandidates;
    for(auto cand = packedCands->cbegin(); cand != packedCands->cend(); ++cand){
        packedCandidates.push_back(*cand);
    }
    if(packedCandidates.size() == 0){
        _TrackMult_pt0[_nL] = 0;
        _TrackMult_pt1[_nL] = 0;
        _TrackMult_pt2[_nL] = 0;
        _TrackMult_pt3[_nL] = 0;
        _TrackMult_pt4[_nL] = 0;
	_TrackMult_pt5[_nL] = 0;
        _TrackMult_noIP_pt0[_nL] = 0;
        _TrackMult_noIP_pt1[_nL] = 0;
        _TrackMult_noIP_pt2[_nL] = 0;
        _TrackMult_noIP_pt3[_nL] = 0;
        _TrackMult_noIP_pt4[_nL] = 0;
	_TrackMult_noIP_pt5[_nL] = 0;
    } else {
        _TrackMult_pt0[_nL] = 0;
        _TrackMult_pt1[_nL] = 0;
        _TrackMult_pt2[_nL] = 0;
        _TrackMult_pt3[_nL] = 0;
        _TrackMult_pt4[_nL] = 0;
	_TrackMult_pt5[_nL] = 0;
        _TrackMult_noIP_pt0[_nL] = 0;
        _TrackMult_noIP_pt1[_nL] = 0;
        _TrackMult_noIP_pt2[_nL] = 0;
        _TrackMult_noIP_pt3[_nL] = 0;
        _TrackMult_noIP_pt4[_nL] = 0;
	_TrackMult_noIP_pt5[_nL] = 0;
	auto lep = lepton.p4();
        TLorentzVector lepVec(lep.px(), lep.py(), lep.pz(), lep.E());
	//std::cout << std::endl << "lep pt, phi and eta: " << lep.Pt() << ", " << lep.Phi() << ", " << lep.Eta() << std::endl;
	//std::cout << std::endl << "lepVec pt, phi and eta: " << lepVec.Pt() << ", " << lepVec.Phi() << ", " << lepVec.Eta() << std::endl;
	//std::cout << "packed Candidates: " << packedCandidates.size() << std::endl;
 	for(unsigned j = 0; j < packedCandidates.size(); ++j){
	    try {//do same try catch thing as in fillLeptonJetVariables just to be safe
		const reco::Track& candTrack = packedCandidates[j].pseudoTrack();
		TLorentzVector trackVec(candTrack.px(), candTrack.py(), candTrack.pz(), candTrack.p());
		double leptondeltaR = trackVec.DeltaR(lepVec);
		bool goodTrack		= leptondeltaR < 0.4 &&
					  candTrack.normalizedChi2() < 5 &&
					  candTrack.charge() != 0;
		bool goodIP		= candTrack.hitPattern().numberOfValidHits() > 7 &&
					  candTrack.hitPattern().numberOfValidPixelHits() > 1 &&
					  fabs(candTrack.dz(vertex.position())) < 17 &&
					  fabs(candTrack.dxy(vertex.position())) < 17;
					  //packedCandidates[j].fromPV() > 1;
		bool is_lepton		= fabs(lepVec.Pt() - trackVec.Pt()) < 0.1 &&
					  leptondeltaR < 0.05;
		if(!(lepton.isElectron() || lepton.isMuon())) is_lepton = false;
		if(goodTrack && goodIP && !is_lepton && candTrack.pt() > 0) ++_TrackMult_pt0[_nL];
		if(goodTrack && goodIP && !is_lepton && candTrack.pt() > 1) ++_TrackMult_pt1[_nL];
		if(goodTrack && goodIP && !is_lepton && candTrack.pt() > 2) ++_TrackMult_pt2[_nL];
		if(goodTrack && goodIP && !is_lepton && candTrack.pt() > 3) ++_TrackMult_pt3[_nL];
		if(goodTrack && goodIP && !is_lepton && candTrack.pt() > 4) ++_TrackMult_pt4[_nL];
		if(goodTrack && goodIP && !is_lepton && candTrack.pt() > 5) ++_TrackMult_pt5[_nL];
		if(goodTrack && !is_lepton && candTrack.pt() > 0) ++_TrackMult_noIP_pt0[_nL];
		if(goodTrack && !is_lepton && candTrack.pt() > 1) ++_TrackMult_noIP_pt1[_nL];
		if(goodTrack && !is_lepton && candTrack.pt() > 2) ++_TrackMult_noIP_pt2[_nL];
		if(goodTrack && !is_lepton && candTrack.pt() > 3) ++_TrackMult_noIP_pt3[_nL];
		if(goodTrack && !is_lepton && candTrack.pt() > 4) ++_TrackMult_noIP_pt4[_nL];
		if(goodTrack && !is_lepton && candTrack.pt() > 5) ++_TrackMult_noIP_pt5[_nL];

		// check pf candidates more explicitly: find lepton!
		/*if( fabs( lepVec.Pt() - packedCandidates[j].pt() ) < 0.2 ){
		    std::cout << "pfcand pt, phi and eta: " << packedCandidates[j].pt() << ", " << packedCandidates[j].phi() << ", " << packedCandidates[j].eta() << std::endl;
		    std::cout << "pfcand track pt, phi and eta: " << candTrack.pt() << ", " << candTrack.phi() << ", " << candTrack.eta() << std::endl;
		    TLorentzVector pfcandVec(packedCandidates[j].px(), packedCandidates[j].py(), packedCandidates[j].pz(), packedCandidates[j].energy());
		    std::cout << "Delta R: " << lepVec.DeltaR(trackVec) << std::endl;
		}*/
	    } catch (...){}
	}
    }
}

void LeptonAnalyzer::fillNutauTrackVariables(const reco::GenParticle& gen_tau, edm::Handle<std::vector<pat::PackedCandidate>>& packedCands, const reco::Vertex& vertex){
    //packedCandidate collection
    std::vector<pat::PackedCandidate> packedCandidates;
    for(auto cand = packedCands->cbegin(); cand != packedCands->cend(); ++cand){
        packedCandidates.push_back(*cand);
    }
    if(packedCandidates.size() == 0){
        _Nutau_TrackMult_pt1[_nL] = 0;
	_Nutau_TrackMult_pt5[_nL] = 0;
    } else {
        _Nutau_TrackMult_pt1[_nL] = 0;
	_Nutau_TrackMult_pt5[_nL] = 0;
	auto lep = gen_tau.p4();
        TLorentzVector lepVec(lep.px(), lep.py(), lep.pz(), lep.E());
	//std::cout << "packed Candidates: " << packedCandidates.size() << std::endl;
 	for(unsigned j = 0; j < packedCandidates.size(); ++j){
	    try {//do same try catch thing as in fillLeptonJetVariables just to be safe
		const reco::Track& candTrack = packedCandidates[j].pseudoTrack();
		TLorentzVector trackVec(candTrack.px(), candTrack.py(), candTrack.pz(), candTrack.p());
		double leptondeltaR = trackVec.DeltaR(lepVec);
		bool goodTrack		= leptondeltaR < 0.4 &&
					  candTrack.charge() != 0 &&
					  candTrack.hitPattern().numberOfValidHits() > 7 &&
					  candTrack.hitPattern().numberOfValidPixelHits() > 1 &&
					  candTrack.normalizedChi2() < 5 &&
					  fabs(candTrack.dz(vertex.position())) < 17 &&
					  fabs(candTrack.dxy(vertex.position())) < 17;
					  //packedCandidates[j].fromPV() > 1;
		if(goodTrack && candTrack.pt() > 1) ++_Nutau_TrackMult_pt1[_nL];
		if(goodTrack && candTrack.pt() > 5) ++_Nutau_TrackMult_pt5[_nL];
	    } catch (...){}
	}
    }
}
