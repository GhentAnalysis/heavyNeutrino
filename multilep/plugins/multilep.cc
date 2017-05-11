//include gen particles for matching
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "multilep.h" // Seems lots of includes are getting imported through this one

multilep::multilep(const edm::ParameterSet& iConfig):
    vtxToken(                         consumes<std::vector<reco::Vertex>>(        iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken(                        consumes<std::vector<pat::Muon>>(           iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                         consumes<std::vector<pat::Electron>>(       iConfig.getParameter<edm::InputTag>("electrons"))),
    eleMvaToken(                      consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMva"))),
    eleMvaHZZToken(                   consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMvaHZZ"))),
    eleCutBasedTightToken(            consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("electronsCutBasedTight"))),
    eleCutBasedMediumToken(           consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("electronsCutBasedMedium"))),
    tauToken(                         consumes<std::vector<pat::Tau>>(            iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken(                      consumes<std::vector<pat::Photon>>(         iConfig.getParameter<edm::InputTag>("photons"))),
    photonCutBasedLooseToken(         consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("photonsCutBasedLoose"))),
    photonCutBasedMediumToken(        consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("photonsCutBasedMedium"))),
    photonCutBasedTightToken(         consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("photonsCutBasedTight"))),
    photonMvaToken(                   consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("photonsMva"))),
    photonChargedIsolationToken(      consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("photonsChargedIsolation"))),
    photonNeutralHadronIsolationToken(consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("photonsNeutralHadronIsolation"))),
    photonPhotonIsolationToken(       consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("photonsPhotonIsolation"))),
    packedCandidatesToken(            consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
    rhoToken(                         consumes<double>(                           iConfig.getParameter<edm::InputTag>("rhoCentralNeutral"))),
    rhoTokenAll(                      consumes<double>(                           iConfig.getParameter<edm::InputTag>("rhoAll"))),
    metToken(                         consumes<std::vector<pat::MET>>(            iConfig.getParameter<edm::InputTag>("met"))),
    jetToken(                         consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jets"))),
  //jecToken(                         consumes<reco::JetCorrector>(               edm::InputTag("ak4PFCHSL1FastL2L3Corrector")))
  //jecToken(                         consumes<reco::JetCorrector>(               edm::InputTag("ak4PFCHSL3Absolute")))
    triggerToken(                     consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("triggers"))),
    recoResultsToken(                 consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("recoResults"))),
    badPFMuonFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
    badChCandFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badChargedCandFilter")))
{
    //usesResource("TFileService");
}

// ------------ method called once each job just before starting event loop  ------------
void multilep::beginJob(){
        //Initialize tree with event info
        //outputTree = new TTree("multiLepTree", "high level event info");
        outputTree = fs->make<TTree>("blackJackAndHookersTree", "blackJackAndHookersTree");

        //Set all branches of the outputTree
        //event labels
        outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
        outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
        outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
        //jet variables
        outputTree->Branch("_nJets",                        &_nJets,                        "_nJets/b");
        outputTree->Branch("_jetPt",                        &_jetPt,                        "_jetPt[_nJets]/D");
        outputTree->Branch("_jetEta",                       &_jetEta,                       "_jetEta[_nJets]/D");
        outputTree->Branch("_jetPhi",                       &_jetPhi,                       "_jetPhi[_nJets]/D");
        outputTree->Branch("_jetE",                         &_jetE,                         "_jetE[_nJets]/D");

        outputTree->Branch("_nPhoton",                      &_nPhoton,                      "_nPhoton/b");
        outputTree->Branch("_photonPt",                     &_photonPt,                     "_photonPt[_nPhoton]/F");
        outputTree->Branch("_photonEta",                    &_photonEta,                    "_photonEta[_nPhoton]/F");
        outputTree->Branch("_photonPhi",                    &_photonPhi,                    "_photonPhi[_nPhoton]/F");
        outputTree->Branch("_photonE",                      &_photonE,                      "_photonE[_nPhoton]/F");
        outputTree->Branch("_photonCutBasedLoose",          &_photonCutBasedLoose,          "_photonCutBasedLoose[_nPhoton]/O");
        outputTree->Branch("_photonCutBasedMedium",         &_photonCutBasedMedium,         "_photonCutBasedMedium[_nPhoton]/O");
        outputTree->Branch("_photonCutBasedLoose",          &_photonCutBasedLoose,          "_photonCutBasedLoose[_nPhoton]/O");
        outputTree->Branch("_photonMva",                    &_photonMva,                    "_photonMva[_nPhoton]/F");
        outputTree->Branch("_photonChargedIsolation",       &_photonChargedIsolation,       "_photonChargedIsolation[_nPhoton]/F");
        outputTree->Branch("_photonNeutralHadronIsolation", &_photonNeutralHadronIsolation, "_photonNeutralHadronIsolation[_nPhoton]/F");
        outputTree->Branch("_photonSigmaIetaIeta",          &_photonSigmaIetaIeta,          "_photonSigmaIetaIeta[_nPhoton]/F");
        outputTree->Branch("_photonHadronicOverEm",         &_photonHadronicOverEm,         "_photonHadronicOverEm[_nPhoton]/F");
        outputTree->Branch("_photonPassElectronVeto",       &_photonPassElectronVeto,       "_photonPassElectronVeto[_nPhoton]/O");
        outputTree->Branch("_photonHasPixelSeed",           &_photonHasPixelSeed,           "_photonHasPixelSeed[_nPhoton]/O");

        //lepton variables
        //number of leptons
        outputTree->Branch("_nL",                           &_nL,                           "_nL/b");
        outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
        outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
        outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");
        outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/b");
        //lepton kinematics
        outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
        outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
        outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
        outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
        //lepton vertex variables
        outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
        outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
        outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
        outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
        //other lepton variables
        outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
        outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
        //MET
        outputTree->Branch("_met",                          &_met,                          "_met/D");
        outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
        //Trigger and MET filter decisions
        outputTree->Branch("_passHnlTrigger",               &_passHnlTrigger,               "_passHnlTrigger[4]/O");
        outputTree->Branch("_metFiltersFlagged",            &_metFiltersFlagged,            "_metFiltersFlagged/O");
        outputTree->Branch("_badMuonFlagged",               &_badMuonFlagged,               "_badMuonFlagged/O");
        outputTree->Branch("_badCloneMuonFlagged",          &_badCloneMuonFlagged,          "_badCloneMuonFlagged/O");
}


// ------------ method called for each event  ------------
void multilep::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

    //determine event number run number and luminosity block
    _runNb     = (unsigned long) iEvent.id().run();
    _lumiBlock = (unsigned long) iEvent.id().luminosityBlock();
    _eventNb   = (unsigned long) iEvent.id().event();

    //Get all objects
    edm::Handle<std::vector<reco::Vertex>> vertices;                 iEvent.getByToken(vtxToken,                          vertices);
    edm::Handle<double> rhoJets;                                     iEvent.getByToken(rhoToken,                          rhoJets);    // For JEC
    edm::Handle<double> rhoJetsAll;                                  iEvent.getByToken(rhoTokenAll,                       rhoJetsAll); // For PUC
    edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(jetToken,                          jets);
  //edm::Handle<reco::JetCorrector> jec;                             iEvent.getByToken(jecToken,                          jec);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(packedCandidatesToken,             packedCands);
    edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(muonToken,                         muons);
    edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(eleToken,                          electrons);
    edm::Handle<edm::ValueMap<float>> electronsMva;                  iEvent.getByToken(eleMvaToken,                       electronsMva);
    edm::Handle<edm::ValueMap<float>> electronsMvaHZZ;               iEvent.getByToken(eleMvaHZZToken,                    electronsMvaHZZ);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedTight;         iEvent.getByToken(eleCutBasedTightToken,             electronsCutBasedTight);
    edm::Handle<edm::ValueMap<bool>> electronsCutBasedMedium;        iEvent.getByToken(eleCutBasedMediumToken,            electronsCutBasedMedium);
    edm::Handle<std::vector<pat::Tau>> taus;                         iEvent.getByToken(tauToken,                          taus);
    edm::Handle<std::vector<pat::Photon>> photons;                   iEvent.getByToken(photonToken,                       photons);
    edm::Handle<edm::ValueMap<bool>> photonsCutBasedLoose;           iEvent.getByToken(photonCutBasedLooseToken,          photonsCutBasedLoose);
    edm::Handle<edm::ValueMap<bool>> photonsCutBasedMedium;          iEvent.getByToken(photonCutBasedMediumToken,         photonsCutBasedMedium);
    edm::Handle<edm::ValueMap<bool>> photonsCutBasedTight;           iEvent.getByToken(photonCutBasedTightToken,          photonsCutBasedTight);
    edm::Handle<edm::ValueMap<float>> photonsMva;                    iEvent.getByToken(photonMvaToken,                    photonsMva);
    edm::Handle<edm::ValueMap<float>> photonsChargedIsolation;       iEvent.getByToken(photonChargedIsolationToken,       photonsChargedIsolation);
    edm::Handle<edm::ValueMap<float>> photonsNeutralHadronIsolation; iEvent.getByToken(photonNeutralHadronIsolationToken, photonsNeutralHadronIsolation);
    edm::Handle<edm::ValueMap<float>> photonsPhotonIsolation;        iEvent.getByToken(photonPhotonIsolationToken,        photonsPhotonIsolation);
    edm::Handle<std::vector<pat::MET>> mets;                         iEvent.getByToken(metToken,                          mets);

    //loop over jets
    _nJets = 0;
    for(const pat::Jet& jet: *jets){
        //double corr = jec->correction(jet);
        //if(jet.pt()< 25) continue;
        if(jet.pt() < 25) continue;
        if(fabs(jet.eta()) > 2.4) continue;
        _jetPt[_nJets] = jet.pt();
        _jetEta[_nJets] = jet.eta();
        _jetPhi[_nJets] = jet.phi();
        _jetE[_nJets] = jet.energy();
        ++_nJets;
    }

    // Loop over photons
    _nPhoton = 0;
    for(auto photon = photons->begin(); photon != photons->end(); ++photon){
      auto photonRef = edm::Ref<std::vector<pat::Photon>>(photons, (photon - photons->begin()));

      _photonPt[_nPhoton]                      = photon->pt();
      _photonEta[_nPhoton]                     = photon->eta();
      _photonPhi[_nPhoton]                     = photon->phi();
      _photonE[_nPhoton]                       = photon->energy();
      _photonCutBasedLoose[_nPhoton]           = (*photonsCutBasedLoose)[photonRef];
      _photonCutBasedMedium[_nPhoton]          = (*photonsCutBasedMedium)[photonRef];
      _photonCutBasedTight[_nPhoton]           = (*photonsCutBasedTight)[photonRef];
      _photonMva[_nPhoton]                     = (*photonsMva)[photonRef];
      _photonChargedIsolation[_nPhoton]        = (*photonsChargedIsolation)[photonRef];
      _photonNeutralHadronIsolation[_nPhoton]  = (*photonsNeutralHadronIsolation)[photonRef];
      _photonPhotonIsolation[_nPhoton]         = (*photonsPhotonIsolation)[photonRef];
      _photonSigmaIetaIeta[_nPhoton]           = photon->full5x5_sigmaIetaIeta();
      _photonHadronicOverEm[_nPhoton]          = photon->hadronicOverEm();
      _photonPassElectronVeto[_nPhoton]        = photon->passElectronVeto();
      _photonHasPixelSeed[_nPhoton]            = photon->hasPixelSeed();

      ++_nPhoton;
    }


    //lepton selection
    _nL   = 0;
    _nMu  = 0;
    _nEle = 0;
    _nTau = 0;
    //loop over muons
    for(const pat::Muon& mu : *muons){
        if(_nL == nL_max)            continue;
        if(mu.innerTrack().isNull()) continue;
        if(mu.pt() < 5)              continue;
        if(fabs(mu.eta()) > 2.4)     continue;
        fillLeptonKinVars(mu);
        fillLeptonGenVars(mu.genParticle());
        _flavor[_nL]  = 1;
        //Vertex variables // better move this too in fillLeptonIdVars
        _dxy[_nL]     = mu.innerTrack()->dxy();
        _dz[_nL]      = mu.innerTrack()->dz();
        _3dIP[_nL]    = mu.dB(pat::Muon::PV3D);
        _3dIPSig[_nL] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);
        //Isolation variables
        _relIso[_nL]  = Tools::getRelIso03(mu, *rhoJets);
        _miniIso[_nL] = Tools::getMiniIso(mu, *packedCands, 0.2, *rhoJets);
//        fillLeptonIdVars(mu);
        ++_nMu;
        ++_nL;
    }

    // Loop over electrons (note: using iterator we can easily get the ref too)
    for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
      auto electronRef = edm::Ref<std::vector<pat::Electron>>(electrons, (ele - electrons->begin()));
      if(_nL == nL_max)            continue;
      if(ele->gsfTrack().isNull()) continue;
      if(ele->pt() < 10)           continue;
      if(fabs(ele->eta()) > 2.5)   continue;
      fillLeptonKinVars(*ele);
      fillLeptonGenVars(ele->genParticle());
      _flavor[_nL]  = 0;
      //Vertex varitables
      _dxy[_nL]     = ele->gsfTrack()->dxy();
      _dz[_nL]      = ele->gsfTrack()->dz();
      _3dIP[_nL]    = ele->dB(pat::Electron::PV3D);
      _3dIPSig[_nL] = ele->dB(pat::Electron::PV3D)/ele->edB(pat::Electron::PV3D);


      //isolation variables
      _relIso[_nL]  = Tools::getRelIso03(*ele, *rhoJets);
      _miniIso[_nL] = Tools::getMiniIso(*ele, *packedCands, 0.2, *rhoJets);

      //id variables TODO: put those in the tree
      float _mvaValue              = (*electronsMva)[electronRef];
/*      float _mvaValue_HZZ          = (*electronsMva)[electronRef];
      bool _passedCutBasedIdTight  = (*electronsCutBasedTight)[electronRef];
      bool _passedCutBasedIdMedium = (*electronsCutBasedMedium)[electronRef];*/
      std::cout << _mvaValue << std::endl;
/*
      _passedMVA_SUSY[leptonCounter][0] = tools::passed_loose_MVA_FR_slidingCut( &*electron, _mvaValue[leptonCounter], _mvaValue_HZZ[leptonCounter]);
      _passedMVA_SUSY[leptonCounter][1] = tools::passed_medium_MVA_FR_slidingCut(&*electron, _mvaValue[leptonCounter]);
      _passedMVA_SUSY[leptonCounter][2] = tools::passed_tight_MVA_FR_slidingCut( &*electron, _mvaValue[leptonCounter]);
*/
//        fillLeptonIdVars(ele);

      ++_nEle;
      ++_nL;
    }
    _nLight = _nEle + _nMu;

    //loop over taus
    for(const pat::Tau& tau : *taus){
        if(_nL == nL_max)         continue;
        if(tau.pt() < 20)         continue; //investigate up to what Pt threshold taus can be properly reconstructed
        if(fabs(tau.eta()) > 2.3) continue;
        fillLeptonKinVars(tau);
        _flavor[_nL]  = 2;
        _dxy[_nL]     = tau.dxy();
        //_dz[_nL]    = tau.dz();
        _3dIP[_nL]    = tau.ip3d();
        _3dIPSig[_nL] = tau.ip3d_Sig();
        ++_nTau;
        ++_nL;
    }
    //Preselect number of leptons here for code efficiency

    //Fill trigger and MET filter decisions
    fillTriggerVars(iEvent);
    fillMetFilterVars(iEvent);

    //determine the met of the event
    const pat::MET& met = (*mets).front();
    _met = met.pt();
    _metPhi = met.phi();
    //store calculated event info in root tree
    outputTree->Fill();

}


// ------------ method called once each job just after ending the event loop  ------------
void multilep::endJob(){

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void multilep::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


//------------- Fill reconstructed lepton variables ---------------------
void multilep::fillLeptonKinVars(const reco::Candidate& lepton){
    //kinematics
    _lPt[_nL]    = lepton.pt();
    _lEta[_nL]   = lepton.eta();
    _lPhi[_nL]   = lepton.phi();
    _lE[_nL]     = lepton.energy();
    _charge[_nL] = lepton.charge();
}


//------------- Fill MC-truth lepton variables -------------------
void multilep::fillLeptonGenVars(const reco::GenParticle* genParticle){
    if(genParticle != nullptr) _isPrompt[_nL] = (genParticle)->isPromptFinalState();
    else                       _isPrompt[_nL] = false;
}
//----------------------------------------------------------

//define this as a plug-in
DEFINE_FWK_MODULE(multilep);
