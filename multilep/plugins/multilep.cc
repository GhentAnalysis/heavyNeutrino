//include gen particles for matching
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "multilep.h" // Seems lots of includes are getting imported through this one
#include "../interface/LeptonAnalyzer.h"

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
    leptonAnalyzer = new LeptonAnalyzer(iConfig, this);
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

        //lepton variables are added by the leptonAnalyzer class
        leptonAnalyzer->beginJob(outputTree);
/*        //number of leptons
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
*/        //MET
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
    edm::Handle<std::vector<pat::Photon>> photons;                   iEvent.getByToken(photonToken,                       photons);
    edm::Handle<edm::ValueMap<bool>> photonsCutBasedLoose;           iEvent.getByToken(photonCutBasedLooseToken,          photonsCutBasedLoose);
    edm::Handle<edm::ValueMap<bool>> photonsCutBasedMedium;          iEvent.getByToken(photonCutBasedMediumToken,         photonsCutBasedMedium);
    edm::Handle<edm::ValueMap<bool>> photonsCutBasedTight;           iEvent.getByToken(photonCutBasedTightToken,          photonsCutBasedTight);
    edm::Handle<edm::ValueMap<float>> photonsMva;                    iEvent.getByToken(photonMvaToken,                    photonsMva);
    edm::Handle<edm::ValueMap<float>> photonsChargedIsolation;       iEvent.getByToken(photonChargedIsolationToken,       photonsChargedIsolation);
    edm::Handle<edm::ValueMap<float>> photonsNeutralHadronIsolation; iEvent.getByToken(photonNeutralHadronIsolationToken, photonsNeutralHadronIsolation);
    edm::Handle<edm::ValueMap<float>> photonsPhotonIsolation;        iEvent.getByToken(photonPhotonIsolationToken,        photonsPhotonIsolation);
    edm::Handle<std::vector<pat::MET>> mets;                         iEvent.getByToken(metToken,                          mets);

    leptonAnalyzer->analyze(iEvent);

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

//define this as a plug-in
DEFINE_FWK_MODULE(multilep);
