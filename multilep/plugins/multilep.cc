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
    photonAnalyzer = new PhotonAnalyzer(iConfig, this);
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

        //lepton and photon branches
        leptonAnalyzer->beginJob(outputTree);
        photonAnalyzer->beginJob(outputTree);

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
    edm::Handle<std::vector<pat::MET>> mets;                         iEvent.getByToken(metToken,                          mets);

    leptonAnalyzer->analyze(iEvent);
    photonAnalyzer->analyze(iEvent);

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
    //Preselect number of leptons here for code efficiency
    // TODO: boolean function in leptonAnalyzer class for the skim

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
