#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


multilep::multilep(const edm::ParameterSet& iConfig):
    vtxToken(                         consumes<std::vector<reco::Vertex>>(        iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken(                        consumes<std::vector<pat::Muon>>(           iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                         consumes<std::vector<pat::Electron>>(       iConfig.getParameter<edm::InputTag>("electrons"))),
    eleMvaToken(                      consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMva"))),
    eleMvaHZZToken(                   consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMvaHZZ"))),
    eleCutBasedLooseToken(            consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("electronsCutBasedLoose"))),
    eleCutBasedMediumToken(           consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("electronsCutBasedMedium"))),
    eleCutBasedTightToken(            consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("electronsCutBasedTight"))),
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
    triggerToken(                     consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("triggers"))),
    recoResultsToken(                 consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("recoResults"))),
    badPFMuonFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
    badChCandFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badChargedCandFilter")))
{
    leptonAnalyzer = new LeptonAnalyzer(iConfig, this);
    photonAnalyzer = new PhotonAnalyzer(iConfig, this);
    jetAnalyzer    = new JetAnalyzer(iConfig, this);
}

// ------------ method called once each job just before starting event loop  ------------
void multilep::beginJob(){
  //Initialize tree with event info
  outputTree = fs->make<TTree>("blackJackAndHookersTree", "blackJackAndHookersTree");

  //Set all branches of the outputTree
  outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
  outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
  outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
  outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/b");

  leptonAnalyzer->beginJob(outputTree);
  photonAnalyzer->beginJob(outputTree);
  jetAnalyzer->beginJob(outputTree);

  outputTree->Branch("_met",                          &_met,                          "_met/D");
  outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");

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
//  edm::Handle<double> rhoJets;                                     iEvent.getByToken(rhoToken,                          rhoJets);    // For JEC
//  edm::Handle<double> rhoJetsAll;                                  iEvent.getByToken(rhoTokenAll,                       rhoJetsAll); // For PUC
  edm::Handle<std::vector<pat::MET>> mets;                         iEvent.getByToken(metToken,                          mets);


  _nVertex = vertices->size();

  leptonAnalyzer->analyze(iEvent, *(vertices->begin()));
  photonAnalyzer->analyze(iEvent);

  //Preselect number of leptons here for code efficiency
  // TODO: boolean function in leptonAnalyzer class for the skim

  //Fill trigger and MET filter decisions
  fillTriggerVars(iEvent);
  fillMetFilterVars(iEvent);

  //determine the met of the event
  const pat::MET& met = (*mets).front();
  _met    = met.pt();
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
