#include "heavyNeutrino/multilep/plugins/multilep.h"


multilep::multilep(const edm::ParameterSet& iConfig):
    vtxToken(                         consumes<std::vector<reco::Vertex>>(        iConfig.getParameter<edm::InputTag>("vertices"))),
    genEventInfoToken(                consumes<GenEventInfoProduct>(              iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    pileUpToken(                      consumes<std::vector<PileupSummaryInfo>>(   iConfig.getParameter<edm::InputTag>("pileUpInfo"))),
    genParticleToken(                 consumes<reco::GenParticleCollection>(      iConfig.getParameter<edm::InputTag>("genParticles"))),
    muonToken(                        consumes<std::vector<pat::Muon>>(           iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                         consumes<std::vector<pat::Electron>>(       iConfig.getParameter<edm::InputTag>("electrons"))),
    eleMvaToken(                      consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMva"))),
    eleMvaHZZToken(                   consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMvaHZZ"))),
    eleCutBasedVetoToken(             consumes<edm::ValueMap<bool>>(              iConfig.getParameter<edm::InputTag>("electronsCutBasedVeto"))),
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
    jetSmearedToken(                  consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jetsSmeared"))),
    jetSmearedUpToken(                consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jetsSmearedUp"))),
    jetSmearedDownToken(              consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jetsSmearedDown"))),
    recoResultsToken(                 consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("recoResults"))),
    triggerToken(                     consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("triggers"))),
    prescalesToken(                   consumes<pat::PackedTriggerPrescales>(      iConfig.getParameter<edm::InputTag>("prescales"))),
    badPFMuonFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
    badChCandFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badChargedCandFilter"))),
    skim(                                                                         iConfig.getUntrackedParameter<std::string>("skim"))
{
    triggerAnalyzer = new TriggerAnalyzer(iConfig, this);
    leptonAnalyzer  = new LeptonAnalyzer(iConfig, this);
    photonAnalyzer  = new PhotonAnalyzer(iConfig, this);
    jetAnalyzer     = new JetAnalyzer(iConfig, this);
    genAnalyzer     = new GenAnalyzer(iConfig, this);
}

// ------------ method called once each job just before starting event loop  ------------
void multilep::beginJob(){
  //Initialize tree with event info
  outputTree = fs->make<TTree>("blackJackAndHookersTree", "blackJackAndHookersTree");

  //Set all branches of the outputTree
  outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
  outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
  outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
  outputTree->Branch("_weight",                       &_weight,                       "_weight/D");
  outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/b");
  outputTree->Branch("_nTrueInt",                     &_nTrueInt,                     "_nTrueInt/F");

  genAnalyzer->beginJob(outputTree);
  triggerAnalyzer->beginJob(outputTree);
  leptonAnalyzer->beginJob(outputTree);
  photonAnalyzer->beginJob(outputTree);
  jetAnalyzer->beginJob(outputTree);

  outputTree->Branch("_met",                          &_met,                          "_met/D");
  outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
}


// ------------ method called for each event  ------------
void multilep::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //determine event number run number and luminosity block
  _runNb     = (unsigned long) iEvent.id().run();
  _lumiBlock = (unsigned long) iEvent.id().luminosityBlock();
  _eventNb   = (unsigned long) iEvent.id().event();

  //Get all objects
  edm::Handle<std::vector<reco::Vertex>> vertices;                 iEvent.getByToken(vtxToken,                          vertices);
  edm::Handle<std::vector<pat::MET>> mets;                         iEvent.getByToken(metToken,                          mets);
  edm::Handle<GenEventInfoProduct> genEventInfo;                   iEvent.getByToken(genEventInfoToken,                 genEventInfo);
  edm::Handle<std::vector<PileupSummaryInfo>>  pileUpInfo;         iEvent.getByToken(pileUpToken,                       pileUpInfo);

  _weight  = genEventInfo->weight();
  _nVertex = vertices->size();
  for(auto puI = pileUpInfo->begin(); puI != pileUpInfo->end(); ++puI){
    if(puI->getBunchCrossing() == 0) _nTrueInt = puI->getTrueNumInteractions(); // getTrueNumInteractions should be the same for all bunch crosssings
  }

  
  if(!leptonAnalyzer->analyze(iEvent, *(vertices->begin()))) return; // returns false if doesn't pass skim condition, so skip event in such case
  if(!photonAnalyzer->analyze(iEvent)) return;
  genAnalyzer->analyze(iEvent);
  triggerAnalyzer->analyze(iEvent);
  jetAnalyzer->analyze(iEvent);

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
