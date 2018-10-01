#include "heavyNeutrino/multilep/plugins/multilep.h"


multilep::multilep(const edm::ParameterSet& iConfig):
    beamSpotToken(		      consumes<reco::BeamSpot>(			  iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
    vtxToken(                         consumes<std::vector<reco::Vertex>>(        iConfig.getParameter<edm::InputTag>("vertices"))),
    genEventInfoToken(                consumes<GenEventInfoProduct>(              iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    genLumiInfoToken(                 consumes<GenLumiInfoHeader, edm::InLumi>(   iConfig.getParameter<edm::InputTag>("genEventInfo"))), //NOT SURE IF THIS WILL WORK, CHECK!
    lheEventInfoToken(                consumes<LHEEventProduct>(                  iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    pileUpToken(                      consumes<std::vector<PileupSummaryInfo>>(   iConfig.getParameter<edm::InputTag>("pileUpInfo"))),
    genParticleToken(                 consumes<reco::GenParticleCollection>(      iConfig.getParameter<edm::InputTag>("genParticles"))),
    packedGenParticleToken(           consumes<std::vector<pat::PackedGenParticle>>(      iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
    muonToken(                        consumes<std::vector<pat::Muon>>(           iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                         consumes<std::vector<pat::Electron>>(       iConfig.getParameter<edm::InputTag>("electrons"))),
    eleMvaToken(                      consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMva"))),
    eleMvaHZZToken(                   consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronsMvaHZZ"))),
    eleMvaFall17IsoToken(             consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronMvaFall17Iso"))),
    eleMvaFall17NoIsoToken(           consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("electronMvaFall17NoIso"))),
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
    photonFull5x5SigmaIEtaIPhiToken(  consumes<edm::ValueMap<float>>(             iConfig.getParameter<edm::InputTag>("photonsFull5x5SigmaIEtaIPhi"))),
    packedCandidatesToken(            consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
    rhoToken(                         consumes<double>(                           iConfig.getParameter<edm::InputTag>("rho"))),
    metToken(                         consumes<std::vector<pat::MET>>(            iConfig.getParameter<edm::InputTag>("met"))),
    jetToken(                         consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jets"))),
    jetSmearedToken(                  consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jetsSmeared"))),
    jetSmearedUpToken(                consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jetsSmearedUp"))),
    jetSmearedDownToken(              consumes<std::vector<pat::Jet>>(            iConfig.getParameter<edm::InputTag>("jetsSmearedDown"))),
    recoResultsPrimaryToken(          consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("recoResultsPrimary"))),
    recoResultsSecondaryToken(        consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("recoResultsSecondary"))),
    triggerToken(                     consumes<edm::TriggerResults>(              iConfig.getParameter<edm::InputTag>("triggers"))),
    prescalesToken(                   consumes<pat::PackedTriggerPrescales>(      iConfig.getParameter<edm::InputTag>("prescales"))),
    badPFMuonFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
    badChCandFilterToken(             consumes<bool>(                             iConfig.getParameter<edm::InputTag>("badChargedCandFilter"))),
    skim(                                                                         iConfig.getUntrackedParameter<std::string>("skim")),
    isData(                                                                       iConfig.getUntrackedParameter<bool>("isData")),
    is2017(                                                                       iConfig.getUntrackedParameter<bool>("is2017")),
    isSUSY(                                                                       iConfig.getUntrackedParameter<bool>("isSUSY"))
    //jecPath(                                                                      iConfig.getParameter<edm::FileInPath>("JECtxtPath").fullPath())
{
    triggerAnalyzer = new TriggerAnalyzer(iConfig, this);
    leptonAnalyzer  = new LeptonAnalyzer(iConfig, this);
    photonAnalyzer  = new PhotonAnalyzer(iConfig, this);
    jetAnalyzer     = new JetAnalyzer(iConfig, this);
    genAnalyzer     = new GenAnalyzer(iConfig, this);
    lheAnalyzer     = new LheAnalyzer(iConfig, this);
    susyMassAnalyzer= new SUSYMassAnalyzer(iConfig, this, lheAnalyzer);

    /*
    //initialize jec txt files
    std::string dirtyHack = "dummy.txt";
    std::string path = jecPath.substr(0, jecPath.size() - dirtyHack.size() );
    jec = new JEC(path, isData, is2017);  //dummy.txt is a dirty hack to give directory parameter in python file
    */
}

multilep::~multilep(){
    delete triggerAnalyzer;
    delete leptonAnalyzer;
    delete photonAnalyzer;
    delete jetAnalyzer;
    delete genAnalyzer;
    delete lheAnalyzer;
    delete susyMassAnalyzer;
    //delete jec;
}

// ------------ method called once each job just before starting event loop  ------------
void multilep::beginJob(){

    //Initialize tree with event info

    outputTree = fs->make<TTree>("blackJackAndHookersTree", "blackJackAndHookersTree");
    nVertices  = fs->make<TH1D>("nVertices", "Number of vertices", 120, 0, 120);

    //Set all branches of the outputTree
    outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
    outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
    outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
    outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/b");

    if(!isData) lheAnalyzer->beginJob(outputTree, fs);
    if(isSUSY)  susyMassAnalyzer->beginJob(outputTree, fs);
    if(!isData) genAnalyzer->beginJob(outputTree);
    triggerAnalyzer->beginJob(outputTree);
    leptonAnalyzer->beginJob(outputTree);
    photonAnalyzer->beginJob(outputTree);
    jetAnalyzer->beginJob(outputTree);

    _runNb = 0;
}

// ------------ method called for each lumi block ---------
void multilep::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){
    if(isSUSY) susyMassAnalyzer->beginLuminosityBlock(iLumi, iSetup);
    _lumiBlock = (unsigned long) iLumi.id().luminosityBlock();
}

//------------- method called for each run -------------
void multilep::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup){

    // HLT results could have different size/order in new run, so look up again de index positions
    triggerAnalyzer->reIndex = true;

    //get Run number
    _runNb = (unsigned long) iRun.id().run();

    //update JEC 
    //jec->updateJEC(_runNb);
}

// ------------ method called for each event  ------------
void multilep::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    edm::Handle<std::vector<reco::Vertex>> vertices; iEvent.getByToken(vtxToken, vertices);
    if(!isData) lheAnalyzer->analyze(iEvent);                          // needs to be run before selection to get correct uncertainties on MC xsection
    if(isSUSY) susyMassAnalyzer->analyze(iEvent);                      // needs to be run after LheAnalyzer, but before all other models

    //extract number of vertices 
    _nVertex = vertices->size();
    nVertices->Fill(_nVertex, lheAnalyzer->getWeight()); 
    if(_nVertex == 0) return;                                      //Don't consider 0 vertex events

    if(!leptonAnalyzer->analyze(iEvent, iSetup, *(vertices->begin()))) return; // returns false if doesn't pass skim condition, so skip event in such case
    if(!isData) genAnalyzer->analyze(iEvent);                          // needs to be run before photonAnalyzer for matching purposes
    if(!photonAnalyzer->analyze(iEvent)) return;
    triggerAnalyzer->analyze(iEvent);
    jetAnalyzer->analyze(iEvent);

    //determine event number run number and luminosity block
    _eventNb   = (unsigned long) iEvent.id().event();

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
