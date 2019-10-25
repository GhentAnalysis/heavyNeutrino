#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "TLorentzVector.h"                               // displaced specific

multilep::multilep(const edm::ParameterSet& iConfig):
    vtxToken(                 consumes<std::vector<reco::Vertex>>(             iConfig.getParameter<edm::InputTag>("vertices"))),
    genEventInfoToken(        consumes<GenEventInfoProduct>(                   iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    genLumiInfoToken(         consumes<GenLumiInfoHeader, edm::InLumi>(        iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    lheEventInfoToken(        consumes<LHEEventProduct>(                       iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    pileUpToken(              consumes<std::vector<PileupSummaryInfo>>(        iConfig.getParameter<edm::InputTag>("pileUpInfo"))),
    genParticleToken(         consumes<reco::GenParticleCollection>(           iConfig.getParameter<edm::InputTag>("genParticles"))),
    particleLevelPhotonsToken(consumes<reco::GenParticleCollection>(           iConfig.getParameter<edm::InputTag>("particleLevelPhotons"))),
    particleLevelLeptonsToken(consumes<reco::GenJetCollection>(                iConfig.getParameter<edm::InputTag>("particleLevelLeptons"))),
    particleLevelJetsToken(   consumes<reco::GenJetCollection>(                iConfig.getParameter<edm::InputTag>("particleLevelJets"))),
    particleLevelMetsToken(   consumes<reco::METCollection>(                   iConfig.getParameter<edm::InputTag>("particleLevelMets"))),
    muonToken(                consumes<std::vector<pat::Muon>>(                iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                 consumes<std::vector<pat::Electron>>(            iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken(                 consumes<std::vector<pat::Tau>>(                 iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken(              consumes<std::vector<pat::Photon>>(              iConfig.getParameter<edm::InputTag>("photons"))),
    packedCandidatesToken(    consumes<std::vector<pat::PackedCandidate>>(     iConfig.getParameter<edm::InputTag>("packedCandidates"))),
    rhoToken(                 consumes<double>(                                iConfig.getParameter<edm::InputTag>("rho"))),
    metToken(                 consumes<std::vector<pat::MET>>(                 iConfig.getParameter<edm::InputTag>("met"))),
    jetToken(                 consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jets"))),
    jetSmearedToken(          consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsSmeared"))),
    jetSmearedUpToken(        consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsSmearedUp"))),
    jetSmearedDownToken(      consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsSmearedDown"))),
    recoResultsPrimaryToken(  consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("recoResultsPrimary"))),
    recoResultsSecondaryToken(consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("recoResultsSecondary"))),
    triggerToken(             consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("triggers"))),
    prescalesToken(           consumes<pat::PackedTriggerPrescales>(           iConfig.getParameter<edm::InputTag>("prescales"))),
    trigObjToken(             consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    prefireWeightToken(       consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProb"))),
    prefireWeightUpToken(     consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
    prefireWeightDownToken(   consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
    skim(                                                                      iConfig.getUntrackedParameter<std::string>("skim")),
    sampleIsData(                                                              iConfig.getUntrackedParameter<bool>("isData")),
    sampleIs2017(                                                              iConfig.getUntrackedParameter<bool>("is2017")),
    sampleIs2018(                                                              iConfig.getUntrackedParameter<bool>("is2018")),
    sampleIsSUSY(                                                              iConfig.getUntrackedParameter<bool>("isSUSY")),
    storeLheParticles(                                                         iConfig.getUntrackedParameter<bool>("storeLheParticles")),
    storeParticleLevel(                                                        iConfig.getUntrackedParameter<bool>("storeParticleLevel")),
    storeAllTauID(                                                                iConfig.getUntrackedParameter<bool>("storeAllTauID"))
{
    if( is2017() || is2018() ) ecalBadCalibFilterToken = consumes<bool>(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));
    triggerAnalyzer       = new TriggerAnalyzer(iConfig, this);
    leptonAnalyzer        = new LeptonAnalyzer(iConfig, this);
    photonAnalyzer        = new PhotonAnalyzer(iConfig, this);
    jetAnalyzer           = new JetAnalyzer(iConfig, this);
    genAnalyzer           = new GenAnalyzer(iConfig, this);
    lheAnalyzer           = new LheAnalyzer(iConfig, this);
    susyMassAnalyzer      = new SUSYMassAnalyzer(iConfig, this, lheAnalyzer);
    particleLevelAnalyzer = new ParticleLevelAnalyzer(iConfig, this);
}

multilep::~multilep(){
    delete triggerAnalyzer;
    delete leptonAnalyzer;
    delete photonAnalyzer;
    delete jetAnalyzer;
    delete genAnalyzer;
    delete particleLevelAnalyzer;
    delete lheAnalyzer;
    delete susyMassAnalyzer;
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
    outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/i");
    outputTree->Branch("_is2017",                       &sampleIs2017,                  "_is2017/O");
    outputTree->Branch("_is2018",                       &sampleIs2018,                  "_is2018/O");

    if( isMC() && !is2018() ){
        outputTree->Branch("_prefireWeight",              &_prefireWeight,                "_prefireWeight/F");
        outputTree->Branch("_prefireWeightUp",            &_prefireWeightUp,              "_prefireWeightUp/F");
        outputTree->Branch("_prefireWeightDown",          &_prefireWeightDown,            "_prefireWeightDown/F");
    }

    if( isMC() ) lheAnalyzer->beginJob(outputTree, fs);
    if( isSUSY() )  susyMassAnalyzer->beginJob(outputTree, fs);
    if( isMC() ) genAnalyzer->beginJob(outputTree);
    if( isMC() && storeParticleLevel) particleLevelAnalyzer->beginJob(outputTree);
    
    triggerAnalyzer->beginJob(outputTree);
    leptonAnalyzer->beginJob(outputTree);
    photonAnalyzer->beginJob(outputTree);
    jetAnalyzer->beginJob(outputTree);

    _runNb = 0;
}

// ------------ method called for each lumi block ---------
void multilep::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){
    if( isSUSY() ) susyMassAnalyzer->beginLuminosityBlock(iLumi, iSetup);
    _lumiBlock = (unsigned long) iLumi.id().luminosityBlock();
}

//------------- method called for each run -------------
void multilep::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup){
    _runNb = (unsigned long) iRun.id().run();
    triggerAnalyzer->reIndex = true;                                   // HLT results could have different size/order in new run, so look up again the index positions
}

// ------------ method called for each event  ------------
void multilep::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    auto vertices = getHandle(iEvent, vtxToken);
    if( isMC() ) lheAnalyzer->analyze(iEvent);                                            // needs to be run before selection to get correct uncertainties on MC xsection
    if( isSUSY() ) susyMassAnalyzer->analyze(iEvent);                                        // needs to be run after LheAnalyzer, but before all other models

    _nVertex = vertices->size();
    nVertices->Fill(_nVertex, lheAnalyzer->getWeight()); 

    bool applySkim; //better not to shadow class variable with name! // Do not skim if event topology is available on particleLevel 
    if( isMC() && storeParticleLevel ) applySkim = !particleLevelAnalyzer->analyze(iEvent);
    else applySkim = true;

    if(_nVertex == 0)                                                        return;          // Don't consider 0 vertex events
    if(!leptonAnalyzer->analyze(iEvent, iSetup, *(vertices->begin())) and applySkim) return;          // returns false if doesn't pass applySkim condition, so skip event in such case
    if(!photonAnalyzer->analyze(iEvent) and applySkim)                       return;
    if(!jetAnalyzer->analyze(iEvent) and applySkim)                          return;
    if( isMC() ) genAnalyzer->analyze(iEvent);
    triggerAnalyzer->analyze(iEvent);

    _eventNb = (unsigned long) iEvent.id().event();

    if(isMC() and !is2018()){
     _prefireWeight = *(getHandle(iEvent, prefireWeightToken));
     _prefireWeightUp = *(getHandle(iEvent, prefireWeightUpToken));
     _prefireWeightDown = *(getHandle(iEvent, prefireWeightDownToken));
    }

    TLorentzVector lepton1;
    TLorentzVector jet1;
    double nJetBackToBack=0;
    double njet=0;
    njet = jetAnalyzer->_nJets;
    lepton1.SetPtEtaPhiE(leptonAnalyzer->_lPt[0],leptonAnalyzer->_lEta[0],leptonAnalyzer->_lPhi[0],leptonAnalyzer->_lE[0]);
    for (int k =0; k < njet; k ++){
        jet1.SetPtEtaPhiE(jetAnalyzer->_jetPt[k],jetAnalyzer->_jetEta[k],jetAnalyzer->_jetPhi[k],jetAnalyzer->_jetE[k]);
        if (jet1.DeltaR(lepton1) > 1 && jetAnalyzer->_jetIsTight[k]  && jetAnalyzer->_jetPt[k] > 30) nJetBackToBack++;
    }
    bool triggerr_FR=triggerAnalyzer->flag["FR_single_lepton"];
    bool triggerr_trilep=triggerAnalyzer->flag["passTrigger_1l"];

    // to be removeeeeeeeeeeeeeeeee --> so the stuff here below is useless?
    if (skim == "FRsM" and leptonAnalyzer->_lFlavor[0] !=1) return;
    if (skim == "FRsM" and leptonAnalyzer->_lPt[0] > 20) return;
    if (skim == "FRsM" and nJetBackToBack == 0) return;
    if (skim == "FRsM" and !triggerr_FR) return;
    
    if(skim == "FR" and nJetBackToBack == 0) return;
    if(skim == "FR" and !triggerr_FR) return;
    if(skim == "trilep" and !triggerr_trilep)return;

    outputTree->Fill();                                                                  //store calculated event info in root tree
}

//define this as a plug-in
DEFINE_FWK_MODULE(multilep);
