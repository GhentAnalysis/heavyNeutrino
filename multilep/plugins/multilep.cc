//include gen particles for matching
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "multilep.h"

multilep::multilep(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
	muonToken_(consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
	eleToken_(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
	tauToken_(consumes<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"))),
	packedCandidatesToken_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
	rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCentralNeutral"))),
	rhoTokenAll_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
	metToken_(consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("met"))),
	jetToken_(consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
	//jecToken_(consumes<reco::JetCorrector>(edm::InputTag("ak4PFCHSL1FastL2L3Corrector")))
	//jecToken_(consumes<reco::JetCorrector>(edm::InputTag("ak4PFCHSL3Absolute")))
	triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggers"))),
	recoResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("recoResults"))),
	badPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
	badChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badChargedCandFilter")))
{
    //usesResource("TFileService");
}


multilep::~multilep()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void multilep::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	//determine event number run number and luminosity block
	_runNb = (unsigned long) iEvent.id().run();
	_lumiBlock = (unsigned long) iEvent.id().luminosityBlock();
	_eventNb = (unsigned long) iEvent.id().event();
    
	//read vertex data
	edm::Handle<std::vector<reco::Vertex> > vertices;
    iEvent.getByToken(vtxToken_, vertices);
	//reco::Vertex::Point primVert = vertices->begin()->position();

	//read energy densities for JEC and PUC
	edm::Handle<double> rhoJets;
	iEvent.getByToken(rhoToken_, rhoJets);
	edm::Handle<double> rhoJetsAll;
	iEvent.getByToken(rhoTokenAll_, rhoJetsAll);


	//read jet data
	edm::Handle<std::vector<pat::Jet>> jets;
	iEvent.getByToken(jetToken_, jets);
	//read JEC data
	//edm::Handle<reco::JetCorrector> jec;
	//iEvent.getByToken(jecToken_, jec);
	//loop over jets
	_nJets = 0;
	for(const pat::Jet& jet: *jets){
		//double corr = jec->correction(jet);
		//if(jet.pt()< 25) continue;
		if(jet.pt() < 25) continue;
		if(fabs(jet.eta()) > 2.4) continue;
		std::cout << "jet pt: " << jet.pt() << "	" << "jet eta: " << jet.eta() << std::endl;
		_jetPt[_nJets] = jet.pt();
		_jetEta[_nJets] = jet.eta();
		_jetPhi[_nJets] = jet.phi();
		_jetE[_nJets] = jet.energy();
		++_nJets;
	}

	//read packed particle collection used in isolation calculations
	edm::Handle<std::vector<pat::PackedCandidate>> packedCands;
	iEvent.getByToken(packedCandidatesToken_, packedCands);
	//lepton selection
	_nL = 0;
	_nMu = 0;
	_nEle = 0;
	_nTau = 0;
	//read muon data
	edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByToken(muonToken_, muons);
	//loop over muons
	for(const pat::Muon& mu : *muons){
		if(_nL == nL_max) continue;
		if(mu.innerTrack().isNull()) continue;
		if(mu.pt() < 5) continue;
		if(fabs(mu.eta()) > 2.4) continue;
		fillLeptonKinVars(mu);
		fillLeptonGenVars(mu.genParticle());
		_flavor[_nL] = 1;
		//Vertex variables
		_dxy[_nL] = mu.innerTrack()->dxy();
		_dz[_nL] = mu.innerTrack()->dz();
		_3dIP[_nL] = mu.dB(pat::Muon::PV3D);
		_3dIPSig[_nL] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);
		//Isolation variables
		_relIso[_nL] = Tools::getRelIso03(mu, *rhoJets);
		_miniIso[_nL] = Tools::getMiniIso(mu, *packedCands, 0.2, *rhoJets);
		++_nMu;
		++_nL;	
	}
	
	//read electron data
	edm::Handle<std::vector<pat::Electron> > electrons;
	iEvent.getByToken(eleToken_, electrons);
	//loop over electrons
	for(const pat::Electron& ele : *electrons){
		if(_nL == nL_max) continue;
		if(ele.gsfTrack().isNull()) continue;
		if(ele.pt() < 10) continue;
		if(fabs(ele.eta()) > 2.5) continue;
		fillLeptonKinVars(ele);
		fillLeptonGenVars(ele.genParticle());
		_flavor[_nL] = 0;
		//Vertex varitables
		_dxy[_nL] = ele.gsfTrack()->dxy();
		_dz[_nL] = ele.gsfTrack()->dz();
		_3dIP[_nL] = ele.dB(pat::Electron::PV3D);
		_3dIPSig[_nL] = ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D);
		//isolation variables
		_relIso[_nL] = Tools::getRelIso03(ele, *rhoJets);
		_miniIso[_nL] = Tools::getMiniIso(ele, *packedCands, 0.2, *rhoJets);
		++_nEle;
		++_nL;
	}
	_nLight = _nEle + _nMu;
	//read tau data
	edm::Handle<std::vector<pat::Tau> > taus;
	iEvent.getByToken(tauToken_, taus);
	//loop over taus
	for(const pat::Tau& tau : *taus){
		if(_nL == nL_max) continue;
		if(tau.pt() < 20) continue; //investigate up to what Pt threshold taus can be properly reconstructed
		if(fabs(tau.eta()) > 2.3) continue;
		fillLeptonKinVars(tau);
		_flavor[_nL] = 2;
		_dxy[_nL] = tau.dxy();
		//_dz[_nL] = tau.dz();
		_3dIP[_nL] = tau.ip3d();
		_3dIPSig[_nL] = tau.ip3d_Sig ();
		++_nTau;
		++_nL;
	}
	//Preselect number of leptons here for code efficiency

	//Fill trigger and MET filter decisions
	fillTriggerVars(iEvent);
	fillMetFilterVars(iEvent);
	//read met data
	edm::Handle<std::vector<pat::MET> > mets;
	iEvent.getByToken(metToken_, mets);
	//determine the met of the event
	const pat::MET& met = (*mets).front();
	_met = met.pt();
	_metPhi = met.phi();
	std::cout << "met  : " << _met << "		" << "met phi : " << _metPhi << std::endl;
	//std::cout << "##################################################" << std::endl;
	//store calculated event info in root tree
	outputTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void multilep::beginJob(){
	//Initialize tree with event info
	//outputTree = new TTree("multiLepTree", "high level event info");
	outputTree = fs->make<TTree>("blackJackAndHookersTree", "blackJackAndHookersTree");
	//Set all branches of the outputTree
	//event labels
	outputTree->Branch("_runNb", &_runNb, "_runNb/l");
	outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
	outputTree->Branch("_eventNb", &_eventNb, "_eventNb/l");
	//jet variables
	outputTree->Branch("_nJets", &_nJets, "_nJets/b");
	outputTree->Branch("_jetPt", &_jetPt, "_jetPt[_nJets]/D");
	outputTree->Branch("_jetEta", &_jetEta, "_jetEta[_nJets]/D");
	outputTree->Branch("_jetPhi", &_jetPhi, "_jetPhi[_nJets]/D");
	outputTree->Branch("_jetE", &_jetE, "_jetE[_nJets]/D");
	//lepton variables
	//number of leptons
	outputTree->Branch("_nL", &_nL, "_nL/b");
	outputTree->Branch("_nMu", &_nMu, "_nMu/b");
	outputTree->Branch("_nEle", &_nEle, "_nEle/b");
	outputTree->Branch("_nLight", &_nLight, "_nLight/b");
	outputTree->Branch("_nTau", &_nTau, "_nTau/b");
	//lepton kinematics
	outputTree->Branch("_lPt", &_lPt, "_lPt[_nL]/D");
	outputTree->Branch("_lEta", &_lEta, "_lEta[_nL]/D");
	outputTree->Branch("_lPhi", &_lPhi, "_lPhi[_nL]/D");
	outputTree->Branch("_lEta", &_lEta, "_lEta[_nL]/D");
	//lepton vertex variables
	outputTree->Branch("_dxy", &_dxy, "_dxy[_nL]/D");
	outputTree->Branch("_dz", &_dz, "_dz[_nL]/D");
	outputTree->Branch("_3dIP", &_3dIP, "_3dIP[_nL]/D");
	outputTree->Branch("_3dIPSig", &_3dIPSig, "_3dIPSig[_nL]/D");
	//other lepton variables
	outputTree->Branch("_relIso", &_relIso, "_relIso[_nLight]/D");
	outputTree->Branch("_miniIso", &_miniIso, "_miniIso[_nLight]/D");
	//MET
	outputTree->Branch("_met", &_met, "_met/D");
	outputTree->Branch("_metPhi", &_metPhi,"_metPhi/D");
	//Trigger and MET filter decisions
	outputTree->Branch("_passHnlTrigger", &_passHnlTrigger, "_passHnlTrigger[4]/O");
	outputTree->Branch("_metFiltersFlagged", &_metFiltersFlagged, "_metFiltersFlagged/O");
	outputTree->Branch("_badMuonFlagged", &_badMuonFlagged, "_badMuonFlagged/O");
	outputTree->Branch("_badCloneMuonFlagged", &_badCloneMuonFlagged, "_badCloneMuonFlagged/O");
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
	_lPt[_nL] = lepton.pt();
	_lEta[_nL] = lepton.eta();
	_lPhi[_nL] = lepton.phi();
	_lE[_nL] = lepton.energy();
	//Charge
	_charge[_nL] = lepton.charge();
}


//------------- Fill MC-truth lepton variables -------------------
void multilep::fillLeptonGenVars(const reco::GenParticle* genParticle){
	if(genParticle != nullptr){
		_isPrompt[_nL] = (genParticle)->isPromptFinalState();
	} else{
		_isPrompt[_nL] = false;
	}
}
//----------------------------------------------------------

//define this as a plug-in
DEFINE_FWK_MODULE(multilep);
