//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//include ROOT classes
#include "TLorentzVector.h"

#include "heavyNeutrino/multilep/interface/BFragAnalyzer.h"

/*
Class storing data for b fragmentation systematics
https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
*/

BFragAnalyzer::BFragAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer){};

void BFragAnalyzer::beginJob(TTree* outputTree){
   
    if( multilepAnalyzer->isData() ) return;
   
    outputTree->Branch("_bf_fragCP5BL",                 &_bf_fragCP5BL,                 "_bf_fragCP5BL/D");
    outputTree->Branch("_bf_fragCP5BLdown",             &_bf_fragCP5BLdown,             "_bf_fragCP5BLdown/D");
    outputTree->Branch("_bf_fragCP5BLup",               &_bf_fragCP5BLup,               "_bf_fragCP5BLup/D");
    outputTree->Branch("_bf_fragCP5Peterson",           &_bf_fragCP5Peterson,           "_bf_fragCP5Peterson/D");
    outputTree->Branch("_bf_fragCP5Petersondown",       &_bf_fragCP5Petersondown,       "_bf_fragCP5Petersondown/D");
    outputTree->Branch("_bf_fragCP5Petersonup",         &_bf_fragCP5Petersonup,         "_bf_fragCP5Petersonup/D");
}

void BFragAnalyzer::analyze(const edm::Event& iEvent){
   
    if( multilepAnalyzer->isData() ) return;
   
    edm::Handle<std::vector<reco::GenJet>> genJets           = getHandle(iEvent, multilepAnalyzer->genJetsToken);
    edm::Handle<edm::ValueMap<float> > fragCP5BL             = getHandle(iEvent, multilepAnalyzer->fragCP5BLToken);
    edm::Handle<edm::ValueMap<float> > fragCP5BLdown         = getHandle(iEvent, multilepAnalyzer->fragCP5BLdownToken);
    edm::Handle<edm::ValueMap<float> > fragCP5BLup           = getHandle(iEvent, multilepAnalyzer->fragCP5BLupToken);
    edm::Handle<edm::ValueMap<float> > fragCP5Peterson       = getHandle(iEvent, multilepAnalyzer->fragCP5PetersonToken);
    edm::Handle<edm::ValueMap<float> > fragCP5Petersondown   = getHandle(iEvent, multilepAnalyzer->fragCP5PetersondownToken);
    edm::Handle<edm::ValueMap<float> > fragCP5Petersonup     = getHandle(iEvent, multilepAnalyzer->fragCP5PetersonupToken);
   
    double wfragCP5BL = 1.;
    double wfragCP5BLdown = 1.;
    double wfragCP5BLup = 1.;
    double wfragCP5Peterson = 1.;
    double wfragCP5Petersondown = 1.;
    double wfragCP5Petersonup = 1.;   
   
    for (auto genJet=genJets->begin(); genJet!=genJets->end(); ++genJet)
     {
	edm::Ref<std::vector<reco::GenJet> > genJetRef(genJets, genJet-genJets->begin());
	
	wfragCP5BL *= (*fragCP5BL)[genJetRef];
	wfragCP5BLdown *= (*fragCP5BLdown)[genJetRef];
	wfragCP5BLup *= (*fragCP5BLup)[genJetRef];
	wfragCP5Peterson *= (*fragCP5Peterson)[genJetRef];
	wfragCP5Petersondown *= (*fragCP5Petersondown)[genJetRef];
	wfragCP5Petersonup *= (*fragCP5Petersonup)[genJetRef];
     }

    _bf_fragCP5BL = wfragCP5BL;
    _bf_fragCP5BLdown = wfragCP5BLdown;
    _bf_fragCP5BLup = wfragCP5BLup;
    _bf_fragCP5Peterson = wfragCP5Peterson;
    _bf_fragCP5Petersondown = wfragCP5Petersondown;
    _bf_fragCP5Petersonup = wfragCP5Petersonup;
}
