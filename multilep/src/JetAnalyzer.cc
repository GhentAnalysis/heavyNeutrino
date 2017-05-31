#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  jecUnc((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath()),
  multilepAnalyzer(multilepAnalyzer)
{};

// Note that here only the uncertainty is saved, the Up and Down variations still need to be calculated later, probably easier to do on python level
// Also b-tagging up/down is easier on python level (see https://github.com/GhentAnalysis/StopsDilepton/blob/leptonSelectionUpdate_80X/tools/python/btagEfficiency.py)
// Storing here jet id variables, id itself also to be implemented at python level https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016, or maybe still here as the loose wp does not change often?
// WARNING: the _nJets is number of stored jets (i.e. including those where JECUp/JERUp passes the cut), do not use as selection
void JetAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/b");
  outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
  outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
  outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
  outputTree->Branch("_jetPt_JERUp",               &_jetPt_JERUp,              "_jetPt_JERUp[_nJets]/D");
  outputTree->Branch("_jetPt_JERDown",             &_jetPt_JERDown,            "_jetPt_JERDown[_nJets]/D");
  outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
  outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
  outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
  outputTree->Branch("_jetBTaggingCSV",            &_jetBTaggingCSV,           "_jetBTaggingCSV[_nJets]/D");
  outputTree->Branch("_jetHadronFlavour",          &_jetHadronFlavour,         "_jetHadronFlavour[_nJets]/D");
  outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
  outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
  outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
  outputTree->Branch("_jetMuonFraction",           &_jetMuonFraction,          "_jetMuonFraction[_nJets]/D");
  outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
  outputTree->Branch("_jetNeutralMult",            &_jetNeutralMult,           "_jetNeutralMult[_nJets]/D");
  outputTree->Branch("_jetChargedMult",            &_jetChargedMult,           "_jetChargedMult[_nJets]/D");
}

void JetAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<pat::Jet>> jets;            iEvent.getByToken(multilepAnalyzer->jetToken,            jets);
  edm::Handle<std::vector<pat::Jet>> jetsSmeared;     iEvent.getByToken(multilepAnalyzer->jetSmearedToken,     jetsSmeared);
  edm::Handle<std::vector<pat::Jet>> jetsSmearedUp;   iEvent.getByToken(multilepAnalyzer->jetSmearedUpToken,   jetsSmearedUp);
  edm::Handle<std::vector<pat::Jet>> jetsSmearedDown; iEvent.getByToken(multilepAnalyzer->jetSmearedDownToken, jetsSmearedDown);

  _nJets = 0;

  auto jet            = jets->begin();
  auto jetSmeared     = jetsSmeared->begin();
  auto jetSmearedUp   = jetsSmearedUp->begin();
  auto jetSmearedDown = jetsSmearedDown->begin();
  for(; jet != jets->end(); ++jet, ++jetSmeared, ++jetSmearedUp, ++jetSmearedDown){
    if(_nJets == nJets_max) break;
    jecUnc.setJetEta(jetSmeared->eta());
    jecUnc.setJetPt(jetSmeared->pt());
    double unc = jecUnc.getUncertainty(true);

    if(std::max((1+unc)*jetSmeared->pt(), std::max(jetSmearedUp->pt(), jetSmearedDown->pt())) < 25) continue;
    if(fabs(jet->eta()) > 2.4)                                                                      continue;

    _jetPt[_nJets]                    = jetSmeared->pt();
    _jetPt_JECDown[_nJets]            = jetSmeared->pt()*(1-unc);
    _jetPt_JECUp[_nJets]              = jetSmeared->pt()*(1+unc);
    _jetPt_JERDown[_nJets]            = jetSmearedDown->pt();
    _jetPt_JERUp[_nJets]              = jetSmearedUp->pt();
    _jetEta[_nJets]                   = jet->eta();
    _jetPhi[_nJets]                   = jet->phi();
    _jetE[_nJets]                     = jet->energy();
    _jetBTaggingCSV[_nJets]           = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    _jetHadronFlavour[_nJets]         = jet->hadronFlavour();                                                                     // Note: these variables need to come from the original collection, not the smeared one
    _jetNeutralHadronFraction[_nJets] = jet->neutralHadronEnergyFraction();
    _jetNeutralEmFraction[_nJets]     = jet->neutralEmEnergyFraction();
    _jetChargedHadronFraction[_nJets] = jet->chargedHadronEnergyFraction();
    _jetMuonFraction[_nJets]          = jet->muonEnergyFraction();
    _jetChargedEmFraction[_nJets]     = jet->chargedEmEnergyFraction();
    _jetNeutralMult[_nJets]           = jet->neutralMultiplicity();
    _jetChargedMult[_nJets]           = jet->chargedMultiplicity();

    ++_nJets;
  }
}
