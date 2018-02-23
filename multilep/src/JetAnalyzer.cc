#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

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
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
    //outputTree->Branch("_jetDeepCsv_cc",             &_jetDeepCsv_cc,            "_jetDeepCsv_cc[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    //outputTree->Branch("_jetId",                     &_jetId,                    "_jetId[_nJets]/i");
    outputTree->Branch("_jetIsLoose",                &_jetIsLoose,               "_jetIsLoose[_nJets]/i");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/i");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/i");
}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Jet>> jets;            iEvent.getByToken(multilepAnalyzer->jetToken,            jets);
    edm::Handle<std::vector<pat::Jet>> jetsSmeared;     iEvent.getByToken(multilepAnalyzer->jetSmearedToken,     jetsSmeared);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedUp;   iEvent.getByToken(multilepAnalyzer->jetSmearedUpToken,   jetsSmearedUp);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedDown; iEvent.getByToken(multilepAnalyzer->jetSmearedDownToken, jetsSmearedDown);

    _nJets = 0;

    //for(auto jetSmeared = jetsSmeared->begin(); jetSmeared != jetsSmeared->end(); ++jetSmeared){
    for(auto& jetSmeared : *jetsSmeared){
        if(_nJets == nJets_max) break;

        //only store loose jets 
        //_jetId[_nJets] = jetId(*jetSmeared, false) + jetId(*jetSmeared, true); // 1: loose, 2: tight
        _jetIsLoose[_nJets] = jetIsLoose(jetSmeared, multilepAnalyzer->is2017);
        //if(!_jetIsLoose[_nJets]) continue;
        _jetIsTight[_nJets] = jetIsTight(jetSmeared, multilepAnalyzer->is2017);
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jetSmeared, multilepAnalyzer->is2017);

        jecUnc.setJetEta(jetSmeared.eta());
        jecUnc.setJetPt(jetSmeared.pt());
        double unc = jecUnc.getUncertainty(true);

        auto jetSmearedUp = jetsSmearedUp->begin();
        for(auto j = jetsSmearedUp->begin() + 1; j != jetsSmearedUp->end(); ++j){
          if(reco::deltaR(jetSmeared, *j) < reco::deltaR(jetSmeared, *jetSmearedUp)) jetSmearedUp = j;
        }

        auto jetSmearedDown = jetsSmearedDown->begin();
        for(auto j = jetsSmearedDown->begin() + 1; j != jetsSmearedDown->end(); ++j){
          if(reco::deltaR(jetSmeared, *j) < reco::deltaR(jetSmeared, *jetSmearedDown)) jetSmearedDown = j;
        }

        double maxpT = std::max( (1+unc)*jetSmeared.pt() ,  std::max(jetSmearedUp->pt(), jetSmearedDown->pt() ) );
        if(maxpT <= 25) continue;
        //if(std::max((1+unc)*jetSmeared->pt(), std::max(jetSmearedUp->pt(), jetSmearedDown->pt())) < 25) continue;

        _jetPt[_nJets]                    = jetSmeared.pt();
        _jetPt_JECDown[_nJets]            = jetSmeared.pt()*(1-unc);
        _jetPt_JECUp[_nJets]              = jetSmeared.pt()*(1+unc);
        _jetPt_JERDown[_nJets]            = jetSmearedDown->pt();
        _jetPt_JERUp[_nJets]              = jetSmearedUp->pt();
        _jetEta[_nJets]                   = jetSmeared.eta();
        _jetPhi[_nJets]                   = jetSmeared.phi();
        _jetE[_nJets]                     = jetSmeared.energy();
        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jetSmeared.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        //new DeepFlavour tagger
        _jetDeepCsv_udsg[_nJets]          = jetSmeared.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jetSmeared.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jetSmeared.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jetSmeared.bDiscriminator("pfDeepCSVJetTags:probbb");
    //  _jetDeepCsv_cc[_nJets]            = jetSmeared->bDiscriminator("pfDeepCSVJetTags:probcc");
        _jetHadronFlavor[_nJets]          = jetSmeared.hadronFlavour();

        ++_nJets;
    }
    if(multilepAnalyzer->skim == "singlejet" and _nJets < 1) return false;
    if(multilepAnalyzer->skim == "FR" and _nJets < 1) return false;
    return true;
}

bool JetAnalyzer::jetIsLoose(const pat::Jet& jet, const bool is2017) const{

    if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >=  0.99) return false;
        if(jet.neutralEmEnergyFraction() >= 0.99) return false;
        if(jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1) return false;
        if(fabs(jet.eta()) <= 2.4){
            if(jet.chargedHadronEnergyFraction() <= 0) return false;
            if(jet.chargedMultiplicity() <= 0) return false;
            if( !is2017 && ( jet.chargedEmEnergyFraction() >= 0.99 ) ) return false;
        }

    } else if(fabs(jet.eta()) <= 3.0){
        if(jet.neutralHadronEnergyFraction() >= 0.98) return false;
        if( !is2017 && (jet.neutralEmEnergyFraction() <= 0.01) ) return false;
        if( is2017 && (jet.neutralEmEnergyFraction() <= 0.02) ) return false;
        if(jet.neutralMultiplicity() <= 2) return false;

    } else {
        if(jet.neutralEmEnergyFraction() >= 0.90) return false;
        if(jet.neutralMultiplicity() <= 10) return false;
        if(is2017 && jet.neutralHadronEnergyFraction() <= 0.02) return false;
    }

    return true;
}

bool JetAnalyzer::jetIsTight(const pat::Jet& jet, const bool is2017) const{
    if( !jetIsLoose(jet, is2017) ) return false;

    if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >= 0.9) return false;
        if(jet.neutralEmEnergyFraction() >= 0.9) return false;

    } else if(fabs(jet.eta()) <= 3.0){
        if(is2017 && jet.neutralEmEnergyFraction() >= 0.99) return false;   
    }
    return true;
}

bool JetAnalyzer::jetIsTightLepVeto(const pat::Jet& jet, const bool is2017) const{
    if( !jetIsTight(jet, is2017) ) return false;
    if(fabs(jet.eta()) <= 2.7){
        if( jet.chargedMuEnergyFraction() >= 0.8 ) return false;
        
        if(fabs(jet.eta()) <= 2.4){
            if(is2017 && (jet.chargedEmEnergyFraction() >= 0.8) ) return false;
            if(!is2017 && (jet.chargedEmEnergyFraction() >= 0.9) ) return false;
        }
    }
    return true;
}
