#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

//include c++ library classes
#include <algorithm>

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{
    std::string jecFile;
    if(multilepAnalyzer->is2018)      jecFile = "jecUncertaintyFile17"; // TODO: update when 2018 JEC become available
    else if(multilepAnalyzer->is2017) jecFile = "jecUncertaintyFile17";
    else                              jecFile = "jecUncertaintyFile16";

    jecUnc = new JetCorrectionUncertainty((iConfig.getParameter<edm::FileInPath>(jecFile)).fullPath());
};

JetAnalyzer::~JetAnalyzer(){
    delete jecUnc;
}

// Note that here only the uncertainty is saved, the Up and Down variations still need to be calculated later, probably easier to do on python level
// Also b-tagging up/down is easier on python level (see https://github.com/GhentAnalysis/StopsDilepton/blob/leptonSelectionUpdate_80X/tools/python/btagEfficiency.py)
// Storing here jet id variables, id itself also to be implemented at python level https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016, or maybe still here as the loose wp does not change often?
// WARNING: the _nJets is number of stored jets (i.e. including those where JECUp/JERUp passes the cut), do not use as selection
void JetAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/i");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetSmearedPt",              &_jetSmearedPt,             "_jetSmearedPt[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECDown",      &_jetSmearedPt_JECDown,     "_jetSmearedPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECUp",        &_jetSmearedPt_JECUp,       "_jetSmearedPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERDown",      &_jetSmearedPt_JERDown,     "_jetSmearedPt_JERDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERUp",        &_jetSmearedPt_JERUp,       "_jetSmearedPt_JERUp[_nJets]/D");
    outputTree->Branch("_jetPt_Uncorrected",         &_jetPt_Uncorrected,        "_jetPt_Uncorrected[_nJets]/D");
    outputTree->Branch("_jetPt_L1",                  &_jetPt_L1,                 "_jetPt_L1[_nJets]/D");
    outputTree->Branch("_jetPt_L2",                  &_jetPt_L2,                 "_jetPt_L2[_nJets]/D");
    outputTree->Branch("_jetPt_L3",                  &_jetPt_L3,                 "_jetPt_L3[_nJets]/D");

    outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetIsLoose",                &_jetIsLoose,               "_jetIsLoose[_nJets]/O");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");

    outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
    outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
    outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
    outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
    outputTree->Branch("_jetHFHadronFraction",       &_jetHFHadronFraction,      "_jetHFHadronFraction[_nJets]/D");
    outputTree->Branch("_jetHFEmFraction",           &_jetHFEmFraction,          "_jetHFEmFraction[_nJets]/D");

    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metRaw",                       &_metRaw,                       "_metRaw/D");
    outputTree->Branch("_metJECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_metJECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");

    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metRawPhi",                    &_metRawPhi,                    "_metRawPhi/D");
    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");

}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Jet>> jets;            iEvent.getByToken(multilepAnalyzer->jetToken,            jets);
    edm::Handle<std::vector<pat::Jet>> jetsSmeared;     iEvent.getByToken(multilepAnalyzer->jetSmearedToken,     jetsSmeared);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedUp;   iEvent.getByToken(multilepAnalyzer->jetSmearedUpToken,   jetsSmearedUp);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedDown; iEvent.getByToken(multilepAnalyzer->jetSmearedDownToken, jetsSmearedDown);
    edm::Handle<std::vector<pat::MET>> mets;            iEvent.getByToken(multilepAnalyzer->metToken, mets);

    //to apply JEC from txt files
    edm::Handle<double> rho;                            iEvent.getByToken(multilepAnalyzer->rhoToken,            rho);

    _nJets = 0;


    for(const auto& jet : *jets){
        if(_nJets == nJets_max) break;

        //only store loose jets
        _jetIsLoose[_nJets]        = jetIsLoose(jet, multilepAnalyzer->is2017 || multilepAnalyzer->is2018);
        if(!_jetIsLoose[_nJets]) continue;
        _jetIsTight[_nJets]        = jetIsTight(jet, multilepAnalyzer->is2017 || multilepAnalyzer->is2018);
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jet, multilepAnalyzer->is2017 || multilepAnalyzer->is2018);

        //find smeared equivalents of nominal jet
        auto jetSmearedIt = jetsSmeared->begin();
        for(auto j = jetsSmeared->cbegin(); j != jetsSmeared->cend(); ++j){
            if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedIt)) jetSmearedIt = j;
        }

        auto jetSmearedUpIt = jetsSmearedUp->begin();
        for(auto j = jetsSmearedUp->cbegin(); j != jetsSmearedUp->cend(); ++j){
            if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedUpIt))  jetSmearedUpIt = j;
        }

        auto jetSmearedDownIt = jetsSmearedDown->begin();
        for(auto j = jetsSmearedDown->cbegin(); j != jetsSmearedDown->cend(); ++j){
            if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedDownIt))  jetSmearedDownIt = j;
        }

        //nominal jet pt and uncertainties
        jecUnc->setJetEta(jet.eta());
        jecUnc->setJetPt(jet.pt());
        double unc = jecUnc->getUncertainty(true);

        _jetPt[_nJets]         = jet.pt();
        _jetPt_JECDown[_nJets] = _jetPt[_nJets]*(1 - unc);
        _jetPt_JECUp[_nJets]   = _jetPt[_nJets]*(1 + unc);

        //smeared jet pt and uncertainties
        jecUnc->setJetEta(jetSmearedIt->eta());
        jecUnc->setJetPt(jetSmearedIt->pt());
        double uncSmeared = jecUnc->getUncertainty(true);
        _jetSmearedPt[_nJets]         = jetSmearedIt->pt();
        _jetSmearedPt_JECDown[_nJets] = _jetPt[_nJets]*( 1. - uncSmeared );
        _jetSmearedPt_JECUp[_nJets]   = _jetPt[_nJets]*( 1. + uncSmeared );
        _jetSmearedPt_JERDown[_nJets] = jetSmearedDownIt->pt();
        _jetSmearedPt_JERUp[_nJets]   = jetSmearedUpIt->pt();

        //find maximum of all pT variations
        std::vector<double> ptVector = {_jetPt[_nJets], _jetPt_JECDown[_nJets], _jetPt_JECUp[_nJets],
            _jetSmearedPt[_nJets], _jetSmearedPt_JECDown[_nJets], _jetSmearedPt_JECUp[_nJets], _jetSmearedPt_JERDown[_nJets],  _jetSmearedPt_JERUp[_nJets]};
        double maxpT = *(std::max_element(ptVector.cbegin(), ptVector.cend()));
        if(maxpT <= 25) continue;

        _jetPt_Uncorrected[_nJets]        = jet.correctedP4("Uncorrected").Pt();
        _jetPt_L1[_nJets]                 = jet.correctedP4("L1FastJet").Pt();
        _jetPt_L2[_nJets]                 = jet.correctedP4("L2Relative").Pt();
        _jetPt_L3[_nJets]                 = jet.correctedP4("L3Absolute").Pt();

        _jetEta[_nJets]                   = jet.eta();
        _jetPhi[_nJets]                   = jet.phi();
        _jetE[_nJets]                     = jet.energy();

        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        //new DeepFlavour tagger
        _jetDeepCsv_udsg[_nJets]          = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _jetHadronFlavor[_nJets]          = jet.hadronFlavour();

        _jetNeutralHadronFraction[_nJets] = jet.neutralHadronEnergyFraction();
        _jetChargedHadronFraction[_nJets] = jet.chargedHadronEnergyFraction();
        _jetNeutralEmFraction[_nJets]     = jet.neutralEmEnergyFraction();
        _jetChargedEmFraction[_nJets]     = jet.chargedEmEnergyFraction();
        _jetHFHadronFraction[_nJets]      = jet.HFHadronEnergyFraction();
        _jetHFEmFraction[_nJets]          = jet.HFEMEnergyFraction();

        ++_nJets;
    }

    //determine the met of the event and its uncertainties
    //nominal MET value
    const pat::MET& met = (*mets).front();
    _met             = met.pt();
    _metPhi          = met.phi();

    //raw met values
    _metRaw          = met.uncorPt();
    _metRawPhi       = met.uncorPhi();
    //met values with uncertainties varied up and down
    _metJECDown      = met.shiftedPt(pat::MET::JetEnDown);
    _metJECUp        = met.shiftedPt(pat::MET::JetEnUp);
    _metUnclDown     = met.shiftedPt(pat::MET::UnclusteredEnDown);
    _metUnclUp       = met.shiftedPt(pat::MET::UnclusteredEnUp);
    _metPhiJECDown   = met.shiftedPhi(pat::MET::JetEnDown);
    _metPhiJECUp     = met.shiftedPhi(pat::MET::JetEnUp);
    _metPhiUnclUp    = met.shiftedPhi(pat::MET::UnclusteredEnUp);
    _metPhiUnclDown  = met.shiftedPhi(pat::MET::UnclusteredEnDown);

    //significance of met
    //note: this is the only one variable which changed between 94X and 102X see https://github.com/cms-sw/cmssw/commit/f7aacfd2ffaac9899ea07d0355afe49bb10a0aeb
    _metSignificance = met.metSignificance();

    if(multilepAnalyzer->skim == "singlejet" and _nJets < 1) return false;
//  if(multilepAnalyzer->skim == "FR" and _nJets < 1)        return false;              // this is in the master branch
    if(multilepAnalyzer->skim == "FR" and   (_nJets < 1 || _jetPt[0] < 30)) return false; // the _jetPt[0] is only ok if you do not use systematics for this skim
    if(multilepAnalyzer->skim == "FRsM" and (_nJets < 1 || _jetPt[0] < 30)) return false; // the _jetPt[0] is only ok if you do not use systematics for this skim

    return true;
}

/*
 * JetID implementations, references:
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
 * 2018, currently takes the same as 2017
 */

// There should be no loose Jet ID for 2017, not sure where the cuts for this below originate, so use at own risk
bool JetAnalyzer::jetIsLoose(const pat::Jet& jet, const bool is2017) const{
    if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >=  0.99)               return false;
        if(jet.neutralEmEnergyFraction() >= 0.99)                    return false;
        if(jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1) return false;
        if(fabs(jet.eta()) <= 2.4){
            if(jet.chargedHadronEnergyFraction() <= 0)               return false;
            if(jet.chargedMultiplicity() <= 0)                       return false;
            if(!is2017 && jet.chargedEmEnergyFraction()>= 0.99)      return false;
        }

    } else if(fabs(jet.eta()) <= 3.0){
        if(jet.neutralHadronEnergyFraction() >= 0.98)                return false;
        if(jet.neutralEmEnergyFraction() <= (is2017 ? 0.02 : 0.01))  return false;
        if(jet.neutralMultiplicity() <= 2)                           return false;

    } else {
        if(jet.neutralEmEnergyFraction() >= 0.90)                    return false;
        if(jet.neutralMultiplicity() <= 10)                          return false;
        if(is2017 && jet.neutralHadronEnergyFraction() <= 0.02)      return false;
    }

    return true;
}

bool JetAnalyzer::jetIsTight(const pat::Jet& jet, const bool is2017) const{
    if(!jetIsLoose(jet, is2017))                            return false;

    if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >= 0.9)        return false;
        if(jet.neutralEmEnergyFraction() >= 0.9)            return false;

    } else if(fabs(jet.eta()) <= 3.0){
        if(is2017 && jet.neutralEmEnergyFraction() >= 0.99) return false;
    }
    return true;
}

bool JetAnalyzer::jetIsTightLepVeto(const pat::Jet& jet, const bool is2017) const{
    if(!jetIsTight(jet, is2017))                                      return false;

    if(fabs(jet.eta()) <= 2.7){
        if(jet.chargedMuEnergyFraction() >= 0.8 )                     return false;

        if(fabs(jet.eta()) <= 2.4){
            if(jet.chargedEmEnergyFraction() >= (is2017 ? 0.8 : 0.9)) return false;
        }
    }
    return true;
}
