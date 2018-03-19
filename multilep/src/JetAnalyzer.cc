#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    jecUnc((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath()),
    multilepAnalyzer(multilepAnalyzer)
{
    if(multilepAnalyzer->isData){
        jecLevel = "L2L3Residual";
    } else {
        jecLevel = "L3Absolute";
    }
};

// Note that here only the uncertainty is saved, the Up and Down variations still need to be calculated later, probably easier to do on python level
// Also b-tagging up/down is easier on python level (see https://github.com/GhentAnalysis/StopsDilepton/blob/leptonSelectionUpdate_80X/tools/python/btagEfficiency.py)
// Storing here jet id variables, id itself also to be implemented at python level https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016, or maybe still here as the loose wp does not change often?
// WARNING: the _nJets is number of stored jets (i.e. including those where JECUp/JERUp passes the cut), do not use as selection
void JetAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/b");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetIsLoose",                &_jetIsLoose,               "_jetIsLoose[_nJets]/i");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/i");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/i");

    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metJECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_metJECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_metJetResDown",                &_metJetResDown,                "_metJetResDown/D");
    outputTree->Branch("_metJetResUp",                  &_metJetResUp,                  "_metJetResUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");

    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhiJetResDown",             &_metPhiJetResDown,             "_metPhiJetResDown/D");
    outputTree->Branch("_metPhiJetResUp",               &_metPhiJetResUp,               "_metPhiJetResUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");

    outputTree->Branch("_rawmetPhi",                    &_rawmetPhi,                    "_rawmetPhi/D");
    outputTree->Branch("_rawmetPhiJECDown",             &_rawmetPhiJECDown,             "_rawmetPhiJECDown/D");
    outputTree->Branch("_rawmetPhiJECUp",               &_rawmetPhiJECUp,               "_rawmetPhiJECUp/D");
    outputTree->Branch("_rawmetPhiJetResDown",          &_rawmetPhiJetResDown,          "_rawmetPhiJetResDown/D");
    outputTree->Branch("_rawmetPhiJetResUp",            &_rawmetPhiJetResUp,            "_rawmetPhiJetResUp/D");
    outputTree->Branch("_rawmetPhiUnclDown",            &_rawmetPhiUnclDown,            "_rawmetPhiUnclDown/D");
    outputTree->Branch("_rawmetPhiUnclUp",              &_rawmetPhiUnclUp,              "_rawmetPhiUnclUp/D");

    outputTree->Branch("_rawmet",                       &_rawmet,                       "_rawmet/D");
    outputTree->Branch("_rawmetJECDown",                &_rawmetJECDown,                "_rawmetJECDown/D");
    outputTree->Branch("_rawmetJECUp",                  &_rawmetJECUp,                  "_rawmetJECUp/D");
    outputTree->Branch("_rawmetJetResDown",             &_rawmetJetResDown,             "_rawmetJetResDown/D");
    outputTree->Branch("_rawmetJetResUp",               &_rawmetJetResUp,               "_rawmetJetResUp/D");
    outputTree->Branch("_rawmetUnclDown",               &_rawmetUnclDown,               "_rawmetUnclDown/D");
    outputTree->Branch("_rawmetUnclUp",                 &_rawmetUnclUp,                 "_rawmetUnclUp/D");
    
    outputTree->Branch("_met_sm",                       &_met_sm,                       "_met_sm/D");
    outputTree->Branch("_metJECDown_sm",                &_metJECDown_sm,                "_metJECDown_sm/D");
    outputTree->Branch("_metJECUp_sm",                  &_metJECUp_sm,                  "_metJECUp_sm/D");
    outputTree->Branch("_metJetResDown_sm",             &_metJetResDown_sm,             "_metJetResDown_sm/D");
    outputTree->Branch("_metJetResUp_sm",               &_metJetResUp_sm,               "_metJetResUp_sm/D");
    outputTree->Branch("_metUnclDown_sm",               &_metUnclDown_sm,               "_metUnclDown_sm/D");
    outputTree->Branch("_metUnclUp_sm",                 &_metUnclUp_sm,                 "_metUnclUp_sm/D");

    outputTree->Branch("_metPhi_sm",                    &_metPhi_sm,                    "_metPhi_sm/D");
    outputTree->Branch("_metPhiJECDown_sm",             &_metPhiJECDown_sm,             "_metPhiJECDown_sm/D");
    outputTree->Branch("_metPhiJECUp_sm",               &_metPhiJECUp_sm,               "_metPhiJECUp_sm/D");
    outputTree->Branch("_metPhiJetResDown_sm",          &_metPhiJetResDown_sm,          "_metPhiJetResDown_sm/D");
    outputTree->Branch("_metPhiJetResUp_sm",            &_metPhiJetResUp_sm,            "_metPhiJetResUp_sm/D");
    outputTree->Branch("_metPhiUnclDown_sm",            &_metPhiUnclDown_sm,            "_metPhiUnclDown_sm/D");
    outputTree->Branch("_metPhiUnclUp_sm",              &_metPhiUnclUp_sm,              "_metPhiUnclUp_sm/D");



}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Jet>> jets;            iEvent.getByToken(multilepAnalyzer->jetToken,            jets);
    edm::Handle<std::vector<pat::MET>> mets;            iEvent.getByToken(multilepAnalyzer->metToken, mets);
    //to apply JEC from txt files
    edm::Handle<double> rho;                            iEvent.getByToken(multilepAnalyzer->rhoToken,            rho);

    _nJets = 0;

    //for(auto jetSmeared = jetsSmeared->begin(); jetSmeared != jetsSmeared->end(); ++jetSmeared){
    for(auto& jet : *jets){
        if(_nJets == nJets_max) break;

        //only store loose jets 
        _jetIsLoose[_nJets] = jetIsLoose(jet, multilepAnalyzer->is2017);
        if(!_jetIsLoose[_nJets]) continue;
        _jetIsTight[_nJets] = jetIsTight(jet, multilepAnalyzer->is2017);
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jet, multilepAnalyzer->is2017);
    
        /*
        jecUnc.setJetEta(jet.eta());
        jecUnc.setJetPt(jet.pt());
        double unc = jecUnc.getUncertainty(true);
        */

        //txt based JEC
        double corrTxt = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), *rho, jet.jetArea(), jecLevel);
        _jetPt[_nJets] = jet.correctedP4("Uncorrected").Pt()*corrTxt;
        /*
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "txt based JetPt = " << _jetPt[_nJets] << std::endl;
        std::cout << "CMSSSW JetPt = " << jet.pt() << std::endl;
        */
        //extract uncertainty based on corrected jet pt!
        double unc = multilepAnalyzer->jec->jetUncertainty(_jetPt[_nJets], jet.eta());
        double maxpT = (1+unc)*_jetPt[_nJets];
        if(maxpT <= 25) continue;

        _jetPt_JECDown[_nJets]            = _jetPt[_nJets]*(1-unc);
        _jetPt_JECUp[_nJets]              = _jetPt[_nJets]*(1+unc);
        _jetEta[_nJets]                   = jet.eta();
        _jetPhi[_nJets]                   = jet.phi();
        _jetE[_nJets]                     = jet.correctedP4("Uncorrected").E()*corrTxt;
        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        //new DeepFlavour tagger
        _jetDeepCsv_udsg[_nJets]          = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _jetHadronFlavor[_nJets]          = jet.hadronFlavour();

        ++_nJets;
    }

    //determine the met of the event and its uncertainties
    //nominal MET value
    const pat::MET& met = (*mets).front();
    _met             = met.pt();
    _metPhi          = met.phi();
    /*
    std::cout << "~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "met = " << _met << std::endl;
    std::cout << "txt corrected met = " << multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho).first << std::endl;
    */
    //_met = multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho).first;
    //_metPhi = multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho).second;

    //raw met values
    _rawmet= met.uncorPt();
    _rawmetJECDown   = met.shiftedPt(pat::MET::JetEnDown, pat::MET::Raw);
    _rawmetJECUp     = met.shiftedPt(pat::MET::JetEnUp, pat::MET::Raw);
    _rawmetJetResDown= met.shiftedPt(pat::MET::JetResDown, pat::MET::Raw);
    _rawmetJetResUp  = met.shiftedPt(pat::MET::JetResUp, pat::MET::Raw);
    _rawmetUnclDown  = met.shiftedPt(pat::MET::UnclusteredEnDown, pat::MET::Raw);
    _rawmetUnclUp    = met.shiftedPt(pat::MET::UnclusteredEnUp, pat::MET::Raw);

    _rawmetPhi= met.uncorPhi();
    _rawmetPhiJECDown   = met.shiftedPhi(pat::MET::JetEnDown, pat::MET::Raw);
    _rawmetPhiJECUp     = met.shiftedPhi(pat::MET::JetEnUp, pat::MET::Raw);
    _rawmetPhiJetResDown= met.shiftedPhi(pat::MET::JetResDown, pat::MET::Raw);
    _rawmetPhiJetResUp  = met.shiftedPhi(pat::MET::JetResUp, pat::MET::Raw);
    _rawmetPhiUnclDown  = met.shiftedPhi(pat::MET::UnclusteredEnDown, pat::MET::Raw);
    _rawmetPhiUnclUp    = met.shiftedPhi(pat::MET::UnclusteredEnUp, pat::MET::Raw);
    //met values with uncertainties varied up and down
    _metJECDown      = met.shiftedPt(pat::MET::JetEnDown);
    _metJECUp        = met.shiftedPt(pat::MET::JetEnUp);
    _metJetResDown   = met.shiftedPt(pat::MET::JetResDown);
    _metJetResUp     = met.shiftedPt(pat::MET::JetResUp);
    _metUnclDown     = met.shiftedPt(pat::MET::UnclusteredEnDown);
    _metUnclUp       = met.shiftedPt(pat::MET::UnclusteredEnUp);

    _metPhiJECDown   = met.shiftedPhi(pat::MET::JetEnDown);
    _metPhiJECUp     = met.shiftedPhi(pat::MET::JetEnUp);
    _metPhiJetResDown= met.shiftedPhi(pat::MET::JetResDown);
    _metPhiJetResUp  = met.shiftedPhi(pat::MET::JetResUp);
    _metPhiUnclDown  = met.shiftedPhi(pat::MET::UnclusteredEnDown);
    _metPhiUnclUp    = met.shiftedPhi(pat::MET::UnclusteredEnUp);
    //smeared met values
    _met_sm          = met.corPt(pat::MET::Type1Smear);
    _metJECDown_sm   = met.shiftedPt(pat::MET::JetEnDown, pat::MET::Type1Smear);
    _metJECUp_sm     = met.shiftedPt(pat::MET::JetEnUp, pat::MET::Type1Smear);
    _metJetResDown_sm= met.shiftedPt(pat::MET::JetResDown, pat::MET::Type1Smear);
    _metJetResUp_sm  = met.shiftedPt(pat::MET::JetResUp, pat::MET::Type1Smear);
    _metUnclDown_sm  = met.shiftedPt(pat::MET::UnclusteredEnDown, pat::MET::Type1Smear);
    _metUnclUp_sm    = met.shiftedPt(pat::MET::UnclusteredEnUp, pat::MET::Type1Smear);

    _metPhi_sm          = met.corPhi(pat::MET::Type1Smear);
    _metPhiJECDown_sm   = met.shiftedPhi(pat::MET::JetEnDown, pat::MET::Type1Smear);
    _metPhiJECUp_sm     = met.shiftedPhi(pat::MET::JetEnUp, pat::MET::Type1Smear);
    _metPhiJetResDown_sm= met.shiftedPhi(pat::MET::JetResDown, pat::MET::Type1Smear);
    _metPhiJetResUp_sm  = met.shiftedPhi(pat::MET::JetResUp, pat::MET::Type1Smear);
    _metPhiUnclDown_sm  = met.shiftedPhi(pat::MET::UnclusteredEnDown, pat::MET::Type1Smear);
    _metPhiUnclUp_sm    = met.shiftedPhi(pat::MET::UnclusteredEnUp, pat::MET::Type1Smear);

    /*
    std::cout << "~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "met = " << _met << " " << _metJECUp << " " << _metJECDown << " " << _metJetResUp << " " << _metJetResDown  << " " << _metUnclUp << " " << _metUnclDown << std::endl;
    std::cout << "raw met = " << _rawmet << " " << _rawmetJECUp << " " << _rawmetJECDown  << " " << _rawmetJetResUp << " " << _rawmetJetResDown << " " << _rawmetUnclUp << " " << _rawmetUnclDown << std::endl;
    std::cout << "smeared met = " << _met_sm << " " << _metJECUp_sm << " " << _metJECDown_sm  << " " << _metJetResUp_sm << " " << _metJetResDown_sm  << " " << _metUnclUp_sm << " " << _metUnclDown_sm << std::endl;
    */
    //std::cout << "txt corrected met = " << multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho).first << std::endl;
    //significance of met
    _metSignificance = met.metSignificance(); 


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
