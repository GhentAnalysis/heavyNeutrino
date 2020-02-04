#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

//include c++ library classes
#include <algorithm>

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{
    if(multilepAnalyzer->isData()){
        jecLevel = "L2L3residual";
    } else {
        jecLevel = "L3Absolute";
    }
    jecUnc = new JetCorrectionUncertainty((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath());
    jecUncPuppi = new JetCorrectionUncertainty((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFilePuppi")).fullPath());
};

JetAnalyzer::~JetAnalyzer(){
    delete jecUnc;
    delete jecUncPuppi;
}

// Note: the JEC are already applied through the GT, if you need back the old way (JEC.cc) check the code before c3564f71a2e7dca3cb963ef69481894cb04bbf53
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
    outputTree->Branch("_jetDeepCsv",                &_jetDeepCsv,               "_jetDeepCsv[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");

    outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
    outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
    outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
    outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
    outputTree->Branch("_jetHFHadronFraction",       &_jetHFHadronFraction,      "_jetHFHadronFraction[_nJets]/D");
    outputTree->Branch("_jetHFEmFraction",           &_jetHFEmFraction,          "_jetHFEmFraction[_nJets]/D");

    outputTree->Branch("_nJetsPuppi",                &_nJetsPuppi,               "_nJetsPuppi/i");
    outputTree->Branch("_jetPuppiPt",                &_jetPuppiPt,               "_jetPuppiPt[_nJetsPuppi]/D");
    outputTree->Branch("_jetPuppiEta",               &_jetPuppiEta,              "_jetPuppiEta[_nJetsPuppi]/D");
    outputTree->Branch("_jetPuppiPhi",               &_jetPuppiPhi,              "_jetPuppiPhi[_nJetsPuppi]/D");
    outputTree->Branch("_jetPuppiPt_JECDown",        &_jetPuppiPt_JECDown,       "_jetPuppiPt_JECDown[_nJetsPuppi]/D");
    outputTree->Branch("_jetPuppiPt_JECUp",          &_jetPuppiPt_JECUp,         "_jetPuppiPt_JECUp[_nJetsPuppi]/D");

    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metType1",                     &_metType1,                     "_metType1/D");
    outputTree->Branch("_metRaw",                       &_metRaw,                       "_metRaw/D");
    outputTree->Branch("_metJECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_metJECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");
    outputTree->Branch("_metResDown",                   &_metResDown,                   "_metResDown/D");
    outputTree->Branch("_metResUp",                     &_metResUp,                     "_metResUp/D");

    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metType1Phi",                  &_metType1Phi,                  "_metType1Phi/D");
    outputTree->Branch("_metRawPhi",                    &_metRawPhi,                    "_metRawPhi/D");
    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    outputTree->Branch("_metPhiResDown",                &_metPhiResDown,                "_metPhiResDown/D");
    outputTree->Branch("_metPhiResUp",                  &_metPhiResUp,                  "_metPhiResUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");

    outputTree->Branch("_metPuppi",                          &_metPuppi,                          "_metPuppi/D");
    outputTree->Branch("_metPuppiRaw",                       &_metPuppiRaw,                       "_metPuppiRaw/D");
    outputTree->Branch("_metPuppiJECDown",                   &_metPuppiJECDown,                   "_metPuppiJECDown/D");
    outputTree->Branch("_metPuppiJECUp",                     &_metPuppiJECUp,                     "_metPuppiJECUp/D");
    outputTree->Branch("_metPuppiUnclDown",                  &_metPuppiUnclDown,                  "_metPuppiUnclDown/D");
    outputTree->Branch("_metPuppiUnclUp",                    &_metPuppiUnclUp,                    "_metPuppiUnclUp/D");
    outputTree->Branch("_metPuppiResDown",                   &_metPuppiResDown,                   "_metPuppiResDown/D");
    outputTree->Branch("_metPuppiResUp",                     &_metPuppiResUp,                     "_metPuppiResUp/D");

    outputTree->Branch("_metPuppiPhi",                       &_metPuppiPhi,                       "_metPuppiPhi/D");
    outputTree->Branch("_metPuppiRawPhi",                    &_metPuppiRawPhi,                    "_metPuppiRawPhi/D");
    outputTree->Branch("_metPuppiPhiJECDown",                &_metPuppiPhiJECDown,                "_metPuppiPhiJECDown/D");
    outputTree->Branch("_metPuppiPhiJECUp",                  &_metPuppiPhiJECUp,                  "_metPuppiPhiJECUp/D");
    outputTree->Branch("_metPuppiPhiUnclDown",               &_metPuppiPhiUnclDown,               "_metPuppiPhiUnclDown/D");
    outputTree->Branch("_metPuppiPhiUnclUp",                 &_metPuppiPhiUnclUp,                 "_metPuppiPhiUnclUp/D");
    outputTree->Branch("_metPuppiPhiResDown",                &_metPuppiPhiResDown,                "_metPuppiPhiResDown/D");
    outputTree->Branch("_metPuppiPhiResUp",                  &_metPuppiPhiResUp,                  "_metPuppiPhiResUp/D");

    if(!multilepAnalyzer->is2018() ) outputTree->Branch("_jetIsLoose", _jetIsLoose, "_jetIsLoose[_nJets]/O"); // WARNING, not recommended to be used, only exists for 2016
}

bool JetAnalyzer::analyze(const edm::Event& iEvent, const int npv){
    edm::Handle<std::vector<pat::Jet>> jets            = getHandle(iEvent, multilepAnalyzer->jetToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmeared     = getHandle(iEvent, multilepAnalyzer->jetSmearedToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedUp   = getHandle(iEvent, multilepAnalyzer->jetSmearedUpToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedDown = getHandle(iEvent, multilepAnalyzer->jetSmearedDownToken);
    edm::Handle<std::vector<pat::Jet>> jetsPuppi       = getHandle(iEvent, multilepAnalyzer->jetPuppiToken);
    edm::Handle<std::vector<pat::MET>> mets            = getHandle(iEvent, multilepAnalyzer->metToken);
    edm::Handle<std::vector<pat::MET>> metsPuppi       = getHandle(iEvent, multilepAnalyzer->metPuppiToken);
    //to apply JEC from txt files
    edm::Handle<double> rho                            = getHandle(iEvent, multilepAnalyzer->rhoToken);

    _nJets = 0;


    for(const auto& jet : *jets){
        if(_nJets == nJets_max) break;

        _jetIsLoose[_nJets]        = jetIsLoose(jet, multilepAnalyzer->is2017() || multilepAnalyzer->is2018() );
        _jetIsTight[_nJets]        = jetIsTight(jet, multilepAnalyzer->is2017(), multilepAnalyzer->is2018() );
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jet, multilepAnalyzer->is2017(), multilepAnalyzer->is2018() );

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

        // JEC implementation based on txt, if needed for newer JEC versions
        //double corrTxt = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), *rho, jet.jetArea(), jecLevel);
        //_jetPt[_nJets] = jet.correctedP4("Uncorrected").Pt()*corrTxt;
        //double unc = multilepAnalyzer->jec->jetUncertainty(_jetPt[_nJets], _jetEta[_nJets]);
        //_jetPt_JECDown[_nJets] = _jetPt[_nJets]*(1 - unc);
        //_jetPt_JECUp[_nJets]   = _jetPt[_nJets]*(1 + unc);

        
        //nominal jet pt and uncertainties, through global tag
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
        _jetSmearedPt_JECDown[_nJets] = _jetSmearedPt[_nJets]*( 1. - uncSmeared );
        _jetSmearedPt_JECUp[_nJets]   = _jetSmearedPt[_nJets]*( 1. + uncSmeared );
        _jetSmearedPt_JERDown[_nJets] = jetSmearedDownIt->pt();
        _jetSmearedPt_JERUp[_nJets]   = jetSmearedUpIt->pt();

        //find maximum of all pT variations
        std::vector<double> ptVector = {_jetPt[_nJets], _jetPt_JECDown[_nJets], _jetPt_JECUp[_nJets],
            _jetSmearedPt[_nJets], _jetSmearedPt_JECDown[_nJets], _jetSmearedPt_JECUp[_nJets], _jetSmearedPt_JERDown[_nJets],  _jetSmearedPt_JERUp[_nJets]};
        double maxpT = *(std::max_element(ptVector.cbegin(), ptVector.cend()));
        if(maxpT <= 20) continue;

        _jetPt_Uncorrected[_nJets]        = jet.correctedP4("Uncorrected").Pt();
        _jetPt_L1[_nJets]                 = jet.correctedP4("L1FastJet").Pt();
        _jetPt_L2[_nJets]                 = jet.correctedP4("L2Relative").Pt();
        _jetPt_L3[_nJets]                 = jet.correctedP4("L3Absolute").Pt();

        _jetEta[_nJets]                   = jet.eta();
        _jetPhi[_nJets]                   = jet.phi();
        _jetE[_nJets]                     = jet.energy();

        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        //new DeepCSV tagger
        _jetDeepCsv_udsg[_nJets]          = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _jetDeepCsv[_nJets]               = _jetDeepCsv_b[_nJets] + _jetDeepCsv_bb[_nJets];
        if( std::isnan( _jetDeepCsv[_nJets] ) ) _jetDeepCsv[_nJets] = 0.;
        _jetHadronFlavor[_nJets]          = jet.hadronFlavour();

        _jetNeutralHadronFraction[_nJets] = jet.neutralHadronEnergyFraction();
        _jetChargedHadronFraction[_nJets] = jet.chargedHadronEnergyFraction();
        _jetNeutralEmFraction[_nJets]     = jet.neutralEmEnergyFraction();
        _jetChargedEmFraction[_nJets]     = jet.chargedEmEnergyFraction();
        _jetHFHadronFraction[_nJets]      = jet.HFHadronEnergyFraction();
        _jetHFEmFraction[_nJets]          = jet.HFEMEnergyFraction();

        ++_nJets;
    }

    for(const auto& jet : *jetsPuppi){
        if(_nJetsPuppi == nJets_max) break;

        // Implementation based on global tag
        jecUncPuppi->setJetEta(jet.eta());
        jecUncPuppi->setJetPt(jet.pt());
        double unc = jecUncPuppi->getUncertainty(true);

        _jetPuppiPt[_nJetsPuppi] = jet.pt();
        _jetPuppiEta[_nJetsPuppi] = jet.eta();
        _jetPuppiPhi[_nJetsPuppi] = jet.phi();

        ////JEC uncertainties
        _jetPuppiPt_JECDown[_nJetsPuppi] = _jetPuppiPt[_nJetsPuppi]*(1 - unc);
        _jetPuppiPt_JECUp[_nJetsPuppi]   = _jetPuppiPt[_nJetsPuppi]*(1 + unc);


        // JEC implementation based on txt
        //double corrTxt = multilepAnalyzer->jecPuppi->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), *rho, jet.jetArea(), jecLevel);
        //_jetPuppiPt[_nJetsPuppi] = jet.correctedP4("Uncorrected").Pt()*corrTxt;
        //_jetPuppiEta[_nJetsPuppi] = jet.eta();
        //_jetPuppiPhi[_nJetsPuppi] = jet.phi();

        ////JEC uncertainties
        //double unc = multilepAnalyzer->jecPuppi->jetUncertainty(_jetPt[_nJetsPuppi], _jetEta[_nJetsPuppi]);
        //_jetPuppiPt_JECDown[_nJetsPuppi] = _jetPuppiPt[_nJetsPuppi]*(1 - unc);
        //_jetPuppiPt_JECUp[_nJetsPuppi]   = _jetPuppiPt[_nJetsPuppi]*(1 + unc);

        ++_nJetsPuppi;
    }

    //determine the met of the event and its uncertainties
    //nominal Type 1 MET value
    const pat::MET& met = (*mets).front();
    _metType1             = met.pt();
    _metType1Phi          = met.phi();

    //Type 1 Corrections if JEC from txt file
    //std::pair< double, double> metAndPhi = multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho);
    //_metType1    = metAndPhi.first;
    //_metType1Phi = metAndPhi.second;

    //XY Corrections
    std::pair< double, double> metAndPhi_XYCorr = multilepAnalyzer->jec->METXYCorr_Met_MetPhi(_metType1, _metType1Phi, npv);
    _met    = metAndPhi_XYCorr.first;
    _metPhi = metAndPhi_XYCorr.second;

    //raw met values
    _metRaw          = met.uncorPt();
    _metRawPhi       = met.uncorPhi();
    //met values with uncertainties varied up and down
    _metJECDown      = met.shiftedPt(pat::MET::JetEnDown);
    _metJECUp        = met.shiftedPt(pat::MET::JetEnUp);
    _metUnclDown     = met.shiftedPt(pat::MET::UnclusteredEnDown);
    _metUnclUp       = met.shiftedPt(pat::MET::UnclusteredEnUp);
    _metResDown      = met.shiftedPt(pat::MET::JetResDown);
    _metResUp        = met.shiftedPt(pat::MET::JetResUp);
    _metPhiJECDown   = met.shiftedPhi(pat::MET::JetEnDown);
    _metPhiJECUp     = met.shiftedPhi(pat::MET::JetEnUp);
    _metPhiUnclUp    = met.shiftedPhi(pat::MET::UnclusteredEnUp);
    _metPhiUnclDown  = met.shiftedPhi(pat::MET::UnclusteredEnDown);
    _metPhiResDown   = met.shiftedPhi(pat::MET::JetResDown);
    _metPhiResUp     = met.shiftedPhi(pat::MET::JetResUp);

    //significance of met
    //note: this is the only one variable which changed between 94X and 102X see https://github.com/cms-sw/cmssw/commit/f7aacfd2ffaac9899ea07d0355afe49bb10a0aeb
    _metSignificance = met.metSignificance();


    //PUPPI MET
    const pat::MET& metPuppi = (*metsPuppi).front();

    //Type 1 Corrections based on txt file
    //std::pair< double, double > metAndPhiPuppi = multilepAnalyzer->jecPuppi->correctedMETAndPhi(metPuppi, *jetsPuppi, *rho);

    _metPuppi             = metPuppi.pt();
    _metPuppiPhi          = metPuppi.phi();

    _metPuppiRaw          = metPuppi.uncorPt();
    _metPuppiRawPhi       = metPuppi.uncorPhi();
    _metPuppiJECDown      = metPuppi.shiftedPt(pat::MET::JetEnDown);
    _metPuppiJECUp        = metPuppi.shiftedPt(pat::MET::JetEnUp);
    _metPuppiUnclDown     = metPuppi.shiftedPt(pat::MET::UnclusteredEnDown);
    _metPuppiUnclUp       = metPuppi.shiftedPt(pat::MET::UnclusteredEnUp);
    _metPuppiResDown      = metPuppi.shiftedPt(pat::MET::JetResDown);
    _metPuppiResUp        = metPuppi.shiftedPt(pat::MET::JetResUp);
    _metPuppiPhiJECDown   = metPuppi.shiftedPhi(pat::MET::JetEnDown);
    _metPuppiPhiJECUp     = metPuppi.shiftedPhi(pat::MET::JetEnUp);
    _metPuppiPhiUnclUp    = metPuppi.shiftedPhi(pat::MET::UnclusteredEnUp);
    _metPuppiPhiUnclDown  = metPuppi.shiftedPhi(pat::MET::UnclusteredEnDown);
    _metPuppiPhiResDown   = metPuppi.shiftedPhi(pat::MET::JetResDown);
    _metPuppiPhiResUp     = metPuppi.shiftedPhi(pat::MET::JetResUp);

    if(multilepAnalyzer->skim == "singlejet" and _nJets < 1) return false;
    if(multilepAnalyzer->skim == "FR" and _nJets < 1)        return false;
    return true;
}

/*
 * JetID implementations, references:
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
 */

// WARNING: There should be no loose Jet ID for 2017, not sure where the cuts for this below originate, so use at own risk
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

bool JetAnalyzer::jetIsTight(const pat::Jet& jet, const bool is2017, const bool is2018) const{
    if(is2018){
      if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >=  0.9)                         return false;
        if(jet.neutralEmEnergyFraction() >= 0.9)                              return false;
        if(jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1)          return false;
        if(jet.chargedHadronEnergyFraction() <= 0 and fabs(jet.eta()) <= 2.6) return false; // only for |eta|<2.6
        if(jet.chargedMultiplicity() <= 0)                                    return false;
      } else if(fabs(jet.eta()) <= 3.0){
        if(jet.neutralHadronEnergyFraction() >= 0.99)                         return false;
        if(jet.neutralEmEnergyFraction() <= 0.02)                             return false;
        if(jet.neutralMultiplicity() <= 2)                                    return false;
      } else {
        if(jet.neutralEmEnergyFraction() >= 0.90)                             return false;
        if(jet.neutralMultiplicity() <= 10)                                   return false;
        if(jet.neutralHadronEnergyFraction() <= 0.02)                         return false;
      }
      return true;
    } else {
      if(!jetIsLoose(jet, is2017))                                            return false;

      if(fabs(jet.eta()) <= 2.7){
          if(jet.neutralHadronEnergyFraction() >= 0.9)                        return false;
          if(jet.neutralEmEnergyFraction() >= 0.9)                            return false;

      } else if(fabs(jet.eta()) <= 3.0){
          if(is2017 && jet.neutralEmEnergyFraction() >= 0.99)                 return false;
      }
      return true;
    }
}

bool JetAnalyzer::jetIsTightLepVeto(const pat::Jet& jet, const bool is2017, const bool is2018) const{
    if(!jetIsTight(jet, is2017, is2018))                                                  return false; // Similar to tight ID except with additional cuts:
    if(fabs(jet.eta()) <= 2.7){
      if(jet.chargedMuEnergyFraction() >= 0.8 )                                           return false; // Muon energy fraction cut
      if(is2018 and jet.chargedEmEnergyFraction() >= 0.8)                                 return false; // EM fraction cut 2018
      else if(is2017 and fabs(jet.eta()) <= 2.4 and jet.chargedEmEnergyFraction() >= 0.8) return false; // EM fraction cut 2017
      else if(fabs(jet.eta()) <= 2.4 and jet.chargedEmEnergyFraction() >= 0.9)            return false; // EM fraction cut 2016
    }
    return true;
}
