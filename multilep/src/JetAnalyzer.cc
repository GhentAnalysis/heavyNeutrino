#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

//include c++ library classes
#include <algorithm>

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{
    jecUnc = new JetCorrectionUncertainty((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath());
};

JetAnalyzer::~JetAnalyzer(){
    delete jecUnc;
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

    // variables for Particle Flow HNL Jet tagger MVA based on deep sets theorem
    outputTree->Branch("_nJetConstituents",      &_nJetConstituents,      "_nJetConstituents[_nJets]/i");
    std::string jetConstituentsArraySize = "[_nJets][" + std::to_string(maxJetSize) + "]";
    outputTree->Branch( "_JetConstituentPt",     &_JetConstituentPt,      std::string("_JetConstituentPt" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentEta",    &_JetConstituentEta,     std::string("_JetConstituentEta" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentPhi",    &_JetConstituentPhi,     std::string("_JetConstituentPhi" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentMass",   &_JetConstituentMass,    std::string("_JetConstituentMass" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentPdgId",  &_JetConstituentPdgId,   std::string("_JetConstituentPdgId" + jetConstituentsArraySize + "/I").c_str() );
    outputTree->Branch( "_JetConstituentPdgIdReduced",  &_JetConstituentPdgIdReduced,   std::string("_JetConstituentPdgIdReduced" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentCharge", &_JetConstituentCharge,  std::string("_JetConstituentCharge" + jetConstituentsArraySize + "/I").c_str() );
    outputTree->Branch( "_JetConstituentdxy",    &_JetConstituentdxy,     std::string("_JetConstituentdxy" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentdz",     &_JetConstituentdz,      std::string("_JetConstituentdz" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentdxyErr", &_JetConstituentdxyErr,  std::string("_JetConstituentdxyErr" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentdzErr",  &_JetConstituentdzErr,   std::string("_JetConstituentdzErr" + jetConstituentsArraySize + "/D").c_str() );
    outputTree->Branch( "_JetConstituentNumberOfHits",      &_JetConstituentNumberOfHits, std::string("_JetConstituentNumberOfHits" + jetConstituentsArraySize + "/I").c_str() );
    outputTree->Branch( "_JetConstituentNumberOfPixelHits", &_JetConstituentNumberOfPixelHits, std::string("_JetConstituentNumberOfPixelHits" + jetConstituentsArraySize + "/I").c_str() );
    outputTree->Branch( "_JetConstituentHasTrack",          &_JetConstituentHasTrack, std::string("_JetConstituentHasTrack" + jetConstituentsArraySize + "/O").c_str() );
    
    if(!multilepAnalyzer->is2018() ) outputTree->Branch("_jetIsLoose", _jetIsLoose, "_jetIsLoose[_nJets]/O"); // WARNING, not recommended to be used, only exists for 2016
}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Jet>> jets            = getHandle(iEvent, multilepAnalyzer->jetToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmeared     = getHandle(iEvent, multilepAnalyzer->jetSmearedToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedUp   = getHandle(iEvent, multilepAnalyzer->jetSmearedUpToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedDown = getHandle(iEvent, multilepAnalyzer->jetSmearedDownToken);
    edm::Handle<std::vector<pat::MET>> mets            = getHandle(iEvent, multilepAnalyzer->metToken);

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

        //new DeepFlavour tagger
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
	    
        /*
         * Jet Particle Flow constituent information for HNL Jet tagger MVA
         */
        _nJetConstituents[_nJets]         = std::min((unsigned) jet.numberOfDaughters(), maxJetSize);
        for(unsigned d = 0; d < _nJetConstituents[_nJets]; ++d){
            const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);
            _JetConstituentPt[_nJets][d]     = daughter->pt();
            _JetConstituentEta[_nJets][d]    = daughter->eta();
            _JetConstituentPhi[_nJets][d]    = daughter->phi();
            _JetConstituentMass[_nJets][d]   = daughter->mass();
            _JetConstituentPdgId[_nJets][d]  = daughter->pdgId();
            _JetConstituentPdgIdReduced[_nJets][d]  = reducedPdgId(daughter->pdgId());

            _JetConstituentCharge[_nJets][d] = daughter->charge();
            if( daughter->hasTrackDetails() ){
                _JetConstituentdxy[_nJets][d] = fabs(daughter->dxy());
                _JetConstituentdz[_nJets][d] = fabs(daughter->dz());
                _JetConstituentdxyErr[_nJets][d] = catchNanOrInf(fabs(daughter->dxyError()));
                _JetConstituentdzErr[_nJets][d] = catchNanOrInf(fabs(daughter->dzError()));
                _JetConstituentNumberOfHits[_nJets][d] = daughter->numberOfHits();
                _JetConstituentNumberOfPixelHits[_nJets][d] = daughter->numberOfPixelHits();
                _JetConstituentHasTrack[_nJets][d] = true;
            } else {
                _JetConstituentdxy[_nJets][d] = -1.;
                _JetConstituentdz[_nJets][d] = -1.;
                _JetConstituentdxyErr[_nJets][d] = -1.;
                _JetConstituentdzErr[_nJets][d] = -1.;
                _JetConstituentNumberOfHits[_nJets][d] = -1;
                _JetConstituentNumberOfPixelHits[_nJets][d] = -1;
                _JetConstituentHasTrack[_nJets][d] = false;
            }
        }
        // second loop to put remaining jet constituents to 0 up to maxJetSize
        if( _nJetConstituents[_nJets] < maxJetSize){
            for(unsigned d = _nJetConstituents[_nJets]; d < maxJetSize; ++d){
                _JetConstituentPt[_nJets][d] = 0.;
                _JetConstituentEta[_nJets][d] = 0.;
                _JetConstituentPhi[_nJets][d] = 0.;
                _JetConstituentMass[_nJets][d] = 0.;
                _JetConstituentPdgId[_nJets][d] = 0;
                _JetConstituentCharge[_nJets][d] = 0;
                _JetConstituentdxy[_nJets][d] = 0.;
                _JetConstituentdz[_nJets][d] = 0.;
                _JetConstituentdxyErr[_nJets][d] = 0.;
                _JetConstituentdzErr[_nJets][d] = 0.;
                _JetConstituentNumberOfHits[_nJets][d] = 0;
                _JetConstituentNumberOfPixelHits[_nJets][d] = 0;
                _JetConstituentHasTrack[_nJets][d] = false;
            }
        }
        
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

    //std::cout << "end jetanalyzer" << std::endl;
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
    	//std::cout << "okislooseif1_2" << std::endl;

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

double JetAnalyzer::reducedPdgId( int pdgId ){
    static const std::map< unsigned, double > pdgIdMap = {
        { 0, 0.},
        { 1, 0.125},
        { 2, 0.25},
        { 11, 0.375},
        { 13, 0.5},
        { 22, 0.625},
        { 130, 0.75},
        { 211, 0.875}
    };
    auto entry = pdgIdMap.find( fabs( pdgId ) );
    if( entry != pdgIdMap.cend() ){
        return entry->second;
    } else {
        return 1;
    }
}

double JetAnalyzer::catchNanOrInf( double value ){
    if( std::isnan( value ) ){
        return -1;
    } else if( std::isinf( value ) ){
        return -1;
    } else{
        return value;
    }
}
