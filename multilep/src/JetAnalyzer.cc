#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"

//include c++ library classes
#include <algorithm>

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{
    jecUnc = new JetCorrectionUncertainty((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath());

    std::vector< std::string > jecSources = {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "FlavorZJet", "FlavorPhotonJet", "FlavorPureGluon", "FlavorPureQuark", "FlavorPureCharm", "FlavorPureBottom", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpMuZero", "PileUpEnvelope", "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", "TotalNoFlavor", "TotalNoTime" ,"TotalNoFlavorNoTime", "Total" };
    for( const auto& source : jecSources ){
        jetSourcesCorParameters[ source ] = std::make_shared< JetCorrectorParameters >( (iConfig.getParameter<edm::FileInPath>("jecUncertaintySourcesFile") ).fullPath(), source );

        //make sure map entries exist
        _jetPt_allVariationsDown[ source ];
        _jetPt_allVariationsUp[ source ];
        _jetSmearedPt_allVariationsDown[ source ];
        _jetSmearedPt_allVariationsUp[ source ];
        _corrMETx_allVariationsDown[ source ];
        _corrMETx_allVariationsUp[ source ];
        _corrMETy_allVariationsDown[ source ];
        _corrMETy_allVariationsUp[ source ];
    }   

    std::vector< std::string > jecSourcesGrouped = { "Absolute", "Absolute_2016", "BBEC1", "BBEC1_2016", "EC2", "EC2_2016", "FlavorQCD", "HF", "HF_2016", "RelativeBal", "RelativeSample_2016", "Total" };

    if( multilepAnalyzer->is2018() || multilepAnalyzer->is2017() ){
        std::string replacement = "2017";
        if( multilepAnalyzer->is2018() ) replacement = "2018";
        for( auto& entry: jecSourcesGrouped ){
            size_t index = entry.find( "2016" );
            if( index == std::string::npos ) continue;
            entry.replace(index, 4, replacement );
        }
    }

    for( const auto& source : jecSourcesGrouped ){
        jetGroupedCorParameters[ source ] = std::make_shared< JetCorrectorParameters >( (iConfig.getParameter<edm::FileInPath>("jecUncertaintyRegroupedFile") ).fullPath(), source );
        _jetPt_groupedVariationsDown[ source ];
        _jetPt_groupedVariationsUp[ source ];
        _jetSmearedPt_groupedVariationsDown[ source ];
        _jetSmearedPt_groupedVariationsUp[ source ];
        _corrMETx_groupedVariationsDown[ source ];
        _corrMETx_groupedVariationsUp[ source ];
        _corrMETy_groupedVariationsDown[ source ];
        _corrMETy_groupedVariationsUp[ source ];       
    }

    std::vector< JetCorrectorParameters > JECParameters;
    JECParameters.push_back(JetCorrectorParameters( (iConfig.getParameter<edm::FileInPath>("jecL1FastJetFile") ).fullPath() ));
    JECParameters.push_back(JetCorrectorParameters( (iConfig.getParameter<edm::FileInPath>("jecL2RelativeFile") ).fullPath() ));
    JECParameters.push_back(JetCorrectorParameters( (iConfig.getParameter<edm::FileInPath>("jecL3AbsoluteFile") ).fullPath() ));
    if( multilepAnalyzer->isData() ) JECParameters.push_back(JetCorrectorParameters( (iConfig.getParameter<edm::FileInPath>("jecL2L3ResidualFile") ).fullPath() ));
    jetCorrector.reset(new FactorizedJetCorrector(JECParameters));
};


JetAnalyzer::~JetAnalyzer(){
    delete jecUnc;
}


std::vector<float> JetAnalyzer::getSubCorrections(double rawPt, double eta, double rho, double area)
{   
    jetCorrector->setJetEta(eta);
    jetCorrector->setRho(rho);
    jetCorrector->setJetA(area);
    jetCorrector->setJetPt(rawPt); 
    std::vector< float > corrections = jetCorrector->getSubCorrections();
    return corrections;
}


double px(double pt, double phi) { return pt*cos(phi); };
double py(double pt, double phi) { return pt*sin(phi); };


std::pair<double, double> JetAnalyzer::getMETCorrectionPxPy(double corrPt, double rawPt, double rawEta, double rawMuonSubtractedPt, double phi, double emf, double rho, double area, const std::string& source, unsigned jetIndex, double jecShift)
{      
    std::vector< float > corrections = getSubCorrections(rawPt, rawEta, rho, area);

//    double l1corrptNoMuon   = rawMuonSubtractedPt*corrections.front();
    double fullcorrpt = rawMuonSubtractedPt*corrections.back();
    double PT_L1L2L3 = rawMuonSubtractedPt*corrections.back() + (rawPt-rawMuonSubtractedPt);
    double PT_L1 = rawPt*corrections.front();
   
    float f = PT_L1L2L3*(jecShift/corrPt);
   
    // the corrections for the MET are the difference between l1fastjet and the full corrections on the jet!
    if( emf > 0.90 or fullcorrpt < 15. || ( std::abs(rawEta) > 9.9 ) ) return {0., 0.};

    float ptdiff = (PT_L1 - f);

//    float ptdiffNominal = (l1corrptNoMuon - fullcorrpt);
   
    std::pair<double, double> corr = {px(ptdiff, phi), py(ptdiff, phi)};

    return corr; 
}

void JetAnalyzer::correctedMETAndPhi(const pat::MET& met, const std::vector< pat::Jet >& jets, const double rho)
{
    for (auto it=_corrMETx_groupedVariationsDown.begin(); it!=_corrMETx_groupedVariationsDown.end(); ++it)
    {	
        std::string source = it->first;
        _corrMETx_groupedVariationsDown[ source ] = met.uncorPx();
        _corrMETx_groupedVariationsUp[ source ] = met.uncorPx();
        _corrMETy_groupedVariationsDown[ source ] = met.uncorPy();
        _corrMETy_groupedVariationsUp[ source ] = met.uncorPy();
    }   

    for (auto it=_corrMETx_allVariationsDown.begin(); it!=_corrMETx_allVariationsDown.end(); ++it)
    {	
        std::string source = it->first;
        _corrMETx_allVariationsDown[ source ] = met.uncorPx();
        _corrMETx_allVariationsUp[ source ] = met.uncorPx();
        _corrMETy_allVariationsDown[ source ] = met.uncorPy();
        _corrMETy_allVariationsUp[ source ] = met.uncorPy();
    }
   
//   std::cout << "-----" << std::endl;
    //loop over all jets
    int iJet = 0;
    for(auto& jet : jets)
    {
        //make lorentzVector of raw jet pt
        TLorentzVector jetV;
        jetV.SetPtEtaPhiE(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jet.correctedP4("Uncorrected").Phi(), jet.correctedP4("Uncorrected").E());
       
        //clean jet from muons
        const std::vector<reco::CandidatePtr>& daughters = jet.daughterPtrVector();
        for(auto& daughterPtr : daughters)
        {	     
            const reco::PFCandidate* daughter = dynamic_cast<const reco::PFCandidate* >( daughterPtr.get() );
            const reco::Candidate* muon = (daughter != nullptr ?  
                    (daughter->muonRef().isNonnull() ? daughter->muonRef().get() : nullptr)
                    : daughterPtr.get() );
            if(muon != nullptr && ( muon->isGlobalMuon() || muon->isStandAloneMuon() ) )
            {

                TLorentzVector muonV(muon->px(), muon->py(), muon->pz(), muon->energy());
                jetV -= muonV;
            }
        }

//       std::cout << "jet #" << iJet << ": pt(uncorrected)=" << jet.correctedP4("Uncorrected").Pt() << " pt(corrected)=" << jet.pt() << " up(total)=" << _jetPt_allVariationsUp["Total"][iJet] << " down(total)=" <<
//	 _jetPt_allVariationsDown["Total"][iJet] << " eta=" << jet.eta() << std::endl;

        for (auto it=_corrMETx_groupedVariationsDown.begin(); it!=_corrMETx_groupedVariationsDown.end(); ++it)
        {	
            std::string source = it->first;

            //get JEC on px and py 
            std::pair<double, double> corrUp = getMETCorrectionPxPy(jet.pt(), jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jetV.Pt(),
                    jetV.Phi(), jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction(), rho, jet.jetArea(), source, iJet, _jetPt_groupedVariationsUp[source][iJet]);
            std::pair<double, double> corrDown = getMETCorrectionPxPy(jet.pt(), jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jetV.Pt(),
                    jetV.Phi(), jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction(), rho, jet.jetArea(), source, iJet, _jetPt_groupedVariationsDown[source][iJet]);
	   
            //apply corrections to current met values
            _corrMETx_groupedVariationsDown[ source ] += corrDown.first;
            _corrMETy_groupedVariationsDown[ source ] += corrDown.second;
            _corrMETx_groupedVariationsUp[ source ] += corrUp.first;
            _corrMETy_groupedVariationsUp[ source ] += corrUp.second;
        }

        for (auto it=_corrMETx_allVariationsDown.begin(); it!=_corrMETx_allVariationsDown.end(); ++it)
        {	
            std::string source = it->first;

            //get JEC on px and py
            std::pair<double, double> corrUp = getMETCorrectionPxPy(jet.pt(), jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jetV.Pt(),
                    jetV.Phi(), jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction(), rho, jet.jetArea(), source, iJet, _jetPt_allVariationsUp[source][iJet]);
            std::pair<double, double> corrDown = getMETCorrectionPxPy(jet.pt(), jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jetV.Pt(),
                    jetV.Phi(), jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction(), rho, jet.jetArea(), source, iJet, _jetPt_allVariationsDown[source][iJet]);
	   
//	   if( source == "Total" )
//	     std::cout << "METx contribution from jet #" << iJet << " Up=" << corrUp.first << " Down=" << corrDown.first << std::endl;
	   
            //apply corrections to current met values
            _corrMETx_allVariationsDown[ source ] += corrDown.first;
            _corrMETy_allVariationsDown[ source ] += corrDown.second;
            _corrMETx_allVariationsUp[ source ] += corrUp.first;
            _corrMETy_allVariationsUp[ source ] += corrUp.second;
        }
       
        iJet ++;
    }

//   float metxUp = _corrMETx_allVariationsUp["Total"];
//   float metyUp = _corrMETy_allVariationsUp["Total"];
//   float metxDown = _corrMETx_allVariationsDown["Total"];
//   float metyDown = _corrMETy_allVariationsDown["Total"];
//   float metPtUp = sqrt(metxUp*metxUp + metyUp*metyUp);
//   float metPtDown = sqrt(metxDown*metxDown + metyDown*metyDown);

//   std::cout << "met uncorrected: " << met.uncorPt() << std::endl;
//   std::cout << "met corrected: " << _met << std::endl;
//   std::cout << "met up: " << _metJECUp << " " << metPtUp << std::endl;
//   std::cout << "met down: " << _metJECDown << " " << metPtDown << std::endl;

//   double correctedMET = sqrt(corrMETx*corrMETx + corrMETy*corrMETy);
//   double correctedMETPhi = atan2(corrMETy, corrMETx);
//   return {correctedMET, correctedMETPhi};
}

//helper function to set branches for jet uncertainty source splitting
template< typename T > void setSourceBranches( TTree* outputTree, T& sourceMap, const std::string& namePrefix, const std::string& namePostfix ){
    for( auto& entry : sourceMap ){
        outputTree->Branch( ( namePrefix + entry.first + namePostfix ).c_str(), &entry.second, ( namePrefix + entry.first + namePostfix + "[_nJets]/D" ).c_str() );
    }
}

template< typename T > void setMETSourceBranches( TTree* outputTree, T& sourceMap, const std::string& namePrefix, const std::string& namePostfix ){
    for( auto& entry : sourceMap ){
        outputTree->Branch( ( namePrefix + entry.first + namePostfix ).c_str(), &entry.second, ( namePrefix + entry.first + namePostfix + "/D" ).c_str() );
    }
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
    outputTree->Branch("_jetDeepFlavor_b",           &_jetDeepFlavor_b,          "_jetDeepFlavor_b[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_bb",          &_jetDeepFlavor_bb,         "_jetDeepFlavor_bb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_lepb",        &_jetDeepFlavor_lepb,       "_jetDeepFlavor_lepb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor",             &_jetDeepFlavor,            "_jetDeepFlavor[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_c",           &_jetDeepFlavor_c,          "_jetDeepFlavor_c[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_uds",         &_jetDeepFlavor_uds,        "_jetDeepFlavor_uds[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_g",           &_jetDeepFlavor_g,          "_jetDeepFlavor_g[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetPartonFlavor",           &_jetPartonFlavor,          "_jetPartonFlavor[_nJets]/i");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");
    outputTree->Branch("_jetPileupIdFullDisc",       &_jetPileupIdFullDisc,      "_jetPileupIdFullDisc[_nJets]/D");
    outputTree->Branch("_jetPileupIdFullId",         &_jetPileupIdFullId,        "_jetPileupIdFullId[_nJets]/I");
   
    if( ! multilepAnalyzer->isData() ) {
       
       outputTree->Branch("_jetHasGen",                 &_jetHasGen,                "_jetHasGen[_nJets]/O");
       outputTree->Branch("_jetGenPt",                  &_jetGenPt,                 "_jetGenPt[_nJets]/D");
       outputTree->Branch("_jetGenEta",                 &_jetGenEta,                "_jetGenEta[_nJets]/D");
       outputTree->Branch("_jetGenPhi",                 &_jetGenPhi,                "_jetGenPhi[_nJets]/D");
       outputTree->Branch("_jetGenE",                   &_jetGenE,                  "_jetGenE[_nJets]/D");
       
    }   

    outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
    outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
    outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
    outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
    outputTree->Branch("_jetHFHadronFraction",       &_jetHFHadronFraction,      "_jetHFHadronFraction[_nJets]/D");
    outputTree->Branch("_jetHFEmFraction",           &_jetHFEmFraction,          "_jetHFEmFraction[_nJets]/D");

    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metRaw",                       &_metRaw,                       "_metRaw/D");
    outputTree->Branch("_met_JECDown",                  &_metJECDown,                   "_met_JECDown/D");
    outputTree->Branch("_met_JECUp",                    &_metJECUp,                     "_met_JECUp/D");
    outputTree->Branch("_met_UnclDown",                 &_metUnclDown,                  "_met_UnclDown/D");
    outputTree->Branch("_met_UnclUp",                   &_metUnclUp,                    "_met_UnclUp/D");

    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metRawPhi",                    &_metRawPhi,                    "_metRawPhi/D");
    outputTree->Branch("_metPhi_JECDown",               &_metPhiJECDown,                "_metPhi_JECDown/D");
    outputTree->Branch("_metPhi_JECUp",                 &_metPhiJECUp,                  "_metPhi_JECUp/D");
    outputTree->Branch("_metPhi_UnclDown",              &_metPhiUnclDown,               "_metPhi_UnclDown/D");
    outputTree->Branch("_metPhi_UnclUp",                &_metPhiUnclUp,                 "_metPhi_UnclUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");

    if(!multilepAnalyzer->is2018() ) outputTree->Branch("_jetIsLoose", _jetIsLoose, "_jetIsLoose[_nJets]/O"); // WARNING, not recommended to be used, only exists for 2016

    // Store splitting of JEC unc into sources
    if( multilepAnalyzer->storeJecSources ){
        setSourceBranches( outputTree, _jetPt_allVariationsDown, "_jetPt_", "_JECSourcesDown" );
        setSourceBranches( outputTree, _jetPt_allVariationsUp, "_jetPt_", "_JECSourcesUp" );
        setSourceBranches( outputTree, _jetPt_groupedVariationsDown, "_jetPt_", "_JECGroupedDown" );
        setSourceBranches( outputTree, _jetPt_groupedVariationsUp, "_jetPt_", "_JECGroupedUp" );

        setSourceBranches( outputTree, _jetSmearedPt_allVariationsDown, "_jetSmearedPt_", "_JECSourcesDown" );
        setSourceBranches( outputTree, _jetSmearedPt_allVariationsUp, "_jetSmearedPt_", "_JECSourcesUp" );
        setSourceBranches( outputTree, _jetSmearedPt_groupedVariationsDown, "_jetSmearedPt_", "_JECGroupedDown" );
        setSourceBranches( outputTree, _jetSmearedPt_groupedVariationsUp, "_jetSmearedPt_", "_JECGroupedUp" );

        setMETSourceBranches( outputTree, _corrMETx_allVariationsDown, "_corrMETx_", "_JECSourcesDown" );
        setMETSourceBranches( outputTree, _corrMETx_allVariationsUp, "_corrMETx_", "_JECSourcesUp" );
        setMETSourceBranches( outputTree, _corrMETx_groupedVariationsDown, "_corrMETx_", "_JECGroupedDown" );
        setMETSourceBranches( outputTree, _corrMETx_groupedVariationsUp, "_corrMETx_", "_JECGroupedUp" );

        setMETSourceBranches( outputTree, _corrMETy_allVariationsDown, "_corrMETy_", "_JECSourcesDown" );
        setMETSourceBranches( outputTree, _corrMETy_allVariationsUp, "_corrMETy_", "_JECSourcesUp" );
        setMETSourceBranches( outputTree, _corrMETy_groupedVariationsDown, "_corrMETy_", "_JECGroupedDown" );
        setMETSourceBranches( outputTree, _corrMETy_groupedVariationsUp, "_corrMETy_", "_JECGroupedUp" );
    }
}


//helper function to fill maps with split JEC uncertainties
template< typename T > void fillJetUncertaintySources( const std::map< std::string, std::shared_ptr< JetCorrectorParameters> >& jetCorrParameters,  T& sourcesMapDown, T& sourcesMapUp, const pat::Jet& jet, const unsigned jetIndex ){
    for( const auto& source : jetCorrParameters ){
        JetCorrectionUncertainty* jetCorUnc = new JetCorrectionUncertainty( *( source.second ) );
        jetCorUnc->setJetPt( jet.pt() );
        jetCorUnc->setJetEta( jet.eta() );

        double uncJec = jetCorUnc->getUncertainty( true );

        sourcesMapDown[ source.first ][ jetIndex ] = jet.pt()*( 1 - uncJec );
        sourcesMapUp[ source.first ][ jetIndex ] = jet.pt()*( 1 + uncJec );

        delete jetCorUnc;
    }
}


bool JetAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Jet>> jets            = getHandle(iEvent, multilepAnalyzer->jetToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmeared     = getHandle(iEvent, multilepAnalyzer->jetSmearedToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedUp   = getHandle(iEvent, multilepAnalyzer->jetSmearedUpToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedDown = getHandle(iEvent, multilepAnalyzer->jetSmearedDownToken);
    edm::Handle<std::vector<pat::MET>> mets            = getHandle(iEvent, multilepAnalyzer->metToken);

    edm::Handle<double> rho                            = getHandle(iEvent, multilepAnalyzer->rhoToken);

    _nJets = 0;

    for(const auto& jet : *jets){
        if(_nJets == nJets_max) break;

        _jetIsLoose[_nJets]        = jetIsLoose(jet, multilepAnalyzer->is2017() || multilepAnalyzer->is2018() );
        _jetIsTight[_nJets]        = jetIsTight(jet, multilepAnalyzer->is2017(), multilepAnalyzer->is2018() );
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jet, multilepAnalyzer->is2017(), multilepAnalyzer->is2018() );

        if( jet.hasUserFloat("pileupJetId:fullDiscriminant") ) _jetPileupIdFullDisc[_nJets] = jet.userFloat("pileupJetId:fullDiscriminant");
        if( jet.hasUserInt("pileupJetIdUpdated:fullId") ) _jetPileupIdFullId[_nJets] = jet.userInt("pileupJetIdUpdated:fullId");

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
        _jetSmearedPt_JECDown[_nJets] = _jetSmearedPt[_nJets]*( 1. - uncSmeared );
        _jetSmearedPt_JECUp[_nJets]   = _jetSmearedPt[_nJets]*( 1. + uncSmeared );
        _jetSmearedPt_JERDown[_nJets] = jetSmearedDownIt->pt();
        _jetSmearedPt_JERUp[_nJets]   = jetSmearedUpIt->pt();


        //split of JEC uncertainties into separate sources
        if( multilepAnalyzer->storeJecSources )
        {	   
            fillJetUncertaintySources( jetGroupedCorParameters, _jetPt_groupedVariationsDown, _jetPt_groupedVariationsUp, jet, _nJets );
            fillJetUncertaintySources( jetSourcesCorParameters, _jetPt_allVariationsDown, _jetPt_allVariationsUp, jet, _nJets );

            fillJetUncertaintySources( jetGroupedCorParameters, _jetSmearedPt_groupedVariationsDown, _jetSmearedPt_groupedVariationsUp, *jetSmearedIt, _nJets );
            fillJetUncertaintySources( jetSourcesCorParameters, _jetSmearedPt_allVariationsDown, _jetSmearedPt_allVariationsUp, *jetSmearedIt, _nJets );
        }       

        _jetPt_Uncorrected[_nJets]        = jet.correctedP4("Uncorrected").Pt();
        _jetPt_L1[_nJets]                 = jet.correctedP4("L1FastJet").Pt();
        _jetPt_L2[_nJets]                 = jet.correctedP4("L2Relative").Pt();
        _jetPt_L3[_nJets]                 = jet.correctedP4("L3Absolute").Pt();

        _jetEta[_nJets]                   = jet.eta();
        _jetPhi[_nJets]                   = jet.phi();
        _jetE[_nJets]                     = jet.energy();

        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        //DeepCSV tagger
        _jetDeepCsv_udsg[_nJets]          = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _jetDeepCsv[_nJets]               = _jetDeepCsv_b[_nJets] + _jetDeepCsv_bb[_nJets];
        if( std::isnan( _jetDeepCsv[_nJets] ) ) _jetDeepCsv[_nJets] = 0.;

        //DeepFlavor taggeer 
        _jetDeepFlavor_b[_nJets]          = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
        _jetDeepFlavor_bb[_nJets]         = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
        _jetDeepFlavor_lepb[_nJets]       = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        _jetDeepFlavor[_nJets]            = _jetDeepFlavor_b[_nJets] + _jetDeepFlavor_bb[_nJets] + _jetDeepFlavor_lepb[_nJets];
        if( std::isnan( _jetDeepFlavor[_nJets] ) ) _jetDeepFlavor[_nJets] = 0.;
        _jetDeepFlavor_c[_nJets]          = jet.bDiscriminator("pfDeepFlavourJetTags:probc");
        _jetDeepFlavor_uds[_nJets]        = jet.bDiscriminator("pfDeepFlavourJetTags:probuds");
        _jetDeepFlavor_g[_nJets]          = jet.bDiscriminator("pfDeepFlavourJetTags:probg");

        _jetHadronFlavor[_nJets]          = jet.hadronFlavour();
        _jetPartonFlavor[_nJets]          = jet.partonFlavour();

        _jetNeutralHadronFraction[_nJets] = jet.neutralHadronEnergyFraction();
        _jetChargedHadronFraction[_nJets] = jet.chargedHadronEnergyFraction();
        _jetNeutralEmFraction[_nJets]     = jet.neutralEmEnergyFraction();
        _jetChargedEmFraction[_nJets]     = jet.chargedEmEnergyFraction();
        _jetHFHadronFraction[_nJets]      = jet.HFHadronEnergyFraction();
        _jetHFEmFraction[_nJets]          = jet.HFEMEnergyFraction();

        _jetHasGen[_nJets] = 0;
        _jetGenPt[_nJets] = 0;
        _jetGenEta[_nJets] = 0;
        _jetGenPhi[_nJets] = 0;
        _jetGenE[_nJets] = 0;
       
        if( ! multilepAnalyzer->isData() ) {	    
	   if( jet.genJet() != NULL ) {
	      
	      _jetHasGen[_nJets] = 1;
	      _jetGenPt[_nJets] = jet.genJet()->pt();
	      _jetGenEta[_nJets] = jet.genJet()->eta();
	      _jetGenPhi[_nJets] = jet.genJet()->phi();
	      _jetGenE[_nJets] = jet.genJet()->energy();
	      
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

    // propagate JEC to MET
    correctedMETAndPhi(met, *jets, *rho);
   
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
