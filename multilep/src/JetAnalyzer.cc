#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    jecUnc((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath()),
    multilepAnalyzer(multilepAnalyzer)
{
    /*
    if(multilepAnalyzer->isData){
        jecLevel = "L2L3Residual";
    } else {
        jecLevel = "L3Absolute";
    }
    */
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

    /*outputTree->Branch("_nDaughters",		     &_nDaughters,	         "_nDaughters/I");
    outputTree->Branch("_jet_tag_for_daughters",     &_jet_tag_for_daughters,    "_jet_tag_for_daughters[_nDaughters]/I");
    outputTree->Branch("_jet_daughter_pdgid",	     &_jet_daughter_pdgid,       "_jet_daughter_pdgid[_nDaughters]/I");
    outputTree->Branch("_jet_daughter_pt",	     &_jet_daughter_pt,          "_jet_daughter_pt[_nDaughters]/D");
    outputTree->Branch("_jet_daughter_eta",	     &_jet_daughter_eta,         "_jet_daughter_eta[_nDaughters]/D");
    outputTree->Branch("_jet_daughter_phi",	     &_jet_daughter_phi,         "_jet_daughter_phi[_nDaughters]/D");
    outputTree->Branch("_jet_daughter_energy",	     &_jet_daughter_energy,      "_jet_daughter_energy[_nDaughters]/D");*/

}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
    //std::cout << "begin jetanalyzer" << std::endl;
    edm::Handle<std::vector<pat::Jet>> jets;            iEvent.getByToken(multilepAnalyzer->jetToken,            jets);
    edm::Handle<std::vector<pat::MET>> mets;            iEvent.getByToken(multilepAnalyzer->metToken, mets);
    //to apply JEC from txt files
    //edm::Handle<double> rho;                            iEvent.getByToken(multilepAnalyzer->rhoToken,            rho);

    _nJets = 0;
    //_nDaughters = 0;
    //std::cout << "ok1" << std::endl;
    //for(auto jetSmeared = jetsSmeared->begin(); jetSmeared != jetsSmeared->end(); ++jetSmeared){
    //std::cout << "ok1.5" << std::endl;
    for(auto& jet : *jets){
        if(_nJets == nJets_max) break;// or _nDaughters == nDaughters_max) break;

    	//std::cout << "ok1.55" << std::endl;
        //only store loose jets 
    	//std::cout << "jet pt: " << jet.pt() << std::endl;
    	//std::cout << "jet eta: " << jet.eta() << std::endl;
    	//std::cout << "jet phi: " << jet.phi() << std::endl;
    	//std::cout << "jet energy: " << jet.energy() << std::endl;
    	//std::cout << "_nJets: " << _nJets << std::endl;
	//if(multilepAnalyzer->is2017 or !multilepAnalyzer->is2017) std::cout << "is2017 ok" << std::endl;
        _jetIsLoose[_nJets] = jetIsLoose(jet, multilepAnalyzer->is2017);
    	//std::cout << "ok1.56" << std::endl;
        if(!_jetIsLoose[_nJets]) continue;
    	//std::cout << "ok1.57" << std::endl;
        _jetIsTight[_nJets] = jetIsTight(jet, multilepAnalyzer->is2017);
    	//std::cout << "ok1.58" << std::endl;
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jet, multilepAnalyzer->is2017);
    	//std::cout << "ok1.6" << std::endl;
    
        jecUnc.setJetEta(jet.eta());
        jecUnc.setJetPt(jet.pt());
        double unc = jecUnc.getUncertainty(true);
    	//std::cout << "ok1.7" << std::endl;

        //txt based JEC
        /*
        double corrTxt = multilepAnalyzer->jec->jetCorrection(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), *rho, jet.jetArea(), jecLevel);
        _jetPt[_nJets] = jet.correctedP4("Uncorrected").Pt()*corrTxt;
        double unc = multilepAnalyzer->jec->jetUncertainty(_jetPt[_nJets], jet.eta());
        */
        //extract uncertainty based on corrected jet pt!

        _jetPt[_nJets] = jet.pt();
        double maxpT = (1+unc)*_jetPt[_nJets];
        if(maxpT <= 25) continue;

        _jetPt_JECDown[_nJets]            = _jetPt[_nJets]*(1-unc);
        _jetPt_JECUp[_nJets]              = _jetPt[_nJets]*(1+unc);
        _jetEta[_nJets]                   = jet.eta();
        _jetPhi[_nJets]                   = jet.phi();
        //_jetE[_nJets]                     = jet.correctedP4("Uncorrected").E()*corrTxt;
        _jetE[_nJets]                     = jet.energy();
        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        //new DeepFlavour tagger
        _jetDeepCsv_udsg[_nJets]          = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _jetHadronFlavor[_nJets]          = jet.hadronFlavour();

    	//std::cout << "ok3" << std::endl;
    
	/*for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
      	    const pat::PackedCandidate* daughter  = (const pat::PackedCandidate*) jet.daughter(d);
            _jet_tag_for_daughters[_nDaughters]   = _nJets;
            _jet_daughter_pdgid[_nDaughters] 	  = daughter->pdgId();
            _jet_daughter_pt[_nDaughters] 	  = daughter->pt();
            _jet_daughter_eta[_nDaughters] 	  = daughter->eta();
            _jet_daughter_phi[_nDaughters] 	  = daughter->phi();
            _jet_daughter_energy[_nDaughters]     = daughter->energy();
            ++_nDaughters;
    	}*/
    	//std::cout << "ok4" << std::endl;
        ++_nJets;
    }

    //determine the met of the event and its uncertainties
    //nominal MET value
    const pat::MET& met = (*mets).front();
    _met             = met.pt();
    _metPhi          = met.phi();
    /*
    _met = multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho).first;
    _metPhi = multilepAnalyzer->jec->correctedMETAndPhi(met, *jets, *rho).second;
    */
    //raw met values
    _metRaw = met.uncorPt();
    _metRawPhi = met.uncorPhi();
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
    _metSignificance = met.metSignificance(); 


    //std::cout << "end jetanalyzer" << std::endl;
    if(multilepAnalyzer->skim == "singlejet" and _nJets < 1) return false;
    if(multilepAnalyzer->skim == "FR" and _nJets < 1) return false;
    return true;
}

bool JetAnalyzer::jetIsLoose(const pat::Jet& jet, const bool is2017) const{

    //std::cout << "okisloose1" << std::endl;
    if(fabs(jet.eta()) <= 2.7){
    	//std::cout << "okislooseif1" << std::endl;
        if(jet.neutralHadronEnergyFraction() >=  0.99) return false;
        if(jet.neutralEmEnergyFraction() >= 0.99) return false;
        if(jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1) return false;
        if(fabs(jet.eta()) <= 2.4){
            if(jet.chargedHadronEnergyFraction() <= 0) return false;
            if(jet.chargedMultiplicity() <= 0) return false;
            if( !is2017 && ( jet.chargedEmEnergyFraction() >= 0.99 ) ) return false;
        }
    	//std::cout << "okislooseif1_2" << std::endl;

    } else if(fabs(jet.eta()) <= 3.0){
    	//std::cout << "okislooseif2" << std::endl;
        if(jet.neutralHadronEnergyFraction() >= 0.98) return false;
        if( !is2017 && (jet.neutralEmEnergyFraction() <= 0.01) ) return false;
        if( is2017 && (jet.neutralEmEnergyFraction() <= 0.02) ) return false;
        if(jet.neutralMultiplicity() <= 2) return false;
    	//std::cout << "okislooseif2_2" << std::endl;

    } else {
    	//std::cout << "okislooseif3" << std::endl;
        if(jet.neutralEmEnergyFraction() >= 0.90) return false;
        if(jet.neutralMultiplicity() <= 10) return false;
        if(is2017 && jet.neutralHadronEnergyFraction() <= 0.02) return false;
    	//std::cout << "okislooseif3_2" << std::endl;
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
