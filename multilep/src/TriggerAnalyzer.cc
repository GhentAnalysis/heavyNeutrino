#include "heavyNeutrino/multilep/interface/TriggerAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"

/*
 * Triggers and MET filters
 * Simply add your triggers to the list below, the key in allFlags[key] takes the OR of the following triggers
 * Currently not only the combined but also the individual triggers are stored, keeping the possibility for trigger studies
 * Might add a flag to switch all the storage of all those individual paths off
 * Note: use "pass" in the combined flag if you want to use it recursively
 */

TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer){
  allFlags["hlt_pfmet"] = {"HLT_PFMET120_PFMHT120_IDTight"};
  //2017 data and MC   
  if(multilepAnalyzer->is2017 || multilepAnalyzer->is2018){
    allFlags["passMETFilters"] = {"Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",                // MET filters
                                  "Flag_goodVertices", "Flag_globalTightHalo2016Filter", "Flag_ecalBadCalibFilter",                             //for Moriond 2018 superTightHalo -> TightHalo
                                  "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter"}; //, "Flag_badMuons", "Flag_duplicateMuons"};      // Duplicate muons still missing in mAOD, not sure how to get those in
    //met filters only to be applied on fullsim and data
    if(!multilepAnalyzer->isSUSY){
        allFlags["passMETFilters"].push_back("Flag_globalTightHalo2016Filter");
    }
    //extra filters only to be applied on data:
    if(multilepAnalyzer->isData){
        allFlags["passMETFilters"].push_back("Flag_eeBadScFilter");
    }

    allFlags["passTrigger_m"]         = {"HLT_IsoMu24", "HLT_IsoMu24_eta2p1", "HLT_IsoMu27", "HLT_IsoMu30", "HLT_Mu50", "HLT_Mu55"};

    allFlags["passTrigger_e"]         = {"HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf", "HLT_Ele38_WPTight_Gsf", "HLT_Ele40_WPTight_Gsf"};

    allFlags["passTrigger_mm"]        = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                                  "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"};

    allFlags["passTrigger_em"]        = {"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"};

    allFlags["passTrigger_ee"]        = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};

                                    //Is it useful to also store MediumChargedIso and TightChargedIso versions of these tau triggers?
    allFlags["passTrigger_et"]        = {"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1", "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1"};

    allFlags["passTrigger_mt"]        = {"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1", 
                                  "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1", "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1", "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1",
                                  "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1"}; 

    allFlags["passTrigger_t"]         = {"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr"};

    allFlags["passTrigger_mmm"]       = {"HLT_TripleMu_10_5_5_DZ", "HLT_TripleMu_5_3_3_Mass3p8to60_DZ", "TripleMu_12_10_5"};
    allFlags["passTrigger_emm"]       = {"HLT_DiMu9_Ele9_CaloIdL_TrackIdL", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ"};
    allFlags["passTrigger_eem"]       = {"HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ"};
    allFlags["passTrigger_eee"]       = {"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"};                                                                    //Bullshit trigger because L1 seeds are higher than HLT, be careful using it

  //2016 data and MC
  } else {
    allFlags["passMETFilters"] = {"Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",                // MET filters (if legacy mAOD vbecomes available, copy the filters listed for 2017)
                                  "Flag_goodVertices", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter"};
    //met filters only to be applied on fullsim and data
    if(!multilepAnalyzer->isSUSY){
        allFlags["passMETFilters"].push_back("Flag_globalTightHalo2016Filter");
    }
    //extra filters only to be applied on data:
    if(multilepAnalyzer->isData){
        allFlags["passMETFilters"].push_back("Flag_eeBadScFilter");
    }

    allFlags["passTrigger_e"]      = {"HLT_Ele27_WPTight_Gsf", "HLT_Ele105_CaloIdVT_GsfTrkIdT", "HLT_Ele115_CaloIdVT_GsfTrkIdT"};
    allFlags["passTrigger_m"]      = {"HLT_IsoMu24", "HLT_IsoTkMu24", "HLT_Mu50", "HLT_TkMu50", "HLT_Mu45_eta2p1"};
    allFlags["passTrigger_ee"]     = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
                               "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};
    allFlags["passTrigger_em"]     = {"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                               "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
                               "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL"};
    allFlags["passTrigger_mm"]     = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                               "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                               "HLT_Mu30_TkMu11"};
    allFlags["passTrigger_et"]     = {"HLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1", "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30", "HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
                               "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20"};
    allFlags["passTrigger_mt"]     = {"HLT_IsoMu19_eta2p1_LooseIsoPFTau20"};
    allFlags["passTrigger_eee"]    = {"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"};
    allFlags["passTrigger_eem"]    = {"HLT_Mu8_DiEle12_CaloIdL_TrackIdL"};
    allFlags["passTrigger_emm"]    = {"HLT_DiMu9_Ele9_CaloIdL_TrackIdL"};
    allFlags["passTrigger_mmm"]    = {"HLT_TripleMu_12_10_5"};
  }
}

void TriggerAnalyzer::beginJob(TTree* outputTree){
  reIndex = true;
  initList(triggersToSave, "HLT");
  initList(filtersToSave, "Flag");

  for(TString f : getAllFlags()){
    outputTree->Branch("_" + f, &flag[f], "_" + f + "/O");
    if(f.Contains("HLT")) outputTree->Branch("_" + f + "_prescale", &prescale[f], "_" + f + "_prescale/I");
  }
}


/*
 * Filters the HLT and MET flags from allFlags (maybe better to change to set...)
 */
void TriggerAnalyzer::initList(std::vector<TString>& list, TString identifier){
  list.clear();
  for(auto& v : allFlags){
    for(auto& t : v.second){
      if(t.Contains(identifier)){
        if(std::find(list.begin(), list.end(), t) == list.end()){
          list.push_back(t);
        }
      }
    }
  }
}


std::vector<TString> TriggerAnalyzer::getAllFlags(){
  std::vector<TString> list;
  for(auto& v : allFlags){
    list.push_back(v.first);
    for(auto& t : v.second){
      if(std::find(list.begin(), list.end(), t) == list.end()){
        list.push_back(t);
      }
    }
  }
  return list;
}

bool TriggerAnalyzer::passCombinedFlagOR(TString combinedFlag){
  for(auto& f : allFlags[combinedFlag]){
    if(f.Contains("pass") and passCombinedFlagOR(f)) return true;
    else if(flag[f])                                 return true;
  }
  return false;
}

bool TriggerAnalyzer::passCombinedFlagAND(TString combinedFlag){
  for(auto& f : allFlags[combinedFlag]){
    if(f.Contains("pass") and not passCombinedFlagAND(f)) return false;
    else if(not flag[f])                                  return false;
  }
  return true;
}

void TriggerAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<edm::TriggerResults> recoResultsPrimary;   iEvent.getByToken(multilepAnalyzer->recoResultsPrimaryToken,   recoResultsPrimary);
  edm::Handle<edm::TriggerResults> recoResultsSecondary; iEvent.getByToken(multilepAnalyzer->recoResultsSecondaryToken, recoResultsSecondary);
  edm::Handle<edm::TriggerResults> triggerResults;       iEvent.getByToken(multilepAnalyzer->triggerToken,              triggerResults);
  edm::Handle<bool> badPFMuonFilter;
  edm::Handle<bool> badChCandFilter;

  if(   !(multilepAnalyzer->is2017 || multilepAnalyzer->is2018)   ){
    iEvent.getByToken(multilepAnalyzer->badPFMuonFilterToken,  badPFMuonFilter);
    iEvent.getByToken(multilepAnalyzer->badChCandFilterToken,  badChCandFilter);
  }

  // Get all flags
  getResults(iEvent, triggerResults,                                                               triggersToSave, true);
  getResults(iEvent, recoResultsPrimary.failedToGet() ? recoResultsSecondary : recoResultsPrimary, filtersToSave,  false);

  reIndex = false;
 
  // These met filters are to be applied on-the-fly in 2016 samples
  if(    !(multilepAnalyzer->is2017 || multilepAnalyzer->is2018)   ){
    flag["Flag_BadPFMuonFilter"] = *badPFMuonFilter;
    flag["Flag_BadChargedCandidateFilter"] = *badChCandFilter;
  }

  for(auto& combinedFlag : allFlags){
    if(combinedFlag.first.Contains("MET")) flag[combinedFlag.first] = passCombinedFlagAND(combinedFlag.first);
    else                                   flag[combinedFlag.first] = passCombinedFlagOR(combinedFlag.first);
  }
}

/*
 * Call this at the first event, checks for available triggers and warns for missing triggers
 * Stores indexes of wanted triggers, which minimizes string comparisons for the next events
 */
void TriggerAnalyzer::indexFlags(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& results, std::vector<TString>& toSave){
  for(TString t : toSave) index[t] = -1;

  std::cout << "Available triggers:" << std::endl;
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*results);
  for (unsigned int i = 0; i < results->size(); ++i){
    std::cout << "  " << triggerNames.triggerName(i);
    for(TString t : toSave){
      TString tt = (t.Contains("HLT") ? t + "_v" : t);
      if(TString(triggerNames.triggerName(i)).Contains(tt)){
        index[t] = i;
        std::cout << "     --> saving to tree";
      }
    }
    std::cout << std::endl;
  }

  for(TString t : toSave){
    if(index[t] == -1) std::cout << "WARNING: " << t << " not found in triggerresult, please check!" << std::endl;
  }
}


/*
 * Saving triggers and prescales
 */
void TriggerAnalyzer::getResults(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& results, std::vector<TString>& toSave, const bool savePrescales){
  if(results.failedToGet()) return;

  if(reIndex) indexFlags(iEvent, results, toSave);

  for(TString t : toSave){
    if(index[t] == -1) flag[t] = false; 
    else               flag[t] = results->accept(index[t]);
  }

  if(savePrescales){
    edm::Handle<pat::PackedTriggerPrescales> prescales;
    iEvent.getByToken(multilepAnalyzer->prescalesToken, prescales);
    for(TString t : toSave){
      if(index[t] == -1) prescale[t] = -1;
      else               prescale[t] = prescales->getPrescaleForIndex(index[t]);
    }
  }
}
