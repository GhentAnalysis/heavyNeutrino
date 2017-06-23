#include "heavyNeutrino/multilep/interface/TriggerAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"

/*
 * Triggers and MET filters
 * Simply add your triggers to the list below, the key in allFlags[key] takes the OR of the following triggers
 * Currently not only the combined but also the individual triggers are stored, keeping the possibility for trigger studies
 * Might add a flag to switch all the storage of all those individual paths off
 */

TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer){
  
  allFlags["passMETFilters"] = {"Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",                // MET filters
                                "Flag_goodVertices", "Flag_eeBadScFilter", "Flag_globalTightHalo2016Filter",                                 // Not sure if eeBadScFilter is still needed
                                "flag_badPFMuonFilter","flag_badChCandFilter"};                                                              // TODO: Check if those are in the miniAOD, and where is the duplicateMuons one?

  allFlags["passHN_1l"]      = {"HLT_Ele27_WPTight_Gsf", "HLT_IsoMu24", "HLT_IsoTkMu24"};                                                    // HN 1l triggers
  allFlags["passHN_eee"]     = {"passHN_1l","HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};           // HN eee
  allFlags["passHN eem"]     = {"passHN_1l","HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",         // HN emm
                                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", 
                                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};
  allFlags["passHN_emm"]     = {"passHN_1l","HLT_DiMu9_Ele9_CaloIdL_TrackIdL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                        // HN emm
                                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", 
                                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", 
                                "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL",  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 
                                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", 
                                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL"};
  allFlags["passHN_mmm"]     = {"passHN_1l","HLT_TripleMu_12_10_5", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                                   // HN mmm
                                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",  
                                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", 
                                "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL"};
  allFlags["passMET"]        = {"HLT_PFMET120_PFMHT120_IDTight", "HLT_PFMET110_PFMHT110_IDTight", "HLT_PFMET170_BeamHaloCleaned",            // To trigger in MET dataset for trigger efficiencies
                                "HLT_PFMET170_HBHECleaned", "HLT_PFMET170_HBHE_BeamHaloCleaned", "HLT_PFMET170_NotCleaned",
                                "HLT_MET200","HLT_MET250","HLT_MET300"};
  allFlags["passTTG_ee"]     = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",                    // TTG e (same as in SUS-17-001)
                                "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW"};
  allFlags["passTTG_e"]      = {"HLT_Ele105_CaloIdVT_GsfTrkIdT", "HLT_Ele115_CaloIdVT_GsfTrkIdT"};                                           // TTG e
  allFlags["passTTG_mm"]     = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu30_TkMu11",                // TTG mm
                                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"};
  allFlags["passTTG_m"]      = {"HLT_Mu50", "HLT_TkMu50", "HLT_Mu45_eta2p1"};                                                                // TTG m
  allFlags["passTTG_em"]     = {"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",       // TTG em
                                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                                "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL"};
};

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

bool TriggerAnalyzer::passCombinedFlag(TString combinedFlag){
  for(auto& f : allFlags[combinedFlag]){
    if(f.Contains("pass") and passCombinedFlag(f)) return true;
    else if(flag[f])                               return true;
  }
  return false;
}

void TriggerAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<edm::TriggerResults> recoResults;       iEvent.getByToken(multilepAnalyzer->recoResultsToken,      recoResults);
  edm::Handle<edm::TriggerResults> triggerResults;    iEvent.getByToken(multilepAnalyzer->triggerToken,          triggerResults);
  edm::Handle<bool> badPFMuonFilter;                  iEvent.getByToken(multilepAnalyzer->badPFMuonFilterToken,  badPFMuonFilter);
  edm::Handle<bool> badChCandFilter;                  iEvent.getByToken(multilepAnalyzer->badChCandFilterToken,  badChCandFilter);

  // Get all flags
  getResults(iEvent, triggerResults, triggersToSave, true);
  getResults(iEvent, recoResults,    filtersToSave,  false);
  flag["passBadPFMuonFilter"] = *badPFMuonFilter;
  flag["passBadChCandFilter"] = *badChCandFilter;

  for(auto& combinedFlag : allFlags) flag[combinedFlag.first] = passCombinedFlag(combinedFlag.first);
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
  reIndex = false;
}


/*
 * Saving triggers and prescales
 */
void TriggerAnalyzer::getResults(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& results, std::vector<TString>& toSave, bool savePrescales){
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
