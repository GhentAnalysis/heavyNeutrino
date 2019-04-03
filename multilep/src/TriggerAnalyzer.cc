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

  // MET Filters: first add common ones for 2016, 2017, 2018
  // MET filter are taken in AND (based on the occurence of capitalized 'MET' in the allFlags key) and always start with "Flag"
  // References: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
  allFlags["passMETFilters"] = {"Flag_goodVertices", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",  "Flag_EcalDeadCellTriggerPrimitiveFilter",
                                "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter"};

  if(multilepAnalyzer->isData){ // This one is only to be applied on data
    allFlags["passMETFilters"].push_back("Flag_eeBadScFilter");
  }

  if(multilepAnalyzer->is2017 || multilepAnalyzer->is2018){ // This one is only for 2017 and 2018 data
    allFlags["passMETFilters"].push_back("updated_ecalBadCalibFilter"); // improved version over the Flag_ecalBadCalibFilter, implementation manually
  }

  // Triggers, grouped per year
  // Triggers are taken in OR and always start with "HLT"
  // WARNING/TODO: the 2018 year is currently a placeholder using the 2017 triggers, please check before starting analysis
  // Also for the other years the triggers might be checked and optimized in detail, could use https://tomc.web.cern.ch/tomc/triggerPrescales to check prescales
  if (multilepAnalyzer->skim == "FR"){
    allFlags["FR_single_lepton"]= {"HLT_Mu8","HLT_Mu17","HLT_Mu3_PFJet40","HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele8_CaloIdM_TrackIdM_PFJet30"};
  }
  if (multilepAnalyzer->skim == "FRsM"){
    allFlags["FR_single_lepton"]= {"HLT_Mu3_PFJet40"};
  }
  if(multilepAnalyzer->is2018){
    allFlags["passTrigger_1l"]   = {"HLT_IsoMu24", "HLT_IsoMu27","HLT_Ele32_WPTight_Gsf"};
  } else if(multilepAnalyzer->is2017){
    allFlags["passTrigger_1l"]   = {"HLT_IsoMu27", "HLT_IsoMu24", "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG"};
  } else {
    allFlags["passTrigger_1l"]   = {"HLT_Ele27_WPTight_Gsf","HLT_IsoMu24", "HLT_IsoTkMu24"};
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

  if(multilepAnalyzer->is2017 || multilepAnalyzer->is2018){ // The updated ecalBadCalibFilter
    edm::Handle<bool> passEcalBadCalibFilterUpdate; iEvent.getByToken(multilepAnalyzer->ecalBadCalibFilterToken, passEcalBadCalibFilterUpdate);
    flag["updated_ecalBadCalibFilter"] = (*passEcalBadCalibFilterUpdate);
  }

  // Get all flags
  getResults(iEvent, triggerResults,                                                               triggersToSave, true);
  getResults(iEvent, recoResultsPrimary.failedToGet() ? recoResultsSecondary : recoResultsPrimary, filtersToSave,  false);

  reIndex = false;

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
