#include "heavyNeutrino/multilep/plugins/multilep.h"

void multilep::fillTriggerVars(const edm::Event& iEvent){
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByToken(triggerToken_, triggerResults);
	for(unsigned f = 0; f < 4; ++f){
		_passHnlTrigger[f] = trigPass(f, triggerResults, iEvent);
	}
}

bool multilep::trigPass(unsigned flavorComb, edm::Handle<edm::TriggerResults>& triggerResults, const edm::Event& iEvent){
	static const char* triggerPath_all[3] =  {"HLT_Ele27_WPTight", "HLT_IsoMu24_v", "HLT_IsoTkMu24_v"};
  	static const char* triggerPath_eee[2] = {"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
	static const char* triggerPath_eem[6] = {"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
	static const char* triggerPath_emm[11]	= {"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"};
	static const char* triggerPath_mmm[7] = {"HLT_TripleMu_12_10_5_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"};
	if(triggerResults.failedToGet()) return false;
	const edm::TriggerNames& triggerN = iEvent.triggerNames(*triggerResults);
	for(unsigned p = 0; p < triggerResults->size(); ++p){
	   	if(triggerResults.product()->accept(p)) {
		    std::string trigPath = triggerN.triggerName(p); 
			for(unsigned t = 0; t < 3; ++t){
				if(trigPath.find(triggerPath_all[t]) != std::string::npos) return true;
			}
			switch(flavorComb){
				case 0:
					for(unsigned t = 0; t < 2; ++t){
						if(trigPath.find(triggerPath_eee[t]) != std::string::npos) return true;
					}
				case 1:
					for(unsigned t = 0; t < 6; ++t){
						if(trigPath.find(triggerPath_eem[t]) != std::string::npos) return true;
					}
				case 2:
					for(unsigned t = 0; t < 11; ++t){
						if(trigPath.find(triggerPath_emm[t]) != std::string::npos) return true;
					}
				case 3:
					for(unsigned t = 0; t < 7; ++t){
						if(trigPath.find(triggerPath_mmm[t]) != std::string::npos) return true;
					}
				default: return false;
			}
		}
	}
	return false;
}

void multilep::fillMetFilterVars(const edm::Event& iEvent){
	_metFiltersFlagged = metFiltersFlagged(iEvent);
}

bool multilep::metFiltersFlagged(const edm::Event& iEvent){
	//MET filters included in MiniAOD
	static const char* metFilters[6] = {"Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter", "Flag_globalTightHalo2016Filter"};
	edm::Handle<edm::TriggerResults> recoResults;
	iEvent.getByToken(recoResultsToken_, recoResults);

    if(recoResults.failedToGet()) return false;
    const edm::TriggerNames& filterNames = iEvent.triggerNames(*recoResults);
	for(unsigned p = 0; p < recoResults->size(); ++p) {
    	if(!recoResults.product()->accept(p)){
			std::string filterPath = filterNames.triggerName(p);
			for(unsigned f = 0; f < 6; ++f){
				if(filterPath.find(metFilters[6]) != std::string::npos) return true;
			}
		}
	}
	//MET filters to be read separately
	//bad PF muon filter
	edm::Handle<bool> badPFMuonFilter;
	iEvent.getByToken(badPFMuonFilterToken_, badPFMuonFilter);
	if(!(*badPFMuonFilter)) return true;
	//bad charged candidate filter
	edm::Handle<bool> badChCandFilter;
	iEvent.getByToken(badChCandFilterToken_, badChCandFilter);
	if(!(*badChCandFilter)) return true;

	return false;
}

