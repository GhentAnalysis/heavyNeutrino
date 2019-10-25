#ifndef TRIGGER_ANALYZER_H
#define TRIGGER_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class TriggerAnalyzer {
  private:
    friend class multilep;
  
    std::map<TString, std::vector<TString>> allFlags;
    std::vector<TString> triggersToSave;
    std::vector<TString> filtersToSave;

    std::map<TString, bool> flag;
    std::map<TString, int>  prescale;
    std::map<TString, int>  index;

    multilep* multilepAnalyzer;

    void indexFlags(const edm::Event&, edm::Handle<edm::TriggerResults>&, std::vector<TString>&);
    void getResults(const edm::Event&, edm::Handle<edm::TriggerResults>&, std::vector<TString>&, const bool);

    void initList(std::vector<TString>&, TString);
    std::vector<TString> getAllFlags();

    bool passCombinedFlagAND(TString combinedFlag);
    bool passCombinedFlagOR(TString combinedFlag);

    bool passEle32WPTight(const edm::Event&, edm::Handle<edm::TriggerResults>&);
    std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjects(const pat::Electron&, const std::vector<pat::TriggerObjectStandAlone>&,const float);

  public:
    TriggerAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~TriggerAnalyzer(){};

    bool reIndex;
    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
};

#endif
