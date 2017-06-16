#ifndef LHE_ANALYZER_H
#define LHE_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class LheAnalyzer {
  private:
    double _lheHTIncoming;
    double _ctauHN;
    double _lheWeight[110];
    TH1D* lheCounter;

    multilep* multilepAnalyzer;

  public:
    LheAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LheAnalyzer();

    void beginJob(TTree* outputTree, edm::Service<TFileService>& fs);
    void analyze(const edm::Event&);
};

#endif
