#ifndef LHE_ANALYZER_H
#define LHE_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;
class SUSYMassAnalyzer;

class LheAnalyzer {
  public:
    double getWeight() const;
  private:
    float  _nTrueInt;
    double _weight;
    double _lheHTIncoming;
    double _ctauHN;

    TH1D*  hCounter;
    TH1D*  lheCounter;
    TH1D*  psCounter;
    TH1D*  nTrueInteractions;

    unsigned _nLheWeights;
    double _lheWeight[110];

    unsigned _nPsWeights;
    double _psWeight[14];

    multilep* multilepAnalyzer;

  public:
    LheAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LheAnalyzer(){};

    void beginJob(TTree* outputTree, edm::Service<TFileService>& fs);
    void analyze(const edm::Event&);
};
#endif
