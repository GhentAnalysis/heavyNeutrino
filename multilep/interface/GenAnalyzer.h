#ifndef GEN_ANALYZER_H
#define GEN_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class GenAnalyzer {
  private:
    static const unsigned nGen_max = 200;

    unsigned _nGen;

    double _genPt[nGen_max];
    double _genEta[nGen_max];
    double _genPhi[nGen_max];
    double _genMass[nGen_max];
    int    _genCharge[nGen_max];
    int    _genPdgId[nGen_max];
    int    _genStatus[nGen_max];
    bool   _genFromHardProcess[nGen_max];
    bool   _genIsPrompt[nGen_max];
    int    _genMotherIndex1[nGen_max];
    int    _genMotherIndex2[nGen_max];
    int    _genDaughterIndex1[nGen_max];
    int    _genDaughterIndex2[nGen_max];

    multilep* multilepAnalyzer;

  public:
    GenAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~GenAnalyzer();

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
};

#endif
