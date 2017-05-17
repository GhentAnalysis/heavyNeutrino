#ifndef JET_ANALYZER_H
#define JET_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class JetAnalyzer {
  private:
    static const unsigned nJets_max = 20;

    unsigned _nJets;
    double _jetPt[nJets_max];
    double _jetEta[nJets_max];
    double _jetPhi[nJets_max];
    double _jetE[nJets_max];

    multilep* multilepAnalyzer;

  public:
    JetAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~JetAnalyzer();

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
};

#endif
