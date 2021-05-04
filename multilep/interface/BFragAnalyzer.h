#ifndef B_FRAG_ANALYZER_H
#define B_FRAG_ANALYZER_H

/*
Class storing data for b fragmentation systematics
*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;
class BFragAnalyzer {

  private:

    double _bf_fragCP5BL = 0;
    double _bf_fragCP5BLdown = 0;
    double _bf_fragCP5BLup = 0;
    double _bf_fragCP5Peterson = 0;
    double _bf_fragCP5Petersondown = 0;
    double _bf_fragCP5Petersonup = 0;

    multilep* multilepAnalyzer;

  public:
    BFragAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~BFragAnalyzer(){};

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
};
#endif
