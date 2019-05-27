#ifndef PARTICLE_LEVEL_ANALYZER_H
#define PARTICLE_LEVEL_ANALYZER_H

/*
Class storing data for unfolding to particle-level in differential measurement
*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;
class ParticleLevelAnalyzer {

  private:
    static const unsigned pl_nL_max = 10;
    static const unsigned pl_nPh_max = 5;
    static const unsigned pl_nJet_max = 10;
   
    //particle level MET
    double   _pl_met;
    double   _pl_metPhi;

    //particle level photons
    unsigned _pl_nPh = 0;
    double   _pl_phPt[pl_nPh_max];
    double   _pl_phEta[pl_nPh_max];
    double   _pl_phPhi[pl_nPh_max];
    double   _pl_phE[pl_nPh_max];

    //particle level leptons
    unsigned _pl_nL = 0;
    double   _pl_lPt[pl_nL_max];
    double   _pl_lEta[pl_nL_max];
    double   _pl_lPhi[pl_nL_max];
    double   _pl_lE[pl_nL_max];
    unsigned _pl_lFlavor[pl_nL_max];
    int      _pl_lCharge[pl_nL_max];

    //particle level jets
    unsigned _pl_nJets = 0;
    double   _pl_jetPt[pl_nJet_max];
    double   _pl_jetEta[pl_nJet_max];
    double   _pl_jetPhi[pl_nJet_max];
    double   _pl_jetE[pl_nJet_max];
    unsigned _pl_jetHadronFlavor[pl_nJet_max];

    multilep* multilepAnalyzer;

  public:
    ParticleLevelAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~ParticleLevelAnalyzer(){};

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&);
};
#endif
