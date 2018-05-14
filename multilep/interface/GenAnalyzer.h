#ifndef GEN_ANALYZER_H
#define GEN_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"

#include "TTree.h"

class multilep;
class PhotonAnalyzer;
class GenMatching;
class GenAnalyzer {
  //class friends
  friend PhotonAnalyzer;
  friend GenMatching;

  private:
    static const unsigned gen_nL_max = 20;
    static const unsigned gen_nPh_max = 10;
   
    unsigned    _ttgEventType;
    unsigned    _zgEventType;

    //generator level MET
    double   _gen_met;
    double   _gen_metPhi;

    //Generator photons
    unsigned _gen_nPh;
    unsigned _gen_phStatus[gen_nPh_max];
    double   _gen_phPt[gen_nPh_max];
    double   _gen_phEta[gen_nPh_max];
    double   _gen_phPhi[gen_nPh_max];
    double   _gen_phE[gen_nPh_max];
    int      _gen_phMomPdg[gen_nPh_max];
    bool     _gen_phIsPrompt[gen_nPh_max];
    bool     _gen_phPassParentage[gen_nPh_max];
    double   _gen_phMinDeltaR[gen_nPh_max];

    //Generator leptons
    unsigned _gen_nL;
    double   _gen_lPt[gen_nL_max];
    double   _gen_lEta[gen_nL_max];
    double   _gen_lPhi[gen_nL_max];
    double   _gen_lE[gen_nL_max];
    unsigned _gen_lFlavor[gen_nL_max];
    int      _gen_lCharge[gen_nL_max];
    int      _gen_lMomPdg[gen_nL_max];
    double   _gen_vertex_x[gen_nL_max];
    double   _gen_vertex_y[gen_nL_max];
    double   _gen_vertex_z[gen_nL_max];
    bool     _gen_lIsPrompt[gen_nL_max];
    bool     _gen_lPassParentage[gen_nL_max];
    double   _gen_lMinDeltaR[gen_nL_max];

    // Array of pointers to genLeptons (NOT saved in the tree!)
    // (only charged leptons for now, no photons)
    //const reco::GenParticle* _gen_refs[gen_nL_max+gen_nPh_max];
    const reco::GenParticle* _gen_lRefs[gen_nL_max];
    // const reco::GenParticle* _gen_phRefs[gen_nPh_max];

   //Generator HT (needed when merging HT binned sample with inclusive one)
    double _gen_HT;

    //Functions to find the mother of a gen particle
    unsigned ttgEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut) const;
    double   getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles) const;

    multilep* multilepAnalyzer;

  public:
    GenAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~GenAnalyzer(){};

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
};
#endif
