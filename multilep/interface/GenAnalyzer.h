#ifndef GEN_ANALYZER_H
#define GEN_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;
class GenAnalyzer {

  private:
    static const unsigned gen_nL_max = 50;
    static const unsigned gen_nPh_max = 10;
   
    unsigned    _ttgEventType;
    unsigned    _zgEventType;
    unsigned    _zgOldEventType;
    unsigned    _hasInternalConversion;

    //generator level MET
    double   _gen_met;
    double   _gen_metPhi;

    //Generator photons
    unsigned _gen_nPh = 0;
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
    unsigned _gen_nL = 0;
    double   _gen_pdgID[gen_nL_max]; // displaced specific, possibly copy from flavor but in strange double format FIXME: get rid of this
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
    bool     _gen_lDecayedHadr[gen_nL_max];
    bool     _gen_lIsPrompt[gen_nL_max];
    bool     _gen_lPassParentage[gen_nL_max];
    double   _gen_lMinDeltaR[gen_nL_max];

    unsigned overlapEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut, double genCone) const;
    bool photonToInternalConversion(const reco::GenParticle& photon, const std::vector<reco::GenParticle>& genParticles) const;

    // Array of pointers to genLeptons (NOT saved in the tree!)
    // (only charged leptons for now, no photons)
    const reco::GenParticle* _gen_lRefs[gen_nL_max];

    //Functions to find the mother of a gen particle
    double   getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles) const;

    multilep* multilepAnalyzer;

  public:
    GenAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~GenAnalyzer(){};

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
    unsigned getGenLeptonIndex(const reco::GenParticle* match);
};
#endif
