#ifndef GEN_ANALYZER_H
#define GEN_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;
class GenAnalyzer {

  private:
    static const unsigned gen_nL_max = 20;
    static const unsigned gen_nPh_max = 10;
    static const unsigned gen_n_max = 1000;

    unsigned    _ttgEventType;
    unsigned    _zgEventType;

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
    double   _gen_lPt[gen_nL_max];
    double   _gen_lEta[gen_nL_max];
    double   _gen_lPhi[gen_nL_max];
    double   _gen_lE[gen_nL_max];
    unsigned _gen_lFlavor[gen_nL_max];
    int      _gen_lCharge[gen_nL_max];
    int      _gen_lMomPdg[gen_nL_max];
    bool     _gen_lDecayedHadr[gen_nL_max];
    bool     _gen_lIsPrompt[gen_nL_max];
    bool     _gen_lPassParentage[gen_nL_max];
    double   _gen_lMinDeltaR[gen_nL_max];
    double   _gen_lVisPt[gen_nL_max];
    double   _gen_lVisEta[gen_nL_max];
    double   _gen_lVisPhi[gen_nL_max];
    double   _gen_lVisE[gen_nL_max];
   
    //Generator particles (all)
    unsigned _gen_n = 0;
    double   _gen_pt[gen_n_max];
    double   _gen_eta[gen_n_max];
    double   _gen_phi[gen_n_max];
    double   _gen_E[gen_n_max];
    int      _gen_pdgId[gen_n_max];
    int      _gen_charge[gen_n_max];
    int      _gen_status[gen_n_max];
    bool     _gen_isPromptFinalState[gen_n_max];
    bool     _gen_isDirectPromptTauDecayProductFinalState[gen_n_max];
    bool     _gen_isLastCopy[gen_n_max];
    int      _gen_index[gen_n_max];
    int      _gen_motherIndex[gen_n_max];
    int      _gen_daughter_n[gen_n_max];
    int      _gen_daughterIndex[gen_n_max][10];
   
    unsigned overlapEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut, double genCone) const;
    double   getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles) const;

    TH1D*    hCounterDirac;
    bool     _gen_isDiracType;
    int      _gen_lProvenanceHNL[gen_nL_max];

    multilep* multilepAnalyzer;

  public:
    GenAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~GenAnalyzer(){};

    void beginJob(TTree* outputTree, edm::Service<TFileService>& fs);
    void analyze(const edm::Event&);
};
#endif
