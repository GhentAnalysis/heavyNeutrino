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
    static const unsigned gen_n_max = 20;
    static const unsigned gen_ndtr_max = 50;
   
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
    double   _gen_vertex_x[gen_nL_max];
    double   _gen_vertex_y[gen_nL_max];
    double   _gen_vertex_z[gen_nL_max];
    bool     _gen_lDecayedHadr[gen_nL_max];
    bool     _gen_lIsPrompt[gen_nL_max];
    bool     _gen_lPassParentage[gen_nL_max];
    double   _gen_lMinDeltaR[gen_nL_max];


    //Generator quarks
    unsigned _gen_nN = 0;
    double   _gen_NPt;
    double   _gen_NEta;
    double   _gen_NPhi;
    double   _gen_NE;
    double   _gen_Nvertex_x;
    double   _gen_Nvertex_y;
    double   _gen_Nvertex_z;

    unsigned _gen_nNPackedDtrs = 0;
    double   _gen_NPackedDtrsPt[gen_ndtr_max];
    double   _gen_NPackedDtrsEta[gen_ndtr_max];
    double   _gen_NPackedDtrsPhi[gen_ndtr_max];
    double   _gen_NPackedDtrsE[gen_ndtr_max];
    int      _gen_NPackedDtrsPdgId[gen_ndtr_max];
    int      _gen_NPackedDtrsCharge[gen_ndtr_max];

    unsigned _gen_nNdaughters = 0;
    int      _gen_Ndaughters_pdg[gen_n_max];
    double   _gen_Ndaughters_Pt[gen_n_max];
    double   _gen_Ndaughters_Eta[gen_n_max];
    double   _gen_Ndaughters_Phi[gen_n_max];
    double   _gen_Ndaughters_E[gen_n_max];
    int      _gen_Ndaughters_Charge[gen_n_max];
    
    unsigned _gen_nq = 0;
    double   _gen_qPt[gen_n_max];
    double   _gen_qEta[gen_n_max];
    double   _gen_qPhi[gen_n_max];
    double   _gen_qE[gen_n_max];
    
    //Functions to find the mother of a gen particle
    const reco::GenParticle* getMother(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    const int getMotherPdgId(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    void     getMotherList(const reco::GenParticle&, const std::vector<reco::GenParticle>&, std::vector<int>&);
    unsigned ttgEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut) const;
    
    unsigned overlapEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut) const;
    double   getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles) const;
    bool     isAncestor(const reco::Candidate*, const reco::Candidate*);

    multilep* multilepAnalyzer;

  public:
    GenAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~GenAnalyzer(){};

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);
};
#endif
