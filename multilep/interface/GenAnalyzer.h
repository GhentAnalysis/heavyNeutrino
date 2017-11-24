#ifndef GEN_ANALYZER_H
#define GEN_ANALYZER_H
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"

#include "TTree.h"

class multilep;
class PhotonAnalyzer;
class GenAnalyzer {
  //class friends
  friend PhotonAnalyzer;

  private:
    static const unsigned gen_nL_max = 20;
    static const unsigned gen_nPh_max = 10;
    static const unsigned gen_n_max = 20;
    static const unsigned gen_ndtr_max = 100;
   
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
    bool     _gen_lIsPrompt[gen_nL_max];
    bool     _gen_lPassParentage[gen_nPh_max];
    double   _gen_lMinDeltaR[gen_nPh_max];

    //Generator HT (needed when merging HT binned sample with inclusive one)
    double _gen_HT;

    //Functions to find the mother of a gen particle
    const reco::GenParticle* getMother(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    const int                getMotherPdgId(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    void                     getMotherList(const reco::GenParticle&, const std::vector<reco::GenParticle>&, std::vector<int>&);
    void		     getDaughterList(const reco::GenParticle&, const std::vector<reco::GenParticle>&, std::vector<reco::GenParticle>&, std::vector<int>&);
    void		     removeDoubleCountedDaughters(std::vector<reco::GenParticle>&);
    int 		     check_for_daughter(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    bool                     inMotherList(std::vector<int>& list, int i);

    unsigned                 ttgEventType(const std::vector<reco::GenParticle>& genParticles, double ptCut, double etaCut);
    double                   getMinDeltaR(const reco::GenParticle& p, const std::vector<reco::GenParticle>& genParticles);

    multilep* multilepAnalyzer;

  public:
    GenAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~GenAnalyzer(){};

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);

    unsigned _gen_nL;
    double   _gen_lPt[gen_nL_max];
    double   _gen_lEta[gen_nL_max];
    double   _gen_lPhi[gen_nL_max];

    unsigned _gen_nPh;
    double   _gen_phPt[gen_nPh_max];
    double   _gen_phEta[gen_nPh_max];
    double   _gen_phPhi[gen_nPh_max];

    unsigned _gen_nW;
    unsigned _gen_WMomPdg[gen_n_max];
    unsigned _gen_nWfromN;

    unsigned _gen_nq[6];
    unsigned _gen_nN;
    unsigned _gen_nNdaughters;
    unsigned _gen_Ndaughters_pdg[gen_n_max];
    unsigned _gen_nstatus23;
    unsigned _gen_nstatus23_fromNorW;
    unsigned _gen_nstatus23_fromN;
    unsigned _gen_nstatus23_fromW;
    unsigned _gen_status23_pdg[gen_n_max];
    unsigned _gen_status23_fromNorW_mompdg[gen_n_max];
    unsigned _gen_status23_fromNorW_pdg[gen_n_max];
    unsigned _gen_status23_fromN_pdg[gen_n_max];
    unsigned _gen_status23_fromW_pdg[gen_n_max];
    unsigned _gen_nq23;
    double   _gen_qPt[gen_n_max];
    double   _gen_qEta[gen_n_max];
    double   _gen_qPhi[gen_n_max];
    double   _gen_qE[gen_n_max];
    unsigned _gen_nq1dtr;
    int	     _gen_q1dtr_status[gen_ndtr_max];
    int	     _gen_q1dtr_pdgid[gen_ndtr_max];
    double   _gen_q1dtr_Pt[gen_ndtr_max];
    double   _gen_q1dtr_Eta[gen_ndtr_max];
    double   _gen_q1dtr_Phi[gen_ndtr_max];
    double   _gen_q1dtr_E[gen_ndtr_max];
    unsigned _gen_nq2dtr;
    int	     _gen_q2dtr_status[gen_ndtr_max];
    int	     _gen_q2dtr_pdgid[gen_ndtr_max];
    double   _gen_q2dtr_Pt[gen_ndtr_max];
    double   _gen_q2dtr_Eta[gen_ndtr_max];
    double   _gen_q2dtr_Phi[gen_ndtr_max];
    double   _gen_q2dtr_E[gen_ndtr_max];

};
#endif
