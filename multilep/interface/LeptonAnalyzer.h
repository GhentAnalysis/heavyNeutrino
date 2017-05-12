#ifndef LEPTON_ANALYZER_H
#define LEPTON_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"


#include "../plugins/multilep.h"
#include "TTree.h"

/*
 * Functions for electron identification
 */
class multilep;

class LeptonAnalyzer {
  private:
    EffectiveAreas electronsEffectiveAreas;
    EffectiveAreas muonsEffectiveAreas;

    static const unsigned nL_max      = 20;                                                          //maximum number of particles stored
    unsigned _nL;                                                                                    //number of leptons
    unsigned _nMu;
    unsigned _nEle;
    unsigned _nLight;
    unsigned _nTau;
    double _lPt[nL_max];                                                                             //lepton kinematics
    double _lEta[nL_max];                                                                            //QUESTION: what was the reason again to use arrays instead of vectors?
    double _lPhi[nL_max];
    double _lE[nL_max];
    unsigned _lFlavor[nL_max];                                                                       //other lepton variables
    int _lCharge[nL_max];
    double _relIso[nL_max];
    double _miniIso[nL_max];
    double _dxy[nL_max];
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _3dIPSig[nL_max];
    float _lElectronMva[nL_max];
    bool _lHNLoose[nL_max];                                                                          //lepton selection decisions
    bool _lHNFO[nL_max];
    bool _lHNTight[nL_max];
    bool _isPrompt[nL_max];                                                                          //generator lepton variables
    bool _truthPdg[nL_max];
    bool _truthMomPdg[nL_max];
    unsigned _origin[nL_max];

    multilep* multilepAnalyzer;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    void analyze(const edm::Event&);

    float dEtaInSeed(const pat::Electron*);
    bool isLooseCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool isTightCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool passTriggerEmulationDoubleEG(const pat::Electron*);
    float slidingCut(float, float, float);
    bool passingElectronMvaHZZ(const pat::Electron*, double);
    bool passingElectronMvaLooseSusy(const pat::Electron*, double, double);
    bool passingElectronMvaMediumSusy(const pat::Electron*, double);
    bool passingElectronMvaTightSusy(const pat::Electron*, double);
    bool passingElectronMvaHeavyNeutrinoFO(const pat::Electron*, double);
    void fillLeptonGenVars(const reco::GenParticle* genParticle);                                    //Fill MC-truth lepton variables
    void fillLeptonKinVars(const reco::Candidate&);                                                  //Fill reconstructed lepton kinematics

    double getRelIso03(const pat::Muon&, const double);
    //double getIsoAlt(const pat::Muon&, double);
    double getRelIso03(const pat::Electron&, const double);
    double getRelIso(const reco::RecoCandidate&, const std::vector<pat::PackedCandidate>& , const double, const double);
    double getMiniIso(const reco::RecoCandidate&, const std::vector<pat::PackedCandidate>&, const double, const double);

    bool eleMuOverlap(const pat::Electron& ele);
    bool isHNLoose(const pat::Electron& lepton);
    bool isHNLoose(const pat::Muon& lepton);
    bool isHNFO(const pat::Electron& lepton);
    bool isHNFO(const pat::Muon& lepton);
    bool isHNTight(const pat::Electron& lepton);
    bool isHNTight(const pat::Muon& lepton);
};

#endif
