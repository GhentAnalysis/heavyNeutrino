#ifndef LEPTON_ANALYZER_H
#define LEPTON_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"
#include "TMVA/Reader.h"

//c++ library classes
#include<memory>                                                                                    //for using std::shared_ptr


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
    double _lEta[nL_max];
    double _lEtaSC[nL_max];
    double _lPhi[nL_max];
    double _lE[nL_max];
    unsigned _lFlavor[nL_max];                                                                       //other lepton variables
    int _lCharge[nL_max];
    double _relIso[nL_max];
    double _miniIso[nL_max];
    double _ptRel[nL_max];                                                                          //variables related to closest Jet
    double _ptRatio[nL_max];
    double _closestJetCsv[nL_max];
    double _selectedTrackMult[nL_max];
    double _dxy[nL_max];
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _3dIPSig[nL_max];
    float _lElectronMva[nL_max];
    bool _lElectronPassEmu[nL_max];
    bool _lHNLoose[nL_max];                                                                          //lepton selection decisions
    bool _lHNFO[nL_max];
    bool _lHNTight[nL_max];
    bool _lPOGVeto[nL_max];
    bool _lPOGLoose[nL_max];
    bool _lPOGMedium[nL_max];
    bool _lPOGTight[nL_max];
    bool _lIsPrompt[nL_max];
    int _lMatchPdgId[nL_max];
    double _leptonMva[nL_max];


    multilep* multilepAnalyzer;

    void fillLeptonGenVars(const reco::GenParticle*);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Muon&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&);  
    bool eleMuOverlap(const pat::Electron& ele);
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&);

    // In leptonAnalyzerIso,cc
    double getRelIso03(const pat::Muon&, const double);
    double getRelIso03(const pat::Electron&, const double);
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection> pfcands, double, double, double, double, bool onlyCharged = false);

    // In LeptonAnalyzerId.cc
    float dEtaInSeed(const pat::Electron*);
    bool  isLooseCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool  isTightCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool  passTriggerEmulationDoubleEG(const pat::Electron*);
    float slidingCut(float, float, float);
    bool  passingElectronMvaHZZ(const pat::Electron*, double);
    bool  passingElectronMvaLooseSusy(const pat::Electron*, double, double);
    bool  passingElectronMvaMediumSusy(const pat::Electron*, double);
    bool  passingElectronMvaTightSusy(const pat::Electron*, double);
    bool  passingElectronMvaHeavyNeutrinoFO(const pat::Electron*, double);
  
    bool  isHNLoose(const pat::Electron& lepton);
    bool  isHNLoose(const pat::Muon& lepton);
    bool  isHNFO(const pat::Electron& lepton);
    bool  isHNFO(const pat::Muon& lepton);
    bool  isHNTight(const pat::Electron& lepton);
    bool  isHNTight(const pat::Muon& lepton);

    void   setCommonMvaVars(const reco::Candidate&);
    double leptonMvaVal(const pat::Muon&);
    double leptonMvaVal(const pat::Electron&);

    //Declare lepton MVA
    std::shared_ptr<TMVA::Reader> readerEle;
    std::shared_ptr<TMVA::Reader> readerMu;
    //Variables for booking electroweakino lepton MVA
    float LepGood_pt, 
    LepGood_eta, 
    LepGood_jetNDauChargedMVASel,
    LepGood_miniRelIsoCharged, 
    LepGood_miniRelIsoNeutral,
    LepGood_jetPtRelv2,
    LepGood_jetPtRatio,
    LepGood_jetBTagCSV,
    LepGood_sip3d,
    LepGood_dxy, 
    LepGood_dz,
    LepGood_segmentCompatibility,
    LepGood_mvaIdSpring16GP;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const reco::Vertex&);
};

#endif
