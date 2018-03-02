#ifndef LEPTON_ANALYZER_H
#define LEPTON_ANALYZER_H
//include other parts of CMSSW
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//include other parts of the framework
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"
#include "heavyNeutrino/multilep/interface/LeptonIdHelper.h"
#include "heavyNeutrino/multilep/interface/GenMatching.h"

//include ROOT classes
#include "TTree.h"

//include c++ library classes
#include <memory>                                                                                   //for std::shared_ptr


/*
 * Functions for electron identification
 */
class multilep;
class GenMatching;

class LeptonAnalyzer {
  //Friend classes and functions
  friend class LeptonIdHelper;
  friend class multilep;
  private:
    //this has to come before the effective areas as their initialization depends on it!
    multilep* multilepAnalyzer;

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

    unsigned _lFlavor[nL_max];                                                                       //lepton flavor and charge
    int _lCharge[nL_max];

    double _relIso[nL_max];                                                                          //lepton isolation variables
    double _relIso0p4[nL_max];                                                                          //lepton isolation variables
    double _relIso0p6[nL_max];                                                                          //lepton isolation variables
    double _relIso0p8[nL_max];                                                                          //lepton isolation variables
    double _relIso1p0[nL_max];                                                                          //lepton isolation variables
    double _relIso0p4Mu[nL_max];                                                                          //lepton isolation variables
    double _miniIso[nL_max];
    double _miniIsoCharged[nL_max];                                                              
    
    double _ptRel[nL_max];                                                                           //variables related to closest jet
    double _ptRatio[nL_max];
    double _closestJetCsvV2[nL_max];
    double _closestJetDeepCsv_b[nL_max];
    double _closestJetDeepCsv_bb[nL_max];
    unsigned _selectedTrackMult[nL_max];

    double _dxy[nL_max];                                                                             //pointing variables
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _3dIPSig[nL_max];

    float _lElectronMva[nL_max];                                                                     //electron specific variables
    float _lElectronMvaHZZ[nL_max];
    float _lElectronMvaFall17Iso[nL_max];
    float _lElectronMvaFall17NoIso[nL_max];
    bool _lElectronPassEmu[nL_max];                                                                  
    bool _lElectronPassConvVeto[nL_max];
    bool _lElectronChargeConst[nL_max];
    unsigned _lElectronMissingHits[nL_max];

    double _lMuonSegComp[nL_max];                                                                     //muon speficic variables
    double _lMuonTrackPt[nL_max];
    double _lMuonTrackPtErr[nL_max];

    bool _tauMuonVeto[nL_max];                                                                       //tau specific variables
    bool _tauEleVeto[nL_max];
    bool _decayModeFindingNew[nL_max];                      
    bool _tauVLooseMvaNew[nL_max];                                                                      //"old tau id's will be stored in the POG id definitions (vloose := veto), however very tight is stored separately
    bool _tauLooseMvaNew[nL_max];
    bool _tauMediumMvaNew[nL_max];
    bool _tauTightMvaNew[nL_max];
    bool _tauVTightMvaNew[nL_max];
    bool _tauVTightMvaOld[nL_max];

    double _tauAgainstElectronMVA6Raw[nL_max];
    double _tauCombinedIsoDBRaw3Hits[nL_max];
    double _tauIsoMVAPWdR03oldDMwLT[nL_max];
    double _tauIsoMVADBdR03oldDMwLT[nL_max];
    double _tauIsoMVADBdR03newDMwLT[nL_max];
    double _tauIsoMVAPWnewDMwLT[nL_max];
    double _tauIsoMVAPWoldDMwLT[nL_max];

    double _leptonMvaSUSY[nL_max];                                                                       //lepton MVA used in ewkino analysis
    double _leptonMvaTTH[nL_max];
    double _leptonMvatZqTTV[nL_max];

    bool _lHNLoose[nL_max];                                                                          //analysis specific lepton selection decisions
    bool _lHNFO[nL_max];
    bool _lHNTight[nL_max];
    bool _lEwkLoose[nL_max];
    bool _lEwkFO[nL_max];
    bool _lEwkTight[nL_max];
    bool _lPOGVeto[nL_max];
    bool _lPOGLoose[nL_max];
    bool _lPOGMedium[nL_max];
    bool _lPOGTight[nL_max];

    bool _lPOGLooseWOIso[nL_max];
    bool _lPOGMediumWOIso[nL_max];
    bool _lPOGTightWOIso[nL_max];

    bool _lIsPrompt[nL_max];                                                                          //MC-truth variables
    int _lMatchPdgId[nL_max];
    unsigned _lProvenance[nL_max];                                                                    
    unsigned _lProvenanceCompressed[nL_max];
    double _lPartonPt[nL_max];
    unsigned _lProvenanceConversion[nL_max];

    //void fillLeptonGenVars(const reco::Candidate&, GenMatching*);
    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, GenMatching* genMatcher);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Muon&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&) const;
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose) const;
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose) const;
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&);

    // In leptonAnalyzerIso,cc

    double getRelIso03(const pat::Muon&, const double) const;
    double getRelIso03(const pat::Electron&, const double) const;
    double getRelIso04(const pat::Muon& mu) const;
    double getRelIso(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, const bool onlyCharged = false) const;
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, double, double, bool onlyCharged = false) const;

    // In LeptonAnalyzerId.cc
    float dEtaInSeed(const pat::Electron*) const;
    bool  isLooseCutBasedElectronWithoutIsolation(const pat::Electron*) const;
    bool  isMediumCutBasedElectronWithoutIsolation(const pat::Electron*) const;
    bool  isTightCutBasedElectronWithoutIsolation(const pat::Electron*) const;
    bool  passTriggerEmulationDoubleEG(const pat::Electron*, const bool hOverE = true) const;               //For ewkino id it needs to be possible to check hOverE separately
    float slidingCut(float, float, float) const;
    bool  passingElectronMvaHZZ(const pat::Electron*, double) const;
    bool  passingElectronMvaLooseSusy(const pat::Electron*, double, double) const;
    bool  passingElectronMvaMediumSusy(const pat::Electron*, double) const;
    bool  passingElectronMvaTightSusy(const pat::Electron*, double) const;
    bool  passingElectronMvaHeavyNeutrinoFO(const pat::Electron*, double) const;
    bool  passElectronMvaEwkFO(const pat::Electron* ele, double mvaValue) const;
  
    bool  isHNLoose(const pat::Electron& lepton) const;                                                     //check HNL id definitions
    bool  isHNLoose(const pat::Muon& lepton) const;
    bool  isHNFO(const pat::Electron& lepton) const;
    bool  isHNFO(const pat::Muon& lepton) const;
    bool  isHNTight(const pat::Electron& lepton) const;
    bool  isHNTight(const pat::Muon& lepton) const;
    
    bool isEwkLoose(const pat::Muon&) const;
    bool isEwkLoose(const pat::Electron&) const;
    bool isEwkLoose(const pat::Tau&) const;
    bool isEwkFO(const pat::Muon&) const;
    bool isEwkFO(const pat::Electron&) const;
    bool isEwkFO(const pat::Tau&) const;
    bool isEwkTight(const pat::Muon&) const;
    bool isEwkTight(const pat::Electron&) const;
    bool isEwkTight(const pat::Tau&) const;

    double leptonMvaVal(const pat::Muon&, LeptonMvaHelper*);                                                            //compute ewkino lepton MVA
    double leptonMvaVal(const pat::Electron&, LeptonMvaHelper*);
    
    //for lepton MVA calculation
    LeptonMvaHelper* leptonMvaComputerSUSY;
    LeptonMvaHelper* leptonMvaComputerTTH;
    LeptonMvaHelper* leptonMvaComputertZqTTV;

    //for generator matching
    GenMatching* genMatcher;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const reco::Vertex&);
};
#endif
