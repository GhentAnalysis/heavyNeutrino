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

#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

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
  private:
    EffectiveAreas electronsEffectiveAreas;
    EffectiveAreas muonsEffectiveAreas;

    static const unsigned nL_max      = 20;                                                          //maximum number of particles stored
    static const unsigned nV_max      = 100;   
  
    unsigned _nL;                                                                                    //number of leptons
    unsigned _nMu;
    unsigned _nEle;
    unsigned _nLight;
    unsigned _nTau;          
    unsigned _nVFit;                     // number vertices re-fitted
    unsigned _nGoodLeading;                     // number vertices re-fitted

  
    double _lIndex[nL_max];              // index assigned to leptons to find back the vertices
    double _vertices[nV_max][12];        // array of the vertices: 9 variables+index for each vertex 

    double _lPt[nL_max];                                                                             //lepton kinematics
    double _lEta[nL_max];
    double _lEtaSC[nL_max];
    double _lPhi[nL_max];
    double _lE[nL_max];

    unsigned _lFlavor[nL_max];                                                                       //lepton flavor and charge
    int _lCharge[nL_max];

    double _relIso[nL_max];                                                                          //lepton isolation variables
    double _relIso0p4Mu[nL_max];                                                                          //lepton isolation variables
    double _miniIso[nL_max];
    double _miniIsoCharged[nL_max];    
    double _puCorr[nL_max];     
    double _absIso03[nL_max];   
    double _absIso04[nL_max];                    
    double _sumNeutralHadronEt04[nL_max];          
    double _sumChargedHadronPt04[nL_max];         
    double _sumPhotonEt04[nL_max];                 
    double _sumNeutralHadronEt03[nL_max];         
    double _sumChargedHadronPt03[nL_max];         
    double _sumPhotonEt03[nL_max];  
    double _trackIso[nL_max];
    double _ecalIso[nL_max];
    double _hcalIso[nL_max];
    double _deltaBIso[nL_max];
    double _ecalPFClusterIso[nL_max];
    double _hcalPFClusterIso[nL_max];
   
    double _ptRel[nL_max];                                                                           //variables related to closest jet
    double _ptRatio[nL_max];
    double _closestJetCsvV2[nL_max];
    double _closestJetDeepCsv_b[nL_max];
    double _closestJetDeepCsv_bb[nL_max];
    unsigned _selectedTrackMult[nL_max];

    double _dxy[nL_max];                                                                             //pointing variables
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _2dIP[nL_max];
    double _3dIPSig[nL_max];
    double _2dIPSig[nL_max];
    float _lElectronMva[nL_max];
    float _lElectronMvaHZZ[nL_max];

    bool _lElectronPassEmu[nL_max];
    bool _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[nL_max];
    bool _lPOGVeto[nL_max];
    bool _lPOGLoose[nL_max];
    bool _lPOGMedium[nL_max];
    bool _lPOGTight[nL_max];
    bool _lPFMuon[nL_max];
    bool _lpassConversionVeto[nL_max];
   
  
  
    bool _lPFparticle[nL_max];
    bool _lGlobalMuon[nL_max];
    bool _lTrackerMuon[nL_max];
    double _lInnerTrackValidFraction[nL_max];
    double _lGlobalTrackNormalizeChi2[nL_max];
    double _lCQChi2Position[nL_max];
    double _lCQTrackKink[nL_max];
    double _muonSegComp[nL_max];
    double _lNumberOfMatchedStation[nL_max];
    double _lNumberOfValidPixelHits[nL_max];
    double _muNumberInnerHits[nL_max];
    double _lTrackerLayersWithMeasurement[nL_max];
  
    bool _tauMuonVeto[nL_max];                                                                       //tau specific variables
    bool _tauEleVeto[nL_max];
    bool _decayModeFindingNew[nL_max];                      
    bool _tauVLooseMvaNew[nL_max];                                                                      //"old tau id's will be stored in the POG id definitions (vloose := veto), however very tight is stored separately
    bool _tauLooseMvaNew[nL_max];
    bool _tauMediumMvaNew[nL_max];
    bool _tauTightMvaNew[nL_max];
    bool _tauVTightMvaNew[nL_max];
    bool _tauVTightMvaOld[nL_max];
  
  
  
    double _eleNumberInnerHitsMissing[nL_max];
    
  

    bool _lPOGLooseWOIso[nL_max];
    bool _lPOGMediumWOIso[nL_max];
    bool _lPOGTightWOIso[nL_max];

    bool _lIsPrompt[nL_max];                                                                          //MC-truth variables
    int _lMatchPdgId[nL_max];
  
    unsigned _lProvenance[nL_max];                                                                    



    multilep* multilepAnalyzer;

    TransientVertex dileptonVertex(const reco::Track& tk1, const reco::Track& tk2);
    void fillLeptonGenVars(const reco::Candidate&, GenMatching*);

    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Muon&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    void fillLeptonIsoVars(const pat::Muon& mu, const double rho);
    void fillLeptonIsoVars(const pat::Electron& ele, const double rho);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&);  
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose);
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose);
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&);

    // In leptonAnalyzerIso,cc

    double getRelIso03(const pat::Muon&, const double);
    double getRelIso04(const pat::Muon&, const double);
    double getRelIso03(const pat::Electron&, const double);
    double getRelIso04(const pat::Muon& mu);
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection> pfcands, double, double, double, double, bool onlyCharged = false);

  
  
  
    // In LeptonAnalyzerId.cc
    float dEtaInSeed(const pat::Electron*);
    bool isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(const pat::Electron* ele);
    bool  isLooseCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool  isMediumCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool  isTightCutBasedElectronWithoutIsolation(const pat::Electron*);
    bool  passTriggerEmulationDoubleEG(const pat::Electron*, const bool hOverE = true);               //For ewkino id it needs to be possible to check hOverE separately
    float slidingCut(float, float, float);
    bool  passingElectronMvaHZZ(const pat::Electron*, double);
    bool  passingElectronMvaLooseSusy(const pat::Electron*, double, double);
    bool  passingElectronMvaMediumSusy(const pat::Electron*, double);
    bool  passingElectronMvaTightSusy(const pat::Electron*, double);
    bool  passingElectronMvaHeavyNeutrinoFO(const pat::Electron*, double);
    bool  passElectronMvaEwkFO(const pat::Electron* ele, double mvaValue);
  
    bool  isHNLoose(const pat::Electron& lepton);                                                     //check HNL id definitions
    bool  isHNLoose(const pat::Muon& lepton);
    bool  isHNFO(const pat::Electron& lepton);
    bool  isHNFO(const pat::Muon& lepton);
    bool  isHNTight(const pat::Electron& lepton);
    bool  isHNTight(const pat::Muon& lepton);
    
    bool isEwkLoose(const pat::Muon&);
    bool isEwkLoose(const pat::Electron&);
    bool isEwkLoose(const pat::Tau&);
    bool isEwkFO(const pat::Muon&);
    bool isEwkFO(const pat::Electron&);
    bool isEwkFO(const pat::Tau&);
    bool isEwkTight(const pat::Muon&);
    bool isEwkTight(const pat::Electron&);
    bool isEwkTight(const pat::Tau&);

    double leptonMvaVal(const pat::Muon&, LeptonMvaHelper*);                                                            //compute ewkino lepton MVA
    double leptonMvaVal(const pat::Electron&, LeptonMvaHelper*);
    
    //for lepton MVA calculation
    LeptonMvaHelper* leptonMvaComputerSUSY;
    LeptonMvaHelper* leptonMvaComputerTTH;

    //for generator matching
    GenMatching* genMatcher;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const reco::Vertex&);
};
#endif
