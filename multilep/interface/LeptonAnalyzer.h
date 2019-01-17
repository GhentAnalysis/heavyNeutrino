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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

//include other parts of the framework
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"
#include "heavyNeutrino/multilep/interface/GenMatching.h"  // displaced specific

//include classes for trigger match
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//include ROOT classes
#include "TTree.h"

//include c++ library classes
#include <memory>                                                                                   //for std::shared_ptr


/*
 * Functions for electron identification
 */
class multilep;
class GenMatching; // displaced specific

class LeptonAnalyzer {
  //Friend classes and functions
  friend class multilep;
  friend class GenMatching;

  private:
    //this has to come before the effective areas as their initialization depends on it!
    multilep* multilepAnalyzer;

    EffectiveAreas electronsEffectiveAreas;
    EffectiveAreas muonsEffectiveAreas;

    static const unsigned nL_max = 20;                                                               //maximum number of particles stored
    static const unsigned nV_max = 50;
    unsigned _nL;                                                                                    //number of leptons
    unsigned _nMu;
    unsigned _nEle;
    unsigned _nLight;
    unsigned _nTau;

    double _pvX;
    double _pvY;
    double _pvZ;
    double _pvXErr;
    double _pvYErr;
    double _pvZErr;

    unsigned _nVFit;                     // number vertices re-fitted
    unsigned _nGoodLeading;              // number vertices re-fitted
    unsigned _nGoodDisplaced;

    int    _lSimType[nL_max];
    int    _lSimExtType[nL_max];
    int    _lSimFlavour[nL_max];

    unsigned _lIndex[nL_max];            // index assigned to leptons to find back the vertices
    double _vertices[nV_max][12];        // array of the vertices: 9 variables+index for each vertex
    double _lDisplaced[nV_max][24];      // array of the displaced lepton momenta and positions at the displaced vertex

    unsigned _lHasTrigger[nL_max];                                                                   //trigger matching

    double _lPt[nL_max];                                                                             //lepton kinematics
    double _lPtCorr[nL_max];
    double _lPtScaleUp[nL_max];
    double _lPtScaleDown[nL_max];
    double _lPtResUp[nL_max];
    double _lPtResDown[nL_max];
    double _lEta[nL_max];
    double _lEtaSC[nL_max];
    double _lPhi[nL_max];
    double _lE[nL_max];
    double _lECorr[nL_max];
    double _lEScaleUp[nL_max];
    double _lEScaleDown[nL_max];
    double _lEResUp[nL_max];
    double _lEResDown[nL_max];

    unsigned _lFlavor[nL_max];                                                                       //lepton flavor and charge
    int _lCharge[nL_max];

    double _relIso[nL_max];                                                                          //lepton isolation variables
    double _relIso0p4[nL_max];                                                                       //lepton isolation variables
    double _relIsoOld[nL_max];                                                                       //lepton isolation variables, not stored
    double _relIso0p4Old[nL_max];                                                                    //lepton isolation variables, not stored
    double _relIso0p4MuDeltaBeta[nL_max];                                                            //lepton isolation variables
    double _miniIso[nL_max];
    double _miniIsoCharged[nL_max];
    double _puCorr[nL_max];
    double _absIso03[nL_max]
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

    float _lElectronMvaSummer16GP[nL_max];                                                           // OLD
    float _lElectronMvaSummer16HZZ[nL_max];                                                          // OLD
    float _lElectronMvaFall17v1NoIso[nL_max];                                                        // OLD
    float _lElectronMvaFall17Iso[nL_max];
    float _lElectronMvaFall17NoIso[nL_max];
    bool _lElectronPassEmu[nL_max];
    bool _lElectronPassConvVeto[nL_max];
    bool _lElectronChargeConst[nL_max];
    unsigned _lElectronMissingHits[nL_max];

    double _lMuonSegComp[nL_max];                                                                     //muon speficic variables
    double _lMuonTrackPt[nL_max];
    double _lMuonTrackPtErr[nL_max];

    double _2dIPSig[nL_max];
    float _lElectronMva[nL_max];
    float _lElectronMvaHZZ[nL_max];

    bool _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[nL_max];
    bool _lPFMuon[nL_max];
    bool _lpassConversionVeto[nL_max];


    /////// µ ID variables
    bool _lGlobalMuon[nL_max];
    bool _lTrackerMuon[nL_max];
    double _lInnerTrackValidFraction[nL_max];
    double _lGlobalTrackNormalizeChi2[nL_max];
    double _lCQChi2Position[nL_max];
    double _lCQTrackKink[nL_max];
    double _muonSegComp[nL_max];
    unsigned _lNumberOfMatchedStation[nL_max];
    unsigned _lNumberOfValidPixelHits[nL_max];
    unsigned _muNumberInnerHits[nL_max];
    unsigned _lTrackerLayersWithMeasurement[nL_max];

    int _muDTStationsWithValidHits[nL_max];
    int _muCSCStationsWithValidHits[nL_max];
    int _muRPCStationsWithValidHits[nL_max];
    int _muMuonStationsWithValidHits[nL_max];

    double _lMuTime[nL_max];
    double _lMuTimeErr[nL_max];
    double _lMuRPCTime[nL_max];
    double _lMuRPCTimeErr[nL_max];
    int    _lMuTimenDof[nL_max];
    int    _lMuRPCTimenDof[nL_max];

    /////// ele ID variabels
    bool _lEleIsEB [nL_max];
    bool _lEleIsEE[nL_max];
    double _lEleSuperClusterOverP[nL_max];
    double _lEleEcalEnergy[nL_max];
    double _lElefull5x5SigmaIetaIeta[nL_max];
    double _lEleDEtaInSeed[nL_max];
    double _lEleDeltaPhiSuperClusterTrackAtVtx[nL_max];
    double _lElehadronicOverEm[nL_max];
    double _lEleInvMinusPInv[nL_max];
    double _eleNumberInnerHitsMissing[nL_max];

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

    double _leptonMvaSUSY16[nL_max];                                                                       //lepton MVA used in ewkino analysis
    double _leptonMvaTTH16[nL_max];
    double _leptonMvaSUSY17[nL_max];                                                                       //lepton MVA used in ewkino analysis
    double _leptonMvaTTH17[nL_max];
    double _leptonMvatZqTTV16[nL_max];
    double _leptonMvatZqTTV17[nL_max];

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

    // Index in the gen-particle list
    unsigned _lGenIndex[nL_max];
    unsigned _lMatchType[nL_max];

    bool _lIsPrompt[nL_max];                                                                          //MC-truth variables
    bool _lIsPromptFinalState[nL_max];                                                                          //MC-truth variables
    bool _lIsPromptDecayed[nL_max];                                                                          //MC-truth variables

    int _lMatchPdgId[nL_max];
    int _lMomPdgId[nL_max];
    unsigned _lProvenance[nL_max];
    unsigned _lProvenanceCompressed[nL_max];
    unsigned _lProvenanceConversion[nL_max];
    double _lMatchPt[nL_max];
    double _lMatchEta[nL_max];
    double _lMatchPhi[nL_max];
    double _lMatchVertexX[nL_max];
    double _lMatchVertexY[nL_max];
    double _lMatchVertexZ[nL_max];

    std::vector<std::string> singleEleTrigs, singleMuoTrigs;


    edm::ESHandle<MagneticField> _bField;
    edm::ESHandle<Propagator> _shProp;
    mutable TransverseImpactPointExtrapolator* _gsfProp;
    TransientVertex dileptonVertex(const reco::Track&, const reco::Track&);
    void cleanDileptonVertexArrays(unsigned);

    void fillDileptonVertexArrays(unsigned, unsigned, unsigned,
                  const TransientVertex&,
                  const reco::Track&, const reco::Track&,
                  bool, bool);

    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles);
    template <typename Lepton> void fillLeptonGenVars(GenMatching* genMatcher, const Lepton& lepton, const reco::GenParticle* match = nullptr, unsigned mtchtype = 6);
    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, GenMatching* genMatcher);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Muon&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    void fillLeptonIsoVars(const pat::Muon& mu, const double rho);
    void fillLeptonIsoVars(const pat::Electron& ele, const double rho);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&) const;
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose) const;
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose) const;
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&, const double rho);
    unsigned matchSingleTrigger(bool, double, double, const edm::TriggerNames&, edm::Handle<pat::TriggerObjectStandAloneCollection>);

    // To synchronize lepton selection
    bool passElectronPreselection(const pat::Electron&, const double rho) const;
    bool passMuonPreselection(const pat::Muon&, const double rho) const;
    bool passTauPreselection(const pat::Tau&, const reco::Vertex::Point&) const;

    // In leptonAnalyzerIso.cc
    double getRelIso03(const pat::Muon&, const double) const;
    double getRelIso03(const pat::Electron&, const double) const;
    double getRelIso04(const pat::Muon& mu, const double, const bool DeltaBeta=false) const;
    double getRelIso(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, const bool onlyCharged=false) const;
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, double, double, bool onlyCharged=false) const;

    // In LeptonAnalyzerId.cc
    float dEtaInSeed(const pat::Electron*) const;                                                                          // displaced specific
    bool  isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(const pat::Electron*) const; // displaced specific
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
    LeptonMvaHelper* leptonMvaComputerSUSY16;
    LeptonMvaHelper* leptonMvaComputerTTH16;
    LeptonMvaHelper* leptonMvaComputerSUSY17;
    LeptonMvaHelper* leptonMvaComputerTTH17;
    LeptonMvaHelper* leptonMvaComputertZqTTV16;
    LeptonMvaHelper* leptonMvaComputertZqTTV17;

    //for generator matching (displaced specific)
    GenMatching* genMatcher;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const edm::EventSetup&, const reco::Vertex&);
};
#endif
