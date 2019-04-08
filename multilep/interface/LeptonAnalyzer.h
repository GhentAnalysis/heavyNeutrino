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

//include ROOT classes
#include "TTree.h"

//include c++ library classes
#include <memory>                                                                                   //for std::shared_ptr


/*
 * Functions for electron identification
 */
class multilep;

class LeptonAnalyzer {
  //Friend classes and functions
  friend class multilep;
  private:
    //this has to come before the effective areas as their initialization depends on it!
    multilep* multilepAnalyzer;

    EffectiveAreas electronsEffectiveAreas;
    EffectiveAreas muonsEffectiveAreas;

    //maximum number of leptons to be stored 
    static const unsigned nL_max = 20;                                                               //maximum number of particles stored

    //number of leptons of each type in the event
    unsigned _nL;
    unsigned _nMu;
    unsigned _nEle;
    unsigned _nLight;
    unsigned _nTau;

    //lepton kinematics and systematic variations
    double _lPt[nL_max];
    double _lPtCorr[nL_max];
    double _lPtScaleUp[nL_max];
    double _lPtScaleDown[nL_max];
    double _lPtResUp[nL_max];
    double _lPtResDown[nL_max];
    double _lEta[nL_max];
    double _lPhi[nL_max];
    double _lE[nL_max];
    double _lECorr[nL_max];
    double _lEScaleUp[nL_max]; //probably useless? I guess ECorr/E = pTCorr/pT (same for up- down variations)
    double _lEScaleDown[nL_max]; //probably useless? I guess ECorr/E = pTCorr/pT (same for up- down variations)
    double _lEResUp[nL_max]; //probably useless? I guess ECorr/E = pTCorr/pT (same for up- down variations)
    double _lEResDown[nL_max]; //probably useless? I guess ECorr/E = pTCorr/pT (same for up- down variations)

    //lepton flavor and charge 
    unsigned _lFlavor[nL_max];
    int _lCharge[nL_max];

    //lepton isolation
    double _relIso[nL_max];
    double _relIso0p4[nL_max];
    double _relIso0p4MuDeltaBeta[nL_max];
    double _miniIso[nL_max];
    double _miniIsoCharged[nL_max];

    //variables based on closest jet to lepton (typically containing lepton)
    double _ptRel[nL_max];
    double _ptRatio[nL_max];
    double _closestJetCsvV2[nL_max];
    double _closestJetDeepCsv_b[nL_max];
    double _closestJetDeepCsv_bb[nL_max];
    unsigned _selectedTrackMult[nL_max];

    //pointing variables
    double _dxy[nL_max];
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _3dIPSig[nL_max];

    //electron properties 
    float _lElectronMvaSummer16GP[nL_max];                                                           // OLD
    float _lElectronMvaSummer16HZZ[nL_max];                                                          // OLD
    float _lElectronMvaFall17v1NoIso[nL_max];                                                        // OLD
    float _lElectronMvaFall17Iso[nL_max];
    float _lElectronMvaFall17NoIso[nL_max];
    bool _lElectronPassEmu[nL_max];
    bool _lElectronPassConvVeto[nL_max];
    bool _lElectronChargeConst[nL_max];
    unsigned _lElectronMissingHits[nL_max];
    double _lEtaSC[nL_max];

    //muon properties
    double _lMuonSegComp[nL_max];
    double _lMuonTrackPt[nL_max];
    double _lMuonTrackPtErr[nL_max];

    bool _tauMuonVetoLoose[nL_max];                                                                       //tau specific variables
    bool _tauMuonVetoTight[nL_max];                                                                       
    bool _tauEleVetoVLoose[nL_max];
    bool _tauEleVetoLoose[nL_max];
    bool _tauEleVetoMedium[nL_max];
    bool _tauEleVetoTight[nL_max];
    bool _tauEleVetoVTight[nL_max];
    bool _decayModeFinding[nL_max];                      
    unsigned int _tauDecayMode[nL_max];                                                         // As in https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Decay_Mode_Reconstruction 
    unsigned int _tauGenStatus[nL_max];                                                         //1: prompt ele, 2:prompt mu, 3: ele from leptonic tau, 4:mu from leptonic tau, 5: hadronically decayed tau, 6:rest 
    bool _tauPOGVLoose2015[nL_max];                                                             //version of ID to use in MiniAOD: MC 80X_mcRun2_asymptotic_2016_TrancheIV_v6, Data 03Feb2017
    bool _tauPOGLoose2015[nL_max];                                                              //More info at https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Isolation
    bool _tauPOGMedium2015[nL_max];
    bool _tauPOGTight2015[nL_max];
    bool _tauPOGVTight2015[nL_max];
    
    bool _tauPOGVVLoose2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGVTight2017v2[nL_max];                                                            //Other WPs contained in _lPOG variables (vloose = veto)
    bool _tauPOGVVTight2017v2[nL_max];

    bool _decayModeFindingNew[nL_max];                      
    bool _tauVLooseMvaNew[nL_max];                                                               
    bool _tauVLooseMvaNew2015[nL_max];
    bool _tauLooseMvaNew2015[nL_max];
    bool _tauMediumMvaNew2015[nL_max];
    bool _tauTightMvaNew2015[nL_max];
    bool _tauVTightMvaNew2015[nL_max];
    
    bool _tauVLooseMvaNew2017v2[nL_max];
    bool _tauLooseMvaNew2017v2[nL_max];
    bool _tauMediumMvaNew2017v2[nL_max];
    bool _tauTightMvaNew2017v2[nL_max];
    bool _tauVTightMvaNew2017v2[nL_max];
    
    double _tauAgainstElectronMVA6Raw[nL_max];
    double _tauCombinedIsoDBRaw3Hits[nL_max];
    double _tauIsoMVAPWdR03oldDMwLT[nL_max];
    double _tauIsoMVADBdR03oldDMwLT[nL_max];
    double _tauIsoMVADBdR03newDMwLT[nL_max];
    double _tauIsoMVAPWnewDMwLT[nL_max];
    double _tauIsoMVAPWoldDMwLT[nL_max];

    //lepton MVA definitions for SUSY (ewkino), TTH and tZq 
    double _leptonMvaSUSY[nL_max];
    double _leptonMvaTTH[nL_max];
    double _leptonMvatZq[nL_max];

    //analysis specific lepton selection (please do NOT delete for now, will be deleted when analysis code is updated 
    bool _lEwkLoose[nL_max];
    bool _lEwkFO[nL_max];
    bool _lEwkTight[nL_max];

    //official POG selection definitions
    bool _lPOGVeto[nL_max];
    bool _lPOGLoose[nL_max];
    bool _lPOGMedium[nL_max];
    bool _lPOGTight[nL_max];

    //MC truth information from matching 
    bool _lIsPrompt[nL_max];
    int _lMatchPdgId[nL_max];
    int _lMomPdgId[nL_max];
    unsigned _lProvenance[nL_max];
    unsigned _lProvenanceCompressed[nL_max];
    unsigned _lProvenanceConversion[nL_max];

    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles);
    void fillTauGenVars(const pat::Tau&, const std::vector<reco::GenParticle>& genParticles);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Muon&, const reco::Vertex&);
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&) const;
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose) const;
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose) const;
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&, const double rho);

    // In leptonAnalyzerIso.cc
    double getRelIso03(const pat::Muon&, const double) const;
    double getRelIso03(const pat::Electron&, const double) const;
    double getRelIso04(const pat::Muon& mu, const double, const bool DeltaBeta=false) const;
    double getRelIso(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, const bool onlyCharged=false) const;
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, double, double, bool onlyCharged=false) const;

    // In LeptonAnalyzerId.cc
    bool  passTriggerEmulationDoubleEG(const pat::Electron*, const bool hOverE = true) const;               //For ewkino id it needs to be possible to check hOverE separately
    float slidingCut(float, float, float) const;
    bool  passingElectronMvaHZZ(const pat::Electron*, double) const;
    bool  passingElectronMvaLooseSusy(const pat::Electron*, double, double) const;
    bool  passingElectronMvaMediumSusy(const pat::Electron*, double) const;
    bool  passingElectronMvaTightSusy(const pat::Electron*, double) const;
    bool  passingElectronMvaHeavyNeutrinoFO(const pat::Electron*, double) const;
    bool  passElectronMvaEwkFO(const pat::Electron* ele, double mvaValue) const;

    //analysis specific lepton selection (please do NOT delete for now, will be deleted when analysis code is updated 
    bool isEwkLoose(const pat::Muon&) const;
    bool isEwkLoose(const pat::Electron&) const;
    bool isEwkLoose(const pat::Tau&) const;
    bool isEwkFO(const pat::Muon&) const;
    bool isEwkFO(const pat::Electron&) const;
    bool isEwkFO(const pat::Tau&) const;
    bool isEwkTight(const pat::Muon&) const;
    bool isEwkTight(const pat::Electron&) const;
    bool isEwkTight(const pat::Tau&) const;

    //compute lepton MVA value 
    double leptonMvaVal(const pat::Muon&, LeptonMvaHelper*);
    double leptonMvaVal(const pat::Electron&, LeptonMvaHelper*);

    //for lepton MVA calculation
    LeptonMvaHelper* leptonMvaComputerSUSY;
    LeptonMvaHelper* leptonMvaComputerTTH;
    LeptonMvaHelper* leptonMvaComputertZq;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const reco::Vertex&);
};
#endif
