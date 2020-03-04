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
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


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
    EffectiveAreas electronsEffectiveAreasMiniIso;
    EffectiveAreas muonsEffectiveAreas;

    //maximum number of leptons to be stored 
    static const unsigned nL_max = 20;                                                               //maximum number of particles stored
    static const unsigned nV_max = 50;

    //number of leptons of each type in the event
    unsigned _nL = 0;                                                                                //number of leptons
    unsigned _nMu = 0;
    unsigned _nEle = 0;
    unsigned _nLight = 0;
    unsigned _nTau = 0;

    double _rho;

    double _pvX;
    double _pvY;
    double _pvZ;
    double _pvXErr;
    double _pvYErr;
    double _pvZErr;

    unsigned _nVFit;                     // number vertices re-fitted
    unsigned _nVFit_os;
    unsigned _nGoodLeading;              // number vertices re-fitted
    unsigned _nGoodDisplaced;

    int    _lSimType[nL_max];
    int    _lSimExtType[nL_max];
    int    _lSimFlavour[nL_max];

    double _vertices[nV_max][12];        // array of the vertices: 9 variables+index for each vertex
    double _lDisplaced[nV_max][24];      // array of the displaced lepton momenta and positions at the displaced vertex
    double _vertices_os[nV_max][12];        // array of the vertices: 9 variables+index for each vertex
    double _lDisplaced_os[nV_max][24];
  
    unsigned _lHasTrigger[nL_max];                                                                   //trigger matching

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
    double _lEnergySC[nL_max];
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

    //variables based on closest jet to lepton (typically containing lepton)
    double _ptRel[nL_max];
    double _ptRatio[nL_max];
    double _closestJetCsvV2[nL_max];
    double _closestJetDeepCsv_b[nL_max];
    double _closestJetDeepCsv_bb[nL_max];
    double _closestJetDeepCsv[nL_max];
    double _closestJetDeepFlavor_b[nL_max];
    double _closestJetDeepFlavor_bb[nL_max];
    double _closestJetDeepFlavor_lepb[nL_max];
    double _closestJetDeepFlavor[nL_max];
    unsigned _selectedTrackMult[nL_max];
    double      _closestJEC[nL_max] ;
    double     _closest_lepAwareJetE[nL_max];
    double     _closest_lepAwareJetPx[nL_max];
    double     _closest_lepAwareJetPy[nL_max];
    double     _closest_lepAwareJetPz[nL_max];
    double     _closest_l1JetE[nL_max];
    double     _closest_l1JetPx[nL_max];
    double     _closest_l1JetPy[nL_max];
    double     _closest_l1JetPz[nL_max];
    double     _closest_lJetE [nL_max];
    double     _closest_lJetPx [nL_max] ;
    double     _closest_lJetPy [nL_max] ;
    double     _closest_lJetPz [nL_max]  ;
  
  
  

    //pointing variables
    double _dxy[nL_max];
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _2dIP[nL_max];
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
    bool _lElectronPassMVAFall17NoIsoWP80[nL_max];
    bool _lElectronPassMVAFall17NoIsoWP90[nL_max];
    bool _lElectronPassMVAFall17NoIsoWPLoose[nL_max];
    double _lElectronSigmaIetaIeta[nL_max];
    double _lElectronDeltaPhiSuperClusterTrack[nL_max];
    double _lElectronDeltaEtaSuperClusterTrack[nL_max];
    double _lElectronEInvMinusPInv[nL_max];
    double _lElectronHOverE[nL_max];


    //muon properties
    double _lMuonSegComp[nL_max];
    double _lMuonTrackPt[nL_max];
    double _lMuonTrackPtErr[nL_max];

    double _2dIPSig[nL_max];
    float _lElectronMva[nL_max];
    float _lElectronMvaHZZ[nL_max];

    bool _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[nL_max];
    bool _lPFMuon[nL_max];


    /////// µ ID variables
    bool _lGlobalMuon[nL_max];
    bool _lTrackerMuon[nL_max];
    double _lInnerTrackValidFraction[nL_max];
    double _lGlobalTrackNormalizeChi2[nL_max];
    double _lCQChi2Position[nL_max];
    double _lCQTrackKink[nL_max];
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

    bool _tauMuonVetoLoose[nL_max];                                                                       //tau specific variables
    bool _tauMuonVetoTight[nL_max];                                                                       
    bool _tauEleVetoVLoose[nL_max];
    bool _tauEleVetoLoose[nL_max];
    bool _tauEleVetoMedium[nL_max];
    bool _tauEleVetoTight[nL_max];
    bool _tauEleVetoVTight[nL_max];
    bool _decayModeFinding[nL_max];                      
    unsigned _tauDecayMode[nL_max];                                                             // As in https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Decay_Mode_Reconstruction 
    unsigned _tauGenStatus[nL_max];                                                             //1: prompt ele, 2:prompt mu, 3: ele from leptonic tau, 4:mu from leptonic tau, 5: hadronically decayed tau, 6:rest 
    bool _tauPOGVLoose2015[nL_max];                                                             //version of ID to use in MiniAOD: MC 80X_mcRun2_asymptotic_2016_TrancheIV_v6, Data 03Feb2017
    bool _tauPOGLoose2015[nL_max];                                                              //More info at https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Isolation
    bool _tauPOGMedium2015[nL_max];                                                             // # WARNING # NO LONGER SUPPORTED BY TAU POG, kept for testing reasons, will remove this soon
    bool _tauPOGTight2015[nL_max];
    bool _tauPOGVTight2015[nL_max];
    
    bool _tauPOGVVLoose2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGVTight2017v2[nL_max];                                                            //Other WPs contained in _lPOG variables (vloose = veto)
    bool _tauPOGVVTight2017v2[nL_max];

    bool _decayModeFindingNew[nL_max];                                           
    bool _decayModeFindingDeepTau[nL_max];                                           
    bool _tauVLooseMvaNew[nL_max];                                                              // # WARNING # NO LONGER SUPPORTED BY TAU POG, kept for testing reasons, will remove this soon                 
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

    bool _tauDeepTauVsJetsRaw[nL_max];    
    bool _tauVVVLooseDeepTauVsJets[nL_max];    
    bool _tauVVLooseDeepTauVsJets[nL_max];    
    bool _tauVLooseDeepTauVsJets[nL_max];    
    bool _tauLooseDeepTauVsJets[nL_max];    
    bool _tauMediumDeepTauVsJets[nL_max];    
    bool _tauTightDeepTauVsJets[nL_max];    
    bool _tauVTightDeepTauVsJets[nL_max];    
    bool _tauVVTightDeepTauVsJets[nL_max];    
    
    bool _tauDeepTauVsEleRaw[nL_max];    
    bool _tauVVVLooseDeepTauVsEle[nL_max];    
    bool _tauVVLooseDeepTauVsEle[nL_max];    
    bool _tauVLooseDeepTauVsEle[nL_max];    
    bool _tauLooseDeepTauVsEle[nL_max];    
    bool _tauMediumDeepTauVsEle[nL_max];    
    bool _tauTightDeepTauVsEle[nL_max];    
    bool _tauVTightDeepTauVsEle[nL_max];    
    bool _tauVVTightDeepTauVsEle[nL_max];    
    
    bool _tauDeepTauVsMuRaw[nL_max];    
    bool _tauVLooseDeepTauVsMu[nL_max];    
    bool _tauLooseDeepTauVsMu[nL_max];    
    bool _tauMediumDeepTauVsMu[nL_max];    
    bool _tauTightDeepTauVsMu[nL_max];    
    
    double _tauAgainstElectronMVA6Raw[nL_max];
    double _tauCombinedIsoDBRaw3Hits[nL_max];
    double _tauIsoMVAPWdR03oldDMwLT[nL_max];
    double _tauIsoMVADBdR03oldDMwLT[nL_max];
    double _tauIsoMVADBdR03newDMwLT[nL_max];
    double _tauIsoMVAPWnewDMwLT[nL_max];
    double _tauIsoMVAPWoldDMwLT[nL_max];
    
    //lepton MVA definitions for TTH and tZq 
    double _leptonMvaTTH[nL_max];
    double _leptonMvatZq[nL_max];

    //official POG selection definitions
    bool _lPOGVeto[nL_max];
    bool _lPOGLoose[nL_max];
    bool _lPOGMedium[nL_max];
    bool _lPOGTight[nL_max];

    // Index in the gen-particle list
    unsigned _lGenIndex[nL_max];
    unsigned _lMatchType[nL_max];

    //MC truth information from matching 
    bool _lIsPrompt[nL_max];
    bool _lIsPromptFinalState[nL_max];
    bool _lIsPromptDecayed[nL_max];

    int _lMatchPdgId[nL_max];
    int _lMatchCharge[nL_max];
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
    const reco::Track& getTrack(const reco::RecoCandidate* lep);
    TransientVertex dileptonVertex(const reco::RecoCandidate* lep1, const reco::RecoCandidate* lep2);

    void cleanDileptonVertexArrays(unsigned);
    void fillDileptonVertexArrays(unsigned, unsigned, const reco::RecoCandidate*, const reco::RecoCandidate*, const bool ensureOppositeSign=false);

    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles);
    void fillTauGenVars(const pat::Tau&, const std::vector<reco::GenParticle>& genParticles);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron& );
    void fillLeptonImpactParameters(const pat::Muon& );
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    void fillLeptonIsoVars(const pat::Muon& mu, const double rho);
    void fillLeptonIsoVars(const pat::Electron& ele, const double rho);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&) const;
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose) const;
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose) const;
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&, const double rho, const bool oldMatching = false);
    unsigned matchSingleTrigger(const edm::Event& iEvent, bool isele, double aeta, double aphi);

    // To synchronize lepton selection
    bool passElectronPreselection(const pat::Electron&, const double rho) const;
    bool passMuonPreselection(const pat::Muon&, const double rho) const;
    bool passTauPreselection(const pat::Tau&, const reco::Vertex::Point&) const;

    // In leptonAnalyzerIso.cc
    double getRelIso03(const pat::Muon&, const double) const;
    double getRelIso03(const pat::Electron&, const double) const;
    double getRelIso04(const pat::Muon&, const double, const bool DeltaBeta=false) const;
    double getRelIso04( const pat::Electron&, const double ) const;
    double getRelIso(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, const bool onlyCharged=false) const;
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, double, double, bool onlyCharged=false) const;
    template< typename T > double getMiniIsolation( const T&, const double rho, const bool onlyCharged = false ) const;

    // In LeptonAnalyzerId.cc
    float dEtaInSeed(const pat::Electron*) const;                                                                          // displaced specific
    bool  isLooseCutBasedElectronWithoutIsolationWithoutMissingInnerhitsWithoutConversionVeto(const pat::Electron*) const; // displaced specific
    bool  passTriggerEmulationDoubleEG(const pat::Electron*, const bool hOverE = true) const;               //For ewkino id it needs to be possible to check hOverE separately
    float slidingCut(float, float, float) const;
    bool  passingElectronMvaHZZ(const pat::Electron*, double) const;
    bool  passingElectronMvaLooseSusy(const pat::Electron*, double, double) const;
    bool  passingElectronMvaMediumSusy(const pat::Electron*, double) const;
    bool  passingElectronMvaTightSusy(const pat::Electron*, double) const;

    //compute lepton MVA value 
    double leptonMvaVal(const pat::Muon&, LeptonMvaHelper*);
    double leptonMvaVal(const pat::Electron&, LeptonMvaHelper*);

    //for lepton MVA calculation
    LeptonMvaHelper* leptonMvaComputerTTH;
    LeptonMvaHelper* leptonMvaComputertZq;

    //for generator matching (displaced specific)
    GenMatching* genMatcher;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const edm::EventSetup&, const reco::Vertex&);
};


double etaForEffectiveArea( const pat::Muon& muon );
double etaForEffectiveArea( const pat::Electron& electron );


template< typename T > double LeptonAnalyzer::getMiniIsolation( const T& lepton, const double rho, const bool onlyCharged ) const{
    auto iso = lepton.miniPFIsolation();
    double absIso;
    if( onlyCharged ){
        absIso = iso.chargedHadronIso();
    } else {
        double cone_size = 10.0 / std::min( std::max( lepton.pt(), 50. ), 200. );
        double effective_area = 0;

        if( lepton.isMuon() ){
            effective_area = muonsEffectiveAreas.getEffectiveArea( etaForEffectiveArea( lepton ) );

        } else if( lepton.isElectron() ){
            effective_area = electronsEffectiveAreasMiniIso.getEffectiveArea( etaForEffectiveArea( lepton ) );
        } else {
            throw std::invalid_argument( "getMiniIsolation is only defined for Muon and Electron objects." );
        }
        effective_area *= ( cone_size*cone_size )/ ( 0.3*0.3 );
        double pu_corr = effective_area*rho;
        absIso = iso.chargedHadronIso() + std::max( iso.neutralHadronIso() + iso.photonIso() - pu_corr, 0. ); 
    }
    return ( absIso / lepton.pt() );
}


#endif
