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
#include "heavyNeutrino/multilep/interface/RoccoR.h"

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
    EffectiveAreas electronsEffectiveAreas_Summer16; // lepton MVA's are using old effective areas
    EffectiveAreas electronsEffectiveAreas_Spring15; // lepton MVA's are using old effective areas
    EffectiveAreas muonsEffectiveAreas;
    EffectiveAreas muonsEffectiveAreas_80X;

    //maximum number of leptons to be stored 
    static const unsigned nL_max = 20;                                                               //maximum number of particles stored

    //number of leptons of each type in the event
    unsigned _nL = 0;
    unsigned _nMu = 0;
    unsigned _nEle = 0;
    unsigned _nLight = 0;
    unsigned _nTau = 0;

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
    double _relIsoDeltaBeta[nL_max];
    double _relIso_Summer16[nL_max];
    double _relIso_80X[nL_max];
    double _relIso0p4[nL_max];
    double _relIso0p4_Summer16[nL_max];
    double _relIso0p4MuDeltaBeta[nL_max];
    double _miniIso[nL_max];
    double _miniIsoCharged[nL_max];
    double _miniIso_Spring15[nL_max];
    double _miniIso_80X[nL_max];

    //variables based on closest jet to lepton (typically containing lepton)
    double _ptRel[nL_max];
    double _ptRatio[nL_max];
    double _ptRatio_Summer16[nL_max];
    double _closestJetCsvV2[nL_max];
    double _closestJetDeepCsv_b[nL_max];
    double _closestJetDeepCsv_bb[nL_max];
    double _closestJetDeepCsv[nL_max];
    double _closestJetDeepFlavor_b[nL_max];
    double _closestJetDeepFlavor_bb[nL_max];
    double _closestJetDeepFlavor_lepb[nL_max];
    double _closestJetDeepFlavor[nL_max];
    unsigned _selectedTrackMult[nL_max];

    //pointing variables
    double _dxy[nL_max];
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _3dIPSig[nL_max];
    double _tauDxyLead[nL_max];
    double _tauDzLead[nL_max];

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

    bool _tauMuonVetoLoose[nL_max];                                                                       //tau specific variables
    bool _tauMuonVetoMVALoose[nL_max];                                                                       //tau specific variables
    bool _tauMuonVetoMVATight[nL_max];                                                                       
    bool _tauEleVetoMVAVLoose[nL_max];
    bool _tauEleVetoLoose[nL_max];
    bool _tauEleVetoMVALoose[nL_max];
    bool _tauEleVetoMVAMedium[nL_max];
    bool _tauEleVetoMVATight[nL_max];
    bool _tauEleVetoMVAVTight[nL_max];
    bool _decayModeFinding[nL_max];                      
    unsigned _tauDecayMode[nL_max];                                                             // As in https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Decay_Mode_Reconstruction 
    unsigned _tauGenStatus[nL_max];                                                             //1: prompt ele, 2:prompt mu, 3: ele from leptonic tau, 4:mu from leptonic tau, 5: hadronically decayed tau, 6:rest 
    
    bool _tauPOGVVLoose2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGVLoose2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGLoose2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGMedium2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGTight2017v2[nL_max];                                                           //version of ID to use in 94X and above
    bool _tauPOGVTight2017v2[nL_max];                                                            //Other WPs contained in _lPOG variables (vloose = veto)
    bool _tauPOGVVTight2017v2[nL_max];

    bool _decayModeFindingNew[nL_max];                                           
    bool _decayModeFindingOld[nL_max];                                           
    bool _tauVLooseMvaNew[nL_max];                                                              // # WARNING # NO LONGER SUPPORTED BY TAU POG, kept for testing reasons, will remove this soon                 

    double _tauDeepTauVsJetsRaw[nL_max];    
    bool _tauVVVLooseDeepTauVsJets[nL_max];    
    bool _tauVVLooseDeepTauVsJets[nL_max];    
    bool _tauVLooseDeepTauVsJets[nL_max];    
    bool _tauLooseDeepTauVsJets[nL_max];    
    bool _tauMediumDeepTauVsJets[nL_max];    
    bool _tauTightDeepTauVsJets[nL_max];    
    bool _tauVTightDeepTauVsJets[nL_max];    
    bool _tauVVTightDeepTauVsJets[nL_max];    
    
    double _tauDeepTauVsEleRaw[nL_max];    
    bool _tauVVVLooseDeepTauVsEle[nL_max];    
    bool _tauVVLooseDeepTauVsEle[nL_max];    
    bool _tauVLooseDeepTauVsEle[nL_max];    
    bool _tauLooseDeepTauVsEle[nL_max];    
    bool _tauMediumDeepTauVsEle[nL_max];    
    bool _tauTightDeepTauVsEle[nL_max];    
    bool _tauVTightDeepTauVsEle[nL_max];    
    bool _tauVVTightDeepTauVsEle[nL_max];    
    
    double _tauDeepTauVsMuRaw[nL_max];    
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
    
    //lepton MVA definitions
    double _leptonMvaTTH[nL_max];
    double _leptonMvatZq[nL_max];
    double _leptonMvaTOP[nL_max];

    //official POG selection definitions
    bool _lPOGVeto[nL_max];
    bool _lPOGLoose[nL_max];
    bool _lPOGMedium[nL_max];
    bool _lPOGTight[nL_max];

    //MC truth information from matching 
    bool _lIsPrompt[nL_max];
    int _lMatchPdgId[nL_max];
    int _lMatchCharge[nL_max];
    double _lMatchPt[nL_max];
    bool _lHasMatch[nL_max];
    int _lMomPdgId[nL_max];
    unsigned _lProvenance[nL_max];
    unsigned _lProvenanceCompressed[nL_max];
    unsigned _lProvenanceConversion[nL_max];

    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles);
    void fillTauGenVars(const pat::Tau&, const std::vector<reco::GenParticle>& genParticles);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron& );
    void fillLeptonImpactParameters(const pat::Muon& );
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&) const;
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose) const;
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose) const;
    void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&, const double rho, const bool oldMatching = false);

    // In leptonAnalyzerIso.cc
    double getRelIso03(const pat::Muon&, const double, const EffectiveAreas& effectiveAreas, const bool DeltaBeta=false) const;
    double getRelIso03(const pat::Electron&, const double, const EffectiveAreas& effectiveAreas) const;
    double getRelIso04(const pat::Muon&, const double, const EffectiveAreas& effectiveAreas, const bool DeltaBeta=false) const;
    double getRelIso04( const pat::Electron&, const double, const EffectiveAreas& effectiveAreas) const;
    double getRelIso(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, const EffectiveAreas& effectiveAreas, const bool onlyCharged=false) const;
    double getMiniIsolation(const reco::RecoCandidate&, edm::Handle<pat::PackedCandidateCollection>, double, double, double, double, const EffectiveAreas& effectiveAreas, bool onlyCharged=false) const;
    template< typename T > double getMiniIsolation( const T&, const double rho, const EffectiveAreas& effectiveAreas, const bool onlyCharged = false ) const;

    // In LeptonAnalyzerId.cc
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
    LeptonMvaHelper* leptonMvaComputerTOP;

    //for rochester corrections
    RoccoR rochesterCorrections;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const reco::Vertex&);
};


double etaForEffectiveArea( const pat::Muon& muon );
double etaForEffectiveArea( const pat::Electron& electron );

// Here because of template I assume
template< typename T > double LeptonAnalyzer::getMiniIsolation( const T& lepton, const double rho, const EffectiveAreas& effectiveAreas, const bool onlyCharged ) const{
    auto iso = lepton.miniPFIsolation();
    double absIso;
    if( onlyCharged ){
        absIso = iso.chargedHadronIso();
    } else {
        double cone_size = 10.0 / std::min( std::max( lepton.pt(), 50. ), 200. );
        double effective_area = effectiveAreas.getEffectiveArea( etaForEffectiveArea( lepton ) );
        effective_area *= ( cone_size*cone_size )/ ( 0.3*0.3 );
        double pu_corr = effective_area*rho;
        absIso = iso.chargedHadronIso() + std::max( iso.neutralHadronIso() + iso.photonIso() - pu_corr, 0. ); 
    }
    return ( absIso / lepton.pt() );
}


#endif
