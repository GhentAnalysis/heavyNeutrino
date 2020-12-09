#ifndef JET_ANALYZER_H
#define JET_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class JetAnalyzer {
  friend class multilep;
  private:
    JetCorrectionUncertainty* jecUnc;

    static const unsigned nJets_max = 20;
   
    static const unsigned njecSourcesRegrouped = 12;
    std::string jecSourcesRegrouped[njecSourcesRegrouped];
    static const unsigned njecSources = 45;
    std::string jecSources[njecSources];
    std::map<std::string, JetCorrectorParameters*> jetSourcesCorParameters;
    std::map<std::string, JetCorrectorParameters*> jetRegroupedCorParameters;

    unsigned _nJets = 0;
    double   _jetPt[nJets_max];
    double   _jetPt_JECUp[nJets_max];
    double   _jetPt_JECDown[nJets_max];
    double   _jetSmearedPt[nJets_max];
    double   _jetSmearedPt_JECDown[nJets_max];
    double   _jetSmearedPt_JECUp[nJets_max];
    double   _jetSmearedPt_JERDown[nJets_max];
    double   _jetSmearedPt_JERUp[nJets_max];
    double   _jetPt_Uncorrected[nJets_max];
    double   _jetPt_L1[nJets_max];
    double   _jetPt_L2[nJets_max];
    double   _jetPt_L3[nJets_max];
    double   _jetEta[nJets_max];
    double   _jetPhi[nJets_max];
    double   _jetE[nJets_max];
    double   _jetCsvV2[nJets_max];
    double   _jetDeepCsv_udsg[nJets_max];
    double   _jetDeepCsv_b[nJets_max];
    double   _jetDeepCsv_c[nJets_max];
    double   _jetDeepCsv_bb[nJets_max];
    double   _jetDeepCsv[nJets_max];
    double   _jetDeepFlavor_b[nJets_max];
    double   _jetDeepFlavor_bb[nJets_max];
    double   _jetDeepFlavor_lepb[nJets_max];
    double   _jetDeepFlavor[nJets_max];
    double   _jetDeepFlavor_c[nJets_max];
    double   _jetDeepFlavor_uds[nJets_max];
    double   _jetDeepFlavor_g[nJets_max];
    unsigned _jetHadronFlavor[nJets_max];
    bool     _jetIsLoose[nJets_max];
    bool     _jetIsTight[nJets_max];
    bool     _jetIsTightLepVeto[nJets_max];
    double   _jetChargedHadronFraction[nJets_max];
    double   _jetNeutralHadronFraction[nJets_max];
    double   _jetNeutralEmFraction[nJets_max];
    double   _jetChargedEmFraction[nJets_max];
    double   _jetHFHadronFraction[nJets_max];
    double   _jetHFEmFraction[nJets_max];

    // Regrouped sources
    
    double   _jetPt_Absolute_JECRegroupedDown[nJets_max];
    double   _jetPt_AbsoluteUncorrelated_JECRegroupedDown[nJets_max];
    double   _jetPt_BBEC1_JECRegroupedDown[nJets_max];
    double   _jetPt_BBEC1Uncorrelated_JECRegroupedDown[nJets_max];
    double   _jetPt_EC2_JECRegroupedDown[nJets_max];
    double   _jetPt_EC2Uncorrelated_JECRegroupedDown[nJets_max];
    double   _jetPt_FlavorQCD_JECRegroupedDown[nJets_max];
    double   _jetPt_HF_JECRegroupedDown[nJets_max];
    double   _jetPt_HFUncorrelated_JECRegroupedDown[nJets_max];
    double   _jetPt_RelativeBal_JECRegroupedDown[nJets_max];
    double   _jetPt_RelativeSampleUncorrelated_JECRegroupedDown[nJets_max];
    double   _jetPt_Total_JECRegroupedDown[nJets_max];

    double   _jetPt_Absolute_JECRegroupedUp[nJets_max];
    double   _jetPt_AbsoluteUncorrelated_JECRegroupedUp[nJets_max];
    double   _jetPt_BBEC1_JECRegroupedUp[nJets_max];
    double   _jetPt_BBEC1Uncorrelated_JECRegroupedUp[nJets_max];
    double   _jetPt_EC2_JECRegroupedUp[nJets_max];
    double   _jetPt_EC2Uncorrelated_JECRegroupedUp[nJets_max];
    double   _jetPt_FlavorQCD_JECRegroupedUp[nJets_max];
    double   _jetPt_HF_JECRegroupedUp[nJets_max];
    double   _jetPt_HFUncorrelated_JECRegroupedUp[nJets_max];
    double   _jetPt_RelativeBal_JECRegroupedUp[nJets_max];
    double   _jetPt_RelativeSampleUncorrelated_JECRegroupedUp[nJets_max];
    double   _jetPt_Total_JECRegroupedUp[nJets_max];
   
    double   _jetSmearedPt_Absolute_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_AbsoluteUncorrelated_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_BBEC1_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_BBEC1Uncorrelated_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_EC2_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_EC2Uncorrelated_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_FlavorQCD_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_HF_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_HFUncorrelated_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_RelativeBal_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedDown[nJets_max];
    double   _jetSmearedPt_Total_JECRegroupedDown[nJets_max];

    double   _jetSmearedPt_Absolute_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_AbsoluteUncorrelated_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_BBEC1_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_BBEC1Uncorrelated_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_EC2_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_EC2Uncorrelated_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_FlavorQCD_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_HF_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_HFUncorrelated_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_RelativeBal_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedUp[nJets_max];
    double   _jetSmearedPt_Total_JECRegroupedUp[nJets_max];
   
    // Sources

    double   _jetPt_AbsoluteStat_JECSourcesDown[nJets_max];
    double   _jetPt_AbsoluteScale_JECSourcesDown[nJets_max];
    double   _jetPt_AbsoluteMPFBias_JECSourcesDown[nJets_max];
    double   _jetPt_Fragmentation_JECSourcesDown[nJets_max];
    double   _jetPt_SinglePionECAL_JECSourcesDown[nJets_max];
    double   _jetPt_SinglePionHCAL_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorQCD_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorZJet_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorPhotonJet_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorPureGluon_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorPureQuark_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorPureCharm_JECSourcesDown[nJets_max];
    double   _jetPt_FlavorPureBottom_JECSourcesDown[nJets_max];
    double   _jetPt_TimePtEta_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeJEREC1_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeJEREC2_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeJERHF_JECSourcesDown[nJets_max];
    double   _jetPt_RelativePtBB_JECSourcesDown[nJets_max];
    double   _jetPt_RelativePtEC1_JECSourcesDown[nJets_max];
    double   _jetPt_RelativePtEC2_JECSourcesDown[nJets_max];
    double   _jetPt_RelativePtHF_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeBal_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeSample_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeFSR_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeStatFSR_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeStatEC_JECSourcesDown[nJets_max];
    double   _jetPt_RelativeStatHF_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpDataMC_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpPtRef_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpPtBB_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpPtEC1_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpPtEC2_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpPtHF_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpMuZero_JECSourcesDown[nJets_max];
    double   _jetPt_PileUpEnvelope_JECSourcesDown[nJets_max];
    double   _jetPt_SubTotalPileUp_JECSourcesDown[nJets_max];
    double   _jetPt_SubTotalRelative_JECSourcesDown[nJets_max];
    double   _jetPt_SubTotalPt_JECSourcesDown[nJets_max];
    double   _jetPt_SubTotalScale_JECSourcesDown[nJets_max];
    double   _jetPt_SubTotalAbsolute_JECSourcesDown[nJets_max];
    double   _jetPt_SubTotalMC_JECSourcesDown[nJets_max];
    double   _jetPt_TotalNoFlavor_JECSourcesDown[nJets_max];
    double   _jetPt_TotalNoTime_JECSourcesDown[nJets_max];
    double   _jetPt_TotalNoFlavorNoTime_JECSourcesDown[nJets_max];
    double   _jetPt_Total_JECSourcesDown[nJets_max];

    double   _jetPt_AbsoluteStat_JECSourcesUp[nJets_max];
    double   _jetPt_AbsoluteScale_JECSourcesUp[nJets_max];
    double   _jetPt_AbsoluteMPFBias_JECSourcesUp[nJets_max];
    double   _jetPt_Fragmentation_JECSourcesUp[nJets_max];
    double   _jetPt_SinglePionECAL_JECSourcesUp[nJets_max];
    double   _jetPt_SinglePionHCAL_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorQCD_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorZJet_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorPhotonJet_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorPureGluon_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorPureQuark_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorPureCharm_JECSourcesUp[nJets_max];
    double   _jetPt_FlavorPureBottom_JECSourcesUp[nJets_max];
    double   _jetPt_TimePtEta_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeJEREC1_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeJEREC2_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeJERHF_JECSourcesUp[nJets_max];
    double   _jetPt_RelativePtBB_JECSourcesUp[nJets_max];
    double   _jetPt_RelativePtEC1_JECSourcesUp[nJets_max];
    double   _jetPt_RelativePtEC2_JECSourcesUp[nJets_max];
    double   _jetPt_RelativePtHF_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeBal_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeSample_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeFSR_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeStatFSR_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeStatEC_JECSourcesUp[nJets_max];
    double   _jetPt_RelativeStatHF_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpDataMC_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpPtRef_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpPtBB_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpPtEC1_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpPtEC2_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpPtHF_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpMuZero_JECSourcesUp[nJets_max];
    double   _jetPt_PileUpEnvelope_JECSourcesUp[nJets_max];
    double   _jetPt_SubTotalPileUp_JECSourcesUp[nJets_max];
    double   _jetPt_SubTotalRelative_JECSourcesUp[nJets_max];
    double   _jetPt_SubTotalPt_JECSourcesUp[nJets_max];
    double   _jetPt_SubTotalScale_JECSourcesUp[nJets_max];
    double   _jetPt_SubTotalAbsolute_JECSourcesUp[nJets_max];
    double   _jetPt_SubTotalMC_JECSourcesUp[nJets_max];
    double   _jetPt_TotalNoFlavor_JECSourcesUp[nJets_max];
    double   _jetPt_TotalNoTime_JECSourcesUp[nJets_max];
    double   _jetPt_TotalNoFlavorNoTime_JECSourcesUp[nJets_max];
    double   _jetPt_Total_JECSourcesUp[nJets_max];

    double   _jetSmearedPt_AbsoluteStat_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_AbsoluteScale_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_AbsoluteMPFBias_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_Fragmentation_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SinglePionECAL_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SinglePionHCAL_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorQCD_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorZJet_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorPhotonJet_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorPureGluon_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorPureQuark_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorPureCharm_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_FlavorPureBottom_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_TimePtEta_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeJEREC1_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeJEREC2_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeJERHF_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativePtBB_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativePtEC1_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativePtEC2_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativePtHF_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeBal_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeSample_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeFSR_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeStatFSR_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeStatEC_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_RelativeStatHF_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpDataMC_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpPtRef_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpPtBB_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpPtEC1_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpPtEC2_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpPtHF_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpMuZero_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_PileUpEnvelope_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SubTotalPileUp_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SubTotalRelative_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SubTotalPt_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SubTotalScale_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SubTotalAbsolute_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_SubTotalMC_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_TotalNoFlavor_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_TotalNoTime_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_TotalNoFlavorNoTime_JECSourcesDown[nJets_max];
    double   _jetSmearedPt_Total_JECSourcesDown[nJets_max];

    double   _jetSmearedPt_AbsoluteStat_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_AbsoluteScale_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_AbsoluteMPFBias_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_Fragmentation_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SinglePionECAL_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SinglePionHCAL_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorQCD_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorZJet_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorPhotonJet_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorPureGluon_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorPureQuark_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorPureCharm_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_FlavorPureBottom_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_TimePtEta_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeJEREC1_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeJEREC2_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeJERHF_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativePtBB_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativePtEC1_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativePtEC2_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativePtHF_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeBal_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeSample_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeFSR_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeStatFSR_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeStatEC_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_RelativeStatHF_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpDataMC_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpPtRef_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpPtBB_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpPtEC1_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpPtEC2_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpPtHF_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpMuZero_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_PileUpEnvelope_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SubTotalPileUp_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SubTotalRelative_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SubTotalPt_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SubTotalScale_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SubTotalAbsolute_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_SubTotalMC_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_TotalNoFlavor_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_TotalNoTime_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_TotalNoFlavorNoTime_JECSourcesUp[nJets_max];
    double   _jetSmearedPt_Total_JECSourcesUp[nJets_max];
   
    double   _met;                                                                              //met kinematics
    double   _metPhi;
    double   _metRaw;
    double   _metRawPhi;
    double   _metJECDown;
    double   _metPhiJECDown;
    double   _metJECUp;
    double   _metPhiJECUp;
    double   _metUnclDown;
    double   _metPhiUnclDown;
    double   _metUnclUp;
    double   _metPhiUnclUp;
    double   _metSignificance;

    multilep* multilepAnalyzer;

    bool jetIsLoose(const pat::Jet& jet, const bool is2017) const;
    bool jetIsTight(const pat::Jet& jet, const bool is2017, const bool is2018) const;
    bool jetIsTightLepVeto(const pat::Jet& jet, const bool is2017, const bool is2018) const;

  public:
    JetAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~JetAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&);
};

#endif
