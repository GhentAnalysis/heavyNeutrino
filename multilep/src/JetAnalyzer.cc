#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

//include c++ library classes
#include <algorithm>

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{
    jecUnc = new JetCorrectionUncertainty((iConfig.getParameter<edm::FileInPath>("jecUncertaintyFile")).fullPath());

   jecSources[0] = "AbsoluteStat";
   jecSources[1] = "AbsoluteScale";
   jecSources[2] = "AbsoluteMPFBias";
   jecSources[3] = "Fragmentation";
   jecSources[4] = "SinglePionECAL";
   jecSources[5] = "SinglePionHCAL";
   jecSources[6] = "FlavorQCD";
   jecSources[7] = "FlavorZJet";
   jecSources[8] = "FlavorPhotonJet";
   jecSources[9] = "FlavorPureGluon";
   jecSources[10] = "FlavorPureQuark";
   jecSources[11] = "FlavorPureCharm";
   jecSources[12] = "FlavorPureBottom";
   jecSources[13] = "TimePtEta";
   jecSources[14] = "RelativeJEREC1";
   jecSources[15] = "RelativeJEREC2";
   jecSources[16] = "RelativeJERHF";
   jecSources[17] = "RelativePtBB";
   jecSources[18] = "RelativePtEC1";
   jecSources[19] = "RelativePtEC2";
   jecSources[20] = "RelativePtHF";
   jecSources[21] = "RelativeBal";
   jecSources[22] = "RelativeSample";
   jecSources[23] = "RelativeFSR";
   jecSources[24] = "RelativeStatFSR";
   jecSources[25] = "RelativeStatEC";
   jecSources[26] = "RelativeStatHF";
   jecSources[27] = "PileUpDataMC";
   jecSources[28] = "PileUpPtRef";
   jecSources[29] = "PileUpPtBB";
   jecSources[30] = "PileUpPtEC1";
   jecSources[31] = "PileUpPtEC2";
   jecSources[32] = "PileUpPtHF";
   jecSources[33] = "PileUpMuZero";
   jecSources[34] = "PileUpEnvelope";
   jecSources[35] = "SubTotalPileUp";
   jecSources[36] = "SubTotalRelative";
   jecSources[37] = "SubTotalPt";
   jecSources[38] = "SubTotalScale";
   jecSources[39] = "SubTotalAbsolute";
   jecSources[40] = "SubTotalMC";
   jecSources[41] = "TotalNoFlavor";
   jecSources[42] = "TotalNoTime";
   jecSources[43] = "TotalNoFlavorNoTime";
   jecSources[44] = "Total";
   
   jecSourcesRegrouped[0] = "Absolute";
   jecSourcesRegrouped[1] = "Absolute_2016";
   jecSourcesRegrouped[2] = "BBEC1";
   jecSourcesRegrouped[3] = "BBEC1_2016";
   jecSourcesRegrouped[4] = "EC2";
   jecSourcesRegrouped[5] = "EC2_2016";
   jecSourcesRegrouped[6] = "FlavorQCD";
   jecSourcesRegrouped[7] = "HF";
   jecSourcesRegrouped[8] = "HF_2016";
   jecSourcesRegrouped[9] = "RelativeBal";
   jecSourcesRegrouped[10] = "RelativeSample_2016";
   jecSourcesRegrouped[11] = "Total";
   
    if( multilepAnalyzer->is2018() ) 
     {
	jecSourcesRegrouped[1] = "Absolute_2018";
	jecSourcesRegrouped[3] = "BBEC1_2018";
	jecSourcesRegrouped[5] = "EC2_2018";
	jecSourcesRegrouped[8] = "HF_2018";
	jecSourcesRegrouped[10] = "RelativeSample_2018";
     }
    else if( multilepAnalyzer->is2017() )
     {
	jecSourcesRegrouped[1] = "Absolute_2017";
	jecSourcesRegrouped[3] = "BBEC1_2017";
	jecSourcesRegrouped[5] = "EC2_2017";
	jecSourcesRegrouped[8] = "HF_2017";
	jecSourcesRegrouped[10] = "RelativeSample_2017";
     }

    for( unsigned int isrc=0;isrc<njecSources;isrc++ )
     {
	std::string shiftName = jecSources[isrc];
	JetCorrectorParameters *p = new JetCorrectorParameters((iConfig.getParameter<edm::FileInPath>("jecUncertaintySourcesFile")).fullPath(), shiftName);
	jetSourcesCorParameters.insert(std::make_pair(shiftName, p));
     }   
   
    for( unsigned int isrc=0;isrc<njecSourcesRegrouped;isrc++ )
     {
	std::string shiftName = jecSourcesRegrouped[isrc];
	JetCorrectorParameters *p = new JetCorrectorParameters((iConfig.getParameter<edm::FileInPath>("jecUncertaintyRegroupedFile")).fullPath(), shiftName);
	jetRegroupedCorParameters.insert(std::make_pair(shiftName, p));
     }
};

JetAnalyzer::~JetAnalyzer(){
    delete jecUnc;
    for( std::map<std::string, JetCorrectorParameters*>::iterator it=jetSourcesCorParameters.begin();it!=jetSourcesCorParameters.end();it++ ) delete it->second;
    jetSourcesCorParameters.clear();
    for( std::map<std::string, JetCorrectorParameters*>::iterator it=jetRegroupedCorParameters.begin();it!=jetRegroupedCorParameters.end();it++ ) delete it->second;
    jetRegroupedCorParameters.clear();
}

// Note: the JEC are already applied through the GT, if you need back the old way (JEC.cc) check the code before c3564f71a2e7dca3cb963ef69481894cb04bbf53
// WARNING: the _nJets is number of stored jets (i.e. including those where JECUp/JERUp passes the cut), do not use as selection
void JetAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/i");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetSmearedPt",              &_jetSmearedPt,             "_jetSmearedPt[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECDown",      &_jetSmearedPt_JECDown,     "_jetSmearedPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECUp",        &_jetSmearedPt_JECUp,       "_jetSmearedPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERDown",      &_jetSmearedPt_JERDown,     "_jetSmearedPt_JERDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERUp",        &_jetSmearedPt_JERUp,       "_jetSmearedPt_JERUp[_nJets]/D");
    outputTree->Branch("_jetPt_Uncorrected",         &_jetPt_Uncorrected,        "_jetPt_Uncorrected[_nJets]/D");
    outputTree->Branch("_jetPt_L1",                  &_jetPt_L1,                 "_jetPt_L1[_nJets]/D");
    outputTree->Branch("_jetPt_L2",                  &_jetPt_L2,                 "_jetPt_L2[_nJets]/D");
    outputTree->Branch("_jetPt_L3",                  &_jetPt_L3,                 "_jetPt_L3[_nJets]/D");

    // Sources

    if( multilepAnalyzer->storeJecSources )
     {	
	outputTree->Branch("_jetPt_AbsoluteStat_JECSourcesDown",           &_jetPt_AbsoluteStat_JECSourcesDown,           "_jetPt_AbsoluteStat_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_AbsoluteScale_JECSourcesDown",          &_jetPt_AbsoluteScale_JECSourcesDown,          "_jetPt_AbsoluteScale_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_AbsoluteMPFBias_JECSourcesDown",        &_jetPt_AbsoluteMPFBias_JECSourcesDown,        "_jetPt_AbsoluteMPFBias_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_Fragmentation_JECSourcesDown",          &_jetPt_Fragmentation_JECSourcesDown,          "_jetPt_Fragmentation_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SinglePionECAL_JECSourcesDown",         &_jetPt_SinglePionECAL_JECSourcesDown,         "_jetPt_SinglePionECAL_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SinglePionHCAL_JECSourcesDown",         &_jetPt_SinglePionHCAL_JECSourcesDown,         "_jetPt_SinglePionHCAL_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorQCD_JECSourcesDown",              &_jetPt_FlavorQCD_JECSourcesDown,              "_jetPt_FlavorQCD_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorZJet_JECSourcesDown",             &_jetPt_FlavorZJet_JECSourcesDown,             "_jetPt_FlavorZJet_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPhotonJet_JECSourcesDown",        &_jetPt_FlavorPhotonJet_JECSourcesDown,        "_jetPt_FlavorPhotonJet_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureGluon_JECSourcesDown",        &_jetPt_FlavorPureGluon_JECSourcesDown,        "_jetPt_FlavorPureGluon_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureQuark_JECSourcesDown",        &_jetPt_FlavorPureQuark_JECSourcesDown,        "_jetPt_FlavorPureQuark_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureCharm_JECSourcesDown",        &_jetPt_FlavorPureCharm_JECSourcesDown,        "_jetPt_FlavorPureCharm_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureBottom_JECSourcesDown",       &_jetPt_FlavorPureBottom_JECSourcesDown,       "_jetPt_FlavorPureBottom_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_TimePtEta_JECSourcesDown",              &_jetPt_TimePtEta_JECSourcesDown,              "_jetPt_TimePtEta_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeJEREC1_JECSourcesDown",         &_jetPt_RelativeJEREC1_JECSourcesDown,         "_jetPt_RelativeJEREC1_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeJEREC2_JECSourcesDown",         &_jetPt_RelativeJEREC2_JECSourcesDown,         "_jetPt_RelativeJEREC2_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeJERHF_JECSourcesDown",          &_jetPt_RelativeJERHF_JECSourcesDown,          "_jetPt_RelativeJERHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtBB_JECSourcesDown",           &_jetPt_RelativePtBB_JECSourcesDown,           "_jetPt_RelativePtBB_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtEC1_JECSourcesDown",          &_jetPt_RelativePtEC1_JECSourcesDown,          "_jetPt_RelativePtEC1_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtEC2_JECSourcesDown",          &_jetPt_RelativePtEC2_JECSourcesDown,          "_jetPt_RelativePtEC2_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtHF_JECSourcesDown",           &_jetPt_RelativePtHF_JECSourcesDown,           "_jetPt_RelativePtHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeBal_JECSourcesDown",            &_jetPt_RelativeBal_JECSourcesDown,            "_jetPt_RelativeBal_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeSample_JECSourcesDown",         &_jetPt_RelativeSample_JECSourcesDown,         "_jetPt_RelativeSample_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeFSR_JECSourcesDown",            &_jetPt_RelativeFSR_JECSourcesDown,            "_jetPt_RelativeFSR_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeStatFSR_JECSourcesDown",        &_jetPt_RelativeStatFSR_JECSourcesDown,        "_jetPt_RelativeStatFSR_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeStatEC_JECSourcesDown",         &_jetPt_RelativeStatEC_JECSourcesDown,         "_jetPt_RelativeStatEC_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeStatHF_JECSourcesDown",         &_jetPt_RelativeStatHF_JECSourcesDown,         "_jetPt_RelativeStatHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpDataMC_JECSourcesDown",           &_jetPt_PileUpDataMC_JECSourcesDown,           "_jetPt_PileUpDataMC_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtRef_JECSourcesDown",            &_jetPt_PileUpPtRef_JECSourcesDown,            "_jetPt_PileUpPtRef_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtBB_JECSourcesDown",             &_jetPt_PileUpPtBB_JECSourcesDown,             "_jetPt_PileUpPtBB_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtEC1_JECSourcesDown",            &_jetPt_PileUpPtEC1_JECSourcesDown,            "_jetPt_PileUpPtEC1_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtEC2_JECSourcesDown",            &_jetPt_PileUpPtEC2_JECSourcesDown,            "_jetPt_PileUpPtEC2_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtHF_JECSourcesDown",             &_jetPt_PileUpPtHF_JECSourcesDown,             "_jetPt_PileUpPtHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpMuZero_JECSourcesDown",           &_jetPt_PileUpMuZero_JECSourcesDown,           "_jetPt_PileUpMuZero_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpEnvelope_JECSourcesDown",         &_jetPt_PileUpEnvelope_JECSourcesDown,         "_jetPt_PileUpEnvelope_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalPileUp_JECSourcesDown",         &_jetPt_SubTotalPileUp_JECSourcesDown,         "_jetPt_SubTotalPileUp_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalRelative_JECSourcesDown",       &_jetPt_SubTotalRelative_JECSourcesDown,       "_jetPt_SubTotalRelative_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalPt_JECSourcesDown",             &_jetPt_SubTotalPt_JECSourcesDown,             "_jetPt_SubTotalPt_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalScale_JECSourcesDown",          &_jetPt_SubTotalScale_JECSourcesDown,          "_jetPt_SubTotalScale_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalAbsolute_JECSourcesDown",       &_jetPt_SubTotalAbsolute_JECSourcesDown,       "_jetPt_SubTotalAbsolute_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalMC_JECSourcesDown",             &_jetPt_SubTotalMC_JECSourcesDown,             "_jetPt_SubTotalMC_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_TotalNoFlavor_JECSourcesDown",          &_jetPt_TotalNoFlavor_JECSourcesDown,          "_jetPt_TotalNoFlavor_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_TotalNoTime_JECSourcesDown",            &_jetPt_TotalNoTime_JECSourcesDown,            "_jetPt_TotalNoTime_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_TotalNoFlavorNoTime_JECSourcesDown",    &_jetPt_TotalNoFlavorNoTime_JECSourcesDown,    "_jetPt_TotalNoFlavorNoTime_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetPt_Total_JECSourcesDown",                  &_jetPt_Total_JECSourcesDown,                  "_jetPt_Total_JECSourcesDown[_nJets]/D");
	
	outputTree->Branch("_jetPt_AbsoluteStat_JECSourcesUp",           &_jetPt_AbsoluteStat_JECSourcesUp,           "_jetPt_AbsoluteStat_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_AbsoluteScale_JECSourcesUp",          &_jetPt_AbsoluteScale_JECSourcesUp,          "_jetPt_AbsoluteScale_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_AbsoluteMPFBias_JECSourcesUp",        &_jetPt_AbsoluteMPFBias_JECSourcesUp,        "_jetPt_AbsoluteMPFBias_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_Fragmentation_JECSourcesUp",          &_jetPt_Fragmentation_JECSourcesUp,          "_jetPt_Fragmentation_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SinglePionECAL_JECSourcesUp",         &_jetPt_SinglePionECAL_JECSourcesUp,         "_jetPt_SinglePionECAL_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SinglePionHCAL_JECSourcesUp",         &_jetPt_SinglePionHCAL_JECSourcesUp,         "_jetPt_SinglePionHCAL_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorQCD_JECSourcesUp",              &_jetPt_FlavorQCD_JECSourcesUp,              "_jetPt_FlavorQCD_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorZJet_JECSourcesUp",             &_jetPt_FlavorZJet_JECSourcesUp,             "_jetPt_FlavorZJet_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPhotonJet_JECSourcesUp",        &_jetPt_FlavorPhotonJet_JECSourcesUp,        "_jetPt_FlavorPhotonJet_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureGluon_JECSourcesUp",        &_jetPt_FlavorPureGluon_JECSourcesUp,        "_jetPt_FlavorPureGluon_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureQuark_JECSourcesUp",        &_jetPt_FlavorPureQuark_JECSourcesUp,        "_jetPt_FlavorPureQuark_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureCharm_JECSourcesUp",        &_jetPt_FlavorPureCharm_JECSourcesUp,        "_jetPt_FlavorPureCharm_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorPureBottom_JECSourcesUp",       &_jetPt_FlavorPureBottom_JECSourcesUp,       "_jetPt_FlavorPureBottom_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_TimePtEta_JECSourcesUp",              &_jetPt_TimePtEta_JECSourcesUp,              "_jetPt_TimePtEta_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeJEREC1_JECSourcesUp",         &_jetPt_RelativeJEREC1_JECSourcesUp,         "_jetPt_RelativeJEREC1_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeJEREC2_JECSourcesUp",         &_jetPt_RelativeJEREC2_JECSourcesUp,         "_jetPt_RelativeJEREC2_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeJERHF_JECSourcesUp",          &_jetPt_RelativeJERHF_JECSourcesUp,          "_jetPt_RelativeJERHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtBB_JECSourcesUp",           &_jetPt_RelativePtBB_JECSourcesUp,           "_jetPt_RelativePtBB_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtEC1_JECSourcesUp",          &_jetPt_RelativePtEC1_JECSourcesUp,          "_jetPt_RelativePtEC1_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtEC2_JECSourcesUp",          &_jetPt_RelativePtEC2_JECSourcesUp,          "_jetPt_RelativePtEC2_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativePtHF_JECSourcesUp",           &_jetPt_RelativePtHF_JECSourcesUp,           "_jetPt_RelativePtHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeBal_JECSourcesUp",            &_jetPt_RelativeBal_JECSourcesUp,            "_jetPt_RelativeBal_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeSample_JECSourcesUp",         &_jetPt_RelativeSample_JECSourcesUp,         "_jetPt_RelativeSample_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeFSR_JECSourcesUp",            &_jetPt_RelativeFSR_JECSourcesUp,            "_jetPt_RelativeFSR_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeStatFSR_JECSourcesUp",        &_jetPt_RelativeStatFSR_JECSourcesUp,        "_jetPt_RelativeStatFSR_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeStatEC_JECSourcesUp",         &_jetPt_RelativeStatEC_JECSourcesUp,         "_jetPt_RelativeStatEC_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeStatHF_JECSourcesUp",         &_jetPt_RelativeStatHF_JECSourcesUp,         "_jetPt_RelativeStatHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpDataMC_JECSourcesUp",           &_jetPt_PileUpDataMC_JECSourcesUp,           "_jetPt_PileUpDataMC_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtRef_JECSourcesUp",            &_jetPt_PileUpPtRef_JECSourcesUp,            "_jetPt_PileUpPtRef_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtBB_JECSourcesUp",             &_jetPt_PileUpPtBB_JECSourcesUp,             "_jetPt_PileUpPtBB_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtEC1_JECSourcesUp",            &_jetPt_PileUpPtEC1_JECSourcesUp,            "_jetPt_PileUpPtEC1_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtEC2_JECSourcesUp",            &_jetPt_PileUpPtEC2_JECSourcesUp,            "_jetPt_PileUpPtEC2_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpPtHF_JECSourcesUp",             &_jetPt_PileUpPtHF_JECSourcesUp,             "_jetPt_PileUpPtHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpMuZero_JECSourcesUp",           &_jetPt_PileUpMuZero_JECSourcesUp,           "_jetPt_PileUpMuZero_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_PileUpEnvelope_JECSourcesUp",         &_jetPt_PileUpEnvelope_JECSourcesUp,         "_jetPt_PileUpEnvelope_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalPileUp_JECSourcesUp",         &_jetPt_SubTotalPileUp_JECSourcesUp,         "_jetPt_SubTotalPileUp_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalRelative_JECSourcesUp",       &_jetPt_SubTotalRelative_JECSourcesUp,       "_jetPt_SubTotalRelative_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalPt_JECSourcesUp",             &_jetPt_SubTotalPt_JECSourcesUp,             "_jetPt_SubTotalPt_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalScale_JECSourcesUp",          &_jetPt_SubTotalScale_JECSourcesUp,          "_jetPt_SubTotalScale_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalAbsolute_JECSourcesUp",       &_jetPt_SubTotalAbsolute_JECSourcesUp,       "_jetPt_SubTotalAbsolute_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_SubTotalMC_JECSourcesUp",             &_jetPt_SubTotalMC_JECSourcesUp,             "_jetPt_SubTotalMC_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_TotalNoFlavor_JECSourcesUp",          &_jetPt_TotalNoFlavor_JECSourcesUp,          "_jetPt_TotalNoFlavor_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_TotalNoTime_JECSourcesUp",            &_jetPt_TotalNoTime_JECSourcesUp,            "_jetPt_TotalNoTime_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_TotalNoFlavorNoTime_JECSourcesUp",    &_jetPt_TotalNoFlavorNoTime_JECSourcesUp,    "_jetPt_TotalNoFlavorNoTime_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetPt_Total_JECSourcesUp",                  &_jetPt_Total_JECSourcesUp,                  "_jetPt_Total_JECSourcesUp[_nJets]/D");
	
	outputTree->Branch("_jetSmearedPt_AbsoluteStat_JECSourcesDown",           &_jetSmearedPt_AbsoluteStat_JECSourcesDown,           "_jetSmearedPt_AbsoluteStat_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_AbsoluteScale_JECSourcesDown",          &_jetSmearedPt_AbsoluteScale_JECSourcesDown,          "_jetSmearedPt_AbsoluteScale_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_AbsoluteMPFBias_JECSourcesDown",        &_jetSmearedPt_AbsoluteMPFBias_JECSourcesDown,        "_jetSmearedPt_AbsoluteMPFBias_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_Fragmentation_JECSourcesDown",          &_jetSmearedPt_Fragmentation_JECSourcesDown,          "_jetSmearedPt_Fragmentation_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SinglePionECAL_JECSourcesDown",         &_jetSmearedPt_SinglePionECAL_JECSourcesDown,         "_jetSmearedPt_SinglePionECAL_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SinglePionHCAL_JECSourcesDown",         &_jetSmearedPt_SinglePionHCAL_JECSourcesDown,         "_jetSmearedPt_SinglePionHCAL_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorQCD_JECSourcesDown",              &_jetSmearedPt_FlavorQCD_JECSourcesDown,              "_jetSmearedPt_FlavorQCD_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorZJet_JECSourcesDown",             &_jetSmearedPt_FlavorZJet_JECSourcesDown,             "_jetSmearedPt_FlavorZJet_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPhotonJet_JECSourcesDown",        &_jetSmearedPt_FlavorPhotonJet_JECSourcesDown,        "_jetSmearedPt_FlavorPhotonJet_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureGluon_JECSourcesDown",        &_jetSmearedPt_FlavorPureGluon_JECSourcesDown,        "_jetSmearedPt_FlavorPureGluon_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureQuark_JECSourcesDown",        &_jetSmearedPt_FlavorPureQuark_JECSourcesDown,        "_jetSmearedPt_FlavorPureQuark_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureCharm_JECSourcesDown",        &_jetSmearedPt_FlavorPureCharm_JECSourcesDown,        "_jetSmearedPt_FlavorPureCharm_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureBottom_JECSourcesDown",       &_jetSmearedPt_FlavorPureBottom_JECSourcesDown,       "_jetSmearedPt_FlavorPureBottom_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TimePtEta_JECSourcesDown",              &_jetSmearedPt_TimePtEta_JECSourcesDown,              "_jetSmearedPt_TimePtEta_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeJEREC1_JECSourcesDown",         &_jetSmearedPt_RelativeJEREC1_JECSourcesDown,         "_jetSmearedPt_RelativeJEREC1_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeJEREC2_JECSourcesDown",         &_jetSmearedPt_RelativeJEREC2_JECSourcesDown,         "_jetSmearedPt_RelativeJEREC2_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeJERHF_JECSourcesDown",          &_jetSmearedPt_RelativeJERHF_JECSourcesDown,          "_jetSmearedPt_RelativeJERHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtBB_JECSourcesDown",           &_jetSmearedPt_RelativePtBB_JECSourcesDown,           "_jetSmearedPt_RelativePtBB_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtEC1_JECSourcesDown",          &_jetSmearedPt_RelativePtEC1_JECSourcesDown,          "_jetSmearedPt_RelativePtEC1_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtEC2_JECSourcesDown",          &_jetSmearedPt_RelativePtEC2_JECSourcesDown,          "_jetSmearedPt_RelativePtEC2_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtHF_JECSourcesDown",           &_jetSmearedPt_RelativePtHF_JECSourcesDown,           "_jetSmearedPt_RelativePtHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeBal_JECSourcesDown",            &_jetSmearedPt_RelativeBal_JECSourcesDown,            "_jetSmearedPt_RelativeBal_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeSample_JECSourcesDown",         &_jetSmearedPt_RelativeSample_JECSourcesDown,         "_jetSmearedPt_RelativeSample_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeFSR_JECSourcesDown",            &_jetSmearedPt_RelativeFSR_JECSourcesDown,            "_jetSmearedPt_RelativeFSR_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeStatFSR_JECSourcesDown",        &_jetSmearedPt_RelativeStatFSR_JECSourcesDown,        "_jetSmearedPt_RelativeStatFSR_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeStatEC_JECSourcesDown",         &_jetSmearedPt_RelativeStatEC_JECSourcesDown,         "_jetSmearedPt_RelativeStatEC_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeStatHF_JECSourcesDown",         &_jetSmearedPt_RelativeStatHF_JECSourcesDown,         "_jetSmearedPt_RelativeStatHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpDataMC_JECSourcesDown",           &_jetSmearedPt_PileUpDataMC_JECSourcesDown,           "_jetSmearedPt_PileUpDataMC_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtRef_JECSourcesDown",            &_jetSmearedPt_PileUpPtRef_JECSourcesDown,            "_jetSmearedPt_PileUpPtRef_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtBB_JECSourcesDown",             &_jetSmearedPt_PileUpPtBB_JECSourcesDown,             "_jetSmearedPt_PileUpPtBB_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtEC1_JECSourcesDown",            &_jetSmearedPt_PileUpPtEC1_JECSourcesDown,            "_jetSmearedPt_PileUpPtEC1_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtEC2_JECSourcesDown",            &_jetSmearedPt_PileUpPtEC2_JECSourcesDown,            "_jetSmearedPt_PileUpPtEC2_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtHF_JECSourcesDown",             &_jetSmearedPt_PileUpPtHF_JECSourcesDown,             "_jetSmearedPt_PileUpPtHF_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpMuZero_JECSourcesDown",           &_jetSmearedPt_PileUpMuZero_JECSourcesDown,           "_jetSmearedPt_PileUpMuZero_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpEnvelope_JECSourcesDown",         &_jetSmearedPt_PileUpEnvelope_JECSourcesDown,         "_jetSmearedPt_PileUpEnvelope_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalPileUp_JECSourcesDown",         &_jetSmearedPt_SubTotalPileUp_JECSourcesDown,         "_jetSmearedPt_SubTotalPileUp_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalRelative_JECSourcesDown",       &_jetSmearedPt_SubTotalRelative_JECSourcesDown,       "_jetSmearedPt_SubTotalRelative_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalPt_JECSourcesDown",             &_jetSmearedPt_SubTotalPt_JECSourcesDown,             "_jetSmearedPt_SubTotalPt_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalScale_JECSourcesDown",          &_jetSmearedPt_SubTotalScale_JECSourcesDown,          "_jetSmearedPt_SubTotalScale_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalAbsolute_JECSourcesDown",       &_jetSmearedPt_SubTotalAbsolute_JECSourcesDown,       "_jetSmearedPt_SubTotalAbsolute_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalMC_JECSourcesDown",             &_jetSmearedPt_SubTotalMC_JECSourcesDown,             "_jetSmearedPt_SubTotalMC_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TotalNoFlavor_JECSourcesDown",          &_jetSmearedPt_TotalNoFlavor_JECSourcesDown,          "_jetSmearedPt_TotalNoFlavor_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TotalNoTime_JECSourcesDown",            &_jetSmearedPt_TotalNoTime_JECSourcesDown,            "_jetSmearedPt_TotalNoTime_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TotalNoFlavorNoTime_JECSourcesDown",    &_jetSmearedPt_TotalNoFlavorNoTime_JECSourcesDown,    "_jetSmearedPt_TotalNoFlavorNoTime_JECSourcesDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_Total_JECSourcesDown",                  &_jetSmearedPt_Total_JECSourcesDown,                  "_jetSmearedPt_Total_JECSourcesDown[_nJets]/D");
	
	outputTree->Branch("_jetSmearedPt_AbsoluteStat_JECSourcesUp",           &_jetSmearedPt_AbsoluteStat_JECSourcesUp,           "_jetSmearedPt_AbsoluteStat_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_AbsoluteScale_JECSourcesUp",          &_jetSmearedPt_AbsoluteScale_JECSourcesUp,          "_jetSmearedPt_AbsoluteScale_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_AbsoluteMPFBias_JECSourcesUp",        &_jetSmearedPt_AbsoluteMPFBias_JECSourcesUp,        "_jetSmearedPt_AbsoluteMPFBias_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_Fragmentation_JECSourcesUp",          &_jetSmearedPt_Fragmentation_JECSourcesUp,          "_jetSmearedPt_Fragmentation_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SinglePionECAL_JECSourcesUp",         &_jetSmearedPt_SinglePionECAL_JECSourcesUp,         "_jetSmearedPt_SinglePionECAL_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SinglePionHCAL_JECSourcesUp",         &_jetSmearedPt_SinglePionHCAL_JECSourcesUp,         "_jetSmearedPt_SinglePionHCAL_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorQCD_JECSourcesUp",              &_jetSmearedPt_FlavorQCD_JECSourcesUp,              "_jetSmearedPt_FlavorQCD_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorZJet_JECSourcesUp",             &_jetSmearedPt_FlavorZJet_JECSourcesUp,             "_jetSmearedPt_FlavorZJet_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPhotonJet_JECSourcesUp",        &_jetSmearedPt_FlavorPhotonJet_JECSourcesUp,        "_jetSmearedPt_FlavorPhotonJet_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureGluon_JECSourcesUp",        &_jetSmearedPt_FlavorPureGluon_JECSourcesUp,        "_jetSmearedPt_FlavorPureGluon_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureQuark_JECSourcesUp",        &_jetSmearedPt_FlavorPureQuark_JECSourcesUp,        "_jetSmearedPt_FlavorPureQuark_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureCharm_JECSourcesUp",        &_jetSmearedPt_FlavorPureCharm_JECSourcesUp,        "_jetSmearedPt_FlavorPureCharm_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorPureBottom_JECSourcesUp",       &_jetSmearedPt_FlavorPureBottom_JECSourcesUp,       "_jetSmearedPt_FlavorPureBottom_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TimePtEta_JECSourcesUp",              &_jetSmearedPt_TimePtEta_JECSourcesUp,              "_jetSmearedPt_TimePtEta_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeJEREC1_JECSourcesUp",         &_jetSmearedPt_RelativeJEREC1_JECSourcesUp,         "_jetSmearedPt_RelativeJEREC1_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeJEREC2_JECSourcesUp",         &_jetSmearedPt_RelativeJEREC2_JECSourcesUp,         "_jetSmearedPt_RelativeJEREC2_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeJERHF_JECSourcesUp",          &_jetSmearedPt_RelativeJERHF_JECSourcesUp,          "_jetSmearedPt_RelativeJERHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtBB_JECSourcesUp",           &_jetSmearedPt_RelativePtBB_JECSourcesUp,           "_jetSmearedPt_RelativePtBB_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtEC1_JECSourcesUp",          &_jetSmearedPt_RelativePtEC1_JECSourcesUp,          "_jetSmearedPt_RelativePtEC1_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtEC2_JECSourcesUp",          &_jetSmearedPt_RelativePtEC2_JECSourcesUp,          "_jetSmearedPt_RelativePtEC2_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativePtHF_JECSourcesUp",           &_jetSmearedPt_RelativePtHF_JECSourcesUp,           "_jetSmearedPt_RelativePtHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeBal_JECSourcesUp",            &_jetSmearedPt_RelativeBal_JECSourcesUp,            "_jetSmearedPt_RelativeBal_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeSample_JECSourcesUp",         &_jetSmearedPt_RelativeSample_JECSourcesUp,         "_jetSmearedPt_RelativeSample_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeFSR_JECSourcesUp",            &_jetSmearedPt_RelativeFSR_JECSourcesUp,            "_jetSmearedPt_RelativeFSR_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeStatFSR_JECSourcesUp",        &_jetSmearedPt_RelativeStatFSR_JECSourcesUp,        "_jetSmearedPt_RelativeStatFSR_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeStatEC_JECSourcesUp",         &_jetSmearedPt_RelativeStatEC_JECSourcesUp,         "_jetSmearedPt_RelativeStatEC_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeStatHF_JECSourcesUp",         &_jetSmearedPt_RelativeStatHF_JECSourcesUp,         "_jetSmearedPt_RelativeStatHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpDataMC_JECSourcesUp",           &_jetSmearedPt_PileUpDataMC_JECSourcesUp,           "_jetSmearedPt_PileUpDataMC_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtRef_JECSourcesUp",            &_jetSmearedPt_PileUpPtRef_JECSourcesUp,            "_jetSmearedPt_PileUpPtRef_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtBB_JECSourcesUp",             &_jetSmearedPt_PileUpPtBB_JECSourcesUp,             "_jetSmearedPt_PileUpPtBB_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtEC1_JECSourcesUp",            &_jetSmearedPt_PileUpPtEC1_JECSourcesUp,            "_jetSmearedPt_PileUpPtEC1_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtEC2_JECSourcesUp",            &_jetSmearedPt_PileUpPtEC2_JECSourcesUp,            "_jetSmearedPt_PileUpPtEC2_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpPtHF_JECSourcesUp",             &_jetSmearedPt_PileUpPtHF_JECSourcesUp,             "_jetSmearedPt_PileUpPtHF_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpMuZero_JECSourcesUp",           &_jetSmearedPt_PileUpMuZero_JECSourcesUp,           "_jetSmearedPt_PileUpMuZero_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_PileUpEnvelope_JECSourcesUp",         &_jetSmearedPt_PileUpEnvelope_JECSourcesUp,         "_jetSmearedPt_PileUpEnvelope_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalPileUp_JECSourcesUp",         &_jetSmearedPt_SubTotalPileUp_JECSourcesUp,         "_jetSmearedPt_SubTotalPileUp_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalRelative_JECSourcesUp",       &_jetSmearedPt_SubTotalRelative_JECSourcesUp,       "_jetSmearedPt_SubTotalRelative_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalPt_JECSourcesUp",             &_jetSmearedPt_SubTotalPt_JECSourcesUp,             "_jetSmearedPt_SubTotalPt_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalScale_JECSourcesUp",          &_jetSmearedPt_SubTotalScale_JECSourcesUp,          "_jetSmearedPt_SubTotalScale_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalAbsolute_JECSourcesUp",       &_jetSmearedPt_SubTotalAbsolute_JECSourcesUp,       "_jetSmearedPt_SubTotalAbsolute_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_SubTotalMC_JECSourcesUp",             &_jetSmearedPt_SubTotalMC_JECSourcesUp,             "_jetSmearedPt_SubTotalMC_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TotalNoFlavor_JECSourcesUp",          &_jetSmearedPt_TotalNoFlavor_JECSourcesUp,          "_jetSmearedPt_TotalNoFlavor_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TotalNoTime_JECSourcesUp",            &_jetSmearedPt_TotalNoTime_JECSourcesUp,            "_jetSmearedPt_TotalNoTime_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_TotalNoFlavorNoTime_JECSourcesUp",    &_jetSmearedPt_TotalNoFlavorNoTime_JECSourcesUp,    "_jetSmearedPt_TotalNoFlavorNoTime_JECSourcesUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_Total_JECSourcesUp",                  &_jetSmearedPt_Total_JECSourcesUp,                  "_jetSmearedPt_Total_JECSourcesUp[_nJets]/D");
	
	// Regrouped sources
   
	outputTree->Branch("_jetPt_Absolute_JECRegroupedDown",                      &_jetPt_Absolute_JECRegroupedDown,                    "_jetPt_Absolute_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_AbsoluteUncorrelated_JECRegroupedDown",          &_jetPt_AbsoluteUncorrelated_JECRegroupedDown,        "_jetPt_AbsoluteUncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_BBEC1_JECRegroupedDown",                         &_jetPt_BBEC1_JECRegroupedDown,                       "_jetPt_BBEC1_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_BBEC1Uncorrelated_JECRegroupedDown",             &_jetPt_BBEC1Uncorrelated_JECRegroupedDown,           "_jetPt_BBEC1Uncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_EC2_JECRegroupedDown",                           &_jetPt_EC2_JECRegroupedDown,                         "_jetPt_EC2_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_EC2Uncorrelated_JECRegroupedDown",               &_jetPt_EC2Uncorrelated_JECRegroupedDown,             "_jetPt_EC2Uncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorQCD_JECRegroupedDown",                     &_jetPt_FlavorQCD_JECRegroupedDown,                   "_jetPt_FlavorQCD_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_HF_JECRegroupedDown",                            &_jetPt_HF_JECRegroupedDown,                          "_jetPt_HF_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_HFUncorrelated_JECRegroupedDown",                &_jetPt_HFUncorrelated_JECRegroupedDown,              "_jetPt_HFUncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeBal_JECRegroupedDown",                   &_jetPt_RelativeBal_JECRegroupedDown,                 "_jetPt_RelativeBal_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeSampleUncorrelated_JECRegroupedDown",    &_jetPt_RelativeSampleUncorrelated_JECRegroupedDown,  "_jetPt_RelativeSampleUncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetPt_Total_JECRegroupedDown",                         &_jetPt_Total_JECRegroupedDown,                       "_jetPt_Total_JECRegroupedDown[_nJets]/D");
	
	outputTree->Branch("_jetPt_Absolute_JECRegroupedUp",                      &_jetPt_Absolute_JECRegroupedUp,                    "_jetPt_Absolute_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_AbsoluteUncorrelated_JECRegroupedUp",          &_jetPt_AbsoluteUncorrelated_JECRegroupedUp,        "_jetPt_AbsoluteUncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_BBEC1_JECRegroupedUp",                         &_jetPt_BBEC1_JECRegroupedUp,                       "_jetPt_BBEC1_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_BBEC1Uncorrelated_JECRegroupedUp",             &_jetPt_BBEC1Uncorrelated_JECRegroupedUp,           "_jetPt_BBEC1Uncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_EC2_JECRegroupedUp",                           &_jetPt_EC2_JECRegroupedUp,                         "_jetPt_EC2_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_EC2Uncorrelated_JECRegroupedUp",               &_jetPt_EC2Uncorrelated_JECRegroupedUp,             "_jetPt_EC2Uncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_FlavorQCD_JECRegroupedUp",                     &_jetPt_FlavorQCD_JECRegroupedUp,                   "_jetPt_FlavorQCD_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_HF_JECRegroupedUp",                            &_jetPt_HF_JECRegroupedUp,                          "_jetPt_HF_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_HFUncorrelated_JECRegroupedUp",                &_jetPt_HFUncorrelated_JECRegroupedUp,              "_jetPt_HFUncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeBal_JECRegroupedUp",                   &_jetPt_RelativeBal_JECRegroupedUp,                 "_jetPt_RelativeBal_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_RelativeSampleUncorrelated_JECRegroupedUp",    &_jetPt_RelativeSampleUncorrelated_JECRegroupedUp,  "_jetPt_RelativeSampleUncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetPt_Total_JECRegroupedUp",                         &_jetPt_Total_JECRegroupedUp,                       "_jetPt_Total_JECRegroupedUp[_nJets]/D");
	
	outputTree->Branch("_jetSmearedPt_Absolute_JECRegroupedDown",                      &_jetSmearedPt_Absolute_JECRegroupedDown,                    "_jetSmearedPt_Absolute_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_AbsoluteUncorrelated_JECRegroupedDown",          &_jetSmearedPt_AbsoluteUncorrelated_JECRegroupedDown,        "_jetSmearedPt_AbsoluteUncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_BBEC1_JECRegroupedDown",                         &_jetSmearedPt_BBEC1_JECRegroupedDown,                       "_jetSmearedPt_BBEC1_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_BBEC1Uncorrelated_JECRegroupedDown",             &_jetSmearedPt_BBEC1Uncorrelated_JECRegroupedDown,           "_jetSmearedPt_BBEC1Uncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_EC2_JECRegroupedDown",                           &_jetSmearedPt_EC2_JECRegroupedDown,                         "_jetSmearedPt_EC2_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_EC2Uncorrelated_JECRegroupedDown",               &_jetSmearedPt_EC2Uncorrelated_JECRegroupedDown,             "_jetSmearedPt_EC2Uncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorQCD_JECRegroupedDown",                     &_jetSmearedPt_FlavorQCD_JECRegroupedDown,                   "_jetSmearedPt_FlavorQCD_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_HF_JECRegroupedDown",                            &_jetSmearedPt_HF_JECRegroupedDown,                          "_jetSmearedPt_HF_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_HFUncorrelated_JECRegroupedDown",                &_jetSmearedPt_HFUncorrelated_JECRegroupedDown,              "_jetSmearedPt_HFUncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeBal_JECRegroupedDown",                   &_jetSmearedPt_RelativeBal_JECRegroupedDown,                 "_jetSmearedPt_RelativeBal_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedDown",    &_jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedDown,  "_jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedDown[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_Total_JECRegroupedDown",                         &_jetSmearedPt_Total_JECRegroupedDown,                       "_jetSmearedPt_Total_JECRegroupedDown[_nJets]/D");
	
	outputTree->Branch("_jetSmearedPt_Absolute_JECRegroupedUp",                      &_jetSmearedPt_Absolute_JECRegroupedUp,                    "_jetSmearedPt_Absolute_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_AbsoluteUncorrelated_JECRegroupedUp",          &_jetSmearedPt_AbsoluteUncorrelated_JECRegroupedUp,        "_jetSmearedPt_AbsoluteUncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_BBEC1_JECRegroupedUp",                         &_jetSmearedPt_BBEC1_JECRegroupedUp,                       "_jetSmearedPt_BBEC1_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_BBEC1Uncorrelated_JECRegroupedUp",             &_jetSmearedPt_BBEC1Uncorrelated_JECRegroupedUp,           "_jetSmearedPt_BBEC1Uncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_EC2_JECRegroupedUp",                           &_jetSmearedPt_EC2_JECRegroupedUp,                         "_jetSmearedPt_EC2_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_EC2Uncorrelated_JECRegroupedUp",               &_jetSmearedPt_EC2Uncorrelated_JECRegroupedUp,             "_jetSmearedPt_EC2Uncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_FlavorQCD_JECRegroupedUp",                     &_jetSmearedPt_FlavorQCD_JECRegroupedUp,                   "_jetSmearedPt_FlavorQCD_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_HF_JECRegroupedUp",                            &_jetSmearedPt_HF_JECRegroupedUp,                          "_jetSmearedPt_HF_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_HFUncorrelated_JECRegroupedUp",                &_jetSmearedPt_HFUncorrelated_JECRegroupedUp,              "_jetSmearedPt_HFUncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeBal_JECRegroupedUp",                   &_jetSmearedPt_RelativeBal_JECRegroupedUp,                 "_jetSmearedPt_RelativeBal_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedUp",    &_jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedUp,  "_jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedUp[_nJets]/D");
	outputTree->Branch("_jetSmearedPt_Total_JECRegroupedUp",                         &_jetSmearedPt_Total_JECRegroupedUp,                       "_jetSmearedPt_Total_JECRegroupedUp[_nJets]/D");
     }
   
    outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
    outputTree->Branch("_jetDeepCsv",                &_jetDeepCsv,               "_jetDeepCsv[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_b",           &_jetDeepFlavor_b,          "_jetDeepFlavor_b[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_bb",          &_jetDeepFlavor_bb,         "_jetDeepFlavor_bb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_lepb",        &_jetDeepFlavor_lepb,       "_jetDeepFlavor_lepb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor",             &_jetDeepFlavor,            "_jetDeepFlavor[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_c",           &_jetDeepFlavor_c,          "_jetDeepFlavor_c[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_uds",         &_jetDeepFlavor_uds,        "_jetDeepFlavor_uds[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_g",           &_jetDeepFlavor_g,          "_jetDeepFlavor_g[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");

    outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
    outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
    outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
    outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
    outputTree->Branch("_jetHFHadronFraction",       &_jetHFHadronFraction,      "_jetHFHadronFraction[_nJets]/D");
    outputTree->Branch("_jetHFEmFraction",           &_jetHFEmFraction,          "_jetHFEmFraction[_nJets]/D");

    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metRaw",                       &_metRaw,                       "_metRaw/D");
    outputTree->Branch("_metJECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_metJECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");

    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metRawPhi",                    &_metRawPhi,                    "_metRawPhi/D");
    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");

    if(!multilepAnalyzer->is2018() ) outputTree->Branch("_jetIsLoose", _jetIsLoose, "_jetIsLoose[_nJets]/O"); // WARNING, not recommended to be used, only exists for 2016
}

bool JetAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Jet>> jets            = getHandle(iEvent, multilepAnalyzer->jetToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmeared     = getHandle(iEvent, multilepAnalyzer->jetSmearedToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedUp   = getHandle(iEvent, multilepAnalyzer->jetSmearedUpToken);
    edm::Handle<std::vector<pat::Jet>> jetsSmearedDown = getHandle(iEvent, multilepAnalyzer->jetSmearedDownToken);
    edm::Handle<std::vector<pat::MET>> mets            = getHandle(iEvent, multilepAnalyzer->metToken);

    _nJets = 0;

    for(const auto& jet : *jets){
        if(_nJets == nJets_max) break;

        _jetIsLoose[_nJets]        = jetIsLoose(jet, multilepAnalyzer->is2017() || multilepAnalyzer->is2018() );
        _jetIsTight[_nJets]        = jetIsTight(jet, multilepAnalyzer->is2017(), multilepAnalyzer->is2018() );
        _jetIsTightLepVeto[_nJets] = jetIsTightLepVeto(jet, multilepAnalyzer->is2017(), multilepAnalyzer->is2018() );

        //find smeared equivalents of nominal jet
        auto jetSmearedIt = jetsSmeared->begin();
        for(auto j = jetsSmeared->cbegin(); j != jetsSmeared->cend(); ++j){
            if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedIt)) jetSmearedIt = j;
        }

        auto jetSmearedUpIt = jetsSmearedUp->begin();
        for(auto j = jetsSmearedUp->cbegin(); j != jetsSmearedUp->cend(); ++j){
            if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedUpIt))  jetSmearedUpIt = j;
        }

        auto jetSmearedDownIt = jetsSmearedDown->begin();
        for(auto j = jetsSmearedDown->cbegin(); j != jetsSmearedDown->cend(); ++j){
            if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedDownIt))  jetSmearedDownIt = j;
        }

        //nominal jet pt and uncertainties
        jecUnc->setJetEta(jet.eta());
        jecUnc->setJetPt(jet.pt());
        double unc = jecUnc->getUncertainty(true);
        _jetPt[_nJets]         = jet.pt();
        _jetPt_JECDown[_nJets] = _jetPt[_nJets]*(1 - unc);
        _jetPt_JECUp[_nJets]   = _jetPt[_nJets]*(1 + unc);
       
        //smeared jet pt and uncertainties
        jecUnc->setJetEta(jetSmearedIt->eta());
        jecUnc->setJetPt(jetSmearedIt->pt());
        double uncSmeared = jecUnc->getUncertainty(true);
        _jetSmearedPt[_nJets]         = jetSmearedIt->pt();
        _jetSmearedPt_JECDown[_nJets] = _jetSmearedPt[_nJets]*( 1. - uncSmeared );
        _jetSmearedPt_JECUp[_nJets]   = _jetSmearedPt[_nJets]*( 1. + uncSmeared );
        _jetSmearedPt_JERDown[_nJets] = jetSmearedDownIt->pt();
        _jetSmearedPt_JERUp[_nJets]   = jetSmearedUpIt->pt();

        // Regrouped sources

        if( multilepAnalyzer->storeJecSources )
	 {	    
	    for( unsigned int isrc=0;isrc<njecSourcesRegrouped;isrc++ )
	      {
		 std::string shiftName = jecSourcesRegrouped[isrc];
		 JetCorrectionUncertainty *jetCorUnc = new JetCorrectionUncertainty(*jetRegroupedCorParameters[shiftName]);
		 
		 jetCorUnc->setJetPt(jet.pt());
		 jetCorUnc->setJetEta(jet.eta());
		 double uncJec = jetCorUnc->getUncertainty(true);
		 double jetPtDown = _jetPt[_nJets]*(1 - uncJec);
		 double jetPtUp = _jetPt[_nJets]*(1 + uncJec);
		 
		 if( shiftName == "Absolute" )
		   {		 
		      _jetPt_Absolute_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_Absolute_JECRegroupedUp[_nJets]   = jetPtUp;
		   }	    
		 else if( shiftName == "Absolute_2018" || shiftName == "Absolute_2017" || shiftName == "Absolute_2016" )
		   {		 
		      _jetPt_AbsoluteUncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_AbsoluteUncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "BBEC1" )
		   {		 
		      _jetPt_BBEC1_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_BBEC1_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "BBEC1_2018" || shiftName == "BBEC1_2017" || shiftName == "BBEC1_2016" )
		   {		 
		      _jetPt_BBEC1Uncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_BBEC1Uncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "EC2" )
		   {		 
		      _jetPt_EC2_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_EC2_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "EC2_2018" || shiftName == "EC2_2017" || shiftName == "EC2_2016" )
		   {		 
		      _jetPt_EC2Uncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_EC2Uncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorQCD" )
		   {		 
		      _jetPt_FlavorQCD_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_FlavorQCD_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "HF" )
		   {		 
		      _jetPt_HF_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_HF_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "HF_2018" || shiftName == "HF_2017" || shiftName == "HF_2016" )
		   {		 
		      _jetPt_HFUncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_HFUncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeBal" )
		   {		 
		      _jetPt_RelativeBal_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_RelativeBal_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeSample_2016" || shiftName == "RelativeSample_2017" || shiftName == "RelativeSample_2018" )
		   {		 
		      _jetPt_RelativeSampleUncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_RelativeSampleUncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "Total" )
		   {		 
		      _jetPt_Total_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetPt_Total_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 
		 jetCorUnc->setJetPt(jetSmearedIt->pt());
		 jetCorUnc->setJetEta(jetSmearedIt->eta());
		 uncJec = jetCorUnc->getUncertainty(true);
		 jetPtDown = _jetSmearedPt[_nJets]*(1 - uncJec);
		 jetPtUp = _jetSmearedPt[_nJets]*(1 + uncJec);
		 
		 if( shiftName == "Absolute" )
		   {		 
		      _jetSmearedPt_Absolute_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_Absolute_JECRegroupedUp[_nJets]   = jetPtUp;
		   }	    
		 else if( shiftName == "Absolute_2018" || shiftName == "Absolute_2017" || shiftName == "Absolute_2016" )
		   {		 
		      _jetSmearedPt_AbsoluteUncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_AbsoluteUncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "BBEC1" )
		   {		 
		      _jetSmearedPt_BBEC1_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_BBEC1_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "BBEC1_2018" || shiftName == "BBEC1_2017" || shiftName == "BBEC1_2016" )
		   {		 
		      _jetSmearedPt_BBEC1Uncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_BBEC1Uncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "EC2" )
		   {		 
		      _jetSmearedPt_EC2_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_EC2_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "EC2_2018" || shiftName == "EC2_2017" || shiftName == "EC2_2016" )
		   {		 
		      _jetSmearedPt_EC2Uncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_EC2Uncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorQCD" )
		   {		 
		      _jetSmearedPt_FlavorQCD_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorQCD_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "HF" )
		   {		 
		      _jetSmearedPt_HF_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_HF_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "HF_2018" || shiftName == "HF_2017" || shiftName == "HF_2016" )
		   {		 
		      _jetSmearedPt_HFUncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_HFUncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeBal" )
		   {		 
		      _jetSmearedPt_RelativeBal_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeBal_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeSample_2016" || shiftName == "RelativeSample_2017" || shiftName == "RelativeSample_2018" )
		   {		 
		      _jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeSampleUncorrelated_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "Total" )
		   {		 
		      _jetSmearedPt_Total_JECRegroupedDown[_nJets] = jetPtDown;
		      _jetSmearedPt_Total_JECRegroupedUp[_nJets]   = jetPtUp;
		   }
		 
		 delete jetCorUnc;
	      }

	    // Sources
	    
	    for( unsigned int isrc=0;isrc<njecSources;isrc++ )
	      {
		 std::string shiftName = jecSources[isrc];
		 JetCorrectionUncertainty *jetCorUnc = new JetCorrectionUncertainty(*jetSourcesCorParameters[shiftName]);
		 
		 jetCorUnc->setJetPt(jet.pt());
		 jetCorUnc->setJetEta(jet.eta());
		 double uncJec = jetCorUnc->getUncertainty(true);
		 double jetPtDown = _jetPt[_nJets]*(1 - uncJec);
		 double jetPtUp = _jetPt[_nJets]*(1 + uncJec);
		 
		 if( shiftName == "AbsoluteStat" )
		   {		 
		      _jetPt_AbsoluteStat_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_AbsoluteStat_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "AbsoluteScale" )
		   {		 
		      _jetPt_AbsoluteScale_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_AbsoluteScale_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "AbsoluteMPFBias" )
		   {		 
		      _jetPt_AbsoluteMPFBias_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_AbsoluteMPFBias_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "Fragmentation" )
		   {		 
		      _jetPt_Fragmentation_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_Fragmentation_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SinglePionECAL" )
		   {		 
		      _jetPt_SinglePionECAL_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SinglePionECAL_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SinglePionHCAL" )
		   {		 
		      _jetPt_SinglePionHCAL_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SinglePionHCAL_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorQCD" )
		   {		 
		      _jetPt_FlavorQCD_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorQCD_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorZJet" )
		   {		 
		      _jetPt_FlavorZJet_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorZJet_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPhotonJet" )
		   {		 
		      _jetPt_FlavorPhotonJet_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorPhotonJet_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureGluon" )
		   {		 
		      _jetPt_FlavorPureGluon_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorPureGluon_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureQuark" )
		   {		 
		      _jetPt_FlavorPureQuark_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorPureQuark_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureCharm" )
		   {		 
		      _jetPt_FlavorPureCharm_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorPureCharm_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureBottom" )
		   {		 
		      _jetPt_FlavorPureBottom_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_FlavorPureBottom_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TimePtEta" )
		   {		 
		      _jetPt_TimePtEta_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_TimePtEta_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeJEREC1" )
		   {		 
		      _jetPt_RelativeJEREC1_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeJEREC1_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeJEREC2" )
		   {		 
		      _jetPt_RelativeJEREC2_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeJEREC2_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeJERHF" )
		   {		 
		      _jetPt_RelativeJERHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeJERHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtBB" )
		   {		 
		      _jetPt_RelativePtBB_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativePtBB_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtEC1" )
		   {		 
		      _jetPt_RelativePtEC1_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativePtEC1_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtEC2" )
		   {		 
		      _jetPt_RelativePtEC2_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativePtEC2_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtHF" )
		   {		 
		      _jetPt_RelativePtHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativePtHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeBal" )
		   {		 
		      _jetPt_RelativeBal_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeBal_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeSample" )
		   {		 
		      _jetPt_RelativeSample_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeSample_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeFSR" )
		   {		 
		      _jetPt_RelativeFSR_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeFSR_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeStatFSR" )
		   {		 
		      _jetPt_RelativeStatFSR_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeStatFSR_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeStatEC" )
		   {		 
		      _jetPt_RelativeStatEC_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeStatEC_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeStatHF" )
		   {		 
		      _jetPt_RelativeStatHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_RelativeStatHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpDataMC" )
		   {		 
		      _jetPt_PileUpDataMC_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpDataMC_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtRef" )
		   {		 
		      _jetPt_PileUpPtRef_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpPtRef_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtBB" )
		   {		 
		      _jetPt_PileUpPtBB_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpPtBB_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtEC1" )
		   {		 
		      _jetPt_PileUpPtEC1_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpPtEC1_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtEC2" )
		   {		 
		      _jetPt_PileUpPtEC2_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpPtEC2_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtHF" )
		   {		 
		      _jetPt_PileUpPtHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpPtHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpMuZero" )
		   {		 
		      _jetPt_PileUpMuZero_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpMuZero_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpEnvelope" )
		   {		 
		      _jetPt_PileUpEnvelope_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_PileUpEnvelope_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalPileUp" )
		   {		 
		      _jetPt_SubTotalPileUp_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SubTotalPileUp_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalRelative" )
		   {		 
		      _jetPt_SubTotalRelative_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SubTotalRelative_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalPt" )
		   {		 
		      _jetPt_SubTotalPt_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SubTotalPt_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalScale" )
		   {		 
		      _jetPt_SubTotalScale_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SubTotalScale_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalAbsolute" )
		   {		 
		      _jetPt_SubTotalAbsolute_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SubTotalAbsolute_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalMC" )
		   {		 
		      _jetPt_SubTotalMC_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_SubTotalMC_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TotalNoFlavor" )
		   {		 
		      _jetPt_TotalNoFlavor_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_TotalNoFlavor_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TotalNoTime" )
		   {		 
		      _jetPt_TotalNoTime_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_TotalNoTime_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TotalNoFlavorNoTime" )
		   {		 
		      _jetPt_TotalNoFlavorNoTime_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_TotalNoFlavorNoTime_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "Total" )
		   {		 
		      _jetPt_Total_JECSourcesDown[_nJets] = jetPtDown;
		      _jetPt_Total_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 
		 jetCorUnc->setJetPt(jetSmearedIt->pt());
		 jetCorUnc->setJetEta(jetSmearedIt->eta());
		 uncJec = jetCorUnc->getUncertainty(true);
		 jetPtDown = _jetSmearedPt[_nJets]*(1 - uncJec);
		 jetPtUp = _jetSmearedPt[_nJets]*(1 + uncJec);
		 
		 if( shiftName == "AbsoluteStat" )
		   {		 
		      _jetSmearedPt_AbsoluteStat_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_AbsoluteStat_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "AbsoluteScale" )
		   {		 
		      _jetSmearedPt_AbsoluteScale_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_AbsoluteScale_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "AbsoluteMPFBias" )
		   {		 
		      _jetSmearedPt_AbsoluteMPFBias_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_AbsoluteMPFBias_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "Fragmentation" )
		   {		 
		      _jetSmearedPt_Fragmentation_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_Fragmentation_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SinglePionECAL" )
		   {		 
		      _jetSmearedPt_SinglePionECAL_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SinglePionECAL_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SinglePionHCAL" )
		   {		 
		      _jetSmearedPt_SinglePionHCAL_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SinglePionHCAL_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorQCD" )
		   {		 
		      _jetSmearedPt_FlavorQCD_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorQCD_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorZJet" )
		   {		 
		      _jetSmearedPt_FlavorZJet_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorZJet_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPhotonJet" )
		   {		 
		      _jetSmearedPt_FlavorPhotonJet_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorPhotonJet_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureGluon" )
		   {		 
		      _jetSmearedPt_FlavorPureGluon_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorPureGluon_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureQuark" )
		   {		 
		      _jetSmearedPt_FlavorPureQuark_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorPureQuark_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureCharm" )
		   {		 
		      _jetSmearedPt_FlavorPureCharm_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorPureCharm_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "FlavorPureBottom" )
		   {		 
		      _jetSmearedPt_FlavorPureBottom_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_FlavorPureBottom_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TimePtEta" )
		   {		 
		      _jetSmearedPt_TimePtEta_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_TimePtEta_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeJEREC1" )
		   {		 
		      _jetSmearedPt_RelativeJEREC1_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeJEREC1_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeJEREC2" )
		   {		 
		      _jetSmearedPt_RelativeJEREC2_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeJEREC2_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeJERHF" )
		   {		 
		      _jetSmearedPt_RelativeJERHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeJERHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtBB" )
		   {		 
		      _jetSmearedPt_RelativePtBB_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativePtBB_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtEC1" )
		   {		 
		      _jetSmearedPt_RelativePtEC1_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativePtEC1_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtEC2" )
		   {		 
		      _jetSmearedPt_RelativePtEC2_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativePtEC2_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativePtHF" )
		   {		 
		      _jetSmearedPt_RelativePtHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativePtHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeBal" )
		   {		 
		      _jetSmearedPt_RelativeBal_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeBal_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeSample" )
		   {		 
		      _jetSmearedPt_RelativeSample_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeSample_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeFSR" )
		   {		 
		      _jetSmearedPt_RelativeFSR_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeFSR_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeStatFSR" )
		   {		 
		      _jetSmearedPt_RelativeStatFSR_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeStatFSR_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeStatEC" )
		   {		 
		      _jetSmearedPt_RelativeStatEC_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeStatEC_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "RelativeStatHF" )
		   {		 
		      _jetSmearedPt_RelativeStatHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_RelativeStatHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpDataMC" )
		   {		 
		      _jetSmearedPt_PileUpDataMC_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpDataMC_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtRef" )
		   {		 
		      _jetSmearedPt_PileUpPtRef_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpPtRef_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtBB" )
		   {		 
		      _jetSmearedPt_PileUpPtBB_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpPtBB_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtEC1" )
		   {		 
		      _jetSmearedPt_PileUpPtEC1_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpPtEC1_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtEC2" )
		   {		 
		      _jetSmearedPt_PileUpPtEC2_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpPtEC2_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpPtHF" )
		   {		 
		      _jetSmearedPt_PileUpPtHF_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpPtHF_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpMuZero" )
		   {		 
		      _jetSmearedPt_PileUpMuZero_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpMuZero_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "PileUpEnvelope" )
		   {		 
		      _jetSmearedPt_PileUpEnvelope_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_PileUpEnvelope_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalPileUp" )
		   {		 
		      _jetSmearedPt_SubTotalPileUp_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SubTotalPileUp_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalRelative" )
		   {		 
		      _jetSmearedPt_SubTotalRelative_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SubTotalRelative_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalPt" )
		   {		 
		      _jetSmearedPt_SubTotalPt_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SubTotalPt_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalScale" )
		   {		 
		      _jetSmearedPt_SubTotalScale_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SubTotalScale_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalAbsolute" )
		   {		 
		      _jetSmearedPt_SubTotalAbsolute_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SubTotalAbsolute_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "SubTotalMC" )
		   {		 
		      _jetSmearedPt_SubTotalMC_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_SubTotalMC_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TotalNoFlavor" )
		   {		 
		      _jetSmearedPt_TotalNoFlavor_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_TotalNoFlavor_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TotalNoTime" )
		   {		 
		      _jetSmearedPt_TotalNoTime_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_TotalNoTime_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "TotalNoFlavorNoTime" )
		   {		 
		      _jetSmearedPt_TotalNoFlavorNoTime_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_TotalNoFlavorNoTime_JECSourcesUp[_nJets]   = jetPtUp;
		   }
		 else if( shiftName == "Total" )
		   {		 
		      _jetSmearedPt_Total_JECSourcesDown[_nJets] = jetPtDown;
		      _jetSmearedPt_Total_JECSourcesUp[_nJets]   = jetPtUp;
		   }
	      }
	 }       
       
        //find maximum of all pT variations
        std::vector<double> ptVector = {_jetPt[_nJets], _jetPt_JECDown[_nJets], _jetPt_JECUp[_nJets],
            _jetSmearedPt[_nJets], _jetSmearedPt_JECDown[_nJets], _jetSmearedPt_JECUp[_nJets], _jetSmearedPt_JERDown[_nJets],  _jetSmearedPt_JERUp[_nJets]};
        double maxpT = *(std::max_element(ptVector.cbegin(), ptVector.cend()));
        if(maxpT <= 20) continue;

        _jetPt_Uncorrected[_nJets]        = jet.correctedP4("Uncorrected").Pt();
        _jetPt_L1[_nJets]                 = jet.correctedP4("L1FastJet").Pt();
        _jetPt_L2[_nJets]                 = jet.correctedP4("L2Relative").Pt();
        _jetPt_L3[_nJets]                 = jet.correctedP4("L3Absolute").Pt();

        _jetEta[_nJets]                   = jet.eta();
        _jetPhi[_nJets]                   = jet.phi();
        _jetE[_nJets]                     = jet.energy();

        //Old csvV2 b-tagger
        _jetCsvV2[_nJets]                 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        //DeepCSV tagger
        _jetDeepCsv_udsg[_nJets]          = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
        _jetDeepCsv_b[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _jetDeepCsv_c[_nJets]             = jet.bDiscriminator("pfDeepCSVJetTags:probc");
        _jetDeepCsv_bb[_nJets]            = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _jetDeepCsv[_nJets]               = _jetDeepCsv_b[_nJets] + _jetDeepCsv_bb[_nJets];
        if( std::isnan( _jetDeepCsv[_nJets] ) ) _jetDeepCsv[_nJets] = 0.;

        //DeepFlavor taggeer 
        _jetDeepFlavor_b[_nJets]          = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
        _jetDeepFlavor_bb[_nJets]         = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
        _jetDeepFlavor_lepb[_nJets]       = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        _jetDeepFlavor[_nJets]            = _jetDeepFlavor_b[_nJets] + _jetDeepFlavor_bb[_nJets] + _jetDeepFlavor_lepb[_nJets];
        if( std::isnan( _jetDeepFlavor[_nJets] ) ) _jetDeepFlavor[_nJets] = 0.;
        _jetDeepFlavor_c[_nJets]          = jet.bDiscriminator("pfDeepFlavourJetTags:probc");
        _jetDeepFlavor_uds[_nJets]        = jet.bDiscriminator("pfDeepFlavourJetTags:probuds");
        _jetDeepFlavor_g[_nJets]          = jet.bDiscriminator("pfDeepFlavourJetTags:probg");

        _jetHadronFlavor[_nJets]          = jet.hadronFlavour();

        _jetNeutralHadronFraction[_nJets] = jet.neutralHadronEnergyFraction();
        _jetChargedHadronFraction[_nJets] = jet.chargedHadronEnergyFraction();
        _jetNeutralEmFraction[_nJets]     = jet.neutralEmEnergyFraction();
        _jetChargedEmFraction[_nJets]     = jet.chargedEmEnergyFraction();
        _jetHFHadronFraction[_nJets]      = jet.HFHadronEnergyFraction();
        _jetHFEmFraction[_nJets]          = jet.HFEMEnergyFraction();

        ++_nJets;
    }

    //determine the met of the event and its uncertainties
    //nominal MET value
    const pat::MET& met = (*mets).front();
    _met             = met.pt();
    _metPhi          = met.phi();

    //raw met values
    _metRaw          = met.uncorPt();
    _metRawPhi       = met.uncorPhi();
    //met values with uncertainties varied up and down
    _metJECDown      = met.shiftedPt(pat::MET::JetEnDown);
    _metJECUp        = met.shiftedPt(pat::MET::JetEnUp);
    _metUnclDown     = met.shiftedPt(pat::MET::UnclusteredEnDown);
    _metUnclUp       = met.shiftedPt(pat::MET::UnclusteredEnUp);
    _metPhiJECDown   = met.shiftedPhi(pat::MET::JetEnDown);
    _metPhiJECUp     = met.shiftedPhi(pat::MET::JetEnUp);
    _metPhiUnclUp    = met.shiftedPhi(pat::MET::UnclusteredEnUp);
    _metPhiUnclDown  = met.shiftedPhi(pat::MET::UnclusteredEnDown);

    //significance of met
    //note: this is the only one variable which changed between 94X and 102X see https://github.com/cms-sw/cmssw/commit/f7aacfd2ffaac9899ea07d0355afe49bb10a0aeb
    _metSignificance = met.metSignificance();

    if(multilepAnalyzer->skim == "singlejet" and _nJets < 1) return false;
    if(multilepAnalyzer->skim == "FR" and _nJets < 1)        return false;
    return true;
}

/*
 * JetID implementations, references:
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
 * https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
 */

// WARNING: There should be no loose Jet ID for 2017, not sure where the cuts for this below originate, so use at own risk
bool JetAnalyzer::jetIsLoose(const pat::Jet& jet, const bool is2017) const{
    if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >=  0.99)               return false;
        if(jet.neutralEmEnergyFraction() >= 0.99)                    return false;
        if(jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1) return false;
        if(fabs(jet.eta()) <= 2.4){
            if(jet.chargedHadronEnergyFraction() <= 0)               return false;
            if(jet.chargedMultiplicity() <= 0)                       return false;
            if(!is2017 && jet.chargedEmEnergyFraction()>= 0.99)      return false;
        }

    } else if(fabs(jet.eta()) <= 3.0){
        if(jet.neutralHadronEnergyFraction() >= 0.98)                return false;
        if(jet.neutralEmEnergyFraction() <= (is2017 ? 0.02 : 0.01))  return false;
        if(jet.neutralMultiplicity() <= 2)                           return false;

    } else {
        if(jet.neutralEmEnergyFraction() >= 0.90)                    return false;
        if(jet.neutralMultiplicity() <= 10)                          return false;
        if(is2017 && jet.neutralHadronEnergyFraction() <= 0.02)      return false;
    }

    return true;
}

bool JetAnalyzer::jetIsTight(const pat::Jet& jet, const bool is2017, const bool is2018) const{
    if(is2018){
      if(fabs(jet.eta()) <= 2.7){
        if(jet.neutralHadronEnergyFraction() >=  0.9)                         return false;
        if(jet.neutralEmEnergyFraction() >= 0.9)                              return false;
        if(jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1)          return false;
        if(jet.chargedHadronEnergyFraction() <= 0 and fabs(jet.eta()) <= 2.6) return false; // only for |eta|<2.6
        if(jet.chargedMultiplicity() <= 0)                                    return false;
      } else if(fabs(jet.eta()) <= 3.0){
        if(jet.neutralHadronEnergyFraction() >= 0.99)                         return false;
        if(jet.neutralEmEnergyFraction() <= 0.02)                             return false;
        if(jet.neutralMultiplicity() <= 2)                                    return false;
      } else {
        if(jet.neutralEmEnergyFraction() >= 0.90)                             return false;
        if(jet.neutralMultiplicity() <= 10)                                   return false;
        if(jet.neutralHadronEnergyFraction() <= 0.02)                         return false;
      }
      return true;
    } else {
      if(!jetIsLoose(jet, is2017))                                            return false;

      if(fabs(jet.eta()) <= 2.7){
          if(jet.neutralHadronEnergyFraction() >= 0.9)                        return false;
          if(jet.neutralEmEnergyFraction() >= 0.9)                            return false;

      } else if(fabs(jet.eta()) <= 3.0){
          if(is2017 && jet.neutralEmEnergyFraction() >= 0.99)                 return false;
      }
      return true;
    }
}

bool JetAnalyzer::jetIsTightLepVeto(const pat::Jet& jet, const bool is2017, const bool is2018) const{
    if(!jetIsTight(jet, is2017, is2018))                                                  return false; // Similar to tight ID except with additional cuts:
    if(fabs(jet.eta()) <= 2.7){
      if(jet.chargedMuEnergyFraction() >= 0.8 )                                           return false; // Muon energy fraction cut
      if(is2018 and jet.chargedEmEnergyFraction() >= 0.8)                                 return false; // EM fraction cut 2018
      else if(is2017 and fabs(jet.eta()) <= 2.4 and jet.chargedEmEnergyFraction() >= 0.8) return false; // EM fraction cut 2017
      else if(fabs(jet.eta()) <= 2.4 and jet.chargedEmEnergyFraction() >= 0.9)            return false; // EM fraction cut 2016
    }
    return true;
}
