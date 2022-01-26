//implementation of LeptonMvaHelper class
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <cmath>

// TODO: clean-up of this class, maybe get rid of older trainings

//Default constructor
//This will set up both MVA readers and book the correct variables
LeptonMvaHelper::LeptonMvaHelper(const edm::ParameterSet& iConfig, const std::string taggerName, const int yearTrain):
    tagger( taggerName ), year( yearTrain )
{
    for(unsigned i = 0; i < 2; ++i){

        //Set up Mva reader 
        reader[i] = std::make_shared<TMVA::Reader>( "!Color:!Silent");
    }

    //TTH Lepton MVA
    if( tagger == "TTH" ){
        for(unsigned i = 0; i < 2; ++i){

            //Book Common variables
            reader[i]->AddVariable( "LepGood_pt", &LepGood_pt );
            reader[i]->AddVariable( "LepGood_eta", &LepGood_eta );
            reader[i]->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
            reader[i]->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
            reader[i]->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
            reader[i]->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
            reader[i]->AddVariable( "LepGood_jetDF", &LepGood_jetBTag );
            reader[i]->AddVariable( "LepGood_jetPtRatio", &LepGood_jetPtRatio );
            reader[i]->AddVariable( "LepGood_dxy", &LepGood_dxy );
            reader[i]->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
            reader[i]->AddVariable( "LepGood_dz", &LepGood_dz );
        }

        //Book specific muon variables
        reader[0]->AddVariable( "LepGood_segmentComp", &LepGood_segmentCompatibility );
        reader[1]->AddVariable( "LepGood_mvaFall17V2noIso", &LepGood_mvaIdFall17v2noIso );

        reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMuttH").fullPath());
        reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsElettH").fullPath());

    //tZq lepton MVA
    } else if( tagger == "TZQ" ) {
        for(unsigned i = 0; i < 2; ++i){
            reader[i]->AddVariable( "pt", &LepGood_pt );
            reader[i]->AddVariable( "eta", &LepGood_eta );
            reader[i]->AddVariable( "trackMultClosestJet", &LepGood_jetNDauChargedMVASel );
            reader[i]->AddVariable( "miniIsoCharged", &LepGood_miniRelIsoCharged );
            reader[i]->AddVariable( "miniIsoNeutral", &LepGood_miniRelIsoNeutral );
            reader[i]->AddVariable( "pTRel", &LepGood_jetPtRelv2 );
            reader[i]->AddVariable( "ptRatio", &LepGood_jetPtRatio );
            reader[i]->AddVariable( "relIso", &LepGood_relIso0p3); 
            reader[i]->AddVariable( "deepCsvClosestJet", &LepGood_jetBTag );
            reader[i]->AddVariable( "sip3d", &LepGood_sip3d );
            reader[i]->AddVariable( "dxy", &LepGood_dxy);
            reader[i]->AddVariable( "dz", &LepGood_dz);
        }
        reader[0]->AddVariable("segmentCompatibility", &LepGood_segmentCompatibility);
        if(  year == 2016  ){
            reader[1]->AddVariable("electronMvaSpring16GP", &LepGood_mvaIdSummer16GP);
        } else{
            reader[1]->AddVariable("electronMvaFall17NoIso", &LepGood_mvaIdFall17v1noIso);
        }
        reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMutZqTTV").fullPath());
        reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEletZqTTV").fullPath());
       
    //TOP lepton MVA
    } else if( tagger == "TOP" ) {
        for(unsigned i = 0; i < 2; ++i){
	    reader[i]->AddVariable( "dxylog", &LepGood_dxy);
	    reader[i]->AddVariable( "miniIsoCharged", &LepGood_miniRelIsoCharged );
	    reader[i]->AddVariable( "miniIsoNeutral", &LepGood_miniRelIsoNeutral );
	    reader[i]->AddVariable( "pTRel", &LepGood_jetPtRelv2 );
	    reader[i]->AddVariable( "sip3d", &LepGood_sip3d );
	}
        reader[0]->AddVariable("segmentCompatibility", &LepGood_segmentCompatibility);
        reader[1]->AddVariable("mvaIdFall17v2noIso", &LepGood_mvaIdFall17v2noIso);
        for(unsigned i = 0; i < 2; ++i){
	    reader[i]->AddVariable( "ptRatio", &LepGood_jetPtRatio );
	    reader[i]->AddVariable( "bTagDeepJetClosestJet", &LepGood_jetBTag );
            reader[i]->AddVariable( "pt", &LepGood_pt );
	    reader[i]->AddVariable( "trackMultClosestJet", &LepGood_jetNDauChargedMVASel );
            reader[i]->AddVariable( "etaAbs", &LepGood_eta );
	    reader[i]->AddVariable( "dzlog", &LepGood_dz);
            reader[i]->AddVariable( "relIso", &LepGood_relIso0p3); 
	}       
        reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMuTOP").fullPath());
        reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEleTOP").fullPath());
    }
   //TOP-UL lepton MVA
   else if( tagger == "TOP-UL" ) {
      XGBoosterCreate(NULL, 0, &booster[0]);
      XGBoosterCreate(NULL, 0, &booster[1]);
      XGBoosterLoadModel(booster[0], iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMuTOPUL").fullPath().c_str());
      XGBoosterLoadModel(booster[1], iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEleTOPUL").fullPath().c_str());
   }
   //TOP-UL (v2) lepton MVA
   else if( tagger == "TOPv2-UL" ) {
      XGBoosterCreate(NULL, 0, &booster[0]);
      XGBoosterCreate(NULL, 0, &booster[1]);
      XGBoosterLoadModel(booster[0], iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMuTOPv2UL").fullPath().c_str());
      XGBoosterLoadModel(booster[1], iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEleTOPv2UL").fullPath().c_str());
   }   	
}


void LeptonMvaHelper::bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double miniIsoNeutralVanilla, double ptRel, double ptRatio, double ptRatioVanilla, double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p3Vanilla)
{
    LepGood_pt = pt;
    if( tagger == "TTH" || tagger == "TOP-UL" || tagger == "TOPv2-UL" ) {
        LepGood_eta = eta;
    } else {
        LepGood_eta = fabs(eta);
    }
    LepGood_jetNDauChargedMVASel = selectedTrackMult;
    LepGood_miniRelIsoCharged = miniIsoCharged;
    LepGood_miniRelIsoNeutral = miniIsoNeutral;
    LepGood_miniRelIsoNeutralVanilla = miniIsoNeutralVanilla;
    LepGood_jetPtRelv2 = ptRel;
    LepGood_jetPtRatio = std::min(ptRatio, 1.5);
    LepGood_jetPtRatioVanilla = std::min(ptRatioVanilla, 1.5);
    if( tagger != "TZQ" ){
        LepGood_jetBTag = std::max( ( std::isnan( closestJetDeepFlavor ) ? 0. : closestJetDeepFlavor ), 0. );
    } else {
        LepGood_jetBTag = std::max( ( std::isnan( closestJetDeepCsv ) ? 0. : closestJetDeepCsv ), 0. );
    }
    LepGood_sip3d = sip3d;
    LepGood_dxy = log(fabs(dxy));
    LepGood_dz = log(fabs(dz));
    LepGood_relIso0p3 = relIso0p3;
    LepGood_relIso0p3Vanilla = relIso0p3Vanilla;
   
    //TOP-UL lepton MVA
    if( tagger == "TOP-UL" || tagger == "TOPv2-UL" ) 
     {	
	for(unsigned i = 0; i < 2; ++i)
	  {	     
	     boosterVars[i][0][0] = LepGood_pt;
	     boosterVars[i][0][1] = LepGood_eta;
	     boosterVars[i][0][2] = LepGood_jetNDauChargedMVASel;
	     boosterVars[i][0][3] = LepGood_miniRelIsoCharged;
	     boosterVars[i][0][4] = LepGood_miniRelIsoNeutralVanilla;
	     boosterVars[i][0][5] = LepGood_jetPtRelv2;
	     boosterVars[i][0][6] = LepGood_jetPtRatioVanilla;
	     boosterVars[i][0][7] = LepGood_relIso0p3Vanilla;
	     boosterVars[i][0][8] = LepGood_jetBTag;
	     boosterVars[i][0][9] = LepGood_sip3d;
	     boosterVars[i][0][10] = LepGood_dxy;
	     boosterVars[i][0][11] = LepGood_dz;
	  }	
     }      
}

double LeptonMvaHelper::leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double miniIsoNeutralVanilla, double ptRel, double ptRatio, double ptRatioVanilla, double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p3DB, double relIso0p3Vanilla, double relIso0p3DBVanilla, double segComp)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, miniIsoNeutralVanilla, ptRel, ptRatio, ptRatioVanilla, closestJetDeepCsv, closestJetDeepFlavor, sip3d, dxy, dz, relIso0p3, relIso0p3Vanilla);
    if( tagger == "TOP" ) LepGood_relIso0p3 = relIso0p3DB;
    boosterVars[0][0][7] = relIso0p3DBVanilla;
    boosterVars[0][0][12] = segComp;
    LepGood_segmentCompatibility = segComp;
    if( tagger == "TOP-UL" || tagger == "TOPv2-UL" )
     {	
	DMatrixHandle dtest;
	int nfeat = 13;
	XGDMatrixCreateFromMat(reinterpret_cast<float*>(boosterVars[0]), 1, nfeat, NAN, &dtest);
	bst_ulong out_len;
	const float *f;
	XGBoosterPredict(booster[0], dtest, 0, 0, &out_len, &f);
	XGDMatrixFree(dtest);
	return f[0];
     } else return reader[0]->EvaluateMVA("BDTG method");
}

double LeptonMvaHelper::leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double miniIsoNeutralVanilla, double ptRel, double ptRatio, double ptRatioVanilla, double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p3Vanilla, double eleMvaSummer16, double eleMvaFall17v1, double eleMvaFall17v2, double eleMissingHits)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, miniIsoNeutralVanilla, ptRel, ptRatio, ptRatioVanilla, closestJetDeepCsv, closestJetDeepFlavor, sip3d, dxy, dz, relIso0p3, relIso0p3Vanilla);
    LepGood_mvaIdSummer16GP = eleMvaSummer16;
    LepGood_mvaIdFall17v1noIso = eleMvaFall17v1;
    LepGood_mvaIdFall17v2noIso = eleMvaFall17v2;
    boosterVars[1][0][12] = eleMvaFall17v2;
    boosterVars[1][0][13] = eleMissingHits;
    if( tagger == "TOP-UL" || tagger == "TOPv2-UL" )
     {	
	DMatrixHandle dtest;
	int nfeat = 13;
	if( tagger == "TOPv2-UL" ) nfeat = 14;
	XGDMatrixCreateFromMat(reinterpret_cast<float*>(boosterVars[1]), 1, nfeat, NAN, &dtest);
	bst_ulong out_len;
	const float *f;
	XGBoosterPredict(booster[1], dtest, 0, 0, &out_len, &f);
	XGDMatrixFree(dtest);
	return f[0];
     } else return reader[1]->EvaluateMVA("BDTG method");
}
