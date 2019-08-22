//implementation of LeptonMvaHelper class
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <cmath>

// TODO: clean-up of this class, maybe get rid of older trainings

//Default constructor
//This will set up both MVA readers and book the correct variables
LeptonMvaHelper::LeptonMvaHelper(const edm::ParameterSet& iConfig, const bool isTTHMVA, const bool sampleIs2017or2018): //0 : ttH, 1: tZq/TTV
    isTTH( isTTHMVA ), is2017Or2018( sampleIs2017or2018 )
{
    for(unsigned i = 0; i < 2; ++i){

        //Set up Mva reader 
        reader[i] = std::make_shared<TMVA::Reader>( "!Color:!Silent");
    }

    //TTH Lepton MVA
    if( isTTH ){
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
    } else {
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
        if(  !is2017Or2018  ){
            reader[1]->AddVariable("electronMvaSpring16GP", &LepGood_mvaIdSummer16GP);
        } else{
            reader[1]->AddVariable("electronMvaFall17NoIso", &LepGood_mvaIdFall17v1noIso);
        }
        reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMutZqTTV").fullPath());
        reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEletZqTTV").fullPath());
    }
}


void LeptonMvaHelper::bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
        double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3 )
{
    LepGood_pt = pt;
    if( isTTH ){
        LepGood_eta = eta;
    } else{
        LepGood_eta = fabs(eta);
    }
    LepGood_jetNDauChargedMVASel = selectedTrackMult;
    LepGood_miniRelIsoCharged = miniIsoCharged;
    LepGood_miniRelIsoNeutral = miniIsoNeutral;
    LepGood_jetPtRelv2 = ptRel;
    LepGood_jetPtRatio = std::min(ptRatio, 1.5);
    if( isTTH ){
        LepGood_jetBTag = std::max( ( std::isnan( closestJetDeepFlavor ) ? 0. : closestJetDeepFlavor ), 0. );
    } else {
        LepGood_jetBTag = std::max( ( std::isnan( closestJetDeepCsv ) ? 0. : closestJetDeepCsv ), 0. );
    }
    LepGood_sip3d = sip3d;
    LepGood_dxy = log(fabs(dxy));
    LepGood_dz = log(fabs(dz));
    LepGood_relIso0p3 = relIso0p3;
}

double LeptonMvaHelper::leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double segComp)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, closestJetDeepCsv, closestJetDeepFlavor, sip3d, dxy, dz, relIso0p3);
    LepGood_segmentCompatibility = segComp;
    return reader[0]->EvaluateMVA("BDTG method");
}

double LeptonMvaHelper::leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double eleMvaSummer16, double eleMvaFall17v1, double eleMvaFall17v2)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, closestJetDeepCsv, closestJetDeepFlavor, sip3d, dxy, dz, relIso0p3);
    LepGood_mvaIdSummer16GP = eleMvaSummer16;
    LepGood_mvaIdFall17v1noIso = eleMvaFall17v1;
    LepGood_mvaIdFall17v2noIso = eleMvaFall17v2;
    return reader[1]->EvaluateMVA("BDTG method");
}
