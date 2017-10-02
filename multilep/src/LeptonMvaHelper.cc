//implementation of LeptonMvaHelper class
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"

//Default constructor
//This will set up both MVA readers and book the correct variables
LeptonMvaHelper::LeptonMvaHelper(){
    //Book common variables
    for(unsigned i = 0; i < 2; ++i){
        reader[i]->AddVariable( "LepGood_pt", &LepGood_pt );
        reader[i]->AddVariable( "LepGood_eta", &LepGood_eta );
        reader[i]->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
        reader[i]->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
        reader[i]->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
        reader[i]->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
        reader[i]->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
        reader[i]->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
        reader[i]->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
        reader[i]->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
        reader[i]->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
    }
    //Book specific muon variables
    reader[0]->AddVariable("LepGood_segmentCompatibility", &LepGood_segmentCompatibility);
    //Book specific electron variables
    reader[1]->AddVariable("LepGood_mvaIdSpring16GP", &LepGood_mvaIdSpring16GP);

    //Read Mva weights
    reader[0]->BookMVA("BDTG method", "heavyNeutrino/multilep/data/mvaWeights/mu_BDTG.weights.xml");
    reader[1]->BookMVA("BDTG method", "heavyNeutrino/mutlitlep/data/mvaWeights/el_BDTG.weights.xml");
}

void leptonMvaHelper::bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz){
    LepGood_pt = pt;
    LepGood_eta = eta;
    LepGood_jetNDauChargedMVASel = selectedTrackMult;
    LepGood_miniRelIsoCharged = miniIsoCharged;
    LepGood_miniRelIsoNeutral = miniIsoNeutral;
    LepGood_jetPtRelv2 = ptRel;
    LepGood_jetPtRatio = ptRatio;
    LepGood_jetBTagCSV = closestJetCsv;
    LepGood_sip3d = sip3d;
    LepGood_dxy = dxy;
    LepGood_dz = dz;
}

double leptonMvaHelper::leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz, const double segComp){
    bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz);
    LepGood_segmentCompatibility = segComp;
    return reader[0]->EvaluateMVA("BDTG method");
}

double leptonMvaHelper::leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz, const double eleMva){
    bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz);
    LepGood_mvaIdSpring16GP = eleMva;
    return reader[1]->EvaluateMVA("BDTG method");
}
