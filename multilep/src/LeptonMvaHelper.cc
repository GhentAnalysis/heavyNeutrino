//implementation of LeptonMvaHelper class
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <cmath>

//Default constructor
//This will set up both MVA readers and book the correct variables
LeptonMvaHelper::LeptonMvaHelper(const edm::ParameterSet& iConfig, const unsigned type){ //0 : SUSY , 1: ttH, 2: tZq/TTV
    for(unsigned i = 0; i < 2; ++i){
        //Set up Mva reader 
        reader[i] = std::make_shared<TMVA::Reader>( "!Color:!Silent");
    }
    if(type < 2){
        for(unsigned i = 0; i < 2; ++i){
            //Book Common variables
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

        //Read Mva weights
        if(type == 0){ //SUSY weights used by default
            //Book specific electron variables
            reader[1]->AddVariable("LepGood_mvaIdSpring16GP", &LepGood_mvaIdSpring16GP);
            reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMu").fullPath());
            reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEle").fullPath());
        } else{
            //Book specific electron variables
            reader[1]->AddVariable("LepGood_mvaIdSpring16HZZ", &LepGood_mvaIdSpring16HZZ);
            reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMuttH").fullPath());
            reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsElettH").fullPath());
        }
    } else{
        for(unsigned i = 0; i < 2; ++i){
            reader[i]->AddVariable( "pt", &LepGood_pt );
            reader[i]->AddVariable( "eta", &LepGood_eta );
            reader[i]->AddVariable( "trackMult", &LepGood_jetNDauChargedMVASel );
            reader[i]->AddVariable( "miniIsoCharged", &LepGood_miniRelIsoCharged );
            reader[i]->AddVariable( "miniIsoNeutral", &LepGood_miniRelIsoNeutral );
            reader[i]->AddVariable( "ptrel", &LepGood_jetPtRelv2 );
            reader[i]->AddVariable( "min(ptratio,1.5)", &LepGood_jetPtRatio );
            reader[i]->AddVariable( "relIso0p3", &LepGood_relIso0p3); 
            reader[i]->AddVariable( "max(jetbtagCSV,0)", &LepGood_jetBTagCSV );
            reader[i]->AddVariable( "sip3d", &LepGood_sip3d );
        }
        reader[0]->AddVariable("segmComp", &LepGood_segmentCompatibility);
        reader[1]->AddVariable("eleMVA", &LepGood_mvaIdSpring16GP);
        reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMutZqTTV").fullPath());
        reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEletZqTTV").fullPath());
    }
}
void LeptonMvaHelper::bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetCsv, double sip3d, double dxy, double dz, double relIso)
{
    LepGood_pt = pt;
    LepGood_eta = eta;
    LepGood_jetNDauChargedMVASel = selectedTrackMult;
    LepGood_miniRelIsoCharged = miniIsoCharged;
    LepGood_miniRelIsoNeutral = miniIsoNeutral;
    LepGood_jetPtRelv2 = ptRel;
    LepGood_jetPtRatio = std::min(ptRatio, 1.5);
    LepGood_jetBTagCSV = std::max(closestJetCsv, 0.);
    LepGood_sip3d = sip3d;
    LepGood_dxy = log(fabs(dxy));
    LepGood_dz = log(fabs(dz));
    LepGood_relIso0p3 = relIso;
}

double LeptonMvaHelper::leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetCsv, double sip3d, double dxy, double dz, double relIso, double segComp)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, closestJetCsv, sip3d, dxy, dz, relIso);
    LepGood_segmentCompatibility = segComp;
    return reader[0]->EvaluateMVA("BDTG method");
}

double LeptonMvaHelper::leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetCsv, double sip3d, double dxy, double dz, double relIso, double eleMva, double eleMvaHZZ)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, closestJetCsv, sip3d, dxy, dz, relIso);
    LepGood_mvaIdSpring16GP = eleMva;
    LepGood_mvaIdSpring16HZZ = eleMvaHZZ;
    return reader[1]->EvaluateMVA("BDTG method");
}
