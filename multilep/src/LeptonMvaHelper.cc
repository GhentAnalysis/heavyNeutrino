//implementation of LeptonMvaHelper class
#include "heavyNeutrino/multilep/interface/LeptonMvaHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <cmath>

//Default constructor
//This will set up both MVA readers and book the correct variables
LeptonMvaHelper::LeptonMvaHelper(const edm::ParameterSet& iConfig, const unsigned typeNumber, const bool sampleIs2017): //0 : SUSY , 1: ttH, 2: tZq/TTV
    type(typeNumber), is2017(sampleIs2017)
{
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
            if(!is2017){
                reader[i]->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
                reader[i]->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
            } else{
                reader[i]->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
                reader[i]->AddVariable( "(LepGood_jetBTagCSV>-5)*min(LepGood_jetPtRatiov2,1.5)+(LepGood_jetBTagCSV<-5)/(1+LepGood_relIso04)", &LepGood_jetPtRatio);
            }
            reader[i]->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
            reader[i]->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
            reader[i]->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
        }

        //Book specific muon variables
        reader[0]->AddVariable("LepGood_segmentCompatibility", &LepGood_segmentCompatibility);

        if(!is2017){
            //Read Mva weights
            if(type == 0){ //SUSY weights used by default
                //Book specific electron variables
                reader[1]->AddVariable("LepGood_mvaIdSpring16GP", &LepGood_mvaIdSpring16GP);
            } else{
                //Book specific electron variables
                reader[1]->AddVariable("LepGood_mvaIdSpring16HZZ", &LepGood_mvaIdSpring16HZZ);
            }
        } else {
            reader[1]->AddVariable("LepGood_mvaIdFall17noIso", &LepGood_mvaIdFall17noIso);
        }
        if(type == 0){
            reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>( std::string("leptonMvaWeightsMuSUSY") + (is2017 ? "17" : "16") ).fullPath());
            reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>( std::string("leptonMvaWeightsEleSUSY") + (is2017 ? "17" : "16") ).fullPath());
        } else{
            reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>( std::string("leptonMvaWeightsMuttH") + (is2017 ? "17" : "16") ).fullPath());
            reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>( std::string("leptonMvaWeightsElettH") + (is2017 ? "17" : "16") ).fullPath());
        }
    } else{
        for(unsigned i = 0; i < 2; ++i){
            reader[i]->AddVariable( "pt", &LepGood_pt );
            reader[i]->AddVariable( "eta", &LepGood_eta );
            reader[i]->AddVariable( "trackMultClosestJet", &LepGood_jetNDauChargedMVASel );
            reader[i]->AddVariable( "miniIsoCharged", &LepGood_miniRelIsoCharged );
            reader[i]->AddVariable( "miniIsoNeutral", &LepGood_miniRelIsoNeutral );
            reader[i]->AddVariable( "pTRel", &LepGood_jetPtRelv2 );
            reader[i]->AddVariable( "ptRatio", &LepGood_jetPtRatio );
            reader[i]->AddVariable( "relIso", &LepGood_relIso0p3); 
            reader[i]->AddVariable( "deepCsvClosestJet", &LepGood_jetBTagCSV );
            reader[i]->AddVariable( "sip3d", &LepGood_sip3d );
            reader[i]->AddVariable( "dxy", &LepGood_dxy);
            reader[i]->AddVariable( "dz", &LepGood_dz);
        }
        reader[0]->AddVariable("segmentCompatibility", &LepGood_segmentCompatibility);
        if(!is2017){
            reader[1]->AddVariable("electronMvaSpring16GP", &LepGood_mvaIdSpring16GP);
        } else{
            reader[1]->AddVariable("electronMvaFall17NoIso", &LepGood_mvaIdFall17noIso);
        }
        if(!is2017){
            reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMutZqTTV16").fullPath());
            reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEletZqTTV16").fullPath());
        } else{
            reader[0]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsMutZqTTV17").fullPath());
            reader[1]->BookMVA("BDTG method", iConfig.getParameter<edm::FileInPath>("leptonMvaWeightsEletZqTTV17").fullPath());
        }
    }
}
void LeptonMvaHelper::bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetCsv, double closestJetDeepCsv, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p4)
{
    LepGood_pt = pt;
    if(type >= 2){
        LepGood_eta = fabs(eta);
    } else{
        LepGood_eta = eta;
    }
    LepGood_jetNDauChargedMVASel = selectedTrackMult;
    LepGood_miniRelIsoCharged = miniIsoCharged;
    LepGood_miniRelIsoNeutral = miniIsoNeutral;
    LepGood_jetPtRelv2 = ptRel;
    LepGood_jetPtRatio = std::min(ptRatio, 1.5);
    if(is2017 || type  >= 2){
        LepGood_jetBTagCSV = std::max( (std::isnan(closestJetDeepCsv) ? 0. : closestJetDeepCsv) , 0.);
    } else{
        LepGood_jetBTagCSV = std::max(closestJetCsv, 0.);
    }
    LepGood_relIso0p4 = relIso0p4;
    //use relIso for closest jet when no close jet for 2017 SUSY and ttH mvas
    if(is2017 && type < 2){
        LepGood_jetPtRatio = 1 + relIso0p4;
        bool goodBTag = (closestJetDeepCsv > -5.) || !std::isnan(closestJetDeepCsv);
        LepGood_jetPtRatio = goodBTag*std::min(ptRatio, 1.5) + (!goodBTag)/(1 + relIso0p4);
    }
    LepGood_sip3d = sip3d;
    LepGood_dxy = log(fabs(dxy));
    LepGood_dz = log(fabs(dz));
    LepGood_relIso0p3 = relIso0p3;
}

double LeptonMvaHelper::leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetCsv, double closestJetDeepCsv, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p4, double segComp)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, closestJetCsv, closestJetDeepCsv, sip3d, dxy, dz, relIso0p3, relIso0p4);
    LepGood_segmentCompatibility = segComp;
    return reader[0]->EvaluateMVA("BDTG method");
}

double LeptonMvaHelper::leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, 
    double closestJetCsv, double closestJetDeepCsv, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p4, double eleMvaSpring16, double eleMvaHZZ, double eleMvaFall17)
{
    bookCommonVars(pt, eta, selectedTrackMult, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, closestJetCsv, closestJetDeepCsv, sip3d, dxy, dz, relIso0p3, relIso0p4);
    LepGood_mvaIdSpring16GP = eleMvaSpring16;
    LepGood_mvaIdSpring16HZZ = eleMvaHZZ;
    LepGood_mvaIdFall17noIso = eleMvaFall17;
    return reader[1]->EvaluateMVA("BDTG method");
}
