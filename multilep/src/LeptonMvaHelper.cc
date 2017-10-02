//implementation of LeptonMvaHelper class
#include "heavyNeutrino/multilepton/interface/LeptonMvaHelper.h"

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
