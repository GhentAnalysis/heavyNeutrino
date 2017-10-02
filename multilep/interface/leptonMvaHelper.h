//Helper class for facilitating the computation of the lepton MVA.
#ifndef lepton_mva_helper
#define lepton_mva_helper

#include "TMVA/Reader.h"
#include <memory>
class leptonMvaHelper{
    public:
        leptonMvaHelper();
        double leptonMvaMuon(const double segComp);
        double leptonMVaElectron(const double eleMva);
    private:
        std::shared_ptr<TMVA:Reader> reader[2]; //First entry is for muons, second one for electrons 
        float LepGood_pt,                       //Variables used in MVA computation
        LepGood_eta,
        LepGood_jetNDauChargedMVASel,
        LepGood_miniRelIsoCharged,
        LepGood_miniRelIsoNeutral,
        LepGood_jetPtRelv2,
        LepGood_jetPtRatio,
        LepGood_jetBTagCSV,
        LepGood_sip3d,
        LepGood_dxy,
        LepGood_dz,
        LepGood_segmentCompatibility,
        LepGood_mvaIdSpring16GP;
        void bookCommonVars(double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz); 
};
#endif
