//Helper class for facilitating the computation of the lepton MVA.
#ifndef Lepton_Mva_Helper
#define Lepton_Mva_Helper

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TMVA/Reader.h"
#include <memory>
class LeptonMvaHelper{
    public:
        LeptonMvaHelper(const edm::ParameterSet& iConfig, const unsigned type);
        double leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz, double relIso, double segComp);
        double leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz, double relIso, double eleMva, double eleMvaHZZ);
    private:
        std::shared_ptr<TMVA::Reader> reader[2]; //First entry is for muons, second one for electrons 
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
        LepGood_mvaIdSpring16GP,
        LepGood_mvaIdSpring16HZZ,
        LepGood_relIso0p3;
        void bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double sip3d, double dxy, double dz, double relIso); 
};
#endif
