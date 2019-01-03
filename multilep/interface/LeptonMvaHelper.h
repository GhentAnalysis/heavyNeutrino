//Helper class for facilitating the computation of the lepton MVA.
#ifndef Lepton_Mva_Helper
#define Lepton_Mva_Helper

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TMVA/Reader.h"
#include <memory>
class LeptonMvaHelper{
    public:
        LeptonMvaHelper(const edm::ParameterSet& iConfig, const unsigned type, const bool sampleIs2017);
        double leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double closestJetDeepCsv, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p4, double segComp);
        double leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double closesJetDeepCsv, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p4, double eleMvaSpring16, double eleMvaHZZ, double eleMvaFall17);
    private:
        unsigned type; //0 = SUSY , 1 = ttH , 2 = tZqttV
        bool is2017;
        bool is2018;
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
        LepGood_mvaIdFall17noIso,
        LepGood_mvaIdSpring16HZZ,
        LepGood_relIso0p3,
        LepGood_relIso0p4;
        void bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double ptRel, double ptRatio, double closestJetCsv, double closestJetDeepCsv, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p4); 
};
#endif
