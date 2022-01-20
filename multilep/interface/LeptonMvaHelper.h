//Helper class for facilitating the computation of the lepton MVA.
#ifndef Lepton_Mva_Helper
#define Lepton_Mva_Helper

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TMVA/Reader.h"
#include <xgboost/c_api.h>
#include <memory>
class LeptonMvaHelper{
    public:
        LeptonMvaHelper(const edm::ParameterSet& iConfig, const std::string tagger, const int year);
        double leptonMvaMuon(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double miniIsoNeutralVanilla, double ptRel, double ptRatio, double ptRatioVanilla, double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p3DB, double relIso0p3Vanilla, double relIso0p3DBVanilla, double segComp);
        double leptonMvaElectron(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double miniIsoNeutralVanilla, double ptRel, double ptRatio, double ptRatioVanilla, double closestJetDeepCsv, double closesJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p3Vanilla, double eleMvaSummer16, double eleMvaFall17v1, double eleMvaFall17v2, double eleMissingHits);
    private:
        std::string tagger;
        int year;
        std::shared_ptr<TMVA::Reader> reader[2]; //First entry is for muons, second one for electrons
        float LepGood_pt,                       //Variables used in MVA computation
        LepGood_eta,
        LepGood_jetNDauChargedMVASel,
        LepGood_miniRelIsoCharged,
        LepGood_miniRelIsoNeutral,
        LepGood_miniRelIsoNeutralVanilla,
        LepGood_jetPtRelv2,
        LepGood_jetPtRatio,
        LepGood_jetPtRatioVanilla,
        LepGood_jetBTag,
        LepGood_sip3d,
        LepGood_dxy,
        LepGood_dz,
        LepGood_segmentCompatibility,
        LepGood_mvaIdSummer16GP,
        LepGood_mvaIdFall17v1noIso,
        LepGood_mvaIdFall17v2noIso,
        LepGood_relIso0p3,
        LepGood_relIso0p3Vanilla;
        void bookCommonVars(double pt, double eta, double selectedTrackMult, double miniIsoCharged, double miniIsoNeutral, double miniIsoNeutralVanilla, double ptRel, double ptRatio, double ptRatioVanilla, double closestJetDeepCsv, double closestJetDeepFlavor, double sip3d, double dxy, double dz, double relIso0p3, double relIso0p3Vanilla);
        BoosterHandle booster[2];
        float boosterVars[2][1][15];
};
#endif
