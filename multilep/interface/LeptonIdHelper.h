#ifndef LeptonIdHelper_H
#define LeptonIdHelper_H

#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
class LeptonAnalyzer;
//Class designed to compute analysis specific ID decisions and store them
class LeptonIdHelper{
    public:
        LeptonIdHelper();
        LeptonIdHelper(const LeptonAnalyzer& lepAn, const unsigned index, const unsigned flav);

        void setVars(const LeptonAnalyzer& lepAn, const unsigned index, const unsigned flav);

        bool ewkLoose();           //return ewkino ID decisions
        bool ewkFO();
        bool ewkTight();

        bool hnlLoose();           //return HNL ID decisions
        bool hnlFO();
        bool hnlTight();
    private:
        unsigned flavor;            //lepton Flavor

        double pt;                  //lepton kinematics
        double eta;

        double miniIsoCharged;      //isolation variables
        double miniIso;
        double relIso;
        
        double ptRel;               //closest jet variables
        double ptRatio;
        double closestJetCsv;
        unsigned selectedTrackMult;

        double sip3d;               //Pointing variables 
        double dxy;
        double dz;

        bool mediumMuon;            //Lepton specific id variables
        double eleMva;
        double eleTrigEmu;
        bool eleMuOverlap;
};
#endif
