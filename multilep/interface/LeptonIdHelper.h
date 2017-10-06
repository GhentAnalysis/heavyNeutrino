#ifndef LeptonIdHelper_H
#define LeptonIdHelper_H

//Class designed to compute analysis specific ID decisions and store them
class LeptonIdHelper{
    public:
        

    private:
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
}
#endif
