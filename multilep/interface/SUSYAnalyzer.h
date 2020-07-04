#ifndef SUSYAnalyzer_H
#define SUSYAnalyzer_H
//Class designed to incorportate masses of SMS scans into SUSY signal trees
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"
#include "TH2D.h"



class multilep;
class LheAnalyzer;

class SUSYAnalyzer{
    private:
        TH2D* hCounterSUSY;             //Event counter for every mass point in the scan sample
        TH2D* hCounterSUSYISRWeighted;  //Event counter containing the ISR weights, this is needed because the total XSection after ISR reweighting should not change
        multilep* multilepAnalyzer;
        LheAnalyzer* lheAnalyzer;  //Needed to get event weight
        
        //EWKino masses
        double _mChi1;              //LSP (neutralino1) mass
        double _mChi2;              //Chargino mass = neutralino2 mass
        double _ptISR;              //ISR pT = pT of the initial susy particle pair
        double _ewkISRWeight;       //ISR weight for given SUSY event
    public:
        SUSYAnalyzer(const edm::ParameterSet&, multilep*, LheAnalyzer*);
        ~SUSYAnalyzer(){};
        void beginJob(TTree* outputTree, edm::Service<TFileService>&);
        void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);
        void analyze(const edm::Event&);
};
#endif
