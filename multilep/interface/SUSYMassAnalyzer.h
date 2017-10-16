#ifndef SUSYMassAnalyzer_H
#define SUSYMassAnalyzer_H
//Class designed to incorportate masses of SMS scans into SUSY signal trees
#include "heavyNeutrino/multilep/plugins/multilep.h"
#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"
#include "TH2D.h"



class multilep;
class LheAnalyzer;

class SUSYMassAnalyzer{
    private:
        TH2D* hCounterSUSY;        //Event counter for every mass point in the scan sample
        multilep* multilepAnalyzer;
        LheAnalyzer* lheAnalyzer;  //Needed to get event weight
        
        //EWKino masses
        double _mChi1;              //LSP (neutralino1) mass
        double _mChi2;              //Chargino mass = neutralino2 mass
    public:
        SUSYMassAnalyzer(const edm::ParameterSet&, multilep*, LheAnalyzer*);
        ~SUSYMassAnalyzer(){};
        void beginJob(TTree* outputTree, edm::Service<TFileService>&);
        void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);
        void analyze(const edm::Event&);
};
#endif
