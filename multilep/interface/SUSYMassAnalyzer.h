//Class designed to incorportate masses of SMS scans into SUSY signal trees
#include "multilep.h"


class multilep;

class SUSYMassAnalyzer{
    private:
        TH1D*  hCounterSUSY;        //Event counter for every mass point in the scan sample
        multilep* multilepAnalyzer;
        
        //EWKino masses
        double _mChi0;              //LSP (neutralino1) mass
        double _mChi1;              //Chargino mass = neutralino2 mass
    public:
        SUSYMassAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
        ~SUSYMassAnalyzer();
        void beginJob(TTree* outputTree, edm::Service<TFileService>& fs);
        void analyze(const edm::Event&);
};
