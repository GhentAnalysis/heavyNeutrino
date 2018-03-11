/*
class to apply JEC from txt files
*/

#ifndef JEC_H
#define JEC_H

//include c++ library classes
#include <map>
#include <string>
#include <vector>
#include <memory>

//include CMSSW classes
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class JEC {
    public:
        JEC(const std::string& JECpath, const bool dataSample, const bool fall17Sample);
        ~JEC();

        void updateJEC(const unsigned long);
        void setJEC(const std::string&);

        double jetCorrection(double rawPt, double eta, double rho, double area, const std::string& level = "L3Absolute");     // this function returns, for a given jet the correction factor
        double jetUncertainty(double pt, double eta);
        std::pair <double, double > correctedMETAndPhi(const pat::MET& met, const std::vector< pat::Jet >& jets, const double rho);

    private:
        std::string path;
        bool isData;
        bool is2017;
        std::string currentJEC; 

        std::shared_ptr<FactorizedJetCorrector> jetCorrector;
        std::shared_ptr<JetCorrectionUncertainty> jetUncertainties;

        double px(double pt, double phi){ return pt*cos(phi); };
        double py(double pt, double phi){ return pt*sin(phi); };
        std::vector<float> getSubCorrections(double rawPt, double eta, double rho, double area);
        std::pair<double, double> getMETCorrectionPxPy(double rawPt, double rawEta, double rawMuonSubtractedPt, double phi, double emf, double rho, double area);

        //JEC naming 
        std::string getJECRunName(const unsigned long);
        std::string getJECName(const unsigned long);
        
};
#endif
