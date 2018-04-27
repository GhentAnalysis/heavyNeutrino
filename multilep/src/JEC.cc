/*
   implementation of JEC class
   */

#include "heavyNeutrino/multilep/interface/JEC.h"

JEC::JEC(const std::string& JECPath, bool dataSample, bool fall17Sample):
    path(JECPath), isData(dataSample), is2017(fall17Sample) , currentJEC("none") {}

JEC::~JEC(){}

void JEC::updateJEC(const unsigned long runNumber){
    if(getJECName(runNumber) != currentJEC){
        currentJEC = getJECName(runNumber);
        setJEC(getJECName(runNumber) );
    }
}

void JEC::setJEC(const std::string& JECName){
    std::vector< JetCorrectorParameters > JECParameters;
    JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L1FastJet_AK4PFchs.txt") );
    JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L2Relative_AK4PFchs.txt") );
    JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L3Absolute_AK4PFchs.txt") );
    if(isData) JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L2L3Residual_AK4PFchs.txt") );
    jetCorrector.reset(new FactorizedJetCorrector(JECParameters) );
    jetUncertainties.reset(new JetCorrectionUncertainty(path + JECName + "_Uncertainty_AK4PFchs.txt") );
}
   
std::string JEC::getJECRunName(const unsigned long runNumber){
    ///////////////////////////
    static const std::string version2016 = "V4";
    static const std::string version2017 = "_V6";
    //////////////////////////
    std::string jecName;
    if(isData){
        if(!is2017){
            if (runNumber < 271658){
                jecName = "A";
                std::cerr << "no JEC available for 2016 run A, seems like JSON file is not applied!" << std::endl;
            }
            else if(runNumber < 276812) jecName = "BCD";
            else if(runNumber < 278809) jecName = "EF";
            else if(runNumber < 280385) jecName = "GV";
            else if(runNumber < 294645) jecName = "H";
            jecName += version2016;
        } else {
            if(runNumber < 297020){
                jecName = "A";
                std::cerr << "no JEC available for 2017 run A, seems like JSON file is not applied!" << std::endl;
            } 
            else if(runNumber < 299337) jecName = "B";
            else if(runNumber < 302030) jecName = "C";
            else if(runNumber < 303435) jecName = "D";
            else if(runNumber < 304911) jecName = "E";
            else if(runNumber < 306464) jecName = "F";
            else{
                jecName = "GH";
                std::cerr << "no JEC available for 2017 runs G-H, they are not 13 TeV data! Seems like JSON file is not applied" << std::endl;
            }
            jecName += version2017;
        }
        jecName += "_DATA";
        return jecName;
    } else{
        if(!is2017){
            jecName = version2016;
        } else{
            jecName = version2017;
        }
        return (jecName + "_MC");
    }
}

std::string JEC::getJECName(const unsigned long runNumber){
    if(!is2017){
        return "Summer16_07Aug2017" + getJECRunName(runNumber);
    } else{
        return "Fall17_17Nov2017" + getJECRunName(runNumber);
    }
}

std::vector<float> JEC::getSubCorrections(double rawPt, double eta, double rho, double area){
    jetCorrector->setJetEta(eta);
    jetCorrector->setRho(rho);
    jetCorrector->setJetA(area);
    jetCorrector->setJetPt(rawPt); 
    std::vector< float > corrections = jetCorrector->getSubCorrections();
    return corrections;
}

double JEC::jetCorrection(double rawPt, double eta, double rho, double area, const std::string& level){

    std::vector< float > corrections = getSubCorrections(rawPt, eta, rho, area); 

    if (level == "L1FastJet")    return corrections.front();
    if (level == "L2Relative")   return corrections[1];
    if (level == "L2L3Residual") return corrections.back();
    return corrections[2];
}


double JEC::jetUncertainty(double pt, double eta){
    if(eta> 5.0) eta = 5.0;
    else if (eta<-5.0) eta =-5.0;
    jetUncertainties->setJetPt(pt);
    jetUncertainties->setJetEta(eta);
    return jetUncertainties->getUncertainty(true);
}


std::pair<double, double> JEC::getMETCorrectionPxPy(double rawPt, double rawEta, double rawMuonSubtractedPt, double phi, double emf, double rho, double area){

    std::vector< float > corrections = getSubCorrections(rawPt, rawEta, rho, area); 

    double l1corrpt   = rawMuonSubtractedPt*corrections.front(); // l1fastjet corrections were pushed pack first
    double fullcorrpt = rawMuonSubtractedPt*corrections.back();  // full corrections are the last in the vector
    // the corretions for the MET are the difference between l1fastjet and the full corrections on the jet!
    if(emf > 0.9 or fullcorrpt < 15. || (fabs(rawEta) > 9.9) ) return {0., 0.}; // skip jets with EMF > 0.9
    
    std::pair<double, double> corr = {px(l1corrpt - fullcorrpt, phi), py(l1corrpt - fullcorrpt, phi)};
    return corr; 
}



std::pair<double, double> JEC::correctedMETAndPhi(const pat::MET& met, const std::vector< pat::Jet >& jets, const double rho){
    double corrMETx = met.uncorPx();
    double corrMETy = met.uncorPy();
    
    //loop over all jets
    for(auto& jet : jets){
        //make lorentzVector of raw jet pt
        TLorentzVector jetV;
        jetV.SetPtEtaPhiE(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jet.correctedP4("Uncorrected").Phi(), jet.correctedP4("Uncorrected").E());
        //clean jet from muons
        const std::vector<reco::CandidatePtr>& daughters = jet.daughterPtrVector();
        for(auto& daughterPtr : daughters){
            const reco::PFCandidate* daughter = dynamic_cast<const reco::PFCandidate* >( daughterPtr.get() );
            const reco::Candidate* muon = (daughter != nullptr ?  
                                            (daughter->muonRef().isNonnull() ? daughter->muonRef().get() : nullptr)
                                            : daughterPtr.get() );
            if(muon != nullptr && ( muon->isGlobalMuon() || muon->isStandAloneMuon() ) ){
                TLorentzVector muonV(muon->px(), muon->py(), muon->pz(), muon->energy());
                jetV -= muonV;
            }
        }
        //get JEC on px and py 
        std::pair<double, double> corr = getMETCorrectionPxPy(jet.correctedP4("Uncorrected").Pt(), jet.correctedP4("Uncorrected").Eta(), jetV.Pt(), 
            jetV.Phi(), jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction(), rho, jet.jetArea() );
        
        //apply corrections to current met values
        corrMETx += corr.first;
        corrMETy += corr.second;
    }
    double correctedMET = sqrt(corrMETx*corrMETx + corrMETy*corrMETy);
    double correctedMETPhi = atan2(corrMETx, corrMETy);
    return {correctedMET, correctedMETPhi};
}
