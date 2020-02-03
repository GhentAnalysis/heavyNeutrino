/*
 * implementation of JEC class
 * WARNING: it seems this code is not used
 * TODO: keep an eye on https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Recommended_for_Data for updates
 */

#include "heavyNeutrino/multilep/interface/JEC.h"

enum TheRunEra{y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,y2017B,y2017C,y2017D,y2017E,y2017F,y2018A,y2018B,y2018C,y2018D,y2016MC,y2017MC,y2018MC};

JEC::JEC(const std::string& JECPath, bool dataSample, bool is2017sample, bool is2018sample, bool isPuppiJEC):
    path(JECPath), isData(dataSample), is2017(is2017sample), is2018(is2018sample), isPuppi(isPuppiJEC), currentJEC("none") {}

JEC::~JEC(){}

void JEC::updateJEC(const unsigned long runNumber){
    if(getJECName(runNumber) != currentJEC){
        currentJEC = getJECName(runNumber);
        setJEC(getJECName(runNumber) );
    }
}

void JEC::setJEC(const std::string& JECName){
    std::string jettype = (isPuppi)? "Puppi" : "chs";
    if(!is2018) jettype = "chs"; //currently Puppi JEC corrections are only loaded for 2018!
    std::vector< JetCorrectorParameters > JECParameters;
    JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L1FastJet_AK4PF" + jettype + ".txt") );
    JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L2Relative_AK4PF" + jettype + ".txt") );
    JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L3Absolute_AK4PF" + jettype + ".txt") );
    if(isData) JECParameters.push_back(JetCorrectorParameters( path + JECName + "_L2L3Residual_AK4PF" + jettype + ".txt") );
    jetCorrector.reset(new FactorizedJetCorrector(JECParameters) );
    jetUncertainties.reset(new JetCorrectionUncertainty(path + JECName + "_Uncertainty_AK4PF" + jettype + ".txt") );
}
   
std::string JEC::getJECRunName(const unsigned long runNumber){
    ///////////////////////////
    static const std::string version2016 = "_V11";
    static const std::string version2017 = "_V32";
    static const std::string version2018 = "_V19";
    //////////////////////////
    std::string jecName;
    if(isData){
        if(is2018){
            if(runNumber < 316998)      jecName = "_RunA";
            else if(runNumber < 319313) jecName = "_RunB";
            else if(runNumber < 320394) jecName = "_RunC";
            else if(runNumber < 325274) jecName = "_RunD";
            else                        std::cerr << "no JEC available for 2018 run E, they are not 13 TeV data! Seems like JSON file is not applied" << std::endl;
            jecName += version2018;
        } else if(is2017){
            if(runNumber < 297020)      std::cerr << "no JEC available for 2017 run A, seems like JSON file is not applied!" << std::endl;
            else if(runNumber < 299337) jecName = "B";
            else if(runNumber < 302030) jecName = "C";
            else if(runNumber < 303435) jecName = "D";
            else if(runNumber < 304911) jecName = "E";
            else if(runNumber < 306464) jecName = "F";
            else                        std::cerr << "no JEC available for 2017 runs G-H, they are not 13 TeV data! Seems like JSON file is not applied" << std::endl;
            jecName += version2017;
        } else {
            if (runNumber < 271658)     std::cerr << "no JEC available for 2016 run A, seems like JSON file is not applied!" << std::endl;
            else if(runNumber < 276812) jecName = "BCD";
            else if(runNumber < 278809) jecName = "EF";
            else if(runNumber < 294645) jecName = "GH";
            jecName += version2016;
        }
        return (jecName + "_DATA");
    } else{
        if(is2018)      jecName = version2018;
        else if(is2017) jecName = version2017;
        else            jecName = version2016;
        return (jecName + "_MC");
    }
}

std::string JEC::getJECName(const unsigned long runNumber){
    if(is2018)      return "Autumn18" + getJECRunName(runNumber);
    else if(is2017) return "Fall17_17Nov2017" + getJECRunName(runNumber);
    else            return "Summer16_07Aug2017" + getJECRunName(runNumber);
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
    double correctedMETPhi = atan2(corrMETy, corrMETx);
    return {correctedMET, correctedMETPhi};
}
