/*
 * implementation of JEC class
 * WARNING: it seems this code is not used
 * TODO: keep an eye on https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Recommended_for_Data for updates
 */

#include "heavyNeutrino/multilep/interface/JEC.h"

enum TheRunEra{y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,y2017B,y2017C,y2017D,y2017E,y2017F,y2018A,y2018B,y2018C,y2018D,y2016MC,y2017MC,y2018MC};

JEC::JEC(const std::string& JECPath, bool dataSample, bool is2017sample, bool is2018sample, bool isPuppiJEC):
    path(JECPath), isData(dataSample), is2017(is2017sample), is2018(is2018sample), isPuppi(isPuppiJEC), currentJEC("none"), runera(-1), usemetv2(false) {}

JEC::~JEC(){}

void JEC::updateJEC(const unsigned long runNumber){
    if(getJECName(runNumber) != currentJEC){
        currentJEC = getJECName(runNumber);
        //setJEC(getJECName(runNumber) );
        setRunEra(runNumber);//purely for XY corrections
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

void JEC::setRunEra(const unsigned long runnb){
    usemetv2 = false;
    if(isData){
        if(runnb >=272007 &&runnb<=275376  ) runera = y2016B;
        else if(runnb >=275657 &&runnb<=276283  ) runera = y2016C;
        else if(runnb >=276315 &&runnb<=276811  ) runera = y2016D;
        else if(runnb >=276831 &&runnb<=277420  ) runera = y2016E;
        else if(runnb >=277772 &&runnb<=278808  ) runera = y2016F;
        else if(runnb >=278820 &&runnb<=280385  ) runera = y2016G;
        else if(runnb >=280919 &&runnb<=284044  ) runera = y2016H;

        else if(runnb >=297020 &&runnb<=299329 ){ runera = y2017B; usemetv2 =true;}
        else if(runnb >=299337 &&runnb<=302029 ){ runera = y2017C; usemetv2 =true;}
        else if(runnb >=302030 &&runnb<=303434 ){ runera = y2017D; usemetv2 =true;}
        else if(runnb >=303435 &&runnb<=304826 ){ runera = y2017E; usemetv2 =true;}
        else if(runnb >=304911 &&runnb<=306462 ){ runera = y2017F; usemetv2 =true;}

        else if(runnb >=315252 &&runnb<=316995 ) runera = y2018A;
        else if(runnb >=316998 &&runnb<=319312 ) runera = y2018B;
        else if(runnb >=319313 &&runnb<=320393 ) runera = y2018C;
        else if(runnb >=320394 &&runnb<=325273 ) runera = y2018D;
        else {
            //Couldn't find data era => no correction applied
            runera = -1;
        }
    }else {
        if(is2017) {runera = y2017MC; usemetv2 =true;}
        else if(is2018) runera = y2018MC;
        else runera = y2016MC;
    }
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


std::pair<double,double> JEC::METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int npv){

    std::pair<double,double>  TheXYCorr_Met_MetPhi(uncormet,uncormet_phi);
    if(npv>100) npv=100;
    if(runera == -1) return TheXYCorr_Met_MetPhi;//return uncorrected if runera was not correctly initialized
    double METxcorr(0.),METycorr(0.);

    if(!usemetv2){//Current recommendation for 2016 and 2018
        if(runera==y2016B) METxcorr = -(-0.0478335*npv -0.108032);
        if(runera==y2016B) METycorr = -(0.125148*npv +0.355672);
        if(runera==y2016C) METxcorr = -(-0.0916985*npv +0.393247);
        if(runera==y2016C) METycorr = -(0.151445*npv +0.114491);
        if(runera==y2016D) METxcorr = -(-0.0581169*npv +0.567316);
        if(runera==y2016D) METycorr = -(0.147549*npv +0.403088);
        if(runera==y2016E) METxcorr = -(-0.065622*npv +0.536856);
        if(runera==y2016E) METycorr = -(0.188532*npv +0.495346);
        if(runera==y2016F) METxcorr = -(-0.0313322*npv +0.39866);
        if(runera==y2016F) METycorr = -(0.16081*npv +0.960177);
        if(runera==y2016G) METxcorr = -(0.040803*npv -0.290384);
        if(runera==y2016G) METycorr = -(0.0961935*npv +0.666096);
        if(runera==y2016H) METxcorr = -(0.0330868*npv -0.209534);
        if(runera==y2016H) METycorr = -(0.141513*npv +0.816732);
        if(runera==y2017B) METxcorr = -(-0.259456*npv +1.95372);
        if(runera==y2017B) METycorr = -(0.353928*npv -2.46685);
        if(runera==y2017C) METxcorr = -(-0.232763*npv +1.08318);
        if(runera==y2017C) METycorr = -(0.257719*npv -1.1745);
        if(runera==y2017D) METxcorr = -(-0.238067*npv +1.80541);
        if(runera==y2017D) METycorr = -(0.235989*npv -1.44354);
        if(runera==y2017E) METxcorr = -(-0.212352*npv +1.851);
        if(runera==y2017E) METycorr = -(0.157759*npv -0.478139);
        if(runera==y2017F) METxcorr = -(-0.232733*npv +2.24134);
        if(runera==y2017F) METycorr = -(0.213341*npv +0.684588);
        if(runera==y2018A) METxcorr = -(0.362865*npv -1.94505);
        if(runera==y2018A) METycorr = -(0.0709085*npv -0.307365);
        if(runera==y2018B) METxcorr = -(0.492083*npv -2.93552);
        if(runera==y2018B) METycorr = -(0.17874*npv -0.786844);
        if(runera==y2018C) METxcorr = -(0.521349*npv -1.44544);
        if(runera==y2018C) METycorr = -(0.118956*npv -1.96434);
        if(runera==y2018D) METxcorr = -(0.531151*npv -1.37568);
        if(runera==y2018D) METycorr = -(0.0884639*npv -1.57089);
        if(runera==y2016MC) METxcorr = -(-0.195191*npv -0.170948);
        if(runera==y2016MC) METycorr = -(-0.0311891*npv +0.787627);
        if(runera==y2017MC) METxcorr = -(-0.217714*npv +0.493361);
        if(runera==y2017MC) METycorr = -(0.177058*npv -0.336648);
        if(runera==y2018MC) METxcorr = -(0.296713*npv -0.141506);
        if(runera==y2018MC) METycorr = -(0.115685*npv +0.0128193);
    }
    else {//these are the corrections for v2 MET recipe (currently recommended for 2017)
        if(runera==y2016B) METxcorr = -(-0.0374977*npv +0.00488262);
        if(runera==y2016B) METycorr = -(0.107373*npv +-0.00732239);
        if(runera==y2016C) METxcorr = -(-0.0832562*npv +0.550742);
        if(runera==y2016C) METycorr = -(0.142469*npv +-0.153718);
        if(runera==y2016D) METxcorr = -(-0.0400931*npv +0.753734);
        if(runera==y2016D) METycorr = -(0.127154*npv +0.0175228);
        if(runera==y2016E) METxcorr = -(-0.0409231*npv +0.755128);
        if(runera==y2016E) METycorr = -(0.168407*npv +0.126755);
        if(runera==y2016F) METxcorr = -(-0.0161259*npv +0.516919);
        if(runera==y2016F) METycorr = -(0.141176*npv +0.544062);
        if(runera==y2016G) METxcorr = -(0.0583851*npv +-0.0987447);
        if(runera==y2016G) METycorr = -(0.0641427*npv +0.319112);
        if(runera==y2016H) METxcorr = -(0.0706267*npv +-0.13118);
        if(runera==y2016H) METycorr = -(0.127481*npv +0.370786);
        if(runera==y2017B) METxcorr = -(-0.19563*npv +1.51859);
        if(runera==y2017B) METycorr = -(0.306987*npv +-1.84713);
        if(runera==y2017C) METxcorr = -(-0.161661*npv +0.589933);
        if(runera==y2017C) METycorr = -(0.233569*npv +-0.995546);
        if(runera==y2017D) METxcorr = -(-0.180911*npv +1.23553);
        if(runera==y2017D) METycorr = -(0.240155*npv +-1.27449);
        if(runera==y2017E) METxcorr = -(-0.149494*npv +0.901305);
        if(runera==y2017E) METycorr = -(0.178212*npv +-0.535537);
        if(runera==y2017F) METxcorr = -(-0.165154*npv +1.02018);
        if(runera==y2017F) METycorr = -(0.253794*npv +0.75776);
        if(runera==y2018A) METxcorr = -(0.362642*npv +-1.55094);
        if(runera==y2018A) METycorr = -(0.0737842*npv +-0.677209);
        if(runera==y2018B) METxcorr = -(0.485614*npv +-2.45706);
        if(runera==y2018B) METycorr = -(0.181619*npv +-1.00636);
        if(runera==y2018C) METxcorr = -(0.503638*npv +-1.01281);
        if(runera==y2018C) METycorr = -(0.147811*npv +-1.48941);
        if(runera==y2018D) METxcorr = -(0.520265*npv +-1.20322);
        if(runera==y2018D) METycorr = -(0.143919*npv +-0.979328);
        if(runera==y2016MC) METxcorr = -(-0.159469*npv +-0.407022);
        if(runera==y2016MC) METycorr = -(-0.0405812*npv +0.570415);
        if(runera==y2017MC) METxcorr = -(-0.182569*npv +0.276542);
        if(runera==y2017MC) METycorr = -(0.155652*npv +-0.417633);
        if(runera==y2018MC) METxcorr = -(0.299448*npv +-0.13866);
        if(runera==y2018MC) METycorr = -(0.118785*npv +0.0889588);
    }

    double CorrectedMET_x = uncormet *cos( uncormet_phi)+METxcorr;
    double CorrectedMET_y = uncormet *sin( uncormet_phi)+METycorr;

    double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
    double CorrectedMETPhi;
    if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
    else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
    else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
    else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
    else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
    else CorrectedMETPhi =0;

    TheXYCorr_Met_MetPhi.first= CorrectedMET;
    TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;
    return TheXYCorr_Met_MetPhi;

}
