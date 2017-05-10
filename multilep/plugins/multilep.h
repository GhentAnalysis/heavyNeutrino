#ifndef MULTILEP_H
#define MULTILEP_H
// system include files

#include <memory>
#include <typeinfo> //To check lepton types

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//new include statements
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


//Dataformats
//#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "heavyNeutrino/multilep/interface/Tools.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


// Include for leptonId functions
#include "heavyNeutrino/multilep/interface/electronIdentification.h"
#include "TLorentzVector.h"



//
// class declaration
//

class multilep : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit multilep(const edm::ParameterSet&);
    ~multilep(){};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //Define EDgetTokens to read data from events
    edm::EDGetTokenT<std::vector<reco::Vertex>>         vtxToken;
    edm::EDGetTokenT<std::vector<pat::Muon>>            muonToken;
    edm::EDGetTokenT<std::vector<pat::Electron>>        eleToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              eleMvaToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              eleMvaHZZToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedTightToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedMediumToken;
    edm::EDGetTokenT<std::vector<pat::Tau>>             tauToken;
    edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidatesToken;    //particle collection used to calculate isolation variables
    edm::EDGetTokenT<double>                            rhoToken;                 //neutal energy density in terms of deltaR used for pileup corrections
    edm::EDGetTokenT<double>                            rhoTokenAll;              //energy density used for JEC
    edm::EDGetTokenT<std::vector<pat::MET>>             metToken;                 //missing transverse energy
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetToken;                 //jet collection
  //edm::EDGetTokenT<reco::JetCorrector>                jecToken;                 //JEC
    edm::EDGetTokenT<edm::TriggerResults>               triggerToken;             //Trigger information
    edm::EDGetTokenT<edm::TriggerResults>               recoResultsToken;         //MET filter information
    edm::EDGetTokenT<bool>                              badPFMuonFilterToken;     //MET filter not stored in miniAOD
    edm::EDGetTokenT<bool>                              badChCandFilterToken;     //MET filter not stored in miniAOD

    edm::Service<TFileService> fs;                                                 //Root tree and file for storing event info
    //FILE* outFile;
    TTree* outputTree;
    static const unsigned nL_max = 20;                                             //maximum number of particles stored
    static const unsigned nJets_max = 20;
    unsigned long _runNb;                                                          //event labels
    unsigned long _lumiBlock;
    unsigned long _eventNb;
    unsigned _nL;                                                                  //number of leptons
    unsigned _nMu;
    unsigned _nEle;
    unsigned _nLight;
    unsigned _nTau;
    double _lPt[nL_max];                                                           //lepton kinematics
    double _lEta[nL_max];                                                          //QUESTION: what was the reason again to use arrays instead of vectors?
    double _lPhi[nL_max];
    double _lE[nL_max];
    unsigned _flavor[nL_max];                                                      //other lepton variables
    int _charge[nL_max];
    double _relIso[nL_max];
    double _miniIso[nL_max];
    double _dxy[nL_max];
    double _dz[nL_max];
    double _3dIP[nL_max];
    double _3dIPSig[nL_max];
    bool _lLoose[nL_max];                                                           //lepton selection decisions
    bool _lFO[nL_max];
    bool _lTight[nL_max];
    bool _isPrompt[nL_max];                                                        //generator lepton variables
    bool _truthPdg[nL_max];
    bool _truthMomPdg[nL_max];
    unsigned _origin[nL_max];
    unsigned _nJets;                                                               //jet variables
    double _jetPt[nJets_max];
    double _jetEta[nJets_max];
    double _jetPhi[nJets_max];
    double _jetE[nJets_max];
    double _met;                                                                   //met kinematics
    double _metPhi;
    unsigned _nVertex;                                                             //Event variables
    bool _metFiltersFlagged;
    bool _passHnlTrigger[4];                                                       //0 = eee, 1 = eem, 2 = emm, 3 = mmm
    bool _badMuonFlagged;
    bool _badCloneMuonFlagged;

    //Additional class functions
    void fillTriggerVars(const edm::Event&);
    bool trigPass(unsigned, edm::Handle<edm::TriggerResults>&, const edm::Event&);
    bool metFilterFlagged();                                                       //Events flagged by a met filter should not be used
    void fillLeptonGenVars(const reco::GenParticle* genParticle);                  //Fill MC-truth lepton variables
    void fillLeptonKinVars(const reco::Candidate&);                                //Fill reconstructed lepton kinematics
    void fillMetFilterVars(const edm::Event&);                                     //Make MET filter decision
    bool metFiltersFlagged(const edm::Event&);                                     //Fill MET filter variables

    /*
     * Everything below is specific to the heavyNeutrino analysis, makes the tuples unusable for other analyses
     * Can't put it in a separate namespace because it accesses member variables of this class
     * Maybe better to rewrite some parts here
     */
    /*
    template <class leptonType> void fillLeptonIdVars(const leptonType& lepton){
      _lLoose[_nL] = isLoose(lepton);
      _lFO[_nL]    = isFO(lepton);
      _lTight[_nL] = isTight(lepton);
    }
    // Important: not official-loose like in POG-loose, but own-made loose, never call this a 'loose' lepton in a presentation
    template <class leptonType> bool isLoose(const leptonType& lepton){
      if(fabs(_dxy[_nL]) > 0.05 || fabs(_dz[_nL]) > 0.1) return false;
      if(_relIso[_nL] > 0.6)                             return false;

      if(lepton.isMuon()){
    //		return (lepton.isPFMuon() && (lepton.isTrackerMuon() || lepton.isGlobalMuon()) && (lepton.pt > 5);
        return (lepton.isLooseMuon() && lepton.pt() > 5); // Don't we apply pt cuts already in multilepton.cc? Should remove these
      } else if(lepton.isElectron()){   // This is the loose electron ?????????? Nothing more???
//        if(lepton.gsfTrack()->hitPattern()->numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) return false; // BROKEN
        if(not lepton.passConversionVeto())                                                         return false;
        if(eleMuOverlap(lepton))                                                                    return false; // Always run electrons after muons because of this
        return (lepton.pt() > 10);
      } else if(lepton.isTau()){
        
      }
    }

    template <class leptonType> bool isFO(const leptonType& lepton){
      if(!_lLoose[_nL])          return false; // own-made loose, not POG-loose
      if(fabs(_3dIPSig[_nL]) > 4) return false;

      if(lepton.isMuon()){
        return lepton.isMediumMuon();
      } else if(lepton.isElectron()){
   //     if(lepton.gsfTrack()->hitPattern()->numberOfHits(reco::HitPattern::MISSING_INNER_HITS) != 0) return false; //BROKEN
        if(!electronIdentification::passTriggerEmulationDoubleEG(*lepton))                           return false;
        return true;
      }

    }

    template <class leptonType> bool isTight(const leptonType& lepton){
      if(!_lFO[_nL]) return false;

    }

    //Check if electron overlaps with loose muon
    bool eleMuOverlap(const pat::Electron& ele){
      TLorentzVector eleV, muV;
      eleV.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy());
      for(unsigned m = 0; m < _nMu; ++m){
        if(_lLoose[m]){
          TLorentzVector muV(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
          if(eleV.DeltaR(muV) < 0.05) return true;
        }
      }
      return false;	
    }
*/
};
#endif

