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
    bool _loose[nL_max];                                                           //lepton selection decisions
    bool _FO[nL_max];
    bool _tight[nL_max];
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
    template <class leptonType> void fillLeptonIdVars(const leptonType&);          //Fill basic lepton ID variables
    void fillMetFilterVars(const edm::Event&);                                     //Make MET filter decision
    bool metFiltersFlagged(const edm::Event&);                                     //Fill MET filter variables
};
#endif

