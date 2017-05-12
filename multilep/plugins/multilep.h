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
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


// Include for leptonId functions
#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "TLorentzVector.h"



//
// class declaration
//
class LeptonAnalyzer;

class multilep : public edm::one::EDAnalyzer<edm::one::SharedResources> {

  public:
    explicit multilep(const edm::ParameterSet&);
    ~multilep(){};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    // Define EDgetTokens to read data from events
    // Public such that we can easily access them in the individual object analyzers)
    edm::EDGetTokenT<std::vector<reco::Vertex>>         vtxToken;
    edm::EDGetTokenT<std::vector<pat::Muon>>            muonToken;
    edm::EDGetTokenT<std::vector<pat::Electron>>        eleToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              eleMvaToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              eleMvaHZZToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedTightToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedMediumToken;
    edm::EDGetTokenT<std::vector<pat::Tau>>             tauToken;
    edm::EDGetTokenT<std::vector<pat::Photon>>          photonToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               photonCutBasedLooseToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               photonCutBasedMediumToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               photonCutBasedTightToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              photonMvaToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              photonChargedIsolationToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              photonNeutralHadronIsolationToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              photonPhotonIsolationToken;
    edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidatesToken;                       //particle collection used to calculate isolation variables
    edm::EDGetTokenT<double>                            rhoToken;                                    //neutal energy density in terms of deltaR used for pileup corrections
    edm::EDGetTokenT<double>                            rhoTokenAll;                                 //energy density used for JEC
    edm::EDGetTokenT<std::vector<pat::MET>>             metToken;                                    //missing transverse energy
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetToken;                                    //jet collection
  //edm::EDGetTokenT<reco::JetCorrector>                jecToken;                                    //JEC
    edm::EDGetTokenT<edm::TriggerResults>               triggerToken;                                //Trigger information
    edm::EDGetTokenT<edm::TriggerResults>               recoResultsToken;                            //MET filter information
    edm::EDGetTokenT<bool>                              badPFMuonFilterToken;                        //MET filter not stored in miniAOD
    edm::EDGetTokenT<bool>                              badChCandFilterToken;                        //MET filter not stored in miniAOD

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    LeptonAnalyzer* leptonAnalyzer;

    edm::Service<TFileService> fs;                                                                   //Root tree and file for storing event info
    //FILE* outFile;
    TTree* outputTree;
    static const unsigned nL_max      = 20;                                                          //maximum number of particles stored
    static const unsigned nJets_max   = 20;
    static const unsigned nPhoton_max = 5;

    unsigned long _runNb;                                                                            //event labels
    unsigned long _lumiBlock;
    unsigned long _eventNb;
    unsigned _nJets;                                                                                 //jet variables
    double _jetPt[nJets_max];
    double _jetEta[nJets_max];
    double _jetPhi[nJets_max];
    double _jetE[nJets_max];
    double _met;                                                                                     //met kinematics
    double _metPhi;
    unsigned _nVertex;                                                                               //Event variables
    bool _metFiltersFlagged;
    bool _passHnlTrigger[4];                                                                         //0 = eee, 1 = eem, 2 = emm, 3 = mmm
    bool _badMuonFlagged;
    bool _badCloneMuonFlagged;

    unsigned _nPhoton;                                                                               // photon variables
    float    _photonPt[nPhoton_max];
    float    _photonEta[nPhoton_max];
    float    _photonPhi[nPhoton_max];
    float    _photonE[nPhoton_max];
    bool     _photonCutBasedLoose[nPhoton_max];
    bool     _photonCutBasedMedium[nPhoton_max];
    bool     _photonCutBasedTight[nPhoton_max];
    float    _photonMva[nPhoton_max];
    float    _photonChargedIsolation[nPhoton_max];
    float    _photonNeutralHadronIsolation[nPhoton_max];
    float    _photonPhotonIsolation[nPhoton_max];
    float    _photonSigmaIetaIeta[nPhoton_max];
    float    _photonHadronicOverEm[nPhoton_max];
    bool     _photonPassElectronVeto[nPhoton_max];
    bool     _photonHasPixelSeed[nPhoton_max];

    //Additional class functions
    void fillTriggerVars(const edm::Event&);
    bool trigPass(unsigned, edm::Handle<edm::TriggerResults>&, const edm::Event&);
    bool metFilterFlagged();                                                                         //Events flagged by a met filter should not be used
    void fillMetFilterVars(const edm::Event&);                                                       //Make MET filter decision
    bool metFiltersFlagged(const edm::Event&);                                                       //Fill MET filter variables
};
#endif

