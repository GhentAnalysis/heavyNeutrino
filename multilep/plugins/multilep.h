#ifndef MULTILEP_H
#define MULTILEP_H
//#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"




//
// class declaration
//
class LeptonAnalyzer;
class PhotonAnalyzer;
class JetAnalyzer;

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
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedLooseToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedMediumToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedTightToken;
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
    edm::EDGetTokenT<std::vector<pat::MET>>             metToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedUpToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedDownToken;
    edm::EDGetTokenT<edm::TriggerResults>               triggerToken;
    edm::EDGetTokenT<edm::TriggerResults>               recoResultsToken;                            //MET filter information
    edm::EDGetTokenT<bool>                              badPFMuonFilterToken;                        //MET filter not stored in miniAOD
    edm::EDGetTokenT<bool>                              badChCandFilterToken;                        //MET filter not stored in miniAOD

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    LeptonAnalyzer* leptonAnalyzer;
    PhotonAnalyzer* photonAnalyzer;
    JetAnalyzer*    jetAnalyzer;

    edm::Service<TFileService> fs;                                                                   //Root tree and file for storing event info
    TTree* outputTree;

    unsigned long _runNb;                                                                            //event labels
    unsigned long _lumiBlock;
    unsigned long _eventNb;
    double _met;                                                                                     //met kinematics
    double _metPhi;
    unsigned _nVertex;                                                                               //Event variables
    bool _metFiltersFlagged;
    bool _passHnlTrigger[4];                                                                         //0 = eee, 1 = eem, 2 = emm, 3 = mmm
    bool _badMuonFlagged;
    bool _badCloneMuonFlagged;
    //Additional class functions
    void fillTriggerVars(const edm::Event&);
    bool trigPass(unsigned, edm::Handle<edm::TriggerResults>&, const edm::Event&);
    bool metFilterFlagged();                                                                         //Events flagged by a met filter should not be used
    void fillMetFilterVars(const edm::Event&);                                                       //Make MET filter decision
    bool metFiltersFlagged(const edm::Event&);                                                       //Fill MET filter variables
};
#endif

