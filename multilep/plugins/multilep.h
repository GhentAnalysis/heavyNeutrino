#ifndef MULTILEP_H
#define MULTILEP_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "heavyNeutrino/multilep/interface/TriggerAnalyzer.h"
#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenAnalyzer.h"
#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"




//
// class declaration
//
class TriggerAnalyzer;
class LeptonAnalyzer;
class PhotonAnalyzer;
class JetAnalyzer;
class GenAnalyzer;
class LheAnalyzer;

class multilep : public edm::one::EDAnalyzer<edm::one::SharedResources> {

  public:
    explicit multilep(const edm::ParameterSet&);
    ~multilep(){};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    // Define EDgetTokens to read data from events
    // Public such that we can easily access them in the individual object analyzers)
    edm::EDGetTokenT<std::vector<reco::Vertex>>         vtxToken;
    edm::EDGetTokenT<GenEventInfoProduct>               genEventInfoToken;
    edm::EDGetTokenT<LHEEventProduct>                   lheEventInfoToken;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>>    pileUpToken;
    edm::EDGetTokenT<reco::GenParticleCollection>       genParticleToken; 
    edm::EDGetTokenT<std::vector<pat::Muon>>            muonToken;
    edm::EDGetTokenT<std::vector<pat::Electron>>        eleToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              eleMvaToken;
    edm::EDGetTokenT<edm::ValueMap<float>>              eleMvaHZZToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>               eleCutBasedVetoToken;
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
    edm::EDGetTokenT<edm::ValueMap<float>>              photonFull5x5SigmaIEtaIPhiToken;
    edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidatesToken;                       //particle collection used to calculate isolation variables
    edm::EDGetTokenT<double>                            rhoToken;                                    //neutal energy density in terms of deltaR used for pileup corrections
    edm::EDGetTokenT<double>                            rhoTokenAll;                                 //energy density used for JEC
    edm::EDGetTokenT<std::vector<pat::MET>>             metToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedUpToken;
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedDownToken;
    edm::EDGetTokenT<edm::TriggerResults>               recoResultsToken;                            //MET filter information
    edm::EDGetTokenT<edm::TriggerResults>               triggerToken;
    edm::EDGetTokenT<pat::PackedTriggerPrescales>       prescalesToken;
    edm::EDGetTokenT<bool>                              badPFMuonFilterToken;                        //MET filter not stored in miniAOD
    edm::EDGetTokenT<bool>                              badChCandFilterToken;                        //MET filter not stored in miniAOD
    std::string                                         skim;
    bool                                                isData;

    GenAnalyzer*     genAnalyzer;                                                                    //Public because the photonAnalyzer uses some of its helper functions

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    TriggerAnalyzer* triggerAnalyzer;
    LeptonAnalyzer*  leptonAnalyzer;
    PhotonAnalyzer*  photonAnalyzer;
    JetAnalyzer*     jetAnalyzer;
    LheAnalyzer*     lheAnalyzer;

    edm::Service<TFileService> fs;                                                                   //Root tree and file for storing event info
    TTree* outputTree;

    unsigned long _runNb;                                                                            //event labels
    unsigned long _lumiBlock;
    unsigned long _eventNb;
    unsigned      _nVertex;                                                                          //Event variables
    double        _met;                                                                              //met kinematics
    double        _metPhi;
    double        _metJECDown;
    double        _metPhiJECDown;
    double        _metJECUp;
    double        _metPhiJECUp;
    double        _metUnclDown;
    double        _metPhiUnclDown;
    double        _metUnclUp;
    double        _metPhiUnclUp;
};
#endif

