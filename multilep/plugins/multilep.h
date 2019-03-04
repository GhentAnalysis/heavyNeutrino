#ifndef MULTILEP_H
#define MULTILEP_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
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
//New for SUSY masses
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "heavyNeutrino/multilep/interface/TriggerAnalyzer.h"
#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"
#include "heavyNeutrino/multilep/interface/JetAnalyzer.h"
#include "heavyNeutrino/multilep/interface/GenAnalyzer.h"
#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"
#include "heavyNeutrino/multilep/interface/SUSYMassAnalyzer.h"

//
// class declaration
//
class TriggerAnalyzer;
class LeptonAnalyzer;
class PhotonAnalyzer;
class JetAnalyzer;
class GenAnalyzer;
class LheAnalyzer;
class SUSYMassAnalyzer;

class multilep : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks, edm::one::WatchRuns, edm::one::SharedResources> {
    //Define other analyzers as friends
    friend TriggerAnalyzer;
    friend LeptonAnalyzer;
    friend PhotonAnalyzer;
    friend JetAnalyzer;
    friend GenAnalyzer;
    friend LheAnalyzer;
    friend SUSYMassAnalyzer;
    public:
        explicit multilep(const edm::ParameterSet&);
        ~multilep();

    private:
        // Define EDgetTokens to read data from events
        edm::EDGetTokenT<reco::BeamSpot>		    beamSpotToken;
        edm::EDGetTokenT<std::vector<reco::Vertex>>         vtxToken;
        edm::EDGetTokenT<GenEventInfoProduct>               genEventInfoToken;
        edm::EDGetTokenT<GenLumiInfoHeader>                 genLumiInfoToken;
        edm::EDGetTokenT<LHEEventProduct>                   lheEventInfoToken;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>>    pileUpToken;
        edm::EDGetTokenT<reco::GenParticleCollection>       genParticleToken; 
        edm::EDGetTokenT<std::vector<pat::PackedGenParticle>>packedGenParticleToken;
        edm::EDGetTokenT<std::vector<pat::Muon>>            muonToken;
        edm::EDGetTokenT<std::vector<pat::Electron>>        eleToken;
        edm::EDGetTokenT<std::vector<pat::Tau>>             tauToken;
        edm::EDGetTokenT<std::vector<pat::Photon>>          photonToken;
        edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidatesToken;                       //particle collection used to calculate isolation variables
        edm::EDGetTokenT<double>                            rhoToken;
        edm::EDGetTokenT<std::vector<pat::MET>>             metToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>             jetToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedUpToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>             jetSmearedDownToken;
        edm::EDGetTokenT<edm::TriggerResults>               recoResultsPrimaryToken;                     //MET filter information
        edm::EDGetTokenT<edm::TriggerResults>               recoResultsSecondaryToken;                   //MET filter information (fallback if primary is not available)
        edm::EDGetTokenT<edm::TriggerResults>               triggerToken;
        edm::EDGetTokenT<pat::PackedTriggerPrescales>       prescalesToken;
        edm::EDGetTokenT<std::vector<reco::Vertex>>         secondaryVerticesToken;
        edm::EDGetTokenT<bool>                              ecalBadCalibFilterToken;
        std::string                                         skim;
        bool                                                isData;
        bool                                                is2017;
        bool                                                is2018;
        bool                                                isSUSY;
        bool                                                storeLheParticles;

        virtual void beginJob() override;
        virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
        virtual void beginRun(const edm::Run&, edm::EventSetup const&) override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        virtual void endRun(const edm::Run&, edm::EventSetup const&) override {}                         //Unused functions
        virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {}
        virtual void endJob() override {};

        TriggerAnalyzer*  triggerAnalyzer;
        LeptonAnalyzer*   leptonAnalyzer;
        PhotonAnalyzer*   photonAnalyzer;
        JetAnalyzer*      jetAnalyzer;
        LheAnalyzer*      lheAnalyzer;
        GenAnalyzer*      genAnalyzer;
        SUSYMassAnalyzer* susyMassAnalyzer;

        edm::Service<TFileService> fs;                                                                   //Root tree and file for storing event info
        TTree* outputTree;

        unsigned long _runNb;
        unsigned long _lumiBlock;
        unsigned long _eventNb;
        unsigned      _nVertex;
        double _BS_x;
        double _BS_y;
        double _BS_z;
        double _PV_x;
        double _PV_y;
        double _PV_z;

        TH1D* nVertices;                                                                                 //Histogram with number of vertices
};
#endif

