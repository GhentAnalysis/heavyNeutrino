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

//include vertices
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//include muons
#include "DataFormats/PatCandidates/interface/Muon.h"
//include electrons
#include "DataFormats/PatCandidates/interface/Electron.h"
//include taus
#include "DataFormats/PatCandidates/interface/Tau.h"
//include MET
#include "DataFormats/PatCandidates/interface/MET.h"
//include Jets
#include "DataFormats/PatCandidates/interface/Jet.h"
//include JEC
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
//include Trigger classes
#include "FWCore/Common/interface/TriggerNames.h"
//include other parts of the code
#include "heavyNeutrino/multilep/interface/Tools.h"


//include root classes
#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class multilep : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit multilep(const edm::ParameterSet&);
        ~multilep();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
		
		//Define EDgetTokens to read data from events
		edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;	//vertices
        edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;	//muon
		edm::EDGetTokenT<std::vector<pat::Electron> > eleToken_; //electron
		edm::EDGetTokenT<std::vector<pat::Tau> > tauToken_;	//tau
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidatesToken_;	//particle collection used to calculate isolation variables
		edm::EDGetTokenT<double> rhoToken_; //neutal energy density in terms of deltaR used for pileup corrections
		edm::EDGetTokenT<double> rhoTokenAll_; //jet energy density used for JEC
		edm::EDGetTokenT<std::vector<pat::MET> > metToken_;	//missing transverse energy
		edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_; //jet collection
		//edm::EDGetTokenT<reco::JetCorrector> jecToken_; //JEC
		edm::EDGetTokenT<edm::TriggerResults> triggerToken_; //Trigger information
		edm::EDGetTokenT<edm::TriggerResults> recoResultsToken_; //MET filter information
		edm::EDGetTokenT<bool> badPFMuonFilterToken_; //MET filter not stored in miniAOD
		edm::EDGetTokenT<bool> badChCandFilterToken_;  //MET filter not stored in miniAOD

		//Root tree and file for storing event info
		edm::Service<TFileService> fs;
		//FILE* outFile;
		TTree* outputTree;
		//maximum number of particles stored
		static const unsigned nL_max = 20;
		static const unsigned nJets_max = 20;
		//event labels
		unsigned long _runNb;
		unsigned long _lumiBlock;
		unsigned long _eventNb;
		//number of leptons
		unsigned _nL;
		unsigned _nMu;
		unsigned _nEle;
		unsigned _nLight;
		unsigned _nTau;
		//lepton kinematics
		double _lPt[nL_max];
		double _lEta[nL_max];
		double _lPhi[nL_max];
		double _lE[nL_max];
		//other lepton variables
		unsigned _flavor[nL_max];
		int _charge[nL_max];
		double _relIso[nL_max];
		double _miniIso[nL_max];
		double _dxy[nL_max];
		double _dz[nL_max];
		double _3dIP[nL_max];
		double _3dIPSig[nL_max];
		//lepton selection decisions
		bool _loose[nL_max];
		bool _FO[nL_max];
		bool _tight[nL_max];
		//generator lepton variables
		bool _isPrompt[nL_max];
		bool _truthPdg[nL_max];
		bool _truthMomPdg[nL_max];
		unsigned _origin[nL_max];
		//jet variables
		unsigned _nJets;
		double _jetPt[nJets_max];
		double _jetEta[nJets_max];
		double _jetPhi[nJets_max];
		double _jetE[nJets_max];		
		//met kinematics
		double _met;
		double _metPhi;
		//Event variables
		unsigned _nVertex;
		bool _metFiltersFlagged;
		bool _passHnlTrigger[4]; //0 = eee, 1 = eem, 2 = emm, 3 = mmm
		bool _badMuonFlagged;
		bool _badCloneMuonFlagged;

		//Additional class functions
		void fillTriggerVars(const edm::Event&);
		bool trigPass(unsigned, edm::Handle<edm::TriggerResults>&, const edm::Event&);
		bool metFilterFlagged(); //Events flagged by a met filter should not be used
		void fillLeptonGenVars(const reco::GenParticle* genParticle); //Fill MC-truth lepton variables
		void fillLeptonKinVars(const reco::Candidate&); //Fill reconstructed lepton kinematics
		template <class leptonType> void fillLeptonIdVars(const leptonType&); //Fill basic lepton ID variables 
		void fillMetFilterVars(const edm::Event&); //Make MET filter decision
		bool metFiltersFlagged(const edm::Event&); //Fill MET filter variables
};
#endif

