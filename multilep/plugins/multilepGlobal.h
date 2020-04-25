#ifndef MULTILEPGLOBAL_H
#define MULTILEPGLOBAL_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TLorentzVector.h"

// Allow for easy way to retrieve handles
namespace {
  template<typename T, typename I> edm::Handle<T> getHandle(const I& iEvent,const edm::EDGetTokenT<T>& token){
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);
    return handle;
  }
}


class multilepGlobal : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit multilepGlobal(const edm::ParameterSet&);
        ~multilepGlobal(){};
        
        bool isData() const{ return sampleIsData; }

    private:
        edm::EDGetTokenT<std::vector<reco::Vertex>>              vtxToken;
        edm::EDGetTokenT<GenEventInfoProduct>                    genEventInfoToken;
        edm::EDGetTokenT<LHEEventProduct>                        lheEventInfoToken;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>>         pileUpToken;

        //void beginJob(edm::Service<TFileService>& fs);
        void analyze(const edm::Event&, const edm::EventSetup&);
        
        edm::Service<TFileService> fs;                                                                   //Root tree and file for storing event info
        
        bool   sampleIsData;
        TH1D*  nVertices;                                                                                 //Histogram with number of vertices
        TH1D*  hCounter;
        TH1D*  lheCounter;
        TH1D*  psCounter;
        TH1D*  tauCounter;
        TH1D*  HTCounter;
        TH1D*  nTrueInteractions;
};
#endif

