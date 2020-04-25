#ifndef MULTILEPFILTER_H
#define MULTILEPFILTER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"

// Allow for easy way to retrieve handles
namespace {
  template<typename T, typename I> edm::Handle<T> getHandle(const I& iEvent,const edm::EDGetTokenT<T>& token){
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);
    return handle;
  }
}


class multilepFilter : public edm::one::EDFilter<edm::one::SharedResources> {
    public:
        explicit multilepFilter(const edm::ParameterSet&);
        ~multilepFilter(){};

    private:
        edm::EDGetTokenT<std::vector<reco::Vertex>>              vtxToken;
        edm::EDGetTokenT<std::vector<pat::Muon>>                 muonToken;
        edm::EDGetTokenT<std::vector<pat::Electron>>             eleToken;
        edm::EDGetTokenT<edm::TriggerResults>                    triggerToken;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjToken;
        bool                                                     sampleIs2017;
        bool                                                     sampleIs2018;

        std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjects(const pat::Electron& ele, const std::vector<pat::TriggerObjectStandAlone>& trigObjs, const float maxDeltaR);
        bool passTrigger(const edm::Event& iEvent, std::string name);
        bool filter(edm::Event&, const edm::EventSetup&) override;
};
#endif

