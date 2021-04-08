#include "heavyNeutrino/multilep/plugins/multilepFilter.h"
#include "heavyNeutrino/multilep/interface/Header.h"


multilepFilter::multilepFilter(const edm::ParameterSet& iConfig):
    vtxToken(                 consumes<std::vector<reco::Vertex>>(             iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken(                consumes<std::vector<pat::Muon>>(                iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                 consumes<std::vector<pat::Electron>>(            iConfig.getParameter<edm::InputTag>("electrons"))),
    triggerToken(             consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("triggers"))),
    trigObjToken(             consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    sampleIs2017(                                                              iConfig.getUntrackedParameter<bool>("is2017")),
    sampleIs2018(                                                              iConfig.getUntrackedParameter<bool>("is2018")),
    skim(                                                                      iConfig.getUntrackedParameter<std::string>("skim"))
{
}

// Some quick code duplication form TriggerAnalyzer here...
std::vector<const pat::TriggerObjectStandAlone*> multilepFilter::getMatchedObjects(const pat::Electron& ele, const std::vector<pat::TriggerObjectStandAlone>& trigObjs, const float maxDeltaR){
  std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
  const float maxDR2 = maxDeltaR*maxDeltaR;
  for(auto& trigObj : trigObjs){
    const float dR2 = reco::deltaR2(ele.superCluster()->eta(), ele.superCluster()->phi(), trigObj.eta(), trigObj.phi());
    if(dR2<maxDR2) matchedObjs.push_back(&trigObj);
  }
  return matchedObjs;
}

bool multilepFilter::passTrigger(const edm::Event& iEvent, std::string name){
  auto triggerResults = getHandle(iEvent, triggerToken);

  if(sampleIs2017 and name=="HLT_Ele32_WPTight_Gsf"){
    // the usual Ele32 emulation for 2017
    auto triggerObjects = getHandle(iEvent, trigObjToken);
    std::vector<pat::TriggerObjectStandAlone> unpackedTriggerObjects;
    for(auto& trigObj : *triggerObjects){
      unpackedTriggerObjects.push_back(trigObj);
      unpackedTriggerObjects.back().unpackFilterLabels(iEvent, *triggerResults);
    }

    auto electrons = getHandle(iEvent, eleToken);
    for(auto& ele : *electrons){
      for(const auto trigObj : getMatchedObjects(ele, unpackedTriggerObjects, 0.1)){
        if(!trigObj->hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter")) continue;
        if(!trigObj->hasFilterLabel("hltEGL1SingleEGOrFilter")) continue;
        return true;
      }
    }
    return false;
  } else {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    for (unsigned i = 0; i < triggerResults->size(); ++i){
      if(TString(triggerNames.triggerName(i)).Contains(name + "_v")){
        return triggerResults->accept(i);
      }
    }
    return false;
  }
} 

// ------------ method called for each event  ------------
bool multilepFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
    if(skim != "dilep") return true; //don't apply this dilepton filter when not running dilep skim

    auto electrons      = getHandle(iEvent, eleToken);
    auto muons          = getHandle(iEvent, muonToken);
    auto vertices       = getHandle(iEvent, vtxToken);
    if(vertices->size() == 0) return false;          // Don't consider 0 vertex events



    // Trigger
    std::vector<std::string> triggers;
    if(sampleIs2018 or sampleIs2017) triggers = {"HLT_Ele32_WPTight_Gsf", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_IsoMu24_eta2p1"};
    else                             triggers = {"HLT_Ele27_WPTight_Gsf", "HLT_IsoMu24", "HLT_IsoTkMu24"};

    bool passedTrigger = false;
    for(auto triggerName : triggers){
      passedTrigger = passedTrigger or passTrigger(iEvent, triggerName);
    }
    if(not passedTrigger) return false;



    int _nLight = 0;
    //loop over muons
    for(const pat::Muon& mu : *muons){
        if(mu.innerTrack().isNull())                   continue;
        if(mu.pt() < 5)                                continue;
        if(fabs(mu.eta()) > 2.4)                       continue;
        if(!mu.isPFMuon())                             continue;
        if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue;
        ++_nLight;
    }

    // Loop over electrons
    for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
        if(ele->gsfTrack().isNull())                   continue;
        if(ele->pt() < 7)                              continue;
        if(fabs(ele->eta()) > 2.5)                     continue;
        ++_nLight;
    }

    return _nLight >= 2;
}

//define this as a plug-in
DEFINE_FWK_MODULE(multilepFilter);
