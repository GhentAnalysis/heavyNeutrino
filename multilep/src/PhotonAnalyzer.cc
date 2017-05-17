#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/*
 * Calculating all photon-related variables
 */


PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer)
{};


void PhotonAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nPhoton",                      &_nPhoton,                      "_nPhoton/b");
  outputTree->Branch("_photonPt",                     &_photonPt,                     "_photonPt[_nPhoton]/F");
  outputTree->Branch("_photonEta",                    &_photonEta,                    "_photonEta[_nPhoton]/F");
  outputTree->Branch("_photonPhi",                    &_photonPhi,                    "_photonPhi[_nPhoton]/F");
  outputTree->Branch("_photonE",                      &_photonE,                      "_photonE[_nPhoton]/F");
  outputTree->Branch("_photonCutBasedLoose",          &_photonCutBasedLoose,          "_photonCutBasedLoose[_nPhoton]/O");
  outputTree->Branch("_photonCutBasedMedium",         &_photonCutBasedMedium,         "_photonCutBasedMedium[_nPhoton]/O");
  outputTree->Branch("_photonCutBasedLoose",          &_photonCutBasedLoose,          "_photonCutBasedLoose[_nPhoton]/O");
  outputTree->Branch("_photonMva",                    &_photonMva,                    "_photonMva[_nPhoton]/F");
  outputTree->Branch("_photonChargedIsolation",       &_photonChargedIsolation,       "_photonChargedIsolation[_nPhoton]/F");
  outputTree->Branch("_photonNeutralHadronIsolation", &_photonNeutralHadronIsolation, "_photonNeutralHadronIsolation[_nPhoton]/F");
  outputTree->Branch("_photonSigmaIetaIeta",          &_photonSigmaIetaIeta,          "_photonSigmaIetaIeta[_nPhoton]/F");
  outputTree->Branch("_photonHadronicOverEm",         &_photonHadronicOverEm,         "_photonHadronicOverEm[_nPhoton]/F");
  outputTree->Branch("_photonPassElectronVeto",       &_photonPassElectronVeto,       "_photonPassElectronVeto[_nPhoton]/O");
  outputTree->Branch("_photonHasPixelSeed",           &_photonHasPixelSeed,           "_photonHasPixelSeed[_nPhoton]/O");
}

void PhotonAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<pat::Photon>> photons;                   iEvent.getByToken(multilepAnalyzer->photonToken,                       photons);
  edm::Handle<edm::ValueMap<bool>> photonsCutBasedLoose;           iEvent.getByToken(multilepAnalyzer->photonCutBasedLooseToken,          photonsCutBasedLoose);
  edm::Handle<edm::ValueMap<bool>> photonsCutBasedMedium;          iEvent.getByToken(multilepAnalyzer->photonCutBasedMediumToken,         photonsCutBasedMedium);
  edm::Handle<edm::ValueMap<bool>> photonsCutBasedTight;           iEvent.getByToken(multilepAnalyzer->photonCutBasedTightToken,          photonsCutBasedTight);
  edm::Handle<edm::ValueMap<float>> photonsMva;                    iEvent.getByToken(multilepAnalyzer->photonMvaToken,                    photonsMva);
  edm::Handle<edm::ValueMap<float>> photonsChargedIsolation;       iEvent.getByToken(multilepAnalyzer->photonChargedIsolationToken,       photonsChargedIsolation);
  edm::Handle<edm::ValueMap<float>> photonsNeutralHadronIsolation; iEvent.getByToken(multilepAnalyzer->photonNeutralHadronIsolationToken, photonsNeutralHadronIsolation);
  edm::Handle<edm::ValueMap<float>> photonsPhotonIsolation;        iEvent.getByToken(multilepAnalyzer->photonPhotonIsolationToken,        photonsPhotonIsolation);

  // Loop over photons
  _nPhoton = 0;
  for(auto photon = photons->begin(); photon != photons->end(); ++photon){
    auto photonRef = edm::Ref<std::vector<pat::Photon>>(photons, (photon - photons->begin()));

    _photonPt[_nPhoton]                      = photon->pt();
    _photonEta[_nPhoton]                     = photon->eta();
    _photonPhi[_nPhoton]                     = photon->phi();
    _photonE[_nPhoton]                       = photon->energy();
    _photonCutBasedLoose[_nPhoton]           = (*photonsCutBasedLoose)[photonRef];
    _photonCutBasedMedium[_nPhoton]          = (*photonsCutBasedMedium)[photonRef];
    _photonCutBasedTight[_nPhoton]           = (*photonsCutBasedTight)[photonRef];
    _photonMva[_nPhoton]                     = (*photonsMva)[photonRef];
    _photonChargedIsolation[_nPhoton]        = (*photonsChargedIsolation)[photonRef];
    _photonNeutralHadronIsolation[_nPhoton]  = (*photonsNeutralHadronIsolation)[photonRef];
    _photonPhotonIsolation[_nPhoton]         = (*photonsPhotonIsolation)[photonRef];
    _photonSigmaIetaIeta[_nPhoton]           = photon->full5x5_sigmaIetaIeta();
    _photonHadronicOverEm[_nPhoton]          = photon->hadronicOverEm();
    _photonPassElectronVeto[_nPhoton]        = photon->passElectronVeto();
    _photonHasPixelSeed[_nPhoton]            = photon->hasPixelSeed();

    ++_nPhoton;
  }
}
