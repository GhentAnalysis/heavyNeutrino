#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <TRandom3.h>
/*
 * Calculating all photon-related variables
 */


PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  chargedEffectiveAreas((iConfig.getParameter<edm::FileInPath>("photonsChargedEffectiveAreas")).fullPath()),
  neutralEffectiveAreas((iConfig.getParameter<edm::FileInPath>("photonsNeutralEffectiveAreas")).fullPath()),
  photonsEffectiveAreas((iConfig.getParameter<edm::FileInPath>("photonsPhotonsEffectiveAreas")).fullPath()),
  multilepAnalyzer(multilepAnalyzer)
{};


void PhotonAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_nPh",                                &_nPh,                            "_nPh/b");
  outputTree->Branch("_phPt",                               &_phPt,                           "_phPt[_nPh]/D");
  outputTree->Branch("_phEta",                              &_phEta,                          "_phEta[_nPh]/D");
  outputTree->Branch("_phEtaSC",                            &_phEtaSC,                        "_phEtaSC[_nPh]/D");
  outputTree->Branch("_phPhi",                              &_phPhi,                          "_phPhi[_nPh]/D");
  outputTree->Branch("_phE",                                &_phE,                            "_phE[_nPh]/D");
  outputTree->Branch("_phCutBasedLoose",                    &_phCutBasedLoose,                "_phCutBasedLoose[_nPh]/O");
  outputTree->Branch("_phCutBasedMedium",                   &_phCutBasedMedium,               "_phCutBasedMedium[_nPh]/O");
  outputTree->Branch("_phCutBasedTight",                    &_phCutBasedTight,                "_phCutBasedTight[_nPh]/O");
  outputTree->Branch("_phMva",                              &_phMva,                          "_phMva[_nPh]/D");
  outputTree->Branch("_phRandomConeChargedIsolation",       &_phRandomConeChargedIsolation,   "_phRandomConeChargedIsolation[_nPh]/D");
  outputTree->Branch("_phChargedIsolation",                 &_phChargedIsolation,             "_phChargedIsolation[_nPh]/D");
  outputTree->Branch("_phNeutralHadronIsolation",           &_phNeutralHadronIsolation,       "_phNeutralHadronIsolation[_nPh]/D");
  outputTree->Branch("_phPhotonIsolation",                  &_phPhotonIsolation,              "_phPhoton[_nPh]/D");
  outputTree->Branch("_phSigmaIetaIeta",                    &_phSigmaIetaIeta,                "_phSigmaIetaIeta[_nPh]/D");
  outputTree->Branch("_phSigmaIetaIphi",                    &_phSigmaIetaIphi,                "_phSigmaIetaIphi[_nPh]/D");
  outputTree->Branch("_phHadronicOverEm",                   &_phHadronicOverEm,               "_phHadronicOverEm[_nPh]/D");
  outputTree->Branch("_phPassElectronVeto",                 &_phPassElectronVeto,             "_phPassElectronVeto[_nPh]/O");
  outputTree->Branch("_phHasPixelSeed",                     &_phHasPixelSeed,                 "_phHasPixelSeed[_nPh]/O");
  outputTree->Branch("_phIsPrompt",                         &_phIsPrompt,                     "_phIsPrompt[_nPh]/O");
  outputTree->Branch("_phMatchPdgId",                       &_phMatchPdgId,                   "_phMatchPdgId[_nPh]/I");
}


bool PhotonAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<std::vector<pat::Photon>> photons;                   iEvent.getByToken(multilepAnalyzer->photonToken,                       photons);
  edm::Handle<edm::ValueMap<bool>> photonsCutBasedLoose;           iEvent.getByToken(multilepAnalyzer->photonCutBasedLooseToken,          photonsCutBasedLoose);
  edm::Handle<edm::ValueMap<bool>> photonsCutBasedMedium;          iEvent.getByToken(multilepAnalyzer->photonCutBasedMediumToken,         photonsCutBasedMedium);
  edm::Handle<edm::ValueMap<bool>> photonsCutBasedTight;           iEvent.getByToken(multilepAnalyzer->photonCutBasedTightToken,          photonsCutBasedTight);
  edm::Handle<edm::ValueMap<float>> photonsMva;                    iEvent.getByToken(multilepAnalyzer->photonMvaToken,                    photonsMva);
  edm::Handle<edm::ValueMap<float>> photonsChargedIsolation;       iEvent.getByToken(multilepAnalyzer->photonChargedIsolationToken,       photonsChargedIsolation);
  edm::Handle<edm::ValueMap<float>> photonsNeutralHadronIsolation; iEvent.getByToken(multilepAnalyzer->photonNeutralHadronIsolationToken, photonsNeutralHadronIsolation);
  edm::Handle<edm::ValueMap<float>> photonsPhotonIsolation;        iEvent.getByToken(multilepAnalyzer->photonPhotonIsolationToken,        photonsPhotonIsolation);
  edm::Handle<edm::ValueMap<float>> photonsFull5x5SigmaIEtaIPhi;   iEvent.getByToken(multilepAnalyzer->photonFull5x5SigmaIEtaIPhiToken,   photonsFull5x5SigmaIEtaIPhi);
  edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
  edm::Handle<std::vector<reco::Vertex>> vertices;                 iEvent.getByToken(multilepAnalyzer->vtxToken,                          vertices);
  edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
  edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
  edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);
  edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);

  // Loop over photons
  _nPh = 0;
  for(auto photon = photons->begin(); photon != photons->end(); ++photon){
    if(_nPh == nPhoton_max) break;
    const auto photonRef = edm::Ref<std::vector<pat::Photon>>(photons, (photon - photons->begin()));

    if(photon->pt()  < 10)        continue;
    if(fabs(photon->eta()) > 2.5) continue;

    double rhoCorrCharged = (*rho)*chargedEffectiveAreas.getEffectiveArea(photon->superCluster()->eta());
    double rhoCorrNeutral = (*rho)*neutralEffectiveAreas.getEffectiveArea(photon->superCluster()->eta());
    double rhoCorrPhotons = (*rho)*photonsEffectiveAreas.getEffectiveArea(photon->superCluster()->eta());


    _phPt[_nPh]                         = photon->pt();
    _phEta[_nPh]                        = photon->eta();
    _phEtaSC[_nPh]                      = photon->superCluster()->eta();
    _phPhi[_nPh]                        = photon->phi();
    _phE[_nPh]                          = photon->energy();
    _phCutBasedLoose[_nPh]              = (*photonsCutBasedLoose)[photonRef];
    _phCutBasedMedium[_nPh]             = (*photonsCutBasedMedium)[photonRef];
    _phCutBasedTight[_nPh]              = (*photonsCutBasedTight)[photonRef];
    _phMva[_nPh]                        = (*photonsMva)[photonRef];
    _phRandomConeChargedIsolation[_nPh] = std::max(0., randomConeIso(photon->superCluster()->eta(), packedCands, *(vertices->begin()), electrons, muons, jets, photons) - rhoCorrCharged);
    _phChargedIsolation[_nPh]           = std::max(0., (*photonsChargedIsolation)[photonRef] - rhoCorrCharged);
    _phNeutralHadronIsolation[_nPh]     = std::max(0., (*photonsNeutralHadronIsolation)[photonRef] - rhoCorrNeutral);
    _phPhotonIsolation[_nPh]            = std::max(0., (*photonsPhotonIsolation)[photonRef] - rhoCorrPhotons);
    _phSigmaIetaIeta[_nPh]              = photon->full5x5_sigmaIetaIeta();
    _phSigmaIetaIphi[_nPh]              = (*photonsFull5x5SigmaIEtaIPhi)[photonRef];
    _phHadronicOverEm[_nPh]             = photon->hadronicOverEm();
    _phPassElectronVeto[_nPh]           = photon->passElectronVeto();
    _phHasPixelSeed[_nPh]               = photon->hasPixelSeed();
    fillPhotonGenVars(photon->genParticle());

    ++_nPh;
  }

  if(multilepAnalyzer->skim == "ttg" and _nPh < 1) return false;
  if(multilepAnalyzer->skim == "singlephoton" and _nPh < 1) return false;
  if(multilepAnalyzer->skim == "diphoton" and _nPh < 2) return false;
  return true;
}

void PhotonAnalyzer::fillPhotonGenVars(const reco::GenParticle* genParticle){
    if(genParticle != nullptr){
        _phIsPrompt[_nPh]   = (genParticle)->isPromptFinalState();
        _phMatchPdgId[_nPh] = (genParticle)->pdgId();
    } else{
        _phIsPrompt[_nPh]   = false;
        _phMatchPdgId[_nPh] = 0;
    }
}



double PhotonAnalyzer::randomConeIso(double eta, edm::Handle<std::vector<pat::PackedCandidate>>& pfcands, const reco::Vertex& vertex,
                                     edm::Handle<std::vector<pat::Electron>>& electrons, edm::Handle<std::vector<pat::Muon>>& muons,
                                     edm::Handle<std::vector<pat::Jet>>& jets, edm::Handle<std::vector<pat::Photon>>& photons){

  // First, find random phi direction which does not overlap with jets, photons or leptons
  auto generator = new TRandom3(0);
  bool overlap   = true;
  int attempt    = 0;
  double randomPhi;
  while(overlap and attempt<20){
    randomPhi = generator->Uniform(-TMath::Pi(),TMath::Pi());

    overlap = false;
    for(auto& p : *electrons) if(p.pt() > 10 and deltaR(eta, randomPhi, p.eta(), p.phi()) < 0.6) overlap = true;
    for(auto& p : *muons)     if(p.pt() > 10 and deltaR(eta, randomPhi, p.eta(), p.phi()) < 0.6) overlap = true;
    for(auto& p : *jets)      if(p.pt() > 20 and deltaR(eta, randomPhi, p.eta(), p.phi()) < 0.6) overlap = true;
    for(auto& p : *photons)   if(p.pt() > 10 and deltaR(eta, randomPhi, p.eta(), p.phi()) < 0.6) overlap = true;
    ++attempt;
  }
  if(overlap) return -1.;

  // Calculate chargedIsolation
  float chargedIsoSum = 0;
  for(auto& iCand : *pfcands){

    if(deltaR(eta, randomPhi, iCand.eta(), iCand.phi()) > 0.3) continue;
    if(abs(iCand.pdgId()) != 211) continue;

    float dxy = iCand.pseudoTrack().dxy(vertex.position());
    float dz  = iCand.pseudoTrack().dz(vertex.position());
    if(fabs(dxy) > 0.1) continue;
    if(fabs(dz) > 0.2)  continue;

    chargedIsoSum += iCand.pt();
  }
  return chargedIsoSum;
}
