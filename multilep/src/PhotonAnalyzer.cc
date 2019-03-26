#include "heavyNeutrino/multilep/interface/PhotonAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "heavyNeutrino/multilep/interface/GenTools.h"
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
    outputTree->Branch("_nPh",                                &_nPh,                            "_nPh/i");
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
    outputTree->Branch("_phPhotonIsolation",                  &_phPhotonIsolation,              "_phPhotonIsolation[_nPh]/D");
    outputTree->Branch("_phSigmaIetaIeta",                    &_phSigmaIetaIeta,                "_phSigmaIetaIeta[_nPh]/D");
    outputTree->Branch("_phHadronicOverEm",                   &_phHadronicOverEm,               "_phHadronicOverEm[_nPh]/D");
    outputTree->Branch("_phPassElectronVeto",                 &_phPassElectronVeto,             "_phPassElectronVeto[_nPh]/O");
    outputTree->Branch("_phHasPixelSeed",                     &_phHasPixelSeed,                 "_phHasPixelSeed[_nPh]/O");
    if( multilepAnalyzer->isMC() ){
      outputTree->Branch("_phIsPrompt",                       &_phIsPrompt,                     "_phIsPrompt[_nPh]/O");
      outputTree->Branch("_phTTGMatchCategory",               &_phTTGMatchCategory,             "_phTTGMatchCategory[_nPh]/I");
      outputTree->Branch("_phTTGMatchPt",                     &_phTTGMatchPt,                   "_phTTGMatchPt[_nPh]/D");
      outputTree->Branch("_phTTGMatchEta",                    &_phTTGMatchEta,                  "_phTTGMatchEta[_nPh]/D");
      outputTree->Branch("_phMatchPdgId",                     &_phMatchPdgId,                   "_phMatchPdgId[_nPh]/I");
    }
    if( !multilepAnalyzer->is2018() ){
      outputTree->Branch("_phPtCorr",                         &_phPtCorr,                       "_phPtCorr[_nPh]/D");
      outputTree->Branch("_phPtScaleUp",                      &_phPtScaleUp,                    "_phPtScaleUp[_nPh]/D");
      outputTree->Branch("_phPtScaleDown",                    &_phPtScaleDown,                  "_phPtScaleDown[_nPh]/D");
      outputTree->Branch("_phPtResUp",                        &_phPtResUp,                      "_phPtResUp[_nPh]/D");
      outputTree->Branch("_phPtResDown",                      &_phPtResDown,                    "_phPtResDown[_nPh]/D");
      outputTree->Branch("_phECorr",                          &_phECorr,                        "_phECorr[_nPh]/D");
      outputTree->Branch("_phEScaleUp",                       &_phEScaleUp,                     "_phEScaleUp[_nPh]/D");
      outputTree->Branch("_phEScaleDown",                     &_phEScaleDown,                   "_phEScaleDown[_nPh]/D");
      outputTree->Branch("_phEResUp",                         &_phEResUp,                       "_phEResUp[_nPh]/D");
      outputTree->Branch("_phEResDown",                       &_phEResDown,                     "_phEResDown[_nPh]/D");
    }

    generator = TRandom3(0);
}


bool PhotonAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<pat::Photon>> photons;                   iEvent.getByToken(multilepAnalyzer->photonToken,                       photons);
    edm::Handle<std::vector<pat::PackedCandidate>> packedCands;      iEvent.getByToken(multilepAnalyzer->packedCandidatesToken,             packedCands);
    edm::Handle<std::vector<reco::Vertex>> vertices;                 iEvent.getByToken(multilepAnalyzer->vtxToken,                          vertices);
    edm::Handle<std::vector<pat::Electron>> electrons;               iEvent.getByToken(multilepAnalyzer->eleToken,                          electrons);
    edm::Handle<std::vector<pat::Muon>> muons;                       iEvent.getByToken(multilepAnalyzer->muonToken,                         muons);
    edm::Handle<std::vector<pat::Jet>> jets;                         iEvent.getByToken(multilepAnalyzer->jetToken,                          jets);
    edm::Handle<std::vector<reco::GenParticle>> genParticles;        iEvent.getByToken(multilepAnalyzer->genParticleToken,                  genParticles);
    edm::Handle<double> rho;                                         iEvent.getByToken(multilepAnalyzer->rhoToken,                          rho);

    // Loop over photons
    _nPh = 0;
    for(auto photon = photons->begin(); photon != photons->end(); ++photon){
        if(_nPh == nPhoton_max) break;
        const auto photonRef = edm::Ref<std::vector<pat::Photon>>(photons, (photon - photons->begin()));

        if(photon->pt()  < 10)        continue;
        if(fabs(photon->eta()) > 2.5) continue;
        double rhoCorrCharged      = (*rho)*chargedEffectiveAreas.getEffectiveArea(photon->superCluster()->eta());
        double rhoCorrNeutral      = (*rho)*neutralEffectiveAreas.getEffectiveArea(photon->superCluster()->eta());
        double rhoCorrPhotons      = (*rho)*photonsEffectiveAreas.getEffectiveArea(photon->superCluster()->eta());

        double randomConeIsoUnCorr = randomConeIso(photon->superCluster()->eta(), packedCands, *(vertices->begin()), electrons, muons, jets, photons);

        _phPt[_nPh]                         = photon->pt();
        _phEta[_nPh]                        = photon->eta();
        _phEtaSC[_nPh]                      = photon->superCluster()->eta();
        _phPhi[_nPh]                        = photon->phi();
        _phE[_nPh]                          = photon->energy();
        _phCutBasedLoose[_nPh]              = photon->photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
        _phCutBasedMedium[_nPh]             = photon->photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
        _phCutBasedTight[_nPh]              = photon->photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
        _phMva[_nPh]                        = photon->userFloat("PhotonMVAEstimatorRunIIFall17v2Values");

        _phRandomConeChargedIsolation[_nPh] = randomConeIsoUnCorr < 0 ? -1 : std::max(0., randomConeIsoUnCorr - rhoCorrCharged); // keep -1 when randomConeIso algorithm failed
        _phChargedIsolation[_nPh]           = std::max(0., photon->userFloat("phoChargedIsolation") - rhoCorrCharged);
        _phNeutralHadronIsolation[_nPh]     = std::max(0., photon->userFloat("phoNeutralHadronIsolation") - rhoCorrNeutral);
        _phPhotonIsolation[_nPh]            = std::max(0., photon->userFloat("phoPhotonIsolation") - rhoCorrPhotons);

        _phSigmaIetaIeta[_nPh]              = photon->full5x5_sigmaIetaIeta();
        _phHadronicOverEm[_nPh]             = photon->hadronicOverEm();
        _phPassElectronVeto[_nPh]           = photon->passElectronVeto();
        _phHasPixelSeed[_nPh]               = photon->hasPixelSeed();

        // Note: for the scale and smearing systematics we use the overall values, assuming we are not very sensitive to these systematics
        // In case these systematics turn out to be important, need to add their individual source to the tree (and propagate to their own templates):
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
        // Currently only available for 2016/2017
        if( !multilepAnalyzer->is2018() ){
          _phPtCorr[_nPh]                   = photon->pt()*photon->userFloat("ecalEnergyPostCorr")/photon->energy();
          _phPtScaleUp[_nPh]                = photon->pt()*photon->userFloat("energyScaleUp")/photon->energy();
          _phPtScaleDown[_nPh]              = photon->pt()*photon->userFloat("energyScaleDown")/photon->energy();
          _phPtResUp[_nPh]                  = photon->pt()*photon->userFloat("energySigmaUp")/photon->energy();
          _phPtResDown[_nPh]                = photon->pt()*photon->userFloat("energySigmaDown")/photon->energy();
          _phECorr[_nPh]                    = photon->userFloat("ecalEnergyPostCorr");
          _phEScaleUp[_nPh]                 = photon->userFloat("energyScaleUp");
          _phEScaleDown[_nPh]               = photon->userFloat("energyScaleDown");
          _phEResUp[_nPh]                   = photon->userFloat("energySigmaUp");
          _phEResDown[_nPh]                 = photon->userFloat("energySigmaDown");
        }

        if( multilepAnalyzer->isMC() ){
            fillPhotonGenVars(photon->genParticle());
            matchCategory(*photon, genParticles);
        }
        ++_nPh;
    }

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
    bool overlap   = true;
    int attempt    = 0;
    double randomPhi;
    while(overlap and attempt<20){
        randomPhi = generator.Uniform(-TMath::Pi(),TMath::Pi());

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
        if(iCand.hasTrackDetails()){
            if(deltaR(eta, randomPhi, iCand.eta(), iCand.phi()) > 0.3) continue;
            if(abs(iCand.pdgId()) != 211) continue;

            float dxy = iCand.pseudoTrack().dxy(vertex.position());
            float dz  = iCand.pseudoTrack().dz(vertex.position());
            if(fabs(dxy) > 0.1) continue;
            if(fabs(dz) > 0.2)  continue;

            chargedIsoSum += iCand.pt();
        }
    }
    return chargedIsoSum;
}


// Photon matching as used in TOP-18-010, following https://indico.cern.ch/event/686540/contributions/2816395/attachments/1578345/2493189/Dec20_TTGammaChanges.pdf
void PhotonAnalyzer::matchCategory(const pat::Photon& photon, edm::Handle<std::vector<reco::GenParticle>>& genParticles){
    enum matchCategory {UNDEFINED, GENUINE, MISIDELE, HADRONICPHOTON, HADRONICFAKE};
    _phTTGMatchCategory[_nPh] = UNDEFINED;
    _phTTGMatchPt[_nPh]       = -1.;
    _phTTGMatchEta[_nPh]      = -10.;

    float minDeltaR = 999;
    const reco::GenParticle* matched = nullptr;

    for(auto& p : *genParticles){
      if(p.status()!=1 and p.status()!=71)  continue;
      if(fabs(p.pt()-photon.pt())/p.pt() > 0.5) continue;
      float myDeltaR = deltaR(p.eta(), p.phi(), photon.eta(), photon.phi());
      if(myDeltaR > 0.1 or myDeltaR > minDeltaR) continue;
      minDeltaR  = myDeltaR;
      matched    = &p;
    }

    if(matched){
      _phTTGMatchPt[_nPh]  = matched->pt();
      _phTTGMatchEta[_nPh] = matched->eta();
      bool passParentage   = GenTools::passParentage(*matched, *genParticles);
      float minOtherDeltaR = GenTools::getMinDeltaR(*matched, *genParticles);
      if(matched and matched->pdgId() == 22){
        if(passParentage and minOtherDeltaR > 0.2)       _phTTGMatchCategory[_nPh] = GENUINE;
        else                                             _phTTGMatchCategory[_nPh] = HADRONICPHOTON;
      } else if(matched and abs(matched->pdgId())==11){
        if(passParentage and minOtherDeltaR > 0.2)       _phTTGMatchCategory[_nPh] = MISIDELE;
      } else                                             _phTTGMatchCategory[_nPh] = HADRONICFAKE;
    } else                                               _phTTGMatchCategory[_nPh] = HADRONICFAKE;
}
