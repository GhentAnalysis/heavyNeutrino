#ifndef PHOTON_ANALYZER_H
#define PHOTON_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"
#include <TRandom3.h>

class multilep;

class PhotonAnalyzer {
  private:
    EffectiveAreas chargedEffectiveAreas;
    EffectiveAreas neutralEffectiveAreas;
    EffectiveAreas photonsEffectiveAreas;

    static const unsigned nPhoton_max = 5;

    unsigned _nPh;
    double   _phPt[nPhoton_max];
    double   _phEta[nPhoton_max];
    double   _phEtaSC[nPhoton_max];
    double   _phPhi[nPhoton_max];
    double   _phE[nPhoton_max];
    bool     _phCutBasedLoose[nPhoton_max];
    bool     _phCutBasedMedium[nPhoton_max];
    bool     _phCutBasedTight[nPhoton_max];
    double   _phMva[nPhoton_max];
    double   _phRandomConeChargedIsolation[nPhoton_max];
    double   _phChargedIsolation[nPhoton_max];
    double   _phNeutralHadronIsolation[nPhoton_max];
    double   _phPhotonIsolation[nPhoton_max];
    double   _phSigmaIetaIeta[nPhoton_max];
    double   _phSigmaIetaIphi[nPhoton_max];
    double   _phHadronicOverEm[nPhoton_max];
    bool     _phPassElectronVeto[nPhoton_max];
    bool     _phHasPixelSeed[nPhoton_max];
    bool     _phIsPrompt[nPhoton_max];
    int      _phMatchMCPhotonAN15165[nPhoton_max];
    int      _phMatchMCLeptonAN15165[nPhoton_max];
    int      _phMatchPdgId[nPhoton_max];

    void fillPhotonGenVars(const reco::GenParticle*);
    double randomConeIso(double, edm::Handle<std::vector<pat::PackedCandidate>>&, const reco::Vertex&,
                         edm::Handle<std::vector<pat::Electron>>&, edm::Handle<std::vector<pat::Muon>>&,
                         edm::Handle<std::vector<pat::Jet>>&, edm::Handle<std::vector<pat::Photon>>&);
    void matchAN15165(const pat::Photon&, edm::Handle<std::vector<reco::GenParticle>>&);

    multilep* multilepAnalyzer;
    TRandom3  generator;

  public:
    PhotonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~PhotonAnalyzer(){};

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&);
};
#endif
