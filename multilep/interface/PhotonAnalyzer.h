#ifndef PHOTON_ANALYZER_H
#define PHOTON_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class PhotonAnalyzer {
  private:
    static const unsigned nPhoton_max = 5;

    unsigned _nPhoton;
    double   _photonPt[nPhoton_max];
    double   _photonEta[nPhoton_max];
    double   _photonPhi[nPhoton_max];
    double   _photonE[nPhoton_max];
    bool     _photonCutBasedLoose[nPhoton_max];
    bool     _photonCutBasedMedium[nPhoton_max];
    bool     _photonCutBasedTight[nPhoton_max];
    double   _photonMva[nPhoton_max];
    double   _photonChargedIsolation[nPhoton_max];
    double   _photonNeutralHadronIsolation[nPhoton_max];
    double   _photonPhotonIsolation[nPhoton_max];
    double   _photonSigmaIetaIeta[nPhoton_max];
    double   _photonHadronicOverEm[nPhoton_max];
    bool     _photonPassElectronVeto[nPhoton_max];
    bool     _photonHasPixelSeed[nPhoton_max];

    multilep* multilepAnalyzer;

  public:
    PhotonAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~PhotonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&);
};

#endif
