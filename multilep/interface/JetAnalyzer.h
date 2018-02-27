#ifndef JET_ANALYZER_H
#define JET_ANALYZER_H
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "heavyNeutrino/multilep/plugins/multilep.h"

#include "TTree.h"

class multilep;

class JetAnalyzer {
  friend class multilep;
  private:
    JetCorrectionUncertainty jecUnc;

    static const unsigned nJets_max = 20;
    static const unsigned nDaughters_max = 270;

    unsigned _nJets;
    double   _jetPt[nJets_max];
    double   _jetPt_JECUp[nJets_max];
    double   _jetPt_JECDown[nJets_max];
    double   _jetPt_JERUp[nJets_max];
    double   _jetPt_JERDown[nJets_max];
    double   _jetEta[nJets_max];
    double   _jetPhi[nJets_max];
    double   _jetE[nJets_max];
    double   _jetCsvV2[nJets_max];
    double   _jetDeepCsv_udsg[nJets_max];
    double   _jetDeepCsv_b[nJets_max];
    double   _jetDeepCsv_c[nJets_max];
    double   _jetDeepCsv_bb[nJets_max];
//  double   _jetDeepCsv_cc[nJets_max];
    unsigned _jetHadronFlavor[nJets_max];
    bool    _jetIsLoose[nJets_max];
    bool    _jetIsTight[nJets_max];
    bool    _jetIsTightLepVeto[nJets_max];
    unsigned _nJetswithDaughters;
    Int_t    _nDaughters;
    int      _jet_tag_for_daughters[nDaughters_max];
    int      _jet_daughter_pdgid[nDaughters_max];
    double   _jet_daughter_pt[nDaughters_max];
    double   _jet_daughter_eta[nDaughters_max];
    double   _jet_daughter_phi[nDaughters_max];
    double   _jet_daughter_energy[nDaughters_max];
    
    multilep* multilepAnalyzer;

    bool jetIsLoose(const pat::Jet& jet, const bool is2017) const;
    bool jetIsTight(const pat::Jet& jet, const bool is2017) const;
    bool jetIsTightLepVeto(const pat::Jet& jet, const bool is2017) const;

  public:
    JetAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~JetAnalyzer(){};
    
    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&);
};

#endif
