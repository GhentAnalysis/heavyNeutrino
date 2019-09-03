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
    JetCorrectionUncertainty* jecUnc;

    static const unsigned nJets_max = 20;
    //static const unsigned nDaughters_max = 500;

    bool is2017;

    unsigned _nJets = 0;
    double   _jetPt[nJets_max];
    double   _jetPt_JECUp[nJets_max];
    double   _jetPt_JECDown[nJets_max];
    double   _jetSmearedPt[nJets_max];
    double   _jetSmearedPt_JECDown[nJets_max];
    double   _jetSmearedPt_JECUp[nJets_max];
    double   _jetSmearedPt_JERDown[nJets_max];
    double   _jetSmearedPt_JERUp[nJets_max];
    double   _jetPt_Uncorrected[nJets_max];
    double   _jetPt_L1[nJets_max];
    double   _jetPt_L2[nJets_max];
    double   _jetPt_L3[nJets_max];
    double   _jetEta[nJets_max];
    double   _jetPhi[nJets_max];
    double   _jetE[nJets_max];
    double   _jetCsvV2[nJets_max];
    double   _jetDeepCsv_udsg[nJets_max];
    double   _jetDeepCsv_b[nJets_max];
    double   _jetDeepCsv_c[nJets_max];
    double   _jetDeepCsv_bb[nJets_max];
    double   _jetDeepCsv[nJets_max];
    double   _jetDeepFlavor_b[nJets_max];
    double   _jetDeepFlavor_bb[nJets_max];
    double   _jetDeepFlavor_lepb[nJets_max];
    double   _jetDeepFlavor[nJets_max];
    double   _jetDeepFlavor_c[nJets_max];
    double   _jetDeepFlavor_uds[nJets_max];
    double   _jetDeepFlavor_g[nJets_max];
    unsigned _jetHadronFlavor[nJets_max];
    bool     _jetIsLoose[nJets_max];
    bool     _jetIsTight[nJets_max];
    bool     _jetIsTightLepVeto[nJets_max];
    double   _jetChargedHadronFraction[nJets_max];
    double   _jetNeutralHadronFraction[nJets_max];
    double   _jetNeutralEmFraction[nJets_max];
    double   _jetChargedEmFraction[nJets_max];
    double   _jetHFHadronFraction[nJets_max];
    double   _jetHFEmFraction[nJets_max];

    //to train HNL Jettagger MVA based on the deep sets theorem
    static const unsigned maxJetSize = 50;
    unsigned _nJetConstituents[nJets_max];
    double   _JetConstituentPt[nJets_max][maxJetSize];
    double   _JetConstituentEta[nJets_max][maxJetSize];
    double   _JetConstituentPhi[nJets_max][maxJetSize];
    double   _JetConstituentMass[nJets_max][maxJetSize];
    int      _JetConstituentPdgId[nJets_max][maxJetSize];
    double   _JetConstituentPdgIdReduced[nJets_max][maxJetSize];
    int      _JetConstituentCharge[nJets_max][maxJetSize];
    double   _JetConstituentdxy[nJets_max][maxJetSize];
    double   _JetConstituentdz[nJets_max][maxJetSize];
    double   _JetConstituentdxyErr[nJets_max][maxJetSize];
    double   _JetConstituentdzErr[nJets_max][maxJetSize];
    int      _JetConstituentNumberOfHits[nJets_max][maxJetSize];    //signed integer so particles without tracks can have negative default values
    int      _JetConstituentNumberOfPixelHits[nJets_max][maxJetSize];
    bool     _JetConstituentHasTrack[nJets_max][maxJetSize];
    
    //correction level for JEC
    //std::string jecLevel;

    double   _met;                                                                              //met kinematics
    double   _metPhi;
    double   _metRaw;
    double   _metRawPhi;
    double   _metJECDown;
    double   _metPhiJECDown;
    double   _metJECUp;
    double   _metPhiJECUp;
    double   _metUnclDown;
    double   _metPhiUnclDown;
    double   _metUnclUp;
    double   _metPhiUnclUp;
    double   _metSignificance;

    multilep* multilepAnalyzer;

    bool jetIsLoose(const pat::Jet& jet, const bool is2017) const;
    bool jetIsTight(const pat::Jet& jet, const bool is2017, const bool is2018) const;
    bool jetIsTightLepVeto(const pat::Jet& jet, const bool is2017, const bool is2018) const;
    double reducedPdgId(int pdgId);
    double catchNanOrInf(double value);

  public:
    JetAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~JetAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&);
};

#endif
