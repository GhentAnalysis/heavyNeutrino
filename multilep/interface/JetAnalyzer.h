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
    JetCorrectionUncertainty* jecUncPuppi;

    static const unsigned nJets_max = 20;

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

    unsigned _nJetsPuppi = 0;
    double   _jetPuppiPt[nJets_max];
    double   _jetPuppiPt_JECUp[nJets_max];
    double   _jetPuppiPt_JECDown[nJets_max];
    double   _jetPuppiEta[nJets_max];
    double   _jetPuppiPhi[nJets_max];

    double   _met;                                                                              //met kinematics
    double   _metPhi;
    double   _metType1;
    double   _metType1Phi;
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
    double   _metResDown;
    double   _metPhiResDown;
    double   _metResUp;
    double   _metPhiResUp;
    double   _metSignificance;

    double   _metPuppi;                                                                              //metPuppi kinematics
    double   _metPuppiPhi;
    double   _metPuppiRaw;
    double   _metPuppiRawPhi;
    double   _metPuppiJECDown;
    double   _metPuppiPhiJECDown;
    double   _metPuppiJECUp;
    double   _metPuppiPhiJECUp;
    double   _metPuppiUnclDown;
    double   _metPuppiPhiUnclDown;
    double   _metPuppiUnclUp;
    double   _metPuppiPhiUnclUp;
    double   _metPuppiResDown;
    double   _metPuppiPhiResDown;
    double   _metPuppiResUp;
    double   _metPuppiPhiResUp;

    std::string jecLevel;

    multilep* multilepAnalyzer;

    bool jetIsLoose(const pat::Jet& jet, const bool is2017) const;
    bool jetIsTight(const pat::Jet& jet, const bool is2017, const bool is2018) const;
    bool jetIsTightLepVeto(const pat::Jet& jet, const bool is2017, const bool is2018) const;

  public:
    JetAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
    ~JetAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const int);
};

#endif
