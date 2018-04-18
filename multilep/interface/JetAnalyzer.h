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

    unsigned _nJets;
    double   _jetPt[nJets_max];
    double   _jetPt_JECUp[nJets_max];
    double   _jetPt_JECDown[nJets_max];
    double   _jetPt_Uncorrected[nJets_max];
    double   _jetPt_L1[nJets_max];
    double   _jetPt_L2[nJets_max];
    double   _jetPt_L3[nJets_max];
    double   _jetPt_L2L3[nJets_max];
    double   _jetEta[nJets_max];
    double   _jetPhi[nJets_max];
    double   _jetE[nJets_max];
    double   _jetCsvV2[nJets_max];
    double   _jetDeepCsv_udsg[nJets_max];
    double   _jetDeepCsv_b[nJets_max];
    double   _jetDeepCsv_c[nJets_max];
    double   _jetDeepCsv_bb[nJets_max];
    unsigned _jetHadronFlavor[nJets_max];
    bool    _jetIsLoose[nJets_max];
    bool    _jetIsTight[nJets_max];
    bool    _jetIsTightLepVeto[nJets_max];
    double _jetChargedHadronFraction[nJets_max];
    double _jetNeutralHadronFraction[nJets_max];
    double _jetNeutralEmFraction[nJets_max];
    double _jetChargedEmFraction[nJets_max];
    double _jetHFHadronFraction[nJets_max];
    double _jetHFEmFraction[nJets_max];

    double        _met;                                                                              //met kinematics
    double        _metJECDown;
    double        _metJECUp;
    double        _metJetResDown;
    double        _metJetResUp;
    double        _metUnclDown;
    double        _metUnclUp;

    double        _met_sm;                                                                              //met kinematics
    double        _metJECDown_sm;
    double        _metJECUp_sm;
    double        _metJetResDown_sm;
    double        _metJetResUp_sm;
    double        _metUnclDown_sm;
    double        _metUnclUp_sm;

    double        _metPhi;
    double        _metPhiJECDown;
    double        _metPhiJECUp;
    double        _metPhiJetResDown;
    double        _metPhiJetResUp;
    double        _metPhiUnclDown;
    double        _metPhiUnclUp;

    double        _metPhi_sm;
    double        _metPhiJECDown_sm;
    double        _metPhiJECUp_sm;
    double        _metPhiJetResDown_sm;
    double        _metPhiJetResUp_sm;
    double        _metPhiUnclDown_sm;
    double        _metPhiUnclUp_sm;

    double        _metSignificance;

    double        _rawmet;
    double        _rawmetJECDown;
    double        _rawmetJECUp;
    double        _rawmetJetResDown;
    double        _rawmetJetResUp;
    double        _rawmetUnclDown;
    double        _rawmetUnclUp;

    double        _rawmetPhi;
    double        _rawmetPhiJECDown;
    double        _rawmetPhiJECUp;
    double        _rawmetPhiJetResDown;
    double        _rawmetPhiJetResUp;
    double        _rawmetPhiUnclDown;
    double        _rawmetPhiUnclUp;

    //correction level for JEC
    //std::string jecLevel;


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
