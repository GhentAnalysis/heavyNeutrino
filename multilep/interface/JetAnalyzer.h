#ifndef JET_ANALYZER_H
#define JET_ANALYZER_H

//c++ standard library
#include <memory>


//CMSSW
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

//ROOT
#include "TTree.h"

//other parts of code
#include "heavyNeutrino/multilep/plugins/multilep.h"

class multilep;

class JetAnalyzer {
    friend class multilep;
    private:
        JetCorrectionUncertainty* jecUnc;
        JetCorrectionUncertainty* jecUncPuppi;

        static const unsigned nJets_max = 100;
        static const unsigned nPFCandidates_max = 1000;

        std::map<std::string, std::shared_ptr< JetCorrectorParameters> > jetSourcesCorParameters;
        std::map<std::string, std::shared_ptr< JetCorrectorParameters> > jetGroupedCorParameters;

        std::shared_ptr<FactorizedJetCorrector> jetCorrector;

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
        unsigned _jetPartonFlavor[nJets_max];
        bool     _jetIsLoose[nJets_max];
        bool     _jetIsTight[nJets_max];
        bool     _jetIsTightLepVeto[nJets_max];
        double   _jetChargedHadronFraction[nJets_max];
        double   _jetNeutralHadronFraction[nJets_max];
        double   _jetNeutralEmFraction[nJets_max];
        double   _jetChargedEmFraction[nJets_max];
        double   _jetHFHadronFraction[nJets_max];
        double   _jetHFEmFraction[nJets_max];
        double   _jetPileupIdFullDisc[nJets_max];
        int      _jetPileupIdFullId[nJets_max];
        bool     _jetHasGen[nJets_max];
        double   _jetGenPt[nJets_max];
        double   _jetGenEta[nJets_max];
        double   _jetGenPhi[nJets_max];
        double   _jetGenE[nJets_max];

        // substructure variables for Higgs analysis
        double   _jetNsubTau1[nJets_max];
        double   _jetNsubTau2[nJets_max];
        double   _jetNsubTau3[nJets_max];
        double   _jetQGLikelihood[nJets_max];
        unsigned _jetNPFCandidates[nJets_max];
        unsigned _nPFCandidates = 0;
        unsigned _pfCandidateJetIndex[nPFCandidates_max];
        double   _pfCandidatePt[nPFCandidates_max];
        double   _pfCandidateEta[nPFCandidates_max];
        double   _pfCandidatePhi[nPFCandidates_max];
        double   _pfCandidateE[nPFCandidates_max];
        int      _pfCandidateCharge[nPFCandidates_max];
        int      _pfCandidatePdgId[nPFCandidates_max];
        int      _pfCandidateLostInnerHits[nPFCandidates_max];
        bool     _pfCandidateTrackHighPurity[nPFCandidates_max];
        int      _pfCandidatePvAssociationQuality[nPFCandidates_max];
        double   _pfCandidateDxy[nPFCandidates_max];
        double   _pfCandidateDz[nPFCandidates_max];
        double   _pfCandidateDxyError[nPFCandidates_max];
        double   _pfCandidateDzError[nPFCandidates_max];

        unsigned _nJetsPuppi = 0;
        double   _jetPuppiPt[nJets_max];
        double   _jetPuppiPt_JECUp[nJets_max];
        double   _jetPuppiPt_JECDown[nJets_max];
        double   _jetPuppiEta[nJets_max];
        double   _jetPuppiPhi[nJets_max];

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

        //split JEC in different sources
        std::map< std::string, double[nJets_max] > _jetPt_groupedVariationsDown;
        std::map< std::string, double[nJets_max] > _jetPt_groupedVariationsUp;
        std::map< std::string, double[nJets_max] > _jetPt_allVariationsDown;
        std::map< std::string, double[nJets_max] > _jetPt_allVariationsUp;

        std::map< std::string, double[nJets_max] > _jetSmearedPt_groupedVariationsDown;
        std::map< std::string, double[nJets_max] > _jetSmearedPt_groupedVariationsUp;
        std::map< std::string, double[nJets_max] > _jetSmearedPt_allVariationsDown;
        std::map< std::string, double[nJets_max] > _jetSmearedPt_allVariationsUp;

        // MET with propagated JEC sources and uncertainties
        std::map< std::string, double > _corrMETx_groupedVariationsDown;
        std::map< std::string, double > _corrMETx_groupedVariationsUp;
        std::map< std::string, double > _corrMETy_groupedVariationsDown;
        std::map< std::string, double > _corrMETy_groupedVariationsUp;

        std::map< std::string, double > _corrMETx_allVariationsDown;
        std::map< std::string, double > _corrMETx_allVariationsUp;
        std::map< std::string, double > _corrMETy_allVariationsDown;
        std::map< std::string, double > _corrMETy_allVariationsUp;

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

        std::vector<float> getSubCorrections(double rawPt, double eta, double rho, double area);
        std::pair<double, double> getMETCorrectionPxPy(double corrPt, double rawPt, double rawEta, double rawMuonSubtractedPt, double phi, double emf, double rho, double area, const std::string& source, unsigned jetIndex, double jecShift);

    public:
        JetAnalyzer(const edm::ParameterSet& iConfig, multilep* vars);
        ~JetAnalyzer();

        void beginJob(TTree* outputTree);
        bool analyze(const edm::Event&);

        void correctedMETAndPhi(const pat::MET& met, const std::vector< pat::Jet >& jets, const double rho);
};

#endif
