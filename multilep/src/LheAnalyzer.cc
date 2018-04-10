#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <math.h> 
#include <algorithm>

/*
 * Accessing LHE information
 * lheHTIncoming could be used to get the low HT bin for HT-binned samples, e.g. DY
 * Might consider skimming directly here for such samples
 * Also saving the ctau of the heavy neutrino
 */
LheAnalyzer::LheAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{};


void LheAnalyzer::beginJob(TTree* outputTree, edm::Service<TFileService>& fs){
    //Counter to determine effect of pdf and scale uncertainties on the MC cross section
    if(multilepAnalyzer->isData) return;
    hCounter   = fs->make<TH1D>("hCounter",   "Events counter",    1,  0,1);
    lheCounter = fs->make<TH1D>("lheCounter", "Lhe weights",       110,0,110);
    nTrueInteractions = fs->make<TH1D>("nTrueInteractions",  "nTrueInteractions", 120,0,120);
    outputTree->Branch("_nTrueInt",      &_nTrueInt,      "_nTrueInt/F");
    outputTree->Branch("_weight",        &_weight,        "_weight/D");
    outputTree->Branch("_lheHTIncoming", &_lheHTIncoming, "_lheHTIncoming/D");
    outputTree->Branch("_ctauHN",        &_ctauHN,        "_ctauHN/D");
    outputTree->Branch("_nLheWeights",   &_nLheWeights,   "_nLheWeights/b");
    outputTree->Branch("_lheWeight",     &_lheWeight,     "_lheWeight[_nLheWeights]/D");
}

void LheAnalyzer::analyze(const edm::Event& iEvent){
    if(multilepAnalyzer->isData) return;

    edm::Handle<GenEventInfoProduct> genEventInfo;          iEvent.getByToken(multilepAnalyzer->genEventInfoToken, genEventInfo);
    edm::Handle<LHEEventProduct> lheEventInfo;              iEvent.getByToken(multilepAnalyzer->lheEventInfoToken, lheEventInfo);
    edm::Handle<std::vector<PileupSummaryInfo>> pileUpInfo; iEvent.getByToken(multilepAnalyzer->pileUpToken,       pileUpInfo);

    _nTrueInt = pileUpInfo->begin()->getTrueNumInteractions(); // getTrueNumInteractions is the same for all bunch crossings
    _weight   = genEventInfo->weight();
    hCounter->Fill(0.5,    _weight);
    nTrueInteractions->Fill(_nTrueInt, _weight);

    _lheHTIncoming = 0.;
    _ctauHN = 0.;

    if(!lheEventInfo.isValid()){
        _nLheWeights = 0;
        return;
    }

    // See http://home.thep.lu.se/~leif/LHEF/classLHEF_1_1HEPEUP.html for more detailes
    int nParticles = lheEventInfo->hepeup().NUP;

    for(int i = 0; i < nParticles; ++i){
        int  status            = lheEventInfo->hepeup().ISTUP[i];
        long pdgId             = lheEventInfo->hepeup().IDUP[i];
        double px              = lheEventInfo->hepeup().PUP[i][0];
        double py              = lheEventInfo->hepeup().PUP[i][1];
        double pt              = std::sqrt(px*px+py*py); 
        int mother1            = lheEventInfo->hepeup().MOTHUP[i].first-1;                                                 // MOTHUP starts counting at 1
        int mother2            = lheEventInfo->hepeup().MOTHUP[i].second-1;
        bool hasIncomingMother = lheEventInfo->hepeup().ISTUP[mother1]==-1 and lheEventInfo->hepeup().ISTUP[mother2]==-1;  // Status -1 means mother is incoming
        bool quarkOrGluon      = (pdgId==21 or (abs(pdgId)>0 and abs(pdgId) < 7));

        if(hasIncomingMother and status==1 and quarkOrGluon) _lheHTIncoming += pt;
        if(pdgId==9900012)                                   _ctauHN         = lheEventInfo->hepeup().VTIMUP[i];
    }

    //Store LHE weights to compute pdf and scale uncertainties, as described on https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
    _nLheWeights = std::min(110, (int) lheEventInfo->weights().size()); // 110 for MC@NLO, 254 for powheg, 446(!) for madgraph, 0 for some old samples,... 
    for(unsigned i = 0; i < _nLheWeights; ++i){
        _lheWeight[i] = lheEventInfo->weights()[i].wgt/lheEventInfo->originalXWGTUP(); 
        lheCounter->Fill(i + 0.5, _lheWeight[i]*_weight);
    }
}

double LheAnalyzer::getWeight() const{
    if(multilepAnalyzer->isData) return 1.;
    return _weight;
}

