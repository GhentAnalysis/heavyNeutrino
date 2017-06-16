#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <math.h> 

/*
 * Accessing LHE information
 * lheHTIncoming could be used to get the low HT bin for HT-binned samples, e.g. DY
 * Might consider skimming directly here for such samples
 * Also saving the ctau of the heavy neutrino
 */
LheAnalyzer::LheAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
  multilepAnalyzer(multilepAnalyzer)
{};


void LheAnalyzer::beginJob(TTree* outputTree){
  outputTree->Branch("_lheHTIncoming", &_lheHTIncoming, "_lheHTIncoming/D");
  outputTree->Branch("_ctauHN",        &_ctauHN,        "_ctauHN/D");
}

void LheAnalyzer::analyze(const edm::Event& iEvent){
  edm::Handle<LHEEventProduct> lheEventInfo; iEvent.getByToken(multilepAnalyzer->lheEventInfoToken, lheEventInfo); 

  if(!lheEventInfo.isValid()) return;
  
  // See http://home.thep.lu.se/~leif/LHEF/classLHEF_1_1HEPEUP.html for more detailes
  int nParticles = lheEventInfo->hepeup().NUP;

  _lheHTIncoming = 0;
  _ctauHN = 0.;
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
  for(int i = 0; i < 110; ++i){
    _lheWeight[i] = lheEventInfo->weights()[i].wgt/lheEventInfo->originalXWGTUP(); 
  }

}
