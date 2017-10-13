#include "heavyNeutrino/multilep/interface/SUSYMassAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
SUSYMassAnalyzer::SUSYMassAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
      multilepAnalyzer(multilepAnalyzer)
{};


void SUSYMassAnalyzer::beginJob(TTree* outputTree, edm::Service<TFileService>& fs){
    if(multilepAnalyzer->isSUSY) return;    //only run this module on SUSY samples
    //Counter to determine the amount of events for every SUSY mass point
    //Note, too small binning is used to be sure the binning is smaller than the sample's mass point separation
    //There is no way to access the amount of mass points or there splitting while running over the sample!!!
    hCounterSUSY = fs->make<TH2D>("hCounterSUSY",   "SUSY Events counter", 400,0,2000, 300, 0, 1500);
    //Store SUSY particle masses for event
    outputTree->Branch("_mChi1", &_mChi1, "_mChi1/D");
    outputTree->Branch("_mChi2", &_mChi2, "_mChi2/D");
}

void SUSYMassAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iEventSetup){
    //Extract mass corresponding to this lumi block.
    edm::Handle<GenLumiInfoHeader> genHeader;
    iLumi.getByToken(multilepAnalyzer->genLumiInfoToken, genHeader);
    std::string model = genHeader->configDescription(); 
    std::cout<< model << std::endl;

}

void SUSYMassAnalyzer::analyze(const edm::Event& iEvent){


}
