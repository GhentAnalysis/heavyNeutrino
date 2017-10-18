#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TLegend.h>
#include <TPad.h>
#include "TF1.h"
#include "TF2.h"
#include <TStyle.h>
#include "TLine.h"
#include "TProfile.h"
#include "TAttFill.h"
#include <iostream>
#include <cstring>
#include <string>
#include <TGraphErrors.h>
#include <Riostream.h>
#include "TFile.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <THStack.h>
#include <TPaveText.h>
#include <Analysis_eff.h>
#include "TApplication.h"
#include "TColor.h"
#include <tuple>
#include <set>



//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

//include other parts of the code
#include "tdrstyle.h"
#include "plotCode_new.h"
#include "Selection.h"




using namespace std;




static Double_t pigreco= TMath::ACos(-1);
ClassImp(Analysis_eff)

//_______________________________________________________default constructor_____
Analysis_eff::Analysis_eff():TObject()

{
}
//_______________________________________________________ constructor_____
Analysis_eff::Analysis_eff(string FileNameTree_in):TObject()

{
}
//________________________________________________________________distruttore_____
Analysis_eff::~Analysis_eff()	 {
  // destructor
}


//_______________________________________________________ constructor_____
void Analysis_eff::printProgress(double progress){
  const unsigned barWidth = 100;
  std::cout << "[";
  unsigned pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << unsigned(progress * 100.0) << " %\r" << std::flush;
}



//==================================================================
void Analysis_eff::analisi(int num_histo_kin
			  ){
  
  cout<<"in analisi"<<endl;
  cout<<"---------------------------"<<endl;

  setTDRStyle();
  //gROOT->SetBatch(kTRUE);
 


   
  TGraphAsymmErrors *fakeRate_mu[3];
  TGraphAsymmErrors *fakeRate_mu_04_from_qcd[3];
  TGraphAsymmErrors *fakeRate_e[3];

  TFile hfile2("/Users/Martina/Desktop/CMS/map_FR/fake_rate_mu.root");
  fakeRate_mu[0] = (TGraphAsymmErrors*)hfile2.Get("fakeRate_mu_eta1");
  fakeRate_mu[1] = (TGraphAsymmErrors*)hfile2.Get("fakeRate_mu_eta2");
  fakeRate_mu[2] = (TGraphAsymmErrors*)hfile2.Get("fakeRate_mu_eta3");
 
  TFile hfile1("/Users/Martina/Desktop/CMS/map_FR/fake_rate_e.root");
  fakeRate_e[0] = (TGraphAsymmErrors*)hfile1.Get("fakeRate_e_eta1");
  fakeRate_e[1] = (TGraphAsymmErrors*)hfile1.Get("fakeRate_e_eta2");
  fakeRate_e[2] = (TGraphAsymmErrors*)hfile1.Get("fakeRate_e_eta3");


  const double glugluToZZkFactor = 1;
  const double WZSF = 1;
  const double ZZSF = 1;
  const double XgammaSF = 1;
  const double low_coupling= 0.0001;
  const double low_coupling_2= 0.00001;
  const double high_coupling = 0.001;
  const double coupling_factor_low= low_coupling / 0.0001;
  const double coupling_factor_low_2= low_coupling_2 / 0.0001;

  const double coupling_factor_high= high_coupling ;
  /// 0.0001;
  //low_coupling/ 0.0001
  //TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",             "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",   

  //, "DoubleEG.root " ,"DoubleMuon.root " , "MuonEG.root " ,

  const int nSamples= 28;
  const int nSamples_eff = 18;
  const int nSamples_signal=11;
    
  const TString fileList[nSamples] = { "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", 

				       "1gev.root", "2gev.root","3gev.root","4gev.root","5gev.root","5gev_prompt.root","5_5gev.root","6gev.root","7gev.root","8gev.root","9gev.root",				

				       "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",                  "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",             "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",             "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",   

				       "ZZTo4L_13TeV_powheg_pythia8.root",

				       "WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",      "WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",   "ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",
 
				         "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",

				       "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root", "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",

				       "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",   "WWToLNuQQ_13TeV-powheg.root",            "ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8.root"   };


    
  const TString names[nSamples] = { "total",

				    "trilMaj1", "trilMaj2", "trilMaj3","trilMaj4","trilMaj5","trilMaj5_prompt", "trilMaj5_5", "trilMaj6","trilMaj7", "trilMaj8","trilMaj9",

				    "DY_TTbar",    "DY_TTbar",    "DY_TTbar",    "DY_TTbar",    "DY_TTbar",   
				   
				    "ZZ/H",

				    "triboson", "triboson", "triboson", 
				   
				    "X+gamma",

				    "TT/T + X","TT/T + X","TT/T + X",
				  
				    "diboson", "diboson", "diboson"
  };
    





  const TString eff_names[nSamples_eff +1 ] = {"total",

					       "trilMaj1", "trilMaj2", "trilMaj3","trilMaj4","trilMaj5","trilMaj5_prompt",  "trilMaj5_5", "trilMaj6","trilMaj7", "trilMaj8","trilMaj9",

					       "DY_TTbar",     
				   
					       "ZZ/H",

					       "triboson", 
				   
					       "X+gamma",

					       "TT/T + X",
				  
					       "diboson", "no-prompt"};

  //18610
  // 18610, 1921.8*3, 87.315, 182.175, 182.75,
  const double xSections[nSamples]= { 0,
				     0.5201 * coupling_factor_low,0.5273 * coupling_factor_low, 0.5217 * coupling_factor_low,0.5216 * coupling_factor_low,0.05190 * coupling_factor_low_2, 0.5238 * coupling_factor_low, 0.05220* coupling_factor_low_2, 0.05226* coupling_factor_low_2, 0.05226* coupling_factor_low_2, 0.05193* coupling_factor_low_2, 0.05173* coupling_factor_low_2,


				      18610, 1921.8*3, 87.315, 182.175, 182.75,
				     

				      1.256*ZZSF,

				      0.2086, 0.1651,  0.01398,

				      489*XgammaSF,

				      0.2043, 0.2529, 0.0493, 

				      58.59*WZSF,49.997, 4.0 };




  const TString names_files[nSamples-1] = { 	  "codes_eff/txt_file/trilMaj1.txt", "codes_eff/txt_file/trilMaj2.txt", "codes_eff/txt_file/trilMaj3.txt","codes_eff/txt_file/trilMaj4.txt","codes_eff/txt_file/trilMaj5.txt","codes_eff/txt_file/trilMaj5_prompt.txt",  "codes_eff/txt_file/trilMaj5_5.txt", "codes_eff/txt_file/trilMaj6.txt","codes_eff/txt_file/trilMaj7.txt", "codes_eff/txt_file/trilMaj8.txt","codes_eff/txt_file/trilMaj9.txt",
    
						  "codes_eff/txt_file/DY_TTbar1.txt", "codes_eff/txt_file/DY_TTbar2.txt", "codes_eff/txt_file/DY_TTbar3.txt", "codes_eff/txt_file/DY_TTbar4.txt", "codes_eff/txt_file/DY_TTbar5.txt",
    
						  "codes_eff/txt_file/ZZ_H.txt",
    
						  "codes_eff/txt_file/triboson1.txt","codes_eff/txt_file/triboson2.txt","codes_eff/txt_file/triboson3.txt",
    
						  "codes_eff/txt_file/X_gamma.txt",
    
						  "codes_eff/txt_file/T_X1.txt", "codes_eff/txt_file/T_X2.txt", "codes_eff/txt_file/T_X3.txt",
    
						  "codes_eff/txt_file/diboson1.txt","codes_eff/txt_file/diboson2.txt","codes_eff/txt_file/diboson3.txt"};







				     
  double luminosity = 35.867;

  TFile *hfile[nSamples];
  TTree *inputTree[nSamples];
  //TApplication* rootapp = new TApplication("example",&argc, argv);


 
  // Declaration of leaf types
  ULong64_t       _runNb;
  ULong64_t       _lumiBlock;
  ULong64_t       _eventNb;
  UChar_t         _nVertex;
  Double_t        _weight;
  Double_t        _lheHTIncoming;
  Double_t        _ctauHN;
  UChar_t         _nLheWeights;
  _nLheWeights= 200;
  Double_t        _lheWeight[110];   //[_nLheWeights]
  Float_t         _nTrueInt;
  Int_t           _ttgEventType;
  Int_t           _zgEventType;
  Double_t        _gen_met;
  Double_t        _gen_metPhi;
  UChar_t         _gen_nPh;
  _gen_nPh= 20;
  Double_t        _gen_phPt[10];   //[_gen_nPh]
  Double_t        _gen_phEta[10];   //[_gen_nPh]
  Double_t        _gen_phPhi[10];   //[_gen_nPh]
  Double_t        _gen_phE[10];   //[_gen_nPh]
  Int_t           _gen_phMomPdg[10];   //[_gen_nPh]
  Bool_t          _gen_phIsPrompt[10];   //[_gen_nPh]
  Double_t        _gen_phMinDeltaR[10];   //[_gen_nPh]
  Bool_t          _gen_phPassParentage[10];   //[_gen_nPh]
  UChar_t         _gen_nL;
  _gen_nL= 20;
  Double_t        _gen_lPt[20];   //[_gen_nL]
  Double_t        _gen_lEta[20];   //[_gen_nL]
  Double_t        _gen_lPhi[20];   //[_gen_nL]
  Double_t        _gen_lE[20];   //[_gen_nL]
  Int_t           _gen_lFlavor[20];   //[_gen_nL]
  Int_t           _gen_lCharge[20];   //[_gen_nL]
  Int_t           _gen_lMomPdg[20];   //[_gen_nL]
  Bool_t          _gen_lIsPrompt[20];   //[_gen_nL]
  Bool_t          _passHN_eem;
  Bool_t          _HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  Int_t           _HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale;
  Bool_t          _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _passHN_1l;
  Bool_t          _HLT_Ele27_WPTight_Gsf;
  Int_t           _HLT_Ele27_WPTight_Gsf_prescale;
  Bool_t          _HLT_IsoMu24;
  Int_t           _HLT_IsoMu24_prescale;
  Bool_t          _HLT_IsoTkMu24;
  Int_t           _HLT_IsoTkMu24_prescale;
  Bool_t          _passHN_eee;
  Bool_t          _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
  Int_t           _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;
  Bool_t          _passHN_emm;
  Bool_t          _HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  Int_t           _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Int_t           _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;
  Bool_t          _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Int_t           _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Int_t           _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;
  Bool_t          _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Int_t           _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;
  Bool_t          _passHN_mmm;
  Bool_t          _HLT_TripleMu_12_10_5;
  Int_t           _HLT_TripleMu_12_10_5_prescale;
  Bool_t          _passMET;
  Bool_t          _HLT_PFMET120_PFMHT120_IDTight;
  Int_t           _HLT_PFMET120_PFMHT120_IDTight_prescale;
  Bool_t          _HLT_PFMET110_PFMHT110_IDTight;
  Int_t           _HLT_PFMET110_PFMHT110_IDTight_prescale;
  Bool_t          _HLT_PFMET170_BeamHaloCleaned;
  Int_t           _HLT_PFMET170_BeamHaloCleaned_prescale;
  Bool_t          _HLT_PFMET170_HBHECleaned;
  Int_t           _HLT_PFMET170_HBHECleaned_prescale;
  Bool_t          _HLT_PFMET170_HBHE_BeamHaloCleaned;
  Int_t           _HLT_PFMET170_HBHE_BeamHaloCleaned_prescale;
  Bool_t          _HLT_PFMET170_NotCleaned;
  Int_t           _HLT_PFMET170_NotCleaned_prescale;
  Bool_t          _HLT_MET200;
  Int_t           _HLT_MET200_prescale;
  Bool_t          _HLT_MET250;
  Int_t           _HLT_MET250_prescale;
  Bool_t          _HLT_MET300;
  Int_t           _HLT_MET300_prescale;
  Bool_t          _passMETFilters;
  Bool_t          _Flag_HBHENoiseFilter;
  Bool_t          _Flag_HBHENoiseIsoFilter;
  Bool_t          _Flag_EcalDeadCellTriggerPrimitiveFilter;
  Bool_t          _Flag_goodVertices;
  Bool_t          _Flag_eeBadScFilter;
  Bool_t          _Flag_globalTightHalo2016Filter;
  Bool_t          _flag_badPFMuonFilter;
  Bool_t          _flag_badChCandFilter;
  Bool_t          _passTTG_e;
  Bool_t          _HLT_Ele105_CaloIdVT_GsfTrkIdT;
  Int_t           _HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;
  Bool_t          _HLT_Ele115_CaloIdVT_GsfTrkIdT;
  Int_t           _HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;
  Bool_t          _passTTG_ee;
  Bool_t          _HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  Int_t           _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale;
  Bool_t          _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;
  Int_t           _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale;
  Bool_t          _passTTG_em;
  Bool_t          _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;
  Int_t           _HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale;
  Bool_t          _passTTG_m;
  Bool_t          _HLT_Mu50;
  Int_t           _HLT_Mu50_prescale;
  Bool_t          _HLT_TkMu50;
  Int_t           _HLT_TkMu50_prescale;
  Bool_t          _HLT_Mu45_eta2p1;
  Int_t           _HLT_Mu45_eta2p1_prescale;
  Bool_t          _passTTG_mm;
  Bool_t          _HLT_Mu30_TkMu11;
  Int_t           _HLT_Mu30_TkMu11_prescale;
  UChar_t         _nL;
  _nL= 20;
  UChar_t         _nMu;
  UChar_t         _nEle;
  UChar_t         _nLight;
_nLight= 20;
  UChar_t         _nTau;
  Double_t        _lPt[11];   //[_nL]
  Double_t        _lEta[11];   //[_nL]
  Double_t        _lEtaSC[7];   //[_nLight]
  Double_t        _lPhi[11];   //[_nL]
  Double_t        _lE[11];   //[_nL]
  Int_t           _lFlavor[11];   //[_nL]
  Int_t           _lCharge[11];   //[_nL]
  Double_t        _dxy[11];   //[_nL]
  Double_t        _dz[11];   //[_nL]
  Double_t        _3dIP[11];   //[_nL]
  Double_t        _3dIPSig[11];   //[_nL]
  Float_t         _lElectronMva[7];   //[_nLight]
  Bool_t          _lHNLoose[11];   //[_nL]
  Bool_t          _lHNFO[11];   //[_nL]
  Bool_t          _lHNTight[11];   //[_nL]
  Bool_t          _lPOGVeto[11];   //[_nL]
  Bool_t          _lPOGLoose[11];   //[_nL]
  Bool_t          _lPOGMedium[11];   //[_nL]
  Bool_t          _lPOGTight[11];   //[_nL]
  Bool_t          _lIsPrompt[11];   //[_nL]
  Int_t           _lMatchPdgId[11];   //[_nL]
  Double_t        _relIso[7];   //[_nLight]
  Double_t        _miniIso[7];   //[_nLight]
  Double_t        _ptRel[7];   //[_nLight]
  Double_t        _ptRatio[7];   //[_nLight]
  UChar_t         _nPh;
  _nPh= 20;
  Double_t        _phPt[5];   //[_nPh]
  Double_t        _phEta[5];   //[_nPh]
  Double_t        _phEtaSC[5];   //[_nPh]
  Double_t        _phPhi[5];   //[_nPh]
  Double_t        _phE[5];   //[_nPh]
  Bool_t          _phCutBasedLoose[5];   //[_nPh]
  Bool_t          _phCutBasedMedium[5];   //[_nPh]
  Bool_t          _phCutBasedTight[5];   //[_nPh]
  Double_t        _phMva[5];   //[_nPh]
  Double_t        _phRandomConeChargedIsolation[5];   //[_nPh]
  Double_t        _phChargedIsolation[5];   //[_nPh]
  Double_t        _phNeutralHadronIsolation[5];   //[_nPh]
  Double_t        _phPhotonIsolation[5];   //[_nPh]
  Double_t        _phSigmaIetaIeta[5];   //[_nPh]
  Double_t        _phSigmaIetaIphi[5];   //[_nPh]
  Double_t        _phHadronicOverEm[5];   //[_nPh]
  Bool_t          _phPassElectronVeto[5];   //[_nPh]
  Bool_t          _phHasPixelSeed[5];   //[_nPh]
  Bool_t          _phIsPrompt[5];   //[_nPh]
  Int_t           _phMatchMCPhotonAN15165[5];   //[_nPh]
  Int_t           _phMatchMCLeptonAN15165[5];   //[_nPh]
  Int_t           _phMatchPdgId[5];   //[_nPh]
  UChar_t         _nJets;
  _nJets= 20;
  Double_t        _jetPt[8];   //[_nJets]
  Double_t        _jetPt_JECUp[8];   //[_nJets]
  Double_t        _jetPt_JECDown[8];   //[_nJets]
  Double_t        _jetPt_JERUp[8];   //[_nJets]
  Double_t        _jetPt_JERDown[8];   //[_nJets]
  Double_t        _jetEta[8];   //[_nJets]
  Double_t        _jetPhi[8];   //[_nJets]
  Double_t        _jetE[8];   //[_nJets]
  Double_t        _jetCsvV2[8];   //[_nJets]
  Double_t        _jetDeepCsv_udsg[8];   //[_nJets]
  Double_t        _jetDeepCsv_b[8];   //[_nJets]
  Double_t        _jetDeepCsv_c[8];   //[_nJets]
  Double_t        _jetDeepCsv_bb[8];   //[_nJets]
  Double_t        _jetDeepCsv_cc[8];   //[_nJets]
  Double_t        _jetHadronFlavour[8];   //[_nJets]
  Int_t           _jetId[8];   //[_nJets]
  Double_t        _met;
  Double_t        _metJECDown;
  Double_t        _metJECUp;
  Double_t        _metUnclDown;
  Double_t        _metUnclUp;
  Double_t        _metPhi;
  Double_t        _metPhiJECDown;
  Double_t        _metPhiJECUp;
  Double_t        _metPhiUnclDown;
  Double_t        _metPhiUnclUp;

  // List of branches
  TBranch        *b__runNb;   //!
  TBranch        *b__lumiBlock;   //!
  TBranch        *b__eventNb;   //!
  TBranch        *b__nVertex;   //!
  TBranch        *b__weight;   //!
  TBranch        *b__lheHTIncoming;   //!
  TBranch        *b__ctauHN;   //!
  TBranch        *b__nLheWeights;   //!
  TBranch        *b__lheWeight;   //!
  TBranch        *b__nTrueInt;   //!
  TBranch        *b__ttgEventType;   //!
  TBranch        *b__zgEventType;   //!
  TBranch        *b__gen_met;   //!
  TBranch        *b__gen_metPhi;   //!
  TBranch        *b__gen_nPh;   //!
  TBranch        *b__gen_phPt;   //!
  TBranch        *b__gen_phEta;   //!
  TBranch        *b__gen_phPhi;   //!
  TBranch        *b__gen_phE;   //!
  TBranch        *b__gen_phMomPdg;   //!
  TBranch        *b__gen_phIsPrompt;   //!
  TBranch        *b__gen_phMinDeltaR;   //!
  TBranch        *b__gen_phPassParentage;   //!
  TBranch        *b__gen_nL;   //!
  TBranch        *b__gen_lPt;   //!
  TBranch        *b__gen_lEta;   //!
  TBranch        *b__gen_lPhi;   //!
  TBranch        *b__gen_lE;   //!
  TBranch        *b__gen_lFlavor;   //!
  TBranch        *b__gen_lCharge;   //!
  TBranch        *b__gen_lMomPdg;   //!
  TBranch        *b__gen_lIsPrompt;   //!
  TBranch        *b__passHN_eem;   //!
  TBranch        *b__passHN_1l;   //!
  TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
  TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_Ele27_WPTight_Gsf;   //!
  TBranch        *b__HLT_Ele27_WPTight_Gsf_prescale;   //!
  TBranch        *b__HLT_IsoMu24;   //!
  TBranch        *b__HLT_IsoMu24_prescale;   //!
  TBranch        *b__HLT_IsoTkMu24;   //!
  TBranch        *b__HLT_IsoTkMu24_prescale;   //!
  TBranch        *b__passHN_eee;   //!
  TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
  TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;   //!
  TBranch        *b__passHN_emm;   //!
  TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
  TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;   //!
  TBranch        *b__passHN_mmm;   //!
  TBranch        *b__HLT_TripleMu_12_10_5;   //!
  TBranch        *b__HLT_TripleMu_12_10_5_prescale;   //!
  TBranch        *b__passMET;   //!
  TBranch        *b__HLT_PFMET120_PFMHT120_IDTight;   //!
  TBranch        *b__HLT_PFMET120_PFMHT120_IDTight_prescale;   //!
  TBranch        *b__HLT_PFMET110_PFMHT110_IDTight;   //!
  TBranch        *b__HLT_PFMET110_PFMHT110_IDTight_prescale;   //!
  TBranch        *b__HLT_PFMET170_BeamHaloCleaned;   //!
  TBranch        *b__HLT_PFMET170_BeamHaloCleaned_prescale;   //!
  TBranch        *b__HLT_PFMET170_HBHECleaned;   //!
  TBranch        *b__HLT_PFMET170_HBHECleaned_prescale;   //!
  TBranch        *b__HLT_PFMET170_HBHE_BeamHaloCleaned;   //!
  TBranch        *b__HLT_PFMET170_HBHE_BeamHaloCleaned_prescale;   //!
  TBranch        *b__HLT_PFMET170_NotCleaned;   //!
  TBranch        *b__HLT_PFMET170_NotCleaned_prescale;   //!
  TBranch        *b__HLT_MET200;   //!
  TBranch        *b__HLT_MET200_prescale;   //!
  TBranch        *b__HLT_MET250;   //!
  TBranch        *b__HLT_MET250_prescale;   //!
  TBranch        *b__HLT_MET300;   //!
  TBranch        *b__HLT_MET300_prescale;   //!
  TBranch        *b__passMETFilters;   //!
  TBranch        *b__Flag_HBHENoiseFilter;   //!
  TBranch        *b__Flag_HBHENoiseIsoFilter;   //!
  TBranch        *b__Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
  TBranch        *b__Flag_goodVertices;   //!
  TBranch        *b__Flag_eeBadScFilter;   //!
  TBranch        *b__Flag_globalTightHalo2016Filter;   //!
  TBranch        *b__flag_badPFMuonFilter;   //!
  TBranch        *b__flag_badChCandFilter;   //!
  TBranch        *b__passTTG_e;   //!
  TBranch        *b__HLT_Ele105_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b__HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;   //!
  TBranch        *b__HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b__HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;   //!
  TBranch        *b__passTTG_ee;   //!
  TBranch        *b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale;   //!
  TBranch        *b__passTTG_em;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale;   //!
  TBranch        *b__passTTG_m;   //!
  TBranch        *b__HLT_Mu50;   //!
  TBranch        *b__HLT_Mu50_prescale;   //!
  TBranch        *b__HLT_TkMu50;   //!
  TBranch        *b__HLT_TkMu50_prescale;   //!
  TBranch        *b__HLT_Mu45_eta2p1;   //!
  TBranch        *b__HLT_Mu45_eta2p1_prescale;   //!
  TBranch        *b__passTTG_mm;   //!
  TBranch        *b__HLT_Mu30_TkMu11;   //!
  TBranch        *b__HLT_Mu30_TkMu11_prescale;   //!
  TBranch        *b__nL;   //!
  TBranch        *b__nMu;   //!
  TBranch        *b__nEle;   //!
  TBranch        *b__nLight;   //!
  TBranch        *b__nTau;   //!
  TBranch        *b__lPt;   //!
  TBranch        *b__lEta;   //!
  TBranch        *b__lEtaSC;   //!
  TBranch        *b__lPhi;   //!
  TBranch        *b__lE;   //!
  TBranch        *b__lFlavor;   //!
  TBranch        *b__lCharge;   //!
  TBranch        *b__dxy;   //!
  TBranch        *b__dz;   //!
  TBranch        *b__3dIP;   //!
  TBranch        *b__3dIPSig;   //!
  TBranch        *b__lElectronMva;   //!
  TBranch        *b__lHNLoose;   //!
  TBranch        *b__lHNFO;   //!
  TBranch        *b__lHNTight;   //!
  TBranch        *b__lPOGVeto;   //!
  TBranch        *b__lPOGLoose;   //!
  TBranch        *b__lPOGMedium;   //!
  TBranch        *b__lPOGTight;   //!
  TBranch        *b__lIsPrompt;   //!
  TBranch        *b__lMatchPdgId;   //!
  TBranch        *b__relIso;   //!
  TBranch        *b__miniIso;   //!
  TBranch        *b__ptRel;   //!
  TBranch        *b__ptRatio;   //!
  TBranch        *b__nPh;   //!
  TBranch        *b__phPt;   //!
  TBranch        *b__phEta;   //!
  TBranch        *b__phEtaSC;   //!
  TBranch        *b__phPhi;   //!
  TBranch        *b__phE;   //!
  TBranch        *b__phCutBasedLoose;   //!
  TBranch        *b__phCutBasedMedium;   //!
  TBranch        *b__phCutBasedTight;   //!
  TBranch        *b__phMva;   //!
  TBranch        *b__phRandomConeChargedIsolation;   //!
  TBranch        *b__phChargedIsolation;   //!
  TBranch        *b__phNeutralHadronIsolation;   //!
  TBranch        *b__phPhotonIsolation;   //!
  TBranch        *b__phSigmaIetaIeta;   //!
  TBranch        *b__phSigmaIetaIphi;   //!
  TBranch        *b__phHadronicOverEm;   //!
  TBranch        *b__phPassElectronVeto;   //!
  TBranch        *b__phHasPixelSeed;   //!
  TBranch        *b__phIsPrompt;   //!
  TBranch        *b__phMatchMCPhotonAN15165;   //!
  TBranch        *b__phMatchMCLeptonAN15165;   //!
  TBranch        *b__phMatchPdgId;   //!
  TBranch        *b__nJets;   //!
  TBranch        *b__jetPt;   //!
  TBranch        *b__jetPt_JECUp;   //!
  TBranch        *b__jetPt_JECDown;   //!
  TBranch        *b__jetPt_JERUp;   //!
  TBranch        *b__jetPt_JERDown;   //!
  TBranch        *b__jetEta;   //!
  TBranch        *b__jetPhi;   //!
  TBranch        *b__jetE;   //!
  TBranch        *b__jetCsvV2;   //!
  TBranch        *b__jetDeepCsv_udsg;   //!
  TBranch        *b__jetDeepCsv_b;   //!
  TBranch        *b__jetDeepCsv_c;   //!
  TBranch        *b__jetDeepCsv_bb;   //!
  TBranch        *b__jetDeepCsv_cc;   //!
  TBranch        *b__jetHadronFlavour;   //!
  TBranch        *b__jetId;   //!
  TBranch        *b__met;   //!
  TBranch        *b__metJECDown;   //!
  TBranch        *b__metJECUp;   //!
  TBranch        *b__metUnclDown;   //!
  TBranch        *b__metUnclUp;   //!
  TBranch        *b__metPhi;   //!
  TBranch        *b__metPhiJECDown;   //!
  TBranch        *b__metPhiJECUp;   //!
  TBranch        *b__metPhiUnclDown;   //!
  TBranch        *b__metPhiUnclUp;   //!
    

  double hcounter[nSamples];
    
  for(int sam = 0; sam < nSamples; ++sam){
   
    cout<<"================================= 1"<<endl;
    cout<<"--------> "<< "name " << names[sam] << endl;
    cout<<"--------> "<< "fileList[sam] " << fileList[sam] << endl;
    hfile[sam] = new TFile("/Users/Martina/Desktop/file/new_sample/"+fileList[sam],"read");
    hfile[sam]->cd("blackJackAndHookers");
    inputTree[sam] = static_cast<TTree*>(hfile[sam]->Get("blackJackAndHookers/blackJackAndHookersTree"));
    inputTree[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
    inputTree[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
    inputTree[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
    inputTree[sam]->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);
    inputTree[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
    inputTree[sam]->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
    inputTree[sam]->SetBranchAddress("_ctauHN", &_ctauHN, &b__ctauHN);
    inputTree[sam]->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
    inputTree[sam]->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
    inputTree[sam]->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
    inputTree[sam]->SetBranchAddress("_ttgEventType", &_ttgEventType, &b__ttgEventType);
    inputTree[sam]->SetBranchAddress("_zgEventType", &_zgEventType, &b__zgEventType);
    inputTree[sam]->SetBranchAddress("_gen_met", &_gen_met, &b__gen_met);
    inputTree[sam]->SetBranchAddress("_gen_metPhi", &_gen_metPhi, &b__gen_metPhi);
    inputTree[sam]->SetBranchAddress("_gen_nPh", &_gen_nPh, &b__gen_nPh);
    inputTree[sam]->SetBranchAddress("_gen_phPt", _gen_phPt, &b__gen_phPt);
    inputTree[sam]->SetBranchAddress("_gen_phEta", _gen_phEta, &b__gen_phEta);
    inputTree[sam]->SetBranchAddress("_gen_phPhi", _gen_phPhi, &b__gen_phPhi);
    inputTree[sam]->SetBranchAddress("_gen_phE", _gen_phE, &b__gen_phE);
    inputTree[sam]->SetBranchAddress("_gen_phMomPdg", _gen_phMomPdg, &b__gen_phMomPdg);
    inputTree[sam]->SetBranchAddress("_gen_phIsPrompt", _gen_phIsPrompt, &b__gen_phIsPrompt);
    inputTree[sam]->SetBranchAddress("_gen_phMinDeltaR", _gen_phMinDeltaR, &b__gen_phMinDeltaR);
    inputTree[sam]->SetBranchAddress("_gen_phPassParentage", _gen_phPassParentage, &b__gen_phPassParentage);
    inputTree[sam]->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
    inputTree[sam]->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
    inputTree[sam]->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
    inputTree[sam]->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
    inputTree[sam]->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
    inputTree[sam]->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
    inputTree[sam]->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
    inputTree[sam]->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
    inputTree[sam]->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
    inputTree[sam]->SetBranchAddress("_passHN_eem", &_passHN_eem, &b__passHN_eem);
    inputTree[sam]->SetBranchAddress("_passHN_1l", &_passHN_1l, &b__passHN_1l);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale, &b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    //    inputTree[sam]->SetBranchAddress("_passHN_1l", &_passHN_1l, &b__passHN_1l);
    inputTree[sam]->SetBranchAddress("_HLT_Ele27_WPTight_Gsf", &_HLT_Ele27_WPTight_Gsf, &b__HLT_Ele27_WPTight_Gsf);
    inputTree[sam]->SetBranchAddress("_HLT_Ele27_WPTight_Gsf_prescale", &_HLT_Ele27_WPTight_Gsf_prescale, &b__HLT_Ele27_WPTight_Gsf_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_IsoMu24", &_HLT_IsoMu24, &b__HLT_IsoMu24);
    inputTree[sam]->SetBranchAddress("_HLT_IsoMu24_prescale", &_HLT_IsoMu24_prescale, &b__HLT_IsoMu24_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_IsoTkMu24", &_HLT_IsoTkMu24, &b__HLT_IsoTkMu24);
    inputTree[sam]->SetBranchAddress("_HLT_IsoTkMu24_prescale", &_HLT_IsoTkMu24_prescale, &b__HLT_IsoTkMu24_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_eee", &_passHN_eee, &b__passHN_eee);
    inputTree[sam]->SetBranchAddress("_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
    inputTree[sam]->SetBranchAddress("_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale", &_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale, &b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_emm", &_passHN_emm, &b__passHN_emm);
    inputTree[sam]->SetBranchAddress("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
    inputTree[sam]->SetBranchAddress("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale, &b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_mmm", &_passHN_mmm, &b__passHN_mmm);
    inputTree[sam]->SetBranchAddress("_HLT_TripleMu_12_10_5", &_HLT_TripleMu_12_10_5, &b__HLT_TripleMu_12_10_5);
    inputTree[sam]->SetBranchAddress("_HLT_TripleMu_12_10_5_prescale", &_HLT_TripleMu_12_10_5_prescale, &b__HLT_TripleMu_12_10_5_prescale);
    inputTree[sam]->SetBranchAddress("_passMET", &_passMET, &b__passMET);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET120_PFMHT120_IDTight", &_HLT_PFMET120_PFMHT120_IDTight, &b__HLT_PFMET120_PFMHT120_IDTight);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET120_PFMHT120_IDTight_prescale", &_HLT_PFMET120_PFMHT120_IDTight_prescale, &b__HLT_PFMET120_PFMHT120_IDTight_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET110_PFMHT110_IDTight", &_HLT_PFMET110_PFMHT110_IDTight, &b__HLT_PFMET110_PFMHT110_IDTight);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET110_PFMHT110_IDTight_prescale", &_HLT_PFMET110_PFMHT110_IDTight_prescale, &b__HLT_PFMET110_PFMHT110_IDTight_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_BeamHaloCleaned", &_HLT_PFMET170_BeamHaloCleaned, &b__HLT_PFMET170_BeamHaloCleaned);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_BeamHaloCleaned_prescale", &_HLT_PFMET170_BeamHaloCleaned_prescale, &b__HLT_PFMET170_BeamHaloCleaned_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_HBHECleaned", &_HLT_PFMET170_HBHECleaned, &b__HLT_PFMET170_HBHECleaned);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_HBHECleaned_prescale", &_HLT_PFMET170_HBHECleaned_prescale, &b__HLT_PFMET170_HBHECleaned_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_HBHE_BeamHaloCleaned", &_HLT_PFMET170_HBHE_BeamHaloCleaned, &b__HLT_PFMET170_HBHE_BeamHaloCleaned);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_HBHE_BeamHaloCleaned_prescale", &_HLT_PFMET170_HBHE_BeamHaloCleaned_prescale, &b__HLT_PFMET170_HBHE_BeamHaloCleaned_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_NotCleaned", &_HLT_PFMET170_NotCleaned, &b__HLT_PFMET170_NotCleaned);
    inputTree[sam]->SetBranchAddress("_HLT_PFMET170_NotCleaned_prescale", &_HLT_PFMET170_NotCleaned_prescale, &b__HLT_PFMET170_NotCleaned_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_MET200", &_HLT_MET200, &b__HLT_MET200);
    inputTree[sam]->SetBranchAddress("_HLT_MET200_prescale", &_HLT_MET200_prescale, &b__HLT_MET200_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_MET250", &_HLT_MET250, &b__HLT_MET250);
    inputTree[sam]->SetBranchAddress("_HLT_MET250_prescale", &_HLT_MET250_prescale, &b__HLT_MET250_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_MET300", &_HLT_MET300, &b__HLT_MET300);
    inputTree[sam]->SetBranchAddress("_HLT_MET300_prescale", &_HLT_MET300_prescale, &b__HLT_MET300_prescale);
    inputTree[sam]->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters);
    inputTree[sam]->SetBranchAddress("_Flag_HBHENoiseFilter", &_Flag_HBHENoiseFilter, &b__Flag_HBHENoiseFilter);
    inputTree[sam]->SetBranchAddress("_Flag_HBHENoiseIsoFilter", &_Flag_HBHENoiseIsoFilter, &b__Flag_HBHENoiseIsoFilter);
    inputTree[sam]->SetBranchAddress("_Flag_EcalDeadCellTriggerPrimitiveFilter", &_Flag_EcalDeadCellTriggerPrimitiveFilter, &b__Flag_EcalDeadCellTriggerPrimitiveFilter);
    inputTree[sam]->SetBranchAddress("_Flag_goodVertices", &_Flag_goodVertices, &b__Flag_goodVertices);
    inputTree[sam]->SetBranchAddress("_Flag_eeBadScFilter", &_Flag_eeBadScFilter, &b__Flag_eeBadScFilter);
    inputTree[sam]->SetBranchAddress("_Flag_globalTightHalo2016Filter", &_Flag_globalTightHalo2016Filter, &b__Flag_globalTightHalo2016Filter);
    inputTree[sam]->SetBranchAddress("_flag_badPFMuonFilter", &_flag_badPFMuonFilter, &b__flag_badPFMuonFilter);
    inputTree[sam]->SetBranchAddress("_flag_badChCandFilter", &_flag_badChCandFilter, &b__flag_badChCandFilter);
    inputTree[sam]->SetBranchAddress("_passTTG_e", &_passTTG_e, &b__passTTG_e);
    inputTree[sam]->SetBranchAddress("_HLT_Ele105_CaloIdVT_GsfTrkIdT", &_HLT_Ele105_CaloIdVT_GsfTrkIdT, &b__HLT_Ele105_CaloIdVT_GsfTrkIdT);
    inputTree[sam]->SetBranchAddress("_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale", &_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale, &b__HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Ele115_CaloIdVT_GsfTrkIdT", &_HLT_Ele115_CaloIdVT_GsfTrkIdT, &b__HLT_Ele115_CaloIdVT_GsfTrkIdT);
    inputTree[sam]->SetBranchAddress("_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale", &_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale, &b__HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_ee", &_passTTG_ee, &b__passTTG_ee);
    inputTree[sam]->SetBranchAddress("_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_em", &_passTTG_em, &b__passTTG_em);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", &_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL, &b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale", &_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale, &b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_m", &_passTTG_m, &b__passTTG_m);
    inputTree[sam]->SetBranchAddress("_HLT_Mu50", &_HLT_Mu50, &b__HLT_Mu50);
    inputTree[sam]->SetBranchAddress("_HLT_Mu50_prescale", &_HLT_Mu50_prescale, &b__HLT_Mu50_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu50", &_HLT_TkMu50, &b__HLT_TkMu50);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu50_prescale", &_HLT_TkMu50_prescale, &b__HLT_TkMu50_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu45_eta2p1", &_HLT_Mu45_eta2p1, &b__HLT_Mu45_eta2p1);
    inputTree[sam]->SetBranchAddress("_HLT_Mu45_eta2p1_prescale", &_HLT_Mu45_eta2p1_prescale, &b__HLT_Mu45_eta2p1_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_mm", &_passTTG_mm, &b__passTTG_mm);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_TkMu11", &_HLT_Mu30_TkMu11, &b__HLT_Mu30_TkMu11);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_TkMu11_prescale", &_HLT_Mu30_TkMu11_prescale, &b__HLT_Mu30_TkMu11_prescale);
    inputTree[sam]->SetBranchAddress("_nL", &_nL, &b__nL);
    inputTree[sam]->SetBranchAddress("_nMu", &_nMu, &b__nMu);
    inputTree[sam]->SetBranchAddress("_nEle", &_nEle, &b__nEle);
    inputTree[sam]->SetBranchAddress("_nLight", &_nLight, &b__nLight);
    inputTree[sam]->SetBranchAddress("_nTau", &_nTau, &b__nTau);
    inputTree[sam]->SetBranchAddress("_lPt", _lPt, &b__lPt);
    inputTree[sam]->SetBranchAddress("_lEta", _lEta, &b__lEta);
    inputTree[sam]->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
    inputTree[sam]->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
    inputTree[sam]->SetBranchAddress("_lE", _lE, &b__lE);
    inputTree[sam]->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
    inputTree[sam]->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
    inputTree[sam]->SetBranchAddress("_dxy", _dxy, &b__dxy);
    inputTree[sam]->SetBranchAddress("_dz", _dz, &b__dz);
    inputTree[sam]->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
    inputTree[sam]->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
    inputTree[sam]->SetBranchAddress("_lElectronMva", _lElectronMva, &b__lElectronMva);
    inputTree[sam]->SetBranchAddress("_lHNLoose", _lHNLoose, &b__lHNLoose);
    inputTree[sam]->SetBranchAddress("_lHNFO", _lHNFO, &b__lHNFO);
    inputTree[sam]->SetBranchAddress("_lHNTight", _lHNTight, &b__lHNTight);
    inputTree[sam]->SetBranchAddress("_lPOGVeto", _lPOGVeto, &b__lPOGVeto);
    inputTree[sam]->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
    inputTree[sam]->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
    inputTree[sam]->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
    inputTree[sam]->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
    inputTree[sam]->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
    inputTree[sam]->SetBranchAddress("_relIso", _relIso, &b__relIso);
    inputTree[sam]->SetBranchAddress("_miniIso", _miniIso, &b__miniIso);
    inputTree[sam]->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
    inputTree[sam]->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
    inputTree[sam]->SetBranchAddress("_nPh", &_nPh, &b__nPh);
    inputTree[sam]->SetBranchAddress("_phPt", _phPt, &b__phPt);
    inputTree[sam]->SetBranchAddress("_phEta", _phEta, &b__phEta);
    inputTree[sam]->SetBranchAddress("_phEtaSC", _phEtaSC, &b__phEtaSC);
    inputTree[sam]->SetBranchAddress("_phPhi", _phPhi, &b__phPhi);
    inputTree[sam]->SetBranchAddress("_phE", _phE, &b__phE);
    inputTree[sam]->SetBranchAddress("_phCutBasedLoose", _phCutBasedLoose, &b__phCutBasedLoose);
    inputTree[sam]->SetBranchAddress("_phCutBasedMedium", _phCutBasedMedium, &b__phCutBasedMedium);
    inputTree[sam]->SetBranchAddress("_phCutBasedTight", _phCutBasedTight, &b__phCutBasedTight);
    inputTree[sam]->SetBranchAddress("_phMva", _phMva, &b__phMva);
    inputTree[sam]->SetBranchAddress("_phRandomConeChargedIsolation", _phRandomConeChargedIsolation, &b__phRandomConeChargedIsolation);
    inputTree[sam]->SetBranchAddress("_phChargedIsolation", _phChargedIsolation, &b__phChargedIsolation);
    inputTree[sam]->SetBranchAddress("_phNeutralHadronIsolation", _phNeutralHadronIsolation, &b__phNeutralHadronIsolation);
    inputTree[sam]->SetBranchAddress("_phPhotonIsolation", _phPhotonIsolation, &b__phPhotonIsolation);
    inputTree[sam]->SetBranchAddress("_phSigmaIetaIeta", _phSigmaIetaIeta, &b__phSigmaIetaIeta);
    inputTree[sam]->SetBranchAddress("_phSigmaIetaIphi", _phSigmaIetaIphi, &b__phSigmaIetaIphi);
    inputTree[sam]->SetBranchAddress("_phHadronicOverEm", _phHadronicOverEm, &b__phHadronicOverEm);
    inputTree[sam]->SetBranchAddress("_phPassElectronVeto", _phPassElectronVeto, &b__phPassElectronVeto);
    inputTree[sam]->SetBranchAddress("_phHasPixelSeed", _phHasPixelSeed, &b__phHasPixelSeed);
    inputTree[sam]->SetBranchAddress("_phIsPrompt", _phIsPrompt, &b__phIsPrompt);
    inputTree[sam]->SetBranchAddress("_phMatchMCPhotonAN15165", _phMatchMCPhotonAN15165, &b__phMatchMCPhotonAN15165);
    inputTree[sam]->SetBranchAddress("_phMatchMCLeptonAN15165", _phMatchMCLeptonAN15165, &b__phMatchMCLeptonAN15165);
    inputTree[sam]->SetBranchAddress("_phMatchPdgId", _phMatchPdgId, &b__phMatchPdgId);
    inputTree[sam]->SetBranchAddress("_nJets", &_nJets, &b__nJets);
    inputTree[sam]->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
    inputTree[sam]->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
    inputTree[sam]->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
    inputTree[sam]->SetBranchAddress("_jetPt_JERUp", _jetPt_JERUp, &b__jetPt_JERUp);
    inputTree[sam]->SetBranchAddress("_jetPt_JERDown", _jetPt_JERDown, &b__jetPt_JERDown);
    inputTree[sam]->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
    inputTree[sam]->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
    inputTree[sam]->SetBranchAddress("_jetE", _jetE, &b__jetE);
    inputTree[sam]->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_cc", _jetDeepCsv_cc, &b__jetDeepCsv_cc);
    inputTree[sam]->SetBranchAddress("_jetHadronFlavour", _jetHadronFlavour, &b__jetHadronFlavour);
    inputTree[sam]->SetBranchAddress("_jetId", _jetId, &b__jetId);
    inputTree[sam]->SetBranchAddress("_met", &_met, &b__met);
    inputTree[sam]->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
    inputTree[sam]->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
    inputTree[sam]->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
    inputTree[sam]->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
    inputTree[sam]->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
    inputTree[sam]->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
    inputTree[sam]->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
    inputTree[sam]->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
    inputTree[sam]->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    _hCounter->Read("hCounter");
    hcounter[sam] = _hCounter->GetBinContent(1);
  }//end for on tree
    
  //******************* HISTO **********************
  const int nCat=2;
  const int nDist = 29;  //Number of distributions to plot
  TH1D* Histos[nDist][nCat][nSamples_eff +1];
  const TString catNames[nCat] ={"ossf", "no_ossf"};
  const TString Histnames_ossf[nDist] = {"DeltaPhi_pair","DeltaPhi_lt","DeltaPhi_st", "LeptonPt_le","LeptonPt_subl", "LeptonPt_tr","Sum3Pt","Sum2Pt_lt","Sum2Pt_st","Sum2Pt_ls","Mlll","Mll","Mll_pair OSSF", "MET", "MT", "NJets", "NbJets","HT", "DeltaR_pair","DeltaR_lt","DeltaR_st","RelIso_l", "RelIso_t", "dxy_l", "dxy_t","dz_l", "dz_t","3dIP_l", "3dIP_t"};
    
  //const TString Histnames_n_ossf[nDist] = { "LeptonPt_le","LeptonPt_subl" "LeptonPt_tr","Sum3Pt","Sum2Pt_lt","Sum2Pt_st","Sum2Pt_ls","Mlll","Mll","Mll_pair NOSSF", "MET", "MT", "NJets", "NbJets","HT", "DeltaR","DeltaR_pair","DeltaR_lt","DeltaR_st","RelIso_l", "RelIso_t", "dxy_l", "dxy_t","dz_l", "dz_t","3dIP_l", "3dIP_t"};
    
  const TString Xaxes[nDist] = {"#Delta #Phi (OSSF pair)", "#Delta #Phi (leading,trailing)", "#Delta #Phi (sub-leading, trailing)","P_{T}(leading l) (GeV)","P_{T}(sub-leading l) (GeV)", "P_{T}(trailing l) (GeV)", "SumP_{T}(3leptons) (GeV)","SumP_{T}(leading+trailing) (GeV)","SumP_{T}(sub-leading+trailing) (GeV)","SumP_{T}(leading+sub-leading) (GeV)",
				"M_{lll} (GeV)","M_{ll} (sub-leading+soft) (GeV)", "M_{ll} OSSF pair (GeV)", "MET (GeV)", "M_{T} (GeV)", "number of jets", "number of b-jets","HT (GeV)" , "#Delta R (OSSF pair)", "#Delta R (leading,trailing)", "#Delta R (sub-leading, trailing)", "relIso(leading)", "relIso(trailing)", "|d_{xy}(leading)| (cm)", "|d_{xy}(trailing)| (cm)", "|d_{z}(leading)| (cm)", "|d_{z}(trailing)| (cm)",  "|3DIP(leading)| (cm)","|3DIP(trailing)| (cm)"};
    
  const TString Units[nDist] = {"","","","GeV", "GeV", "GeV", "GeV", "GeV","GeV","GeV","GeV","GeV","GeV","GeV","GeV", "", "", "GeV",  "", "", "",  "cm", "cm","cm", "cm", "cm", "cm","cm","cm"};
    
  const double HistMin[nDist] = {-0.5,-0.5,-0.5, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
  const double HistMax[nDist] = {4,4,4,200, 200, 200, 200, 200,  200,  200,  200,  200,  200, 200, 100, 10,10, 100, 4,4,4,0.99,0.99,5,5, 6, 6, 4, 4};
  const int nBins[nDist] = {20,20,20,50, 50, 50, 50, 50, 50,50,50,50,50,50,20,10,10,20,25,25,25,70,70,10,10,10,10,10,10};

  double data_sum1=0;
  double data_sum2=0;
  double mc_sum1=0;
  double mc_sum2=0;

 
  double tot_1gev=0;
  double tot_2gev=0;
  double fo_1gev=0;
  double fo_2gev=0;
  double tight_1gev=0;
  double tight_2gev=0;
  double tot_mu_1gev=0;
  double tot_mu_2gev=0;
  double fo_mu_1gev=0;
  double fo_mu_2gev=0;
  double tight_mu_1gev=0;
  double tight_mu_2gev=0;
  double tot_e_1gev=0;
  double tot_e_2gev=0;
  double fo_e_1gev=0;
  double fo_e_2gev=0;
  double tight_e_1gev=0;
  double tight_e_2gev=0;
  double fo_3_1gev=0;
  double fo_3_2gev=0;
  double fo_3_mu_1gev=0;
  double fo_3_mu_2gev=0;
  double fo_3_e_1gev=0;
  double fo_3_e_2gev=0;


  double tot_1gev_prompt=0;
  double tot_2gev_prompt=0;
  double fo_1gev_prompt=0;
  double fo_2gev_prompt=0;
  double tight_1gev_prompt=0;
  double tight_2gev_prompt=0;
  double tot_mu_1gev_prompt=0;
  double tot_mu_2gev_prompt=0;
  double fo_mu_1gev_prompt=0;
  double fo_mu_2gev_prompt=0;
  double tight_mu_1gev_prompt=0;
  double tight_mu_2gev_prompt=0;
  double tot_e_1gev_prompt=0;
  double tot_e_2gev_prompt=0;
  double fo_e_1gev_prompt=0;
  double fo_e_2gev_prompt=0;
  double tight_e_1gev_prompt=0;
  double tight_e_2gev_prompt=0;
  double fo_3_1gev_prompt=0;
  double fo_3_2gev_prompt=0;
  double fo_3_mu_1gev_prompt=0;
  double fo_3_mu_2gev_prompt=0;
  double fo_3_e_1gev_prompt=0;
  double fo_3_e_2gev_prompt=0;


  double denominatore_ossf=0;
  double denominatore_n_ossf=0;    
  const int bin_num=100;
  double numerator_ossf[bin_num];
  double numerator_n_ossf[bin_num];


  cout<<"------ 1"<<endl;
  for(int i = 0; i < nDist; ++i){
    // cout<<"-------> numero isto: "<<i<<endl;
    float BinWidth = (HistMax[i] - HistMin[i])/nBins[i];
    std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
    for(int effsam = 0; effsam < nSamples_eff + 1; ++effsam){
      //cout<<"------------> numero sample: "<<effsam<<endl;
      for(int cat = 0; cat < nCat; ++cat){
	//cout<<"--------------------------> numero cat: "<<cat<<endl;


	Histos[i][cat][effsam] = new TH1D(eff_names[effsam] + catNames[cat] + Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + "; events /" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);
	Histos[i][cat][effsam]->Sumw2();
	
	//cout<<i<<" - "<< Histnames_ossf[i] <<")   "<<effsam<<"  -  "<<eff_names[effsam]<<"]   "<<cat<<"  -  "<<catNames[cat]<<"} "<<" --->  "<<Xaxes[i]<<"  "<<Units[i]<<"   "<<HistMin[i]<<"   "<<HistMax[i]<<"   "<<nBins[i]<<endl;
      }
    }
  }
  //Calculate the center of the maximum bin of each histogram
  double maxBinC[nDist];
  for(int i = 0; i < nDist; ++i){
    maxBinC[i] = Histos[i][0][0]->GetBinCenter(Histos[i][0][0]->GetNbinsX());
  }

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARAMETERS AND CUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const int number_veto_leptons=3;
  const double b_jets_wp= 0.5426;
  const double b_jets_pt= 25;
  const double met_cuts = 75;
  const double mlll_cuts = 80;
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  //Loop over all samples
  Double_t Counter[nSamples];
  for(int i = 0; i  < nSamples; ++i){
    Counter[i] = 0;
  }
  Double_t scale[nSamples];
  

  //-------------->  LOOP ON ALL SAMPLES
  for(int sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
    Long64_t nEntries = inputTree[sam]->GetEntries();
    if(sam != 0){
      if(names[sam] == names[sam -1]) --effsam;
    }
    if(sam >= 0){ // 1 data set
      scale[sam] = xSections[sam]*luminosity*1000/(hcounter[sam]);
    }
    //if (effsam == 0) continue;
    std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
    cout << effsam << "   "<<sam<<endl;	
    //if (sam < 26) continue;

    //if (effsam !=1) continue;
	  
        double progress = 0; 	//For printing progress bar


    //--------------> LOOOP ON THE TREE"S ENTRIES
    for (Long64_t it = 0; it < nEntries; ++it){
      inputTree[sam]->GetEntry(it);    
      // if (it%10000 == 0) cout<<'.'<<flush;

       if(it%100 == 0 && it != 0){
	progress += (double) (100./nEntries);
	printProgress(progress);
      } else if(it == nEntries -1){
	progress = 1.;
	printProgress(progress);
      }



      double scal = 0;
      scal = scale[sam]*_weight; 
      /*if(effsam == 0) scal = 0;
      else{
	scal = scale[sam]*_weight; 
	}   */
      //scal=1;
      //if (effsam !=1) continue;
      //number of leptons in tree
      tot_1gev= tot_1gev + 1*scal;

      // if (effsam < 15) continue;
      if (effsam == 0) continue;

      if(_nL < number_veto_leptons) continue;
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTORS AND VARIABLES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      int               nBjets = 0;
      unsigned*         ind = new unsigned[_nL];	//new indices of good leptons
      unsigned          tightC = 0;//Count number of T leptons
      unsigned          promptC = 0;
      double            low_mass_pt_base[1];
      double*           conePt = new double[_nL];
      double            faxtore_FR=0;


      double            prov_index[3]={-1,-1,-1};
      double            prov_number_tight[1]= {-1};
      double            prov_number_prompt[1]= {-1};
      double            prov_fattore[1]= {-1};
      int               skip_event[1]= {-1};
      double            faxtore[1]= {-100};

      TLorentzVector    lepton_reco[3];
      int               flavors_3l[3];
      int               charge_3l[3];
      TLorentzVector    sum_3l_rec;	//M_3l 
      TLorentzVector    pair [3];
      int               kind[1] ={-1}; // 0 no-ossf
      TLorentzVector    sum_2l_rec_pair; 	//M_2l best Z candidate
      int               event_clas[1]={-1}; // 1* = eee   2* = emm   3* = eem  4* = mmm
      int               check_mt= -1;
      Double_t          delta_R_max=-1;
      Double_t          delta_R_min=-1;
      Double_t          _mll_min=50000;  
      TLorentzVector    lepton_transv;
      TLorentzVector    METvec;
      double            m_T=0; 
      double            MET=0;
      double            MET_PHI=0;

      int               nominal=-1;
      int               jec_check=-1; 
      int               unclustered_met=-1;
      int               up=-1;
      int               down=-1;
      double            new_met[1]= {0};
      double            new_met_phi[1]= {0};
      int               new_number_bjets[1]= {0};

      int               new_pt_checks= -1;
      bool              trigger_fired = false;
      bool              low_pt_event=false;
      bool              high_pt_event = false;

      unsigned*         _SR_low= new unsigned[8];
      double             search_region_fill[1]={-1};
      bool              data_control_region=false;

       //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //------------------------------------------------------------ selection class 
      faxtore_FR=1;
      TString prova = "/Users/Martina/Desktop/file/willem/"+fileList[sam];
      Selection right_lepton(prova);
      right_lepton.selezione(_gen_nL,_nL,  _lPt  , _lEta  , _lE  ,  _lPOGMedium  ,  _lPOGLoose ,_lFlavor  , _lCharge  , _relIso  , _dxy  , _dz  , _3dIP  , _3dIPSig  ,  _lHNLoose  ,  _lHNFO  ,  _lHNTight  , _lElectronMva     , _lIsPrompt,  prov_index,prov_number_tight,prov_number_prompt,skip_event, effsam, *&fakeRate_mu, *&fakeRate_e, faxtore, _gen_lPt  , _gen_lEta, _gen_lIsPrompt , _gen_lFlavor);
      tightC = 0;
      if (skip_event[0] == -1) continue;
      // if (effsam== 1 ||effsam== 2 ||effsam== 3 ||effsam== 4 ||effsam== 5 ) cout<< "signal1"<<endl;

      ind[0] = prov_index[0];
      ind[1] = prov_index[1];
      ind[2] = prov_index[2];
      tightC= prov_number_tight[0];
      promptC = prov_number_prompt[0]; 
      faxtore_FR= faxtore[0];
      bool tightFail = (tightC < 3);
      if (tightFail) continue;
      if ((effsam == 1 ||effsam == 2 ||effsam == 3 ||effsam == 4 ||effsam == 5 ||effsam == 6 ||effsam == 7 ||effsam == 8 ||effsam == 9 ||effsam == 10 ||effsam == 11 ) && prov_number_prompt[0] != 3) continue;

      //--------------------------------------------------------------------------------	
      fo_3_1gev= fo_3_1gev + 1*scal;   
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      for (int i =0; i < 3; i ++){
	lepton_reco[i].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]); 
	flavors_3l[i]=_lFlavor[ind[i]];
	charge_3l[i]=_lCharge[ind[i]];
	conePt[i] =  _lPt[ind[i]];
      }
      //M_3l
      sum_3l_rec.SetPtEtaPhiE(0,0,0,0);
      sum_3l_rec= (lepton_reco[0]+ lepton_reco[1]+lepton_reco[2] );
      
      //================== event classification ========================  
      // ---------------- > OSSF or NO_OSSF
      ossf_no_ossf( kind, pair,lepton_reco[0], lepton_reco[1], lepton_reco[2], flavors_3l, charge_3l);
      if (kind[0]  == -1) continue;  
      //M_2l best Z candidate
      sum_2l_rec_pair.SetPtEtaPhiE(0,0,0,0);
      sum_2l_rec_pair= (pair[0]+pair[1] );	
      bool ossf_event= false;
      if (kind[0] == 1) ossf_event = true;
      //if (kind[0] == 0) continue;  		
      // ---------------- > CHANNELS
      class_os( event_clas,  flavors_3l, charge_3l);
      if (event_clas[0] == -1) continue;  
      ////// ===============>  LOW mass selection:   
      //if(lepton_reco[0].Pt() > 55) continue;
      // ---------------- > cut on M_3L > M_W
      //if (sum_3l_rec.M() > mlll_cuts) continue;
      //=============================================================  

      //================== Binning variables ========================  
      // ---------------- > Delta R
      delta_R_max = lepton_reco[0].DeltaR(lepton_reco[1]);
      if (lepton_reco[0].DeltaR(lepton_reco[2]) > delta_R_max) delta_R_max = lepton_reco[0].DeltaR(lepton_reco[2]);
      if (lepton_reco[1].DeltaR(lepton_reco[2]) > delta_R_max) delta_R_max = lepton_reco[1].DeltaR(lepton_reco[2]);
      delta_R_min = lepton_reco[0].DeltaR(lepton_reco[1]);
      if (lepton_reco[0].DeltaR(lepton_reco[2]) < delta_R_min) delta_R_min = lepton_reco[0].DeltaR(lepton_reco[2]);
      if (lepton_reco[1].DeltaR(lepton_reco[2]) < delta_R_min) delta_R_min = lepton_reco[1].DeltaR(lepton_reco[2]);
      // ---------------- > M_min OS
      check_mt=-1;
      if ((_lCharge[ind[0]] != _lCharge[ind[1]] )) {
	_mll_min = (lepton_reco[0]+lepton_reco[1]).M();
	check_mt= 2;
      }
      if ((_lCharge[ind[0]] != _lCharge[ind[2]] ) && (lepton_reco[0]+lepton_reco[2]).M() < _mll_min){
	_mll_min = (lepton_reco[0]+lepton_reco[2]).M();
	check_mt= 1;
      }
      if ((_lCharge[ind[1]] != _lCharge[ind[2]] ) && (lepton_reco[1]+lepton_reco[2]).M() < _mll_min) {
	_mll_min = (lepton_reco[1]+lepton_reco[2]).M();
	check_mt= 0;
      }   
      // ---------------- > to built later M_T_third
      if (check_mt == 0 ) lepton_transv.SetPtEtaPhiE(lepton_reco[0].Pt(),0, lepton_reco[0].Phi(), lepton_reco[0].Pt());
      if (check_mt == 1 ) lepton_transv.SetPtEtaPhiE(lepton_reco[1].Pt(),0, lepton_reco[1].Phi(), lepton_reco[1].Pt());
      if (check_mt == 2 ) lepton_transv.SetPtEtaPhiE(lepton_reco[2].Pt(),0, lepton_reco[2].Phi(), lepton_reco[2].Pt());	
      //=============================================================  

      //=============== PT cocktail  and PT categories ==============
      new_pt_checks= 0;
      if (event_clas[0]== 30){
	new_pt_checks= 1;
        if((_lFlavor[ind[2]]== 0 && lepton_reco[2].Pt() < 15) && lepton_reco[0].Pt() < 23)  new_pt_checks= 0;
	if((_lFlavor[ind[2]]== 1 && lepton_reco[2].Pt()< 8) && (lepton_reco[0].Pt() < 25   || lepton_reco[1].Pt() < 15))  new_pt_checks= 0;
	if((_lFlavor[ind[2]]== 1 && lepton_reco[2].Pt() > 8) && (lepton_reco[0].Pt() < 23   &&  lepton_reco[1].Pt() < 15))  new_pt_checks= 0;
      }
      if (event_clas[0]== 20){
	new_pt_checks= 1;
	if((_lFlavor[ind[2]]== 1 && lepton_reco[2].Pt()< 9) && lepton_reco[0].Pt() < 23)  new_pt_checks= 0;
      }   
      //if (new_pt_checks == 0) continue;
      // ---------------- > 2 PT categories 
      low_pt_event= false;
      if (lepton_reco[0].Pt()<= 30) low_pt_event= true;
      if (lepton_reco[0].Pt()> 30 && lepton_reco[0].Pt() < 55)  high_pt_event = true;
      //=============================================================  
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      TRIGGER     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // 1* = eee
      // 2* = emm
      // 3* = eem
      // 4* = mmm
      trigger_fired = false;
      if ((event_clas[0]== 1 ||  event_clas[0]== 10)   &&  (_passHN_1l   || _passHN_eee)) trigger_fired = true;
      if ((event_clas[0]== 2 ||  event_clas[0]== 20)   &&  (_passHN_1l   || _passHN_emm)) trigger_fired = true;
      if ((event_clas[0]== 3 ||  event_clas[0]== 30)   &&  (_passHN_1l   || _passHN_eem)) trigger_fired = true;
      if ((event_clas[0]== 4 ||  event_clas[0]== 40)   &&  (_passHN_1l   || _passHN_mmm)) trigger_fired = true;     
      //if ( !trigger_fired ) continue;
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      MET and MT     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      nominal=1;
      jec_check=0;
      unclustered_met=0;
      up=0;
      down=0;
      // ---------------- > MET and B JETS according with SR-unc
      //  met_mt_bjets(nominal,  jec_check,  unclustered_met,  up,  down,   _metJECUp,       _met_phiJECUp,       _metOtherUp,       _met_phiOtherUp,       _metJECDown,       _met_phiJECDown,       _metOtherDown,       _met_phiOtherDown,   _jetPt,   _jetPtUp,   _jetPtDown, _met,  _met_phi, _nJets,  _csv, new_met, new_met_phi,  new_number_bjets);
      MET = new_met[0];
      MET_PHI = new_met_phi[0];
      nBjets = new_number_bjets[0];
      // ------------------- cuts according to new met and bjets veto
      //if( tightFail) continue;
      if (nBjets != 0) continue;
      METvec.SetPtEtaPhiE(MET, 0, MET_PHI,MET);        
      m_T= (lepton_transv + METvec).Mag();
      //if (MET > met_cuts) continue;
      tight_1gev= tight_1gev + 1*scal;    
      // ---------------------------------------------------------------------------------------------------------------- 
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      

      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    SR DEFINITION     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
      for (int i =0; i < 8; i ++){
	_SR_low[i]= false;
      }
      //LOW_SR    
      if (low_pt_event){
	if (_mll_min <= 10 )_SR_low[0] = true;
	if (_mll_min > 10  && _mll_min <= 20 )_SR_low[1] = true;
	if (_mll_min > 20  && _mll_min <= 30 )_SR_low[2] = true;
	if (_mll_min >30 )_SR_low[3] = true;
      }
      if (high_pt_event){
	if (_mll_min <= 10 )_SR_low[4] = true;
	if (_mll_min > 10  && _mll_min <= 20 )_SR_low[5] = true;
	if (_mll_min > 20  && _mll_min <= 30 )_SR_low[6] = true;
	if (_mll_min >30 )_SR_low[7] = true;
      }	    
      
      if (_SR_low[0] || _SR_low[4]) search_region_fill[0] = 1;
      if (_SR_low[1] || _SR_low[5]) search_region_fill[0] = 2;
      if (_SR_low[2] || _SR_low[6]) search_region_fill[0] = 3;
      if (_SR_low[3] || _SR_low[7]) search_region_fill[0] = 4;

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //   if (effsam== 1 ||effsam== 2 ||effsam== 3 ||effsam== 4 ||effsam== 5 ) cout<< "signal"<<endl;


      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      NO_PROMPT     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
      //if(effsam == 0 && tightFail)  data_control_region=true;
      unsigned fill = effsam;
      // if( data_control_region ) fill = nSamples_eff;     
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      HISTOGRAMS     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
      double pt_cone_leading=          lepton_reco[0].Pt() ;
      double pt_cone_sub_leading=      lepton_reco[1].Pt();
      double pt_cone_trailing=         lepton_reco[2].Pt();

     

    //      double values_sr[nDist_sr] = {search_region_fill[0], search_region_fill[0]};
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //******************** PT STUDIES *****************************************************
      double numero_bin_pt= 100;

      if (!ossf_event) {
	 denominatore_n_ossf= denominatore_n_ossf +1*scal;
	  for (int i = 0; i< bin_num; i++){

	    double j = i;
	    //   cout<<i<<")  "<<j/200<< "   "<<j/100<< "  "<< _dxy[ind[2]]<< endl;
	    //|| _dz[ind[2]] > j/100
	    //_dxy[ind[2]] < j/200
	    // TMath::Abs(_dxy[ind[2]]) > j/200 && TMath::Abs(_dxy[ind[1]]) > j/200

	    if ((TMath::Abs(_dxy[ind[1]]) > j/200 || TMath::Abs(_dz[ind[1]]) > j/200  || _3dIPSig[ind[1]] > j/2.5)  && (TMath::Abs(_dxy[ind[2]]) > j/200 || TMath::Abs(_dz[ind[2]]) > j/200  || _3dIPSig[ind[2]] > j/2.5) ) numerator_n_ossf[i] = numerator_n_ossf[i] +1*scal;
	  }
      }    
      if (ossf_event) {
	denominatore_ossf= denominatore_ossf +1*scal;
	 for (int  i = 0; i< bin_num; i++){
	   double j = i;
	   if ((TMath::Abs(_dxy[ind[1]]) > j/200 || TMath::Abs(_dz[ind[1]]) > j/200  || _3dIPSig[ind[1]] > j/2.5)  && (TMath::Abs(_dxy[ind[2]]) > j/200 || TMath::Abs(_dz[ind[2]]) > j/200  || _3dIPSig[ind[2]] > j/2.5)  ) numerator_ossf[i] = numerator_ossf[i]+1*scal;
		    
	 }

      }
      //  cout<<"---> "<< sam<<endl;
      for (int i =0; i< nSamples; i++){
	if (i ==0) continue;
	if (sam == i  ) {
	  //	  cout<< sam<< "    "<< names_files[i-1]<<endl;
	  std::ofstream zg_ptcut_max( names_files[i-1]);
	  for (int j =0; j < bin_num ; j++) {

	  
	    zg_ptcut_max<<j<<" "<< numerator_ossf[j]<<" "<<numerator_n_ossf[j]<<" "<<" "<<denominatore_ossf<<" "<<denominatore_n_ossf<<endl;
	  }
	  
	}

      }


    

      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      SF and FR     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>         
      //scal = scal * faxtore_FR;
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  FILLING  HISTOGRAMS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
   
      
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }
  }

  
  









}// end analisi






//==================================================================
double Analysis_eff::maximum(double a, double b){
  double massimo=0;
  massimo = a;
  if (massimo < b) massimo = b; 
  return massimo+1; 
}

//___________________________________________________________________

void Analysis_eff::from_TGraph_to_TH1D (TGraphAsymmErrors *graph, TH1D *histo, int number_point){

  const int numero= number_point;

  double x_graph[numero];
  double y_graph[numero];  
  for (int i =0; i <number_point; i ++){
    x_graph[i]=0;
    y_graph[i]=0;
  }
  for (int i =0; i <number_point; i ++){
    graph -> GetPoint(i, x_graph[i], y_graph[i]);
    histo->SetBinContent (i+1, x_graph[i],  y_graph[i]);

    //cout<<i<<") "<<y_graph[i]<<"  "<< histo->GetBinContent (i+1)<<endl;

  }
}
/*
//==================================================================
double Analysis_mc::fakeWeight(const unsigned ind, const int flavors, const double conePt, const double eta, const bool tight, TH2D* frMap, const unsigned lCount){
unsigned nFO = 0;
double fr[4];
for(unsigned l = 0; l < lCount; ++l){
if(!tight[ind[l]]){
fr[nFO] = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(TMath::Min(conePt[l], 99.), fabs(eta[ind[l]])));
++nFO;
}
}
double weight = 1;
for(unsigned f = 0; f < nFO; ++f){
weight *= fr[f]/(1-fr[f]);
}
if(nFO == 2) weight*= -1;
return weight;
}
*/






//==================================================================
double Analysis_eff::FR_factor(TGraphAsymmErrors *fakeRate_mu[3],
                              TGraphAsymmErrors *fakeRate_e[3],
                              double eta,
                              double flavors,
                              double lptcone
                              ){
    
    
  eta = fabs(eta);
    
    
  TH1D *fakeRate_mu_histo[3];
  TH1D *fakeRate_e_histo[3];
  Double_t newBins_mu1[7] = {5,10, 15, 25, 35, 50, 70};
  Double_t newBins_e1[6] = {10, 15, 25, 35, 50, 70};
  fakeRate_mu_histo[0]= new TH1D("fake_rate_mu_histo_eta1","",6,newBins_mu1);
  fakeRate_mu_histo[1]= new TH1D("fake_rate_mu_histo_eta2","",6,newBins_mu1);
  fakeRate_mu_histo[2]= new TH1D("fake_rate_mu_histo_eta3","",6,newBins_mu1);
  fakeRate_e_histo[0]= new TH1D("fake_rate_e_histo_eta1","",5,newBins_e1);
  fakeRate_e_histo[1]= new TH1D("fake_rate_e_histo_eta2","",5,newBins_e1);
  fakeRate_e_histo[2]= new TH1D("fake_rate_e_histo_eta3","",5,newBins_e1);
    
  for (int i=0; i< 3; i++){
    from_TGraph_to_TH1D(*&fakeRate_mu[i],*&fakeRate_mu_histo[i],6);
    from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],5);
  }
    
  //double momentum = part.Pt() * maximum( 1, iso - 0.1);
  double momentum = lptcone;
    
  double factor=0;
  double factore=0;
  if (flavors == 0){
    if (momentum < 49){
      if (eta < 0.8){
	factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(momentum));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(momentum));
      }//eta1
      else {
	factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(momentum));
      }//eta1
    }// <70
    else {
      if (eta < 0.8){
	factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(45));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(45));
      }//eta1
      else {
	factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(45));
      }//eta1
    }
  }//e
    
  if (flavors == 1){
    if (momentum < 49){
      if (eta < 0.8){
	factore = fakeRate_mu_histo[0]->GetBinContent(fakeRate_mu_histo[0]->FindBin(momentum));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_mu_histo[1]->GetBinContent(fakeRate_mu_histo[1]->FindBin(momentum));
      }//eta1
      else {
	factore = fakeRate_mu_histo[2]->GetBinContent(fakeRate_mu_histo[2]->FindBin(momentum));
      }//eta1
    }// <70
    else {
      if (eta < 0.8){
	factore = fakeRate_mu_histo[0]->GetBinContent(fakeRate_mu_histo[0]->FindBin(45));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_mu_histo[1]->GetBinContent(fakeRate_mu_histo[1]->FindBin(45));
      }//eta1
      else {
	factore = fakeRate_mu_histo[2]->GetBinContent(fakeRate_mu_histo[2]->FindBin(45));
      }//eta1
    }
  }//e
    
  delete fakeRate_mu_histo[0];
  delete fakeRate_mu_histo[1];
  delete fakeRate_mu_histo[2];
  delete fakeRate_e_histo[0];
  delete fakeRate_e_histo[1];
  delete fakeRate_e_histo[2];
    
    
    
  return factore;
    
}





//___________________________________________________________________
void Analysis_eff::class_os(int event_clas[1], int  flavors_3l[3], int  charge_3l[3]){
    
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];


  // 1* = eee
  // 2* = emm
  // 3* = eem
  // 4* = mmm

  


    
  event_clas[0]=-1;


  if (ch_lepton1 == ch_lepton2 && ch_lepton1 == ch_lepton3 && ch_lepton3 == ch_lepton2)   event_clas[0]=-1;


  if (fl_lepton2 == 0 ||  fl_lepton3  == 0 ||  fl_lepton1== 0) {
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 0 ) event_clas[0] = 10; //e e e
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 1 ) event_clas[0] = 30; //e e mu
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 2 ) {
      if (fl_lepton2 == 1 ||  fl_lepton3  == 1 ||  fl_lepton1== 1) event_clas[0]=20; // e mu mu
    }
  }// at least an electron
  else if (fl_lepton2 == 1 &&  fl_lepton3  == 1 &&  fl_lepton1 == 1) {
    event_clas[0] = 40;
  }
  else {
    event_clas[0] =-1;
  }

  if (event_clas[0]  == 30){
    if ((fl_lepton1 == 0 && fl_lepton2 == 0 && fl_lepton3 == 1) ) { // e e mu
      if ((ch_lepton1 + ch_lepton2) == 0) event_clas[0] = 3;
    }
    if (fl_lepton1 == 0 && fl_lepton2 == 1 && fl_lepton3 == 0) { // e mu e
      if ((ch_lepton1 + ch_lepton3) == 0) event_clas[0] = 3;
    }
    if (fl_lepton1 == 1 && fl_lepton2 == 0 && fl_lepton3 == 0 ) { // mu e e
      if ((ch_lepton2 + ch_lepton3) == 0) event_clas[0] = 3;
    }
  }
  
 
  if (event_clas[0]  == 20){
    if ((fl_lepton1 == 1 && fl_lepton2 == 1 && fl_lepton3 == 0) ) { // e e mu
      if ((ch_lepton1 + ch_lepton2) == 0) event_clas[0] = 2;
    }
    if (fl_lepton1 == 1 && fl_lepton2 == 0 && fl_lepton3 == 1 ) { // e mu e
      if ((ch_lepton1 + ch_lepton3) == 0) event_clas[0] = 2;
    }
    if (fl_lepton1 == 0 && fl_lepton2 == 1 && fl_lepton3 == 1 ) { // mu e e
      if  ((ch_lepton2 + ch_lepton3) == 0)event_clas[0] = 2;
    }
  }


  if (event_clas[0]  == 10 ) {
    if ((ch_lepton1 + ch_lepton2 + ch_lepton3) == 1 || (ch_lepton1 + ch_lepton2 + ch_lepton3) == -1)  event_clas[0] =1;
  }

  if (event_clas[0]  == 40 ) {
    if ((ch_lepton1 + ch_lepton2 + ch_lepton3) == 1 || (ch_lepton1 + ch_lepton2 + ch_lepton3) == -1)  event_clas[0] =4;
  }


}



//___________________________________________________________________
void Analysis_eff::ossf_no_ossf(int kind[1],TLorentzVector pair[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int  flavors_3l[3], int  charge_3l[3]){
    
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];
    
    
  kind[0] = -1;

  if (ch_lepton1 == ch_lepton2 && ch_lepton1 == ch_lepton3 && ch_lepton3 == ch_lepton2)   kind[0] = -1;

    
  // OSSF
  if (     ((ch_lepton1 != ch_lepton2)    && (fl_lepton1 == fl_lepton2))  || ((ch_lepton1 != ch_lepton3)   && (fl_lepton1 == fl_lepton3)) || ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) ){ // ossf
    //cout<<"in function where kind is 1: "<<kind[0]<<endl;


    kind[0] = 1;
    double i_m[3]={33333,33333,33333};
    double mass_inv=0;
    int index_inv=100;
    double min_mass=999;
    if ((ch_lepton1 != ch_lepton2)  && (fl_lepton1 == fl_lepton2)) i_m[0]= TMath:: Abs((leep1 + leep2).Mag() - 91.1876);
    if ((ch_lepton1 != ch_lepton3)  && (fl_lepton1 == fl_lepton3)) i_m[1]= TMath:: Abs((leep1 + leep3).Mag() - 91.1876);
    if ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) i_m[2]= TMath:: Abs((leep2 + leep3).Mag() - 91.1876);
    for (int i =0; i < 3; i++){
      if (i_m[i] == 33333) continue;
      mass_inv = i_m[i];
      if (min_mass > mass_inv ){
	min_mass = mass_inv;
	index_inv = i;
      }
    }
    if (index_inv == 0) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[2].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());

    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());

    }
  }// end ossf
  // No_OSSF
  else if (   ((ch_lepton1 + ch_lepton2) == 0  )  || ((ch_lepton1 + ch_lepton3) == 0   ) || ((ch_lepton3 + ch_lepton2) == 0   )   ){
    //cout<<"in function where kind is 0: "<<kind[0]<<endl;
    kind[0] = 0;
    double i_m[3]={33333,33333,33333};
    double mass_inv=0;
    int index_inv=100;
    double min_mass=999;
    if ((ch_lepton1 != ch_lepton2) ) i_m[0]= TMath:: Abs((leep1 + leep2).Mag() - 91.1876);
    if ((ch_lepton1 != ch_lepton3)  ) i_m[1]= TMath:: Abs((leep1 + leep3).Mag() - 91.1876);
    if ((ch_lepton2 != ch_lepton3) ) i_m[2]= TMath:: Abs((leep2 + leep3).Mag() - 91.1876);
    for (int i =0; i < 3; i++){
      if (i_m[i] == 33333) continue;
      mass_inv = i_m[i];
      if (min_mass > mass_inv ){
	min_mass = mass_inv;
	index_inv = i;
      }
    }
    if (index_inv == 0) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[2].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());

    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());

    }
        
        
  }//end no-ossf
    
  /*
    cout<< ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   "<<kind<<endl;
    cout<<"1: "<<ch_lepton1<<"  "<<fl_lepton1<< "  "<< leep1.Pt()<< endl;
    cout<<"2: "<<ch_lepton2<<"  "<<fl_lepton2<< "  "<< leep2.Pt()<< endl;
    cout<<"3: "<<ch_lepton3<<"  "<<fl_lepton3<< "  "<< leep3.Pt()<< endl;
    cout<<"pair 1: "<<pair1.Pt()<< endl;
    cout<<"pair 2: "<<pair2.Pt()<< endl;
  */
  //cout<<"in function "<<kind[0]<<endl;
}





//___________________________________________________________________
void Analysis_eff::fr_selection(int number, TLorentzVector lepton_fake_order[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int index_leptons[3],  int flavor_leptons[3], int origin_leptons[3],int index_3l[3],  int flavor_3l[3], int origin_3l[3]){

  lepton_fake_order[0].SetPtEtaPhiE(0,0,0,0);
  lepton_fake_order[1].SetPtEtaPhiE(0,0,0,0);
  lepton_fake_order[2].SetPtEtaPhiE(0,0,0,0);
  for(int i =0; i< 3; i++){
    index_leptons[i]=  -5;
    flavor_leptons[i]= -5;
    origin_leptons[i]= -5;
  }

  if (number == 3) {
    lepton_fake_order[0].SetPtEtaPhiE(leep1.Pt(),leep1.Eta(),leep1.Phi(),leep1.E());
    lepton_fake_order[1].SetPtEtaPhiE(leep2.Pt(),leep2.Eta(),leep2.Phi(),leep2.E());
    lepton_fake_order[2].SetPtEtaPhiE(leep3.Pt(),leep3.Eta(),leep3.Phi(),leep3.E());
    index_leptons[0]=index_3l[0];
    index_leptons[1]=index_3l[1];
    index_leptons[2]=index_3l[2];
    flavor_leptons[0]=flavor_3l[0];
    flavor_leptons[1]=flavor_3l[1];
    flavor_leptons[2]=flavor_3l[2];
    origin_leptons[0]=origin_3l[0];
    origin_leptons[1]=origin_3l[1];
    origin_leptons[2]=origin_3l[2];
  }
  if (number == 0) {
    lepton_fake_order[0].SetPtEtaPhiE(leep1.Pt(),leep1.Eta(),leep1.Phi(),leep1.E());
    lepton_fake_order[1].SetPtEtaPhiE(leep2.Pt(),leep2.Eta(),leep2.Phi(),leep2.E());
    lepton_fake_order[2].SetPtEtaPhiE(leep3.Pt(),leep3.Eta(),leep3.Phi(),leep3.E());
    index_leptons[0]=index_3l[0];
    index_leptons[1]=index_3l[1];
    index_leptons[2]=index_3l[2];
    flavor_leptons[0]=flavor_3l[0];
    flavor_leptons[1]=flavor_3l[1];
    flavor_leptons[2]=flavor_3l[2];
    origin_leptons[0]=origin_3l[0];
    origin_leptons[1]=origin_3l[1];
    origin_leptons[2]=origin_3l[2];
  }

  

}//end fr
