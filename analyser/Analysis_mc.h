#ifndef ANALYSIS_MC_H
#define ANALYSIS_MC_H
#include "TObject.h"

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




class Analysis_mc : public TObject {
    
 public:
    
    
  Analysis_mc();
  Analysis_mc(string FileNameTree_in);
  virtual ~Analysis_mc();
    
  void analisi(    int num_histo_kin	   
		   );
    
 double maximum(double a, double b);
    void from_TGraph_to_TH1D (TGraphAsymmErrors *graph, TH1D *histo, int number_point);


  double FR_factor(TGraphAsymmErrors *fakeRate_mu[3],
			      TGraphAsymmErrors *fakeRate_e[3],
			      double eta,
			      double flavors,
			      double lptcone
		   );

  void class_os(int event_clas[1], int  flavors_3l[3], int  charge_3l[3]);
  void ossf_no_ossf(int kind[1],TLorentzVector pair1[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int  flavors_3l[3], int  charge_3l[3]);

  void fr_selection (int number, TLorentzVector lepton_fake_order[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int index_leptons[3],  int flavor_leptons[3], int origin_leptons[3],int index_3l[3],  int flavor_3l[3], int origin_3l[3]);


void printProgress(double progress) ;
 private:
    
    
    
  ClassDef(Analysis_mc,1) 
    };

#endif



