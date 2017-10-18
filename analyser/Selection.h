#ifndef SELECTION_H
#define SELECTION_H
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




class Selection : public TObject {
    
public:
    
    
    Selection();
    Selection(TString FileNameTree_in);
    virtual ~Selection();
    


    void selezione( Int_t _gen_nL,
		   Int_t           _nL,
		    Double_t        _lPt[7], 
		    Double_t        _lEta[7],   
		    Double_t        _lE[7],   
		    Bool_t          _lPOGMedium[7],   
		    Bool_t          _lPOGLoose[7],   
		    Int_t           _lFlavor[7],   
		    Int_t        _lCharge[7],   
		    Double_t        _relIso[7],
		    Double_t        _ipPV[7],   
		    Double_t        _ipZPV[7],   
		    Double_t        _3dIP[7],   
		    Double_t        _3dIPsig[7], 
		    Bool_t          _lHNLoose[7],   
		    Bool_t          _lHNFO[7],   
		    Bool_t          _lHNTight[7],  
		    Float_t        _lElectronMva[7],   
		    Bool_t           _lIsPrompt[7],
		    double index[3],
		    double number_tight[1],
		    double number_prompt[1],
		    int skip_event[1],
			      
		    int effsam,
		    TGraphAsymmErrors *fakeRate_mu[3],
		    TGraphAsymmErrors *fakeRate_e[3],
		    double faxtore[1],
		    Double_t        _gen_lPt[7], 
		    Double_t        _gen_lEta[7],
		    Bool_t        _gen_lIsPrompt[22],
		    Int_t            _gen_lFlavor[20]
		    


		    );
    
    double maximum(double a, double b);
    double find_eta(double b);
  void check_low_mass_pt(double check[1], double pt_cone_1,double pt_cone_2,double pt_cone_3,double flavor_3);
    
    
  void from_TGraph_to_TH1D (TGraphAsymmErrors *graph, TH1D *histo, int number_point);


  double FR_factor(TGraphAsymmErrors *fakeRate_mu[3],
		   TGraphAsymmErrors *fakeRate_e[3],
		   double eta,
		   double flavors,
		   double lptcone
		   );
    
 private:
    
    
    
    ClassDef(Selection,1) 
};

#endif



