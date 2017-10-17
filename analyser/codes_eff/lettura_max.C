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
#include <TGaxis.h>


#include <THStack.h>
#include <TPaveText.h>
#include <THStack.h>


using namespace std;
const int num_histo_kin= 3;
const int num_histo= 12;


const Int_t nDati=100;

//__________________________________________________________
void cosmesi(Int_t numero,TGraphErrors *histo_positive,char nameassex[] = "asse x",char nameassey[] = "asse y"){
    
    histo_positive->GetXaxis()->SetLimits(0,0.5);
    histo_positive->GetYaxis()->SetNdivisions(506);
    histo_positive->GetXaxis()->SetNdivisions(506);
    histo_positive->SetMinimum( 0);
    histo_positive->SetMaximum( 1.2);
    histo_positive->SetTitle("");
    histo_positive->GetXaxis()->SetTitle(nameassex);
    histo_positive->GetXaxis()->SetTitleOffset(0.80);
    histo_positive->GetXaxis()->SetTitleSize(0.03);
    histo_positive->GetXaxis() ->SetTitleFont(132);
    histo_positive->GetXaxis()->SetLabelSize(0.03);
    histo_positive->GetXaxis()->SetLabelOffset(0.00);
    histo_positive->GetXaxis()->SetLabelFont(132);
    
    
    
    
    
    histo_positive->GetYaxis()->SetLabelFont(132);
    
    histo_positive->GetYaxis()->SetNdivisions(506);
    histo_positive->GetYaxis()->SetTitle(nameassey);
    histo_positive->	GetYaxis ()->SetTitleOffset(0.75);
    histo_positive->	GetYaxis()->SetTitleSize(0.06);
    histo_positive->	GetYaxis() ->SetTitleFont(132);
    histo_positive->GetYaxis()->SetLabelOffset(0.005);
    histo_positive->GetYaxis()->SetLabelSize(0.04);
    histo_positive->GetYaxis()-> CenterTitle();
    if (numero == 11){
        histo_positive->SetMarkerStyle(25);
        histo_positive->SetMarkerColor(kRed);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 22 ){
        histo_positive->SetMarkerStyle(25);
        histo_positive->SetMarkerColor(kBlue);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 33){
        histo_positive->SetMarkerStyle(25);
        histo_positive->SetMarkerColor(kGreen);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 44 ){
        histo_positive->SetMarkerStyle(25);
        histo_positive->SetMarkerColor(kOrange );
        histo_positive->SetMarkerSize(0.2);
    }
    
    if (numero == 1){
        histo_positive->SetMarkerStyle(21);
        histo_positive->SetMarkerColor(kRed-4);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 2 ){
        histo_positive->SetMarkerStyle(21);
        histo_positive->SetMarkerColor(kBlue);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 3){
        histo_positive->SetMarkerStyle(21);
        histo_positive->SetMarkerColor(8);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 4 ){
        histo_positive->SetMarkerStyle(21);
        histo_positive->SetMarkerColor(kOrange);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 5){
        histo_positive->SetMarkerStyle(20);
        histo_positive->SetMarkerColor(kBlack);
        histo_positive->SetMarkerSize(0.2);
    }
    if (numero == 6){
        histo_positive->SetMarkerStyle(4);
        histo_positive->SetMarkerColor(38);
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 7){
        histo_positive->SetMarkerStyle(4);
        histo_positive->SetMarkerColor(kCyan +2);
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 8){
        histo_positive->SetMarkerStyle(4);
        histo_positive->SetMarkerColor(kCyan);
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 9){
        histo_positive->SetMarkerStyle(4);
        histo_positive->SetMarkerColor(kMagenta);
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 66){
        histo_positive->SetMarkerStyle(3);
        histo_positive->SetMarkerColor(38);
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 77){
        histo_positive->SetMarkerStyle(3);
        histo_positive->SetMarkerColor(kCyan +2);
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 88){
        histo_positive->SetMarkerStyle(3);
        histo_positive->SetMarkerColor(kCyan );
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 99){
        histo_positive->SetMarkerStyle(3);
        histo_positive->SetMarkerColor(kMagenta );
        histo_positive->SetMarkerSize(0.5);
    }
    if (numero == 10){
        histo_positive->SetMarkerStyle(4);
        histo_positive->SetMarkerColor(kMagenta-1);
        histo_positive->SetMarkerSize(0.5);
    }
    
    if (numero == 11){
        histo_positive->SetMarkerStyle(21);
        histo_positive->SetMarkerColor(kMagenta+2);
        histo_positive->SetMarkerSize(0.2);
    }
    
    if (numero == 0){
        histo_positive->SetMarkerStyle(21);
        histo_positive->SetMarkerColor(kBlack);
        histo_positive->SetMarkerSize(0.5);
    }
    
    
    
}//end cosmesi




void lett (string inputRootFile= "file_mva_gen1.root");


// ********************************************************************
void lett(string inputRootFile){
    Double_t pigreco= TMath::ACos(-1);
    cout<<"in analisi"<<endl;
    cout<<"---------------------------"<<endl;
    const int npoint=100;
    static TGraphErrors* g_pt_dy_ossf;
    static TGraphErrors* g_pt_diboson_ossf;
    static TGraphErrors* g_pt_zz_ossf;
    static TGraphErrors* g_pt_triboson_ossf;
    static TGraphErrors* g_pt_tx_ossf;
    static TGraphErrors* g_pt_xg_ossf;
    static TGraphErrors* g_pt_signal_1_ossf;
    static TGraphErrors* g_pt_signal_2_ossf;
    static TGraphErrors* g_pt_signal_5_ossf;
    static TGraphErrors* g_pt_signal_4_ossf;
    static TGraphErrors* g_pt_signal_6_ossf;
    static TGraphErrors* g_pt_signal_7_ossf;
    static TGraphErrors* g_pt_signal_8_ossf;
    static TGraphErrors* g_pt_signal_9_ossf;
    static TGraphErrors* g_pt_signal_3_ossf;
    
    static TGraphErrors* g_pt_tot_ossf;
    
    Double_t num_dy_ossf[npoint];
    Double_t num_diboson_ossf[npoint];
    Double_t num_zz_ossf[npoint];
    double num_triboson_ossf[npoint];
    double num_triboson1_ossf[npoint];
    double num_triboson2_ossf[npoint];
    double num_triboson3_ossf[npoint];
    double num_triboson4_ossf[npoint];
    double num_tx_ossf[npoint];
    double num_tx1_ossf[npoint];
    double num_tx2_ossf[npoint];
    double num_tx3_ossf[npoint];
    
    double num_dy1_ossf[npoint];
    double num_dy2_ossf[npoint];
    double num_dy3_ossf[npoint];
    double num_dy4_ossf[npoint];
    double num_dy5_ossf[npoint];
    
    double num_diboson1_ossf[npoint];
    double num_diboson2_ossf[npoint];
    double num_diboson3_ossf[npoint];
    double num_tx12_ossf[npoint];
    double num_xg_ossf[npoint];
    double num_xg1_ossf[npoint];
    double num_xg2_ossf[npoint];
    double num_xg3_ossf[npoint];
    Double_t num_signal_1_ossf[npoint];
    Double_t num_signal_2_ossf[npoint];
    Double_t num_signal_5_ossf[npoint];
    Double_t num_signal_4_ossf[npoint];
    Double_t num_signal_6_ossf[npoint];
    Double_t num_signal_7_ossf[npoint];
    Double_t num_signal_8_ossf[npoint];
    Double_t num_signal_9_ossf[npoint];
    Double_t num_signal_3_ossf[npoint];
    
    Double_t den_zz_ossf[npoint];
    double den_triboson_ossf[npoint];
    double den_triboson1_ossf[npoint];
    double den_triboson2_ossf[npoint];
    double den_triboson3_ossf[npoint];
    double den_triboson4_ossf[npoint];
    double den_tx_ossf[npoint];
    double den_tx1_ossf[npoint];
    double den_tx2_ossf[npoint];
    double den_tx3_ossf[npoint];
    double den_dy_ossf[npoint];
    double den_dy1_ossf[npoint];
    double den_dy2_ossf[npoint];
    double den_dy3_ossf[npoint];
    double den_dy4_ossf[npoint];
    double den_dy5_ossf[npoint];
    double den_diboson_ossf[npoint];
    double den_diboson1_ossf[npoint];
    double den_diboson2_ossf[npoint];
    double den_diboson3_ossf[npoint];
    double den_tx12_ossf[npoint];
    double den_xg_ossf[npoint];
    double den_xg1_ossf[npoint];
    double den_xg2_ossf[npoint];
    double den_xg3_ossf[npoint];
    Double_t den_signal_1_ossf[npoint];
    Double_t den_signal_2_ossf[npoint];
    Double_t den_signal_5_ossf[npoint];
    Double_t den_signal_4_ossf[npoint];
    Double_t den_signal_6_ossf[npoint];
    Double_t den_signal_7_ossf[npoint];
    Double_t den_signal_8_ossf[npoint];
    Double_t den_signal_9_ossf[npoint];
    Double_t den_signal_3_ossf[npoint];
    
    
    Double_t zz_ossf[npoint];
    double triboson_ossf[npoint];
    double triboson1_ossf[npoint];
    double triboson2_ossf[npoint];
    double triboson3_ossf[npoint];
    double triboson4_ossf[npoint];
    double tx_ossf[npoint];
    double tx1_ossf[npoint];
    double tx2_ossf[npoint];
    double tx3_ossf[npoint];
    double dy_ossf[npoint];
    double dy1_ossf[npoint];
    double dy2_ossf[npoint];
    double dy3_ossf[npoint];
    double dy4_ossf[npoint];
    double dy5_ossf[npoint];
    double diboson_ossf[npoint];
    double diboson1_ossf[npoint];
    double diboson2_ossf[npoint];
    double diboson3_ossf[npoint];
    double tx12_ossf[npoint];
    double xg_ossf[npoint];
    double xg1_ossf[npoint];
    double xg2_ossf[npoint];
    double xg3_ossf[npoint];
    Double_t signal_1_ossf[npoint];
    Double_t signal_2_ossf[npoint];
    Double_t signal_5_ossf[npoint];
    Double_t signal_4_ossf[npoint];
    Double_t signal_6_ossf[npoint];
    Double_t signal_7_ossf[npoint];
    Double_t signal_8_ossf[npoint];
    Double_t signal_9_ossf[npoint];
    Double_t signal_3_ossf[npoint];
    
    
    static TGraphErrors* g_pt_diboson_n_ossf;
    static TGraphErrors* g_pt_dy_n_ossf;
    static TGraphErrors* g_pt_zz_n_ossf;
    static TGraphErrors* g_pt_triboson_n_ossf;
    static TGraphErrors* g_pt_tx_n_ossf;
    static TGraphErrors* g_pt_xg_n_ossf;
    static TGraphErrors* g_pt_signal_1_n_ossf;
    static TGraphErrors* g_pt_signal_2_n_ossf;
    static TGraphErrors* g_pt_signal_5_n_ossf;
    static TGraphErrors* g_pt_signal_4_n_ossf;
    static TGraphErrors* g_pt_signal_6_n_ossf;
    static TGraphErrors* g_pt_signal_7_n_ossf;
    static TGraphErrors* g_pt_signal_8_n_ossf;
    static TGraphErrors* g_pt_signal_9_n_ossf;
    static TGraphErrors* g_pt_signal_3_n_ossf;
    
    static TGraphErrors* g_pt_tot_n_ossf;
    
    Double_t num_zz_n_ossf[npoint];
    double num_triboson_n_ossf[npoint];
    double num_triboson1_n_ossf[npoint];
    double num_triboson2_n_ossf[npoint];
    double num_triboson3_n_ossf[npoint];
    double num_triboson4_n_ossf[npoint];
    double num_tx_n_ossf[npoint];
    double num_tx1_n_ossf[npoint];
    double num_tx2_n_ossf[npoint];
    double num_tx3_n_ossf[npoint];
    double num_dy_n_ossf[npoint];
    double num_dy1_n_ossf[npoint];
    double num_dy2_n_ossf[npoint];
    double num_dy3_n_ossf[npoint];
    double num_dy4_n_ossf[npoint];
    double num_dy5_n_ossf[npoint];
    double num_diboson_n_ossf[npoint];
    double num_diboson1_n_ossf[npoint];
    double num_diboson2_n_ossf[npoint];
    double num_diboson3_n_ossf[npoint];
    double num_tx12_n_ossf[npoint];
    double num_xg_n_ossf[npoint];
    double num_xg1_n_ossf[npoint];
    double num_xg2_n_ossf[npoint];
    double num_xg3_n_ossf[npoint];
    Double_t num_signal_1_n_ossf[npoint];
    Double_t num_signal_2_n_ossf[npoint];
    Double_t num_signal_5_n_ossf[npoint];
    Double_t num_signal_4_n_ossf[npoint];
    Double_t num_signal_6_n_ossf[npoint];
    Double_t num_signal_7_n_ossf[npoint];
    Double_t num_signal_8_n_ossf[npoint];
    Double_t num_signal_9_n_ossf[npoint];
    Double_t num_signal_3_n_ossf[npoint];
    
    Double_t den_zz_n_ossf[npoint];
    double den_triboson_n_ossf[npoint];
    double den_triboson1_n_ossf[npoint];
    double den_triboson2_n_ossf[npoint];
    double den_triboson3_n_ossf[npoint];
    double den_triboson4_n_ossf[npoint];
    double den_tx_n_ossf[npoint];
    double den_tx1_n_ossf[npoint];
    double den_tx2_n_ossf[npoint];
    double den_tx3_n_ossf[npoint];
    double den_dy_n_ossf[npoint];
    double den_dy1_n_ossf[npoint];
    double den_dy2_n_ossf[npoint];
    double den_dy3_n_ossf[npoint];
    double den_dy4_n_ossf[npoint];
    double den_dy5_n_ossf[npoint];
    double den_diboson_n_ossf[npoint];
    double den_diboson1_n_ossf[npoint];
    double den_diboson2_n_ossf[npoint];
    double den_diboson3_n_ossf[npoint];
    double den_tx12_n_ossf[npoint];
    double den_xg_n_ossf[npoint];
    double den_xg1_n_ossf[npoint];
    double den_xg2_n_ossf[npoint];
    double den_xg3_n_ossf[npoint];
    Double_t den_signal_1_n_ossf[npoint];
    Double_t den_signal_2_n_ossf[npoint];
    Double_t den_signal_5_n_ossf[npoint];
    Double_t den_signal_4_n_ossf[npoint];
    Double_t den_signal_6_n_ossf[npoint];
    Double_t den_signal_7_n_ossf[npoint];
    Double_t den_signal_8_n_ossf[npoint];
    Double_t den_signal_9_n_ossf[npoint];
    Double_t den_signal_3_n_ossf[npoint];
    
    
    Double_t zz_n_ossf[npoint];
    double triboson_n_ossf[npoint];
    double triboson1_n_ossf[npoint];
    double triboson2_n_ossf[npoint];
    double triboson3_n_ossf[npoint];
    double triboson4_n_ossf[npoint];
    double tx_n_ossf[npoint];
    double tx1_n_ossf[npoint];
    double tx2_n_ossf[npoint];
    double tx3_n_ossf[npoint];
    double dy_n_ossf[npoint];
    double dy1_n_ossf[npoint];
    double dy2_n_ossf[npoint];
    double dy3_n_ossf[npoint];
    double dy4_n_ossf[npoint];
    double dy5_n_ossf[npoint];
    double diboson_n_ossf[npoint];
    double diboson1_n_ossf[npoint];
    double diboson2_n_ossf[npoint];
    double diboson3_n_ossf[npoint];
    double tx12_n_ossf[npoint];
    double xg_n_ossf[npoint];
    double xg1_n_ossf[npoint];
    double xg2_n_ossf[npoint];
    double xg3_n_ossf[npoint];
    Double_t signal_1_n_ossf[npoint];
    Double_t signal_2_n_ossf[npoint];
    Double_t signal_5_n_ossf[npoint];
    Double_t signal_4_n_ossf[npoint];
    Double_t signal_6_n_ossf[npoint];
    Double_t signal_7_n_ossf[npoint];
    Double_t signal_8_n_ossf[npoint];
    Double_t signal_9_n_ossf[npoint];
    Double_t signal_3_n_ossf[npoint];
    
    
    
    double back_tot_ossf[npoint];
    double back_tot_n_ossf[npoint];
    double back_den_tot_ossf[npoint];
    double back_den_tot_n_ossf[npoint];
    double back_value_tot_ossf[npoint];
    double back_value_tot_n_ossf[npoint];
    
    
    double back_sum_ossf[npoint];
    double back_sum_n_ossf[npoint];
    double back_den_ossf[npoint];
    double back_den_n_ossf[npoint];
    double back_value_ossf[npoint];
    double back_value_n_ossf[npoint];
    
    Double_t signal_bgk_1_n_ossf[npoint];
    Double_t signal_bgk_2_n_ossf[npoint];
    Double_t signal_bgk_5_n_ossf[npoint];
    Double_t signal_bgk_4_n_ossf[npoint];
    Double_t signal_bgk_6_n_ossf[npoint];
    Double_t signal_bgk_7_n_ossf[npoint];
    Double_t signal_bgk_8_n_ossf[npoint];
    Double_t signal_bgk_9_n_ossf[npoint];
    Double_t signal_bgk_3_n_ossf[npoint];
    Double_t signal_bgk_1_ossf[npoint];
    Double_t signal_bgk_2_ossf[npoint];
    Double_t signal_bgk_5_ossf[npoint];
    Double_t signal_bgk_4_ossf[npoint];
    Double_t signal_bgk_6_ossf[npoint];
    Double_t signal_bgk_7_ossf[npoint];
    Double_t signal_bgk_8_ossf[npoint];
    Double_t signal_bgk_9_ossf[npoint];
    Double_t signal_bgk_3_ossf[npoint];
    
    static TGraphErrors* g_pt_s_b_1_ossf;
    static TGraphErrors* g_pt_s_b_2_ossf;
    static TGraphErrors* g_pt_s_b_5_ossf;
    static TGraphErrors* g_pt_s_b_4_ossf;
    static TGraphErrors* g_pt_s_b_6_ossf;
    static TGraphErrors* g_pt_s_b_7_ossf;
    static TGraphErrors* g_pt_s_b_8_ossf;
    static TGraphErrors* g_pt_s_b_9_ossf;
    static TGraphErrors* g_pt_s_b_3_ossf;
    
    static TGraphErrors* g_pt_s_b_1_n_ossf;
    static TGraphErrors* g_pt_s_b_2_n_ossf;
    static TGraphErrors* g_pt_s_b_5_n_ossf;
    static TGraphErrors* g_pt_s_b_4_n_ossf;
    static TGraphErrors* g_pt_s_b_6_n_ossf;
    static TGraphErrors* g_pt_s_b_7_n_ossf;
    static TGraphErrors* g_pt_s_b_8_n_ossf;
    static TGraphErrors* g_pt_s_b_9_n_ossf;
    static TGraphErrors* g_pt_s_b_3_n_ossf;
    
    
    
    
    double x_point[npoint];
    for (int i =0; i< npoint; i++){
        
        back_sum_ossf[i]= 0;
        back_sum_n_ossf[i]= 0;
        back_den_ossf[i]= 0;
        back_den_n_ossf[i]= 0;
        back_value_ossf[i]= 0;
        back_value_n_ossf[i]= 0;
        num_diboson_ossf[i] = 0.;
        num_zz_ossf[i] = 0.;
        num_triboson_ossf[i] = 0.;
        num_triboson1_ossf[i] = 0.;
        num_triboson2_ossf[i] = 0.;
        num_triboson3_ossf[i] = 0.;
        num_triboson4_ossf[i] = 0.;
        num_tx_ossf[i] = 0.;
        num_tx1_ossf[i] = 0.;
        num_tx2_ossf[i] = 0.;
        num_tx3_ossf[i] = 0.;
        num_dy1_ossf[i] = 0.;
        num_dy2_ossf[i] = 0.;
        num_dy3_ossf[i] = 0.;
        num_dy4_ossf[i] = 0.;
        num_dy5_ossf[i] = 0.;
        num_diboson1_ossf[i] = 0.;
        num_diboson2_ossf[i] = 0.;
        num_diboson3_ossf[i] = 0.;
        num_tx12_ossf[i] = 0.;
        num_xg_ossf[i] = 0.;
        num_xg1_ossf[i] = 0.;
        num_xg2_ossf[i] = 0.;
        num_xg3_ossf[i] = 0.;
        num_signal_6_ossf[i] = 0.;
        num_signal_7_ossf[i] = 0.;
        num_signal_8_ossf[i] = 0.;
        num_signal_9_ossf[i] = 0.;
        num_signal_3_ossf[i] = 0.;
        den_diboson_ossf[i] = 0.;
        den_zz_ossf[i] = 0.;
        den_triboson_ossf[i] = 0.;
        den_triboson1_ossf[i] = 0.;
        den_triboson2_ossf[i] = 0.;
        den_triboson3_ossf[i] = 0.;
        den_triboson4_ossf[i] = 0.;
        den_tx_ossf[i] = 0.;
        den_tx1_ossf[i] = 0.;
        den_tx2_ossf[i] = 0.;
        den_tx3_ossf[i] = 0.;
        den_dy1_ossf[i] = 0.;
        den_dy2_ossf[i] = 0.;
        den_dy3_ossf[i] = 0.;
        den_dy4_ossf[i] = 0.;
        den_dy5_ossf[i] = 0.;
        den_diboson1_ossf[i] = 0.;
        den_diboson2_ossf[i] = 0.;
        den_diboson3_ossf[i] = 0.;
        den_tx12_ossf[i] = 0.;
        den_xg_ossf[i] = 0.;
        den_xg1_ossf[i] = 0.;
        den_xg2_ossf[i] = 0.;
        den_xg3_ossf[i] = 0.;
        den_signal_6_ossf[i] = 0.;
        den_signal_7_ossf[i] = 0.;
        den_signal_8_ossf[i] = 0.;
        den_signal_9_ossf[i] = 0.;
        den_signal_3_ossf[i] = 0.;
        num_diboson_n_ossf[i] = 0.;
        num_zz_n_ossf[i] = 0.;
        num_triboson_n_ossf[i] = 0.;
        num_triboson1_n_ossf[i] = 0.;
        num_triboson2_n_ossf[i] = 0.;
        num_triboson3_n_ossf[i] = 0.;
        num_triboson4_n_ossf[i] = 0.;
        num_tx_n_ossf[i] = 0.;
        num_tx1_n_ossf[i] = 0.;
        num_tx2_n_ossf[i] = 0.;
        num_tx3_n_ossf[i] = 0.;
        num_dy1_n_ossf[i] = 0.;
        num_dy2_n_ossf[i] = 0.;
        num_dy3_n_ossf[i] = 0.;
        num_dy4_n_ossf[i] = 0.;
        num_dy5_n_ossf[i] = 0.;
        num_diboson1_n_ossf[i] = 0.;
        num_diboson2_n_ossf[i] = 0.;
        num_diboson3_n_ossf[i] = 0.;
        num_tx12_n_ossf[i] = 0.;
        num_xg_n_ossf[i] = 0.;
        num_xg1_n_ossf[i] = 0.;
        num_xg2_n_ossf[i] = 0.;
        num_xg3_n_ossf[i] = 0.;
        num_signal_6_n_ossf[i] = 0.;
        num_signal_7_n_ossf[i] = 0.;
        num_signal_8_n_ossf[i] = 0.;
        num_signal_9_n_ossf[i] = 0.;
        num_signal_3_n_ossf[i] = 0.;
        den_diboson_n_ossf[i] = 0.;
        den_zz_n_ossf[i] = 0.;
        den_triboson_n_ossf[i] = 0.;
        den_triboson1_n_ossf[i] = 0.;
        den_triboson2_n_ossf[i] = 0.;
        den_triboson3_n_ossf[i] = 0.;
        den_triboson4_n_ossf[i] = 0.;
        den_tx_n_ossf[i] = 0.;
        den_tx1_n_ossf[i] = 0.;
        den_tx2_n_ossf[i] = 0.;
        den_tx3_n_ossf[i] = 0.;
        den_dy1_n_ossf[i] = 0.;
        den_dy2_n_ossf[i] = 0.;
        den_dy3_n_ossf[i] = 0.;
        den_dy4_n_ossf[i] = 0.;
        den_dy5_n_ossf[i] = 0.;
        den_diboson1_n_ossf[i] = 0.;
        den_diboson2_n_ossf[i] = 0.;
        den_diboson3_n_ossf[i] = 0.;
        den_tx12_n_ossf[i] = 0.;
        den_xg_n_ossf[i] = 0.;
        den_xg1_n_ossf[i] = 0.;
        den_xg2_n_ossf[i] = 0.;
        den_xg3_n_ossf[i] = 0.;
        den_signal_6_n_ossf[i] = 0.;
        den_signal_7_n_ossf[i] = 0.;
        den_signal_8_n_ossf[i] = 0.;
        den_signal_9_n_ossf[i] = 0.;
        den_signal_3_n_ossf[i] = 0.;
    }
    
    const int nSamples= 25;
    
    const TString names_files[nSamples] = { "txt_file/trilMaj6.txt","txt_file/trilMaj7.txt", "txt_file/trilMaj8.txt","txt_file/trilMaj9.txt", "txt_file/trilMaj3.txt",
        
        
        "txt_file/ZZ_H.txt","txt_file/X_gamma.txt",
        
        
        "txt_file/triboson1.txt","txt_file/triboson2.txt","txt_file/triboson3.txt",
        
        
        "txt_file/T_X1.txt", "txt_file/T_X2.txt", "txt_file/T_X3.txt",
        
        
        "txt_file/DY_TTbar1.txt", "txt_file/DY_TTbar2.txt", "txt_file/DY_TTbar3.txt", "txt_file/DY_TTbar4.txt", "txt_file/DY_TTbar5.txt",
        
        "txt_file/trilMaj1.txt", "txt_file/trilMaj2.txt","txt_file/trilMaj5.txt","txt_file/trilMaj4.txt",
        
        "txt_file/diboson1.txt","txt_file/diboson2.txt","txt_file/diboson3.txt"};
    
    
    
    
    
    
    
    
    
    ifstream ffsig6(names_files[0]);
    if(!ffsig6){
        cout<<"Il file non esiste 1"<<endl;
        return;
    }
    int sig6 = 0;
    while((ffsig6>>x_point[sig6]>>num_signal_6_ossf[sig6] >>num_signal_6_n_ossf[sig6]>>den_signal_6_ossf[sig6] >>den_signal_6_n_ossf[sig6])){
        sig6++;
    }
    ffsig6.close();
    
    ifstream ffsig7(names_files[1]);
    if(!ffsig7){
        cout<<"Il file non esiste2 "<<endl;
        return;
    }
    int sig7 = 0;
    while((ffsig7>>x_point[sig7]>>num_signal_7_ossf[sig7] >>num_signal_7_n_ossf[sig7]>>den_signal_7_ossf[sig7] >>den_signal_7_n_ossf[sig7])){
        sig7++;
    }
    ffsig7.close();
    
    ifstream ffsig8(names_files[2]);
    if(!ffsig8){
        cout<<"Il file non esiste 3"<<endl;
        return;
    }
    int sig8 = 0;
    while((ffsig8>>x_point[sig8]>>num_signal_8_ossf[sig8] >>num_signal_8_n_ossf[sig8]>>den_signal_8_ossf[sig8] >>den_signal_8_n_ossf[sig8])){
        sig8++;
    }
    ffsig8.close();
    
    ifstream ffsig9(names_files[3]);
    if(!ffsig9){
        cout<<"Il file non esiste 4"<<endl;
        return;
    }
    int sig9 = 0;
    while((ffsig9>>x_point[sig9]>>num_signal_9_ossf[sig9] >>num_signal_9_n_ossf[sig9]>>den_signal_9_ossf[sig9] >>den_signal_9_n_ossf[sig9])){
        sig9++;
    }
    ffsig9.close();
    
    ifstream ffsig3(names_files[4]);
    if(!ffsig3){
        cout<<"Il file non esiste5"<<endl;
        return;
    }
    int sig3 = 0;
    while((ffsig3>>x_point[sig3]>>num_signal_3_ossf[sig3] >>num_signal_3_n_ossf[sig3]>>den_signal_3_ossf[sig3] >>den_signal_3_n_ossf[sig3])){
        sig3++;
    }
    ffsig3.close();
    
    
    
    ifstream ffzz(names_files[5]);
    if(!ffzz){
        cout<<"Il file non esiste6 "<<endl;
        return;
    }
    int zz = 0;
    while((ffzz>>x_point[zz]>>num_zz_ossf[zz] >>num_zz_n_ossf[zz]>>den_zz_ossf[zz] >>den_zz_n_ossf[zz])){
        zz++;
    }
    ffzz.close();
    
    
    ifstream ffxg1(names_files[6]);
    if(!ffxg1){
        cout<<"Il file non esiste7 "<<endl;
        return;
    }
    int xg1 = 0;
    while((ffxg1>>x_point[xg1]>>num_xg1_ossf[xg1] >>num_xg1_n_ossf[xg1]>>den_xg1_ossf[xg1] >>den_xg1_n_ossf[xg1])){
        xg1++;
    }
    ffxg1.close();
    
    
    
    
    ifstream fftriboson1(names_files[7]);
    if(!fftriboson1){
        cout<<"Il file non esiste 8"<<endl;
        return;
    }
    int triboson1 = 0;
    while((fftriboson1>>x_point[triboson1]>>num_triboson1_ossf[triboson1] >>num_triboson1_n_ossf[triboson1]>>den_triboson1_ossf[triboson1] >>den_triboson1_n_ossf[triboson1])){
        triboson1++;
    }
    fftriboson1.close();
    
    ifstream fftriboson2(names_files[8]);
    if(!fftriboson2){
        cout<<"Il file non esiste 9"<<endl;
        return;
    }
    int triboson2 = 0;
    while((fftriboson2>>x_point[triboson2]>>num_triboson2_ossf[triboson2] >>num_triboson2_n_ossf[triboson2]>>den_triboson2_ossf[triboson2] >>den_triboson2_n_ossf[triboson2])){
        triboson2++;
    }
    fftriboson2.close();
    
    ifstream fftriboson3(names_files[9]);
    if(!fftriboson3){
        cout<<"Il file non esiste 10"<<endl;
        return;
    }
    int triboson3 = 0;
    while((fftriboson3>>x_point[triboson3]>>num_triboson3_ossf[triboson3] >>num_triboson3_n_ossf[triboson3]>>den_triboson3_ossf[triboson3] >>den_triboson3_n_ossf[triboson3])){
        triboson3++;
    }
    fftriboson3.close();
    
    
    
    
    ifstream fftx1(names_files[10]);
    if(!fftx1){
        cout<<"Il file non esiste 11"<<endl;
        return;
    }
    int tx1 = 0;
    while((fftx1>>x_point[tx1]>>num_tx1_ossf[tx1] >>num_tx1_n_ossf[tx1]>>den_tx1_ossf[tx1] >>den_tx1_n_ossf[tx1])){
        tx1++;
    }
    fftx1.close();
    
    ifstream fftx2(names_files[11]);
    if(!fftx2){
        cout<<"Il file non esiste 12"<<endl;
        return;
    }
    int tx2 = 0;
    while((fftx2>>x_point[tx2]>>num_tx2_ossf[tx2] >>num_tx2_n_ossf[tx2]>>den_tx2_ossf[tx2] >>den_tx2_n_ossf[tx2])){
        tx2++;
    }
    fftx2.close();
    
    ifstream fftx3(names_files[12]);
    if(!fftx3){
        cout<<"Il file non esiste 13"<<endl;
        return;
    }
    int tx3 = 0;
    while((fftx3>>x_point[tx3]>>num_tx3_ossf[tx3] >>num_tx3_n_ossf[tx3]>>den_tx3_ossf[tx3] >>den_tx3_n_ossf[tx3])){
        tx3++;
    }
    fftx3.close();
    
    ifstream ffdy1(names_files[13]);
    if(!ffdy1){
        cout<<"Il file non esiste 14"<<endl;
        return;
    }
    int dy1 = 0;
    while((ffdy1>>x_point[dy1]>>num_dy1_ossf[dy1] >>num_dy1_n_ossf[dy1]>>den_dy1_ossf[dy1] >>den_dy1_n_ossf[dy1])){
        dy1++;
    }
    ffdy1.close();
    
    ifstream ffdy2(names_files[14]);
    if(!ffdy2){
        cout<<"Il file non esiste 15"<<endl;
        return;
    }
    int dy2 = 0;
    while((ffdy2>>x_point[dy2]>>num_dy2_ossf[dy2] >>num_dy2_n_ossf[dy2]>>den_dy2_ossf[dy2] >>den_dy2_n_ossf[dy2])){
        dy2++;
    }
    ffdy2.close();
    
    ifstream ffdy3(names_files[15]);
    if(!ffdy3){
        cout<<"Il file non esiste 16"<<endl;
        return;
    }
    int dy3 = 0;
    while((ffdy3>>x_point[dy3]>>num_dy3_ossf[dy3] >>num_dy3_n_ossf[dy3]>>den_dy3_ossf[dy3] >>den_dy3_n_ossf[dy3])){
        dy3++;
    }
    ffdy3.close();
    
    ifstream ffdy4(names_files[16]);
    if(!ffdy4){
        cout<<"Il file non esiste 17"<<endl;
        return;
    }
    int dy4 = 0;
    while((ffdy4>>x_point[dy4]>>num_dy4_ossf[dy4] >>num_dy4_n_ossf[dy4]>>den_dy4_ossf[dy4] >>den_dy4_n_ossf[dy4])){
        dy4++;
    }
    ffdy4.close();
    
    
    ifstream ffdy5(names_files[17]);
    if(!ffdy5){
        cout<<"Il file non esiste 18"<<endl;
        return;
    }
    int dy5 = 0;
    while((ffdy5>>x_point[dy5]>>num_dy5_ossf[dy5] >>num_dy5_n_ossf[dy5]>>den_dy5_ossf[dy5] >>den_dy5_n_ossf[dy5])){
        dy5++;
    }
    ffdy5.close();
    
    
    
    
    ifstream ffsig1(names_files[18]);
    if(!ffsig1){
        cout<<"Il file non esiste 19"<<endl;
        return;
    }
    int sig1 = 0;
    while((ffsig1>>x_point[sig1]>>num_signal_1_ossf[sig1] >>num_signal_1_n_ossf[sig1]>>den_signal_1_ossf[sig1] >>den_signal_1_n_ossf[sig1])){
        sig1++;
    }
    ffsig1.close();
    ifstream ffsig2(names_files[19]);
    if(!ffsig2){
        cout<<"Il file non esiste 20"<<endl;
        return;
    }
    int sig2 = 0;
    while((ffsig2>>x_point[sig2]>>num_signal_2_ossf[sig2] >>num_signal_2_n_ossf[sig2]>>den_signal_2_ossf[sig2] >>den_signal_2_n_ossf[sig2])){
        sig2++;
    }
    ffsig2.close();
    
    ifstream ffsig5(names_files[20]);
    if(!ffsig5){
        cout<<"Il file non esiste 21"<<endl;
        return;
    }
    int sig5 = 0;
    while((ffsig5>>x_point[sig5]>>num_signal_5_ossf[sig5] >>num_signal_5_n_ossf[sig5]>>den_signal_5_ossf[sig5] >>den_signal_5_n_ossf[sig5])){
        sig5++;
    }
    ffsig5.close();
    
    ifstream ffsig4(names_files[21]);
    if(!ffsig4){
        cout<<"Il file non esiste 22"<<endl;
        return;
    }
    int sig4 = 0;
    while((ffsig4>>x_point[sig4]>>num_signal_4_ossf[sig4] >>num_signal_4_n_ossf[sig4]>>den_signal_4_ossf[sig4] >>den_signal_4_n_ossf[sig4])){
        sig4++;
    }
    ffsig4.close();
    
    
    
    ifstream ffdiboson1(names_files[22]);
    if(!ffdiboson1){
        cout<<"Il file non esiste 23"<<endl;
        return;
    }
    int diboson1 = 0;
    while((ffdiboson1>>x_point[diboson1]>>num_diboson1_ossf[diboson1] >>num_diboson1_n_ossf[diboson1]>>den_diboson1_ossf[diboson1] >>den_diboson1_n_ossf[diboson1])){
        diboson1++;
    }
    ffdiboson1.close();
    
    ifstream ffdiboson2(names_files[23]);
    if(!ffdiboson2){
        cout<<"Il file non esiste 24"<<endl;
        return;
    }
    int diboson2 = 0;
    while((ffdiboson2>>x_point[diboson2]>>num_diboson2_ossf[diboson2] >>num_diboson2_n_ossf[diboson2]>>den_diboson2_ossf[diboson2] >>den_diboson2_n_ossf[diboson2])){
        diboson2++;
    }
    ffdiboson2.close();
    
    ifstream ffdiboson3(names_files[24]);
    if(!ffdiboson3){
        cout<<"Il file non esiste 25"<<endl;
        return;
    }
    int diboson3 = 0;
    while((ffdiboson3>>x_point[diboson3]>>num_diboson3_ossf[diboson3] >>num_diboson3_n_ossf[diboson3]>>den_diboson3_ossf[diboson3] >>den_diboson3_n_ossf[diboson3])){
        diboson3++;
    }
    ffdiboson3.close();
    
    
    
    
    
    
    // =============== BGK curves =============
    for (int i =0; i< npoint; i++){
        cout<<"x;  "<<x_point[i]<<endl;
        num_triboson_ossf[i] = num_triboson1_ossf[i] + num_triboson2_ossf[i] + num_triboson3_ossf[i] ;
        num_triboson_n_ossf[i] = num_triboson1_n_ossf[i] + num_triboson2_n_ossf[i] + num_triboson3_n_ossf[i] ;
        den_triboson_ossf[i] = den_triboson1_ossf[i] + den_triboson2_ossf[i] + den_triboson3_ossf[i] ;
        den_triboson_n_ossf[i] = den_triboson1_n_ossf[i] + den_triboson2_n_ossf[i] + den_triboson3_n_ossf[i] ;
        
        num_xg_ossf[i] = num_xg1_ossf[i] ;
        num_xg_n_ossf[i] = num_xg1_n_ossf[i] ;
        den_xg_ossf[i] = den_xg1_ossf[i] ;
        den_xg_n_ossf[i] = den_xg1_n_ossf[i]  ;
        
        num_tx_ossf[i] = num_tx1_ossf[i] + num_tx2_ossf[i] + num_tx3_ossf[i] ;
        num_tx_n_ossf[i] = num_tx1_n_ossf[i] + num_tx2_n_ossf[i] + num_tx3_n_ossf[i] ;
        den_tx_ossf[i] = den_tx1_ossf[i] + den_tx2_ossf[i] + den_tx3_ossf[i] ;
        den_tx_n_ossf[i] = den_tx1_n_ossf[i] + den_tx2_n_ossf[i] + den_tx3_n_ossf[i] ;
        
        
        num_dy_ossf[i] =  num_dy1_ossf[i] + num_dy2_ossf[i] + num_dy3_ossf[i] + num_dy4_ossf[i] + num_dy5_ossf[i];
        num_dy_n_ossf[i] =  num_dy1_n_ossf[i] + num_dy2_n_ossf[i] + num_dy3_n_ossf[i] + num_dy4_n_ossf[i] + num_dy5_n_ossf[i] ;
        den_dy_ossf[i] =  den_dy1_ossf[i] + den_dy2_ossf[i] + den_dy3_ossf[i] + den_dy4_ossf[i] + den_dy5_ossf[i] ;
        den_dy_n_ossf[i] =  den_dy1_n_ossf[i] + den_dy2_n_ossf[i] + den_dy3_n_ossf[i] + den_dy4_n_ossf[i] + den_dy5_n_ossf[i];
        
        
        
        num_diboson_ossf[i] =  num_diboson1_ossf[i] + num_diboson2_ossf[i] + num_diboson3_ossf[i] ;
        num_diboson_n_ossf[i] = num_diboson1_n_ossf[i] + num_diboson2_n_ossf[i] + num_diboson3_n_ossf[i] ;
        den_diboson_ossf[i] =  den_diboson1_ossf[i] + den_diboson2_ossf[i] + den_diboson3_ossf[i] ;
        den_diboson_n_ossf[i] =  den_diboson1_n_ossf[i] + den_diboson2_n_ossf[i] + den_diboson3_n_ossf[i] ;
        
        
        
        
        
        back_sum_ossf[i]=  num_triboson_ossf[i] + num_xg_ossf[i]+ num_tx_ossf[i]  + num_zz_ossf[i]+ num_dy_ossf[i] + num_diboson_ossf[i] ;
        back_sum_n_ossf[i]=  num_triboson_n_ossf[i] + num_xg_n_ossf[i]+ num_tx_n_ossf[i]  + num_zz_n_ossf[i]+ num_dy_n_ossf[i] + num_diboson_n_ossf[i] ;
        
        back_tot_ossf[i]=  num_triboson_ossf[i] + num_xg_ossf[i]+ num_tx_ossf[i]  + num_zz_ossf[i]+ num_dy_ossf[i] + num_diboson_ossf[i] ;
        back_tot_n_ossf[i]=  num_triboson_n_ossf[i] + num_xg_n_ossf[i]+ num_tx_n_ossf[i]  + num_zz_n_ossf[i]+ num_dy_n_ossf[i] + num_diboson_n_ossf[i] ;
        
        back_den_tot_ossf[i]=  den_triboson_ossf[i] + den_xg_ossf[i]+ den_tx_ossf[i]  + den_zz_ossf[i]+ den_dy_ossf[i] + den_diboson_ossf[i] ;
        back_den_tot_n_ossf[i]=  den_triboson_n_ossf[i] + den_xg_n_ossf[i]+ den_tx_n_ossf[i]  + den_zz_n_ossf[i]+ den_dy_n_ossf[i] + den_diboson_n_ossf[i] ;
    }
    
    
    
    
    
    
    for (int i =0; i< npoint; i++){
        double j = i;
        x_point[i] = j/200;
        
        
        signal_bgk_1_ossf[i] = num_signal_1_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_1_ossf[i]);
        signal_bgk_2_ossf[i] = num_signal_2_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_2_ossf[i]);
        signal_bgk_5_ossf[i] = num_signal_5_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_5_ossf[i]);
        signal_bgk_4_ossf[i] = num_signal_4_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_4_ossf[i]);
        signal_bgk_6_ossf[i] = num_signal_6_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_6_ossf[i]);
        signal_bgk_7_ossf[i] = num_signal_7_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_7_ossf[i]);
        signal_bgk_8_ossf[i] = num_signal_8_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_8_ossf[i]);
        signal_bgk_9_ossf[i] = num_signal_9_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_9_ossf[i]);
        signal_bgk_3_ossf[i] = num_signal_3_ossf[i]/TMath::Sqrt (back_sum_ossf[i] + num_signal_3_ossf[i]);
        
        signal_bgk_1_n_ossf[i] = num_signal_1_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_1_n_ossf[i]);
        signal_bgk_2_n_ossf[i] = num_signal_2_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_2_n_ossf[i]);
        signal_bgk_5_n_ossf[i] = num_signal_5_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_5_n_ossf[i]);
        signal_bgk_4_n_ossf[i] = num_signal_4_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_4_n_ossf[i]);
        signal_bgk_6_n_ossf[i] = num_signal_6_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_6_n_ossf[i]);
        signal_bgk_7_n_ossf[i] = num_signal_7_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_7_n_ossf[i]);
        signal_bgk_8_n_ossf[i] = num_signal_8_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_8_n_ossf[i]);
        signal_bgk_9_n_ossf[i] = num_signal_9_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_9_n_ossf[i]);
        signal_bgk_3_n_ossf[i] = num_signal_3_n_ossf[i]/TMath::Sqrt (back_sum_n_ossf[i] + num_signal_3_n_ossf[i]);
        
        
        if (den_dy_ossf[i] == 0 ) dy_ossf[i]=0;
        else {
            dy_ossf[i] = num_dy_ossf[i]/den_dy_ossf[i];
        }
        if (den_dy_n_ossf[i] == 0 ) dy_n_ossf[i]=0;
        else {
            dy_n_ossf[i] = num_dy_n_ossf[i]/den_dy_n_ossf[i];
        }
        
        
        if (den_diboson_ossf[i] == 0 ) diboson_ossf[i]=0;
        else {
            diboson_ossf[i] = num_diboson_ossf[i]/den_diboson_ossf[i];
        }
        
        
        
        if (den_zz_ossf[i] == 0 ) zz_ossf[i]=0;
        else {
            zz_ossf[i] = num_zz_ossf[i]/den_zz_ossf[i];
        }
        
        if (den_triboson_ossf[i] == 0 ) triboson_ossf[i]=0;
        else {
            triboson_ossf[i] = num_triboson_ossf[i]/den_triboson_ossf[i];
        }
        
        if (den_xg_ossf[i] == 0 ) xg_ossf[i]=0;
        else {
            xg_ossf[i] = num_xg_ossf[i]/den_xg_ossf[i];
        }
        if (den_tx_ossf[i] == 0 ) tx_ossf[i]=0;
        else {
            tx_ossf[i] = num_tx_ossf[i]/den_tx_ossf[i];
        }
        
        if (den_signal_1_ossf[i] == 0 ) signal_1_ossf[i]=0;
        else {
            signal_1_ossf[i] = num_signal_1_ossf[i]/den_signal_1_ossf[i];
        }
        if (den_signal_2_ossf[i] == 0 ) signal_2_ossf[i]=0;
        else {
            signal_2_ossf[i] = num_signal_2_ossf[i]/den_signal_2_ossf[i];
        }
        if (den_signal_5_ossf[i] == 0 ) signal_5_ossf[i]=0;
        else {
            signal_5_ossf[i] = num_signal_5_ossf[i]/den_signal_5_ossf[i];
        }
        if (den_signal_4_ossf[i] == 0 ) signal_4_ossf[i]=0;
        else {
            signal_4_ossf[i] = num_signal_4_ossf[i]/den_signal_4_ossf[i];
        }
        
        
        if (den_signal_6_ossf[i] == 0 ) signal_6_ossf[i]=0;
        else {
            signal_6_ossf[i] = num_signal_6_ossf[i]/den_signal_6_ossf[i];
        }
        if (den_signal_7_ossf[i] == 0 ) signal_7_ossf[i]=0;
        else {
            signal_7_ossf[i] = num_signal_7_ossf[i]/den_signal_7_ossf[i];
        }
        if (den_signal_8_ossf[i] == 0 ) signal_8_ossf[i]=0;
        else {
            signal_8_ossf[i] = num_signal_8_ossf[i]/den_signal_8_ossf[i];
        }
        if (den_signal_9_ossf[i] == 0 ) signal_9_ossf[i]=0;
        else {
            signal_9_ossf[i] = num_signal_9_ossf[i]/den_signal_9_ossf[i];
        }
        if (den_signal_3_ossf[i] == 0 ) signal_3_ossf[i]=0;
        else {
            signal_3_ossf[i] = num_signal_3_ossf[i]/den_signal_3_ossf[i];
        }
        if (den_diboson_n_ossf[i] == 0 ) diboson_n_ossf[i]=0;
        else {
            diboson_n_ossf[i] = num_diboson_n_ossf[i]/den_diboson_n_ossf[i];
        }
        if (den_zz_n_ossf[i] == 0 ) zz_n_ossf[i]=0;
        else {
            zz_n_ossf[i] = num_zz_n_ossf[i]/den_zz_n_ossf[i];
        }
        if (den_triboson_n_ossf[i] == 0 ) triboson_n_ossf[i]=0;
        else {
            triboson_n_ossf[i] = num_triboson_n_ossf[i]/den_triboson_n_ossf[i];
        }
        if (den_xg_n_ossf[i] == 0 ) xg_n_ossf[i]=0;
        else {
            xg_n_ossf[i] = num_xg_n_ossf[i]/den_xg_n_ossf[i];
        }
        if (den_tx_n_ossf[i] == 0 ) tx_n_ossf[i]=0;
        else {
            tx_n_ossf[i] = num_tx_n_ossf[i]/den_tx_n_ossf[i];
        }
        
        
        if (den_signal_1_n_ossf[i] == 0 ) signal_1_n_ossf[i]=0;
        else {
            signal_1_n_ossf[i] = num_signal_1_n_ossf[i]/den_signal_1_n_ossf[i];
        }
        if (den_signal_2_n_ossf[i] == 0 ) signal_2_n_ossf[i]=0;
        else {
            signal_2_n_ossf[i] = num_signal_2_n_ossf[i]/den_signal_2_n_ossf[i];
        }
        if (den_signal_5_n_ossf[i] == 0 ) signal_5_n_ossf[i]=0;
        else {
            signal_5_n_ossf[i] = num_signal_5_n_ossf[i]/den_signal_5_n_ossf[i];
        }
        if (den_signal_4_n_ossf[i] == 0 ) signal_4_n_ossf[i]=0;
        else {
            signal_4_n_ossf[i] = num_signal_4_n_ossf[i]/den_signal_4_n_ossf[i];
        }
        if (den_signal_6_n_ossf[i] == 0 ) signal_6_n_ossf[i]=0;
        else {
            signal_6_n_ossf[i] = num_signal_6_n_ossf[i]/den_signal_6_n_ossf[i];
        }
        if (den_signal_7_n_ossf[i] == 0 ) signal_7_n_ossf[i]=0;
        else {
            signal_7_n_ossf[i] = num_signal_7_n_ossf[i]/den_signal_7_n_ossf[i];
        }
        if (den_signal_8_n_ossf[i] == 0 ) signal_8_n_ossf[i]=0;
        else {
            signal_8_n_ossf[i] = num_signal_8_n_ossf[i]/den_signal_8_n_ossf[i];
        }
        if (den_signal_9_n_ossf[i] == 0 ) signal_9_n_ossf[i]=0;
        else {
            signal_9_n_ossf[i] = num_signal_9_n_ossf[i]/den_signal_9_n_ossf[i];
        }
        if (den_signal_3_n_ossf[i] == 0 ) signal_3_n_ossf[i]=0;
        else {
            signal_3_n_ossf[i] = num_signal_3_n_ossf[i]/den_signal_3_n_ossf[i];
        }
        
        
        if (back_den_tot_ossf[i] == 0 )  back_value_tot_ossf[i]=0;
        else {
            back_value_tot_ossf[i] = back_tot_ossf[i]/back_den_tot_ossf[i];
        }
        if (back_den_tot_n_ossf[i] == 0 )  back_value_tot_n_ossf[i]=0;
        else {
            back_value_tot_n_ossf[i] = back_tot_n_ossf[i]/back_den_tot_n_ossf[i];
        }
        
        
    }// end for
    
    
    g_pt_s_b_1_ossf = new TGraphErrors(npoint,x_point , signal_bgk_1_ossf,  0, 0);
    g_pt_s_b_2_ossf = new TGraphErrors(npoint,x_point , signal_bgk_2_ossf,  0, 0);
    g_pt_s_b_5_ossf = new TGraphErrors(npoint,x_point , signal_bgk_5_ossf,  0, 0);
    g_pt_s_b_4_ossf = new TGraphErrors(npoint,x_point , signal_bgk_4_ossf,  0, 0);
    g_pt_s_b_6_ossf = new TGraphErrors(npoint,x_point , signal_bgk_6_ossf,  0, 0);
    g_pt_s_b_7_ossf = new TGraphErrors(npoint,x_point , signal_bgk_7_ossf,  0, 0);
    g_pt_s_b_8_ossf = new TGraphErrors(npoint,x_point , signal_bgk_8_ossf,  0, 0);
    g_pt_s_b_9_ossf = new TGraphErrors(npoint,x_point , signal_bgk_9_ossf,  0, 0);
    g_pt_s_b_3_ossf = new TGraphErrors(npoint,x_point , signal_bgk_3_ossf,  0, 0);
    
    g_pt_s_b_1_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_1_n_ossf,  0, 0);
    g_pt_s_b_2_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_2_n_ossf,  0, 0);
    g_pt_s_b_5_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_5_n_ossf,  0, 0);
    g_pt_s_b_4_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_4_n_ossf,  0, 0);
    g_pt_s_b_6_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_6_n_ossf,  0, 0);
    g_pt_s_b_7_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_7_n_ossf,  0, 0);
    g_pt_s_b_8_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_8_n_ossf,  0, 0);
    g_pt_s_b_9_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_9_n_ossf,  0, 0);
    g_pt_s_b_3_n_ossf = new TGraphErrors(npoint,x_point , signal_bgk_3_n_ossf,  0, 0);
    
    
    cosmesi(11,*&g_pt_s_b_1_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(22,*&g_pt_s_b_2_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(33,*&g_pt_s_b_5_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(44,*&g_pt_s_b_4_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(1,*&g_pt_s_b_6_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(2,*&g_pt_s_b_7_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(3,*&g_pt_s_b_8_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(4,*&g_pt_s_b_9_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(5,*&g_pt_s_b_3_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(11,*&g_pt_s_b_1_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(22,*&g_pt_s_b_2_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(33,*&g_pt_s_b_5_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(44,*&g_pt_s_b_4_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(1,*&g_pt_s_b_6_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(2,*&g_pt_s_b_7_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(3,*&g_pt_s_b_8_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(4,*&g_pt_s_b_9_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    cosmesi(5,*&g_pt_s_b_3_n_ossf , "d_{xy} (cm)","S/#sqrt{B}" );
    
    g_pt_tot_ossf = new TGraphErrors(npoint, x_point, back_value_tot_ossf,  0, 0);
    g_pt_dy_ossf = new TGraphErrors(npoint, x_point, dy_ossf,  0, 0);
    g_pt_diboson_ossf = new TGraphErrors(npoint, x_point, diboson_ossf,  0, 0);
    g_pt_zz_ossf = new TGraphErrors(npoint, x_point, zz_ossf,  0, 0);
    g_pt_triboson_ossf = new TGraphErrors(npoint, x_point, triboson_ossf,  0, 0);
    g_pt_tx_ossf = new TGraphErrors(npoint, x_point, tx_ossf,  0, 0);
    g_pt_xg_ossf = new TGraphErrors(npoint, x_point, xg_ossf,  0, 0);
    
    g_pt_signal_1_ossf = new TGraphErrors(npoint, x_point, signal_1_ossf,  0, 0);
    g_pt_signal_2_ossf = new TGraphErrors(npoint, x_point, signal_2_ossf,  0, 0);
    g_pt_signal_5_ossf = new TGraphErrors(npoint, x_point, signal_5_ossf,  0, 0);
    g_pt_signal_4_ossf = new TGraphErrors(npoint, x_point, signal_4_ossf,  0, 0);
    
    g_pt_signal_6_ossf = new TGraphErrors(npoint, x_point, signal_6_ossf,  0, 0);
    g_pt_signal_7_ossf = new TGraphErrors(npoint, x_point, signal_7_ossf,  0, 0);
    g_pt_signal_8_ossf = new TGraphErrors(npoint, x_point, signal_8_ossf,  0, 0);
    g_pt_signal_9_ossf = new TGraphErrors(npoint, x_point, signal_9_ossf,  0, 0);
    g_pt_signal_3_ossf = new TGraphErrors(npoint, x_point, signal_3_ossf,  0, 0);
    
    g_pt_tot_n_ossf = new TGraphErrors(npoint, x_point, back_value_tot_n_ossf,  0, 0);
    g_pt_dy_n_ossf = new TGraphErrors(npoint, x_point, dy_n_ossf,  0, 0);
    g_pt_diboson_n_ossf = new TGraphErrors(npoint, x_point, diboson_n_ossf,  0, 0);
    g_pt_zz_n_ossf = new TGraphErrors(npoint, x_point, zz_n_ossf,  0, 0);
    g_pt_triboson_n_ossf = new TGraphErrors(npoint, x_point, triboson_n_ossf,  0, 0);
    g_pt_tx_n_ossf = new TGraphErrors(npoint, x_point, tx_n_ossf,  0, 0);
    g_pt_xg_n_ossf = new TGraphErrors(npoint, x_point, xg_n_ossf,  0, 0);
    
    g_pt_signal_1_n_ossf = new TGraphErrors(npoint, x_point, signal_1_n_ossf,  0, 0);
    g_pt_signal_2_n_ossf = new TGraphErrors(npoint, x_point, signal_2_n_ossf,  0, 0);
    g_pt_signal_5_n_ossf = new TGraphErrors(npoint, x_point, signal_5_n_ossf,  0, 0);
    g_pt_signal_4_n_ossf = new TGraphErrors(npoint, x_point, signal_4_n_ossf,  0, 0);
    g_pt_signal_6_n_ossf = new TGraphErrors(npoint, x_point, signal_6_n_ossf,  0, 0);
    g_pt_signal_7_n_ossf = new TGraphErrors(npoint, x_point, signal_7_n_ossf,  0, 0);
    g_pt_signal_8_n_ossf = new TGraphErrors(npoint, x_point, signal_8_n_ossf,  0, 0);
    g_pt_signal_9_n_ossf = new TGraphErrors(npoint, x_point, signal_9_n_ossf,  0, 0);
    g_pt_signal_3_n_ossf = new TGraphErrors(npoint, x_point, signal_3_n_ossf,  0, 0);
    
    
    
    cosmesi(0,*&g_pt_tot_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(1,*&g_pt_diboson_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(2,*&g_pt_zz_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(3,*&g_pt_triboson_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(4,*&g_pt_tx_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(5,*&g_pt_xg_ossf, "d_{xy} (cm)","Efficiency" );
    
    
    cosmesi(66,*&g_pt_signal_1_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(77,*&g_pt_signal_2_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(88,*&g_pt_signal_5_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(99,*&g_pt_signal_4_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(6,*&g_pt_signal_6_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(7,*&g_pt_signal_7_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(8,*&g_pt_signal_8_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(9,*&g_pt_signal_9_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(10,*&g_pt_signal_3_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(11,*&g_pt_dy_ossf, "d_{xy} (cm)","Efficiency" );
    
    cosmesi(0,*&g_pt_tot_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(1,*&g_pt_diboson_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(2,*&g_pt_zz_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(3,*&g_pt_triboson_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(4,*&g_pt_tx_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(5,*&g_pt_xg_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(66,*&g_pt_signal_1_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(77,*&g_pt_signal_2_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(88,*&g_pt_signal_5_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(99,*&g_pt_signal_4_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(6,*&g_pt_signal_6_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(7,*&g_pt_signal_7_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(8,*&g_pt_signal_8_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(9,*&g_pt_signal_9_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(10,*&g_pt_signal_3_n_ossf, "d_{xy} (cm)","Efficiency" );
    cosmesi(11,*&g_pt_dy_n_ossf, "d_{xy} (cm)","Efficiency" );
    
    
    
    
    TLegend *leg_n_ossf_pt_eff, *leg_ossf_pt_eff,*leg_n_ossf_pt_eff_2d, *leg_ossf_pt_eff_2d;
    
    TCanvas *c_pt_eff= new TCanvas("c_pt_eff","c_pt_eff");
    c_pt_eff->Divide(2,1);
    c_pt_eff->cd(1);
    gPad->SetBottomMargin(0.2);
    g_pt_signal_6_ossf->Draw("ap");
    g_pt_signal_1_ossf->Draw("psame");
    g_pt_signal_2_ossf->Draw("psame");
    g_pt_signal_5_ossf->Draw("psame");
    g_pt_signal_4_ossf->Draw("psame");
    g_pt_signal_7_ossf->Draw("psame");
    g_pt_signal_8_ossf->Draw("psame");
    g_pt_signal_9_ossf->Draw("psame");
    g_pt_signal_3_ossf->Draw("psame");
    g_pt_diboson_ossf->Draw("psame");
    g_pt_zz_ossf->Draw("psame");
    g_pt_triboson_ossf->Draw("psame");
    g_pt_tx_ossf->Draw("psame");
    g_pt_xg_ossf->Draw("psame");
    g_pt_dy_ossf->Draw("psame");
    g_pt_tot_ossf->Draw("psame");
    leg_ossf_pt_eff  = new TLegend(0.1,0.77,0.56,0.98);
    //leg_ossf_pt_eff->SetTextSize(0.3);
    leg_ossf_pt_eff ->SetHeader("OSSF event "); // option "C" allows to center the header
    leg_ossf_pt_eff ->AddEntry(g_pt_tot_ossf,"Total bgk","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_dy_ossf,"DY and TTbar","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_diboson_ossf,"diboson","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_zz_ossf,"ZZ","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_triboson_ossf,"triboson","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_tx_ossf,"TT/T + X","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_xg_ossf,"X + #gamma","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_1_ossf,"m_{N} = 1 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_2_ossf,"m_{N} = 2 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_3_ossf,"m_{N} = 3 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_4_ossf,"m_{N} = 4 GeV","p");

    leg_ossf_pt_eff ->AddEntry(g_pt_signal_5_ossf,"m_{N} = 5 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_6_ossf,"m_{N} = 6 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_7_ossf,"m_{N} = 7 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_8_ossf,"m_{N} = 8 GeV","p");
    leg_ossf_pt_eff ->AddEntry(g_pt_signal_9_ossf,"m_{N} = 9 GeV","p");
    
    
    leg_ossf_pt_eff ->Draw();
    TF1 *f1=new TF1("f1","x",0,1);
    TGaxis *A1 = new TGaxis(0,-0.12,0.5,-0.12,"f1",506);
    A1->SetTitleOffset(0.80);
    A1->SetTitleSize(0.03);
    A1 ->SetTitleFont(132);
    A1->SetLabelSize(0.03);
    A1->SetLabelOffset(0.00);
    A1->SetLabelFont(132);
    A1->SetTitle("d_{z} (cm)");
    A1->Draw();
    
    TF1 *f2=new TF1("f2","x",0,40);
    TGaxis *A2 = new TGaxis(0,-0.25,0.5,-0.25,"f2",506);
    A2->SetTitle("3DIPsig");
    A2->SetTitleOffset(0.80);
    A2->SetTitleSize(0.03);
    A2 ->SetTitleFont(132);
    A2->SetLabelSize(0.03);
    A2->SetLabelOffset(0.00);
    A2->SetLabelFont(132);
    A2->Draw();
    
    c_pt_eff->cd(2);
    gPad->SetBottomMargin(0.2);
    
    
    g_pt_signal_6_n_ossf->Draw("ap");
    g_pt_signal_1_n_ossf->Draw("psame");
    g_pt_signal_2_n_ossf->Draw("psame");
    g_pt_signal_5_n_ossf->Draw("psame");
    g_pt_signal_4_n_ossf->Draw("psame");
    g_pt_signal_7_n_ossf->Draw("psame");
    g_pt_signal_8_n_ossf->Draw("psame");
    g_pt_signal_9_n_ossf->Draw("psame");
    g_pt_signal_3_n_ossf->Draw("psame");
    g_pt_diboson_n_ossf->Draw("psame");
    g_pt_zz_n_ossf->Draw("psame");
    g_pt_triboson_n_ossf->Draw("psame");
    g_pt_tx_n_ossf->Draw("psame");
    //g_pt_xg_n_ossf->Draw("psame");
    g_pt_dy_n_ossf->Draw("psame");
    g_pt_tot_n_ossf->Draw("psame");
    leg_n_ossf_pt_eff  = new TLegend(0.1,0.77,0.56,0.98);
    //leg_n_ossf_pt_eff->SetTextSize(0.3);
    leg_n_ossf_pt_eff ->SetHeader("NO_OSSF event"); // option "C" allows to center the header
    leg_n_ossf_pt_eff ->AddEntry(g_pt_tot_n_ossf,"Total bgk","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_dy_n_ossf,"DY and TTbar","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_diboson_n_ossf,"diboson","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_zz_n_ossf,"ZZ","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_triboson_n_ossf,"triboson","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_tx_n_ossf,"TT/T + X","p");
    //leg_n_ossf_pt_eff ->AddEntry(g_pt_xg_n_ossf,"X + #gamma","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_1_n_ossf,"m_{N} = 1 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_2_n_ossf,"m_{N} = 2 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_3_n_ossf,"m_{N} = 3 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_4_n_ossf,"m_{N} = 4 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_5_n_ossf,"m_{N} = 5 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_6_n_ossf,"m_{N} = 6 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_7_n_ossf,"m_{N} = 7 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_8_n_ossf,"m_{N} = 8 GeV","p");
    leg_n_ossf_pt_eff ->AddEntry(g_pt_signal_9_n_ossf,"m_{N} = 9 GeV","p");
    
    
    TF1 *f3=new TF1("f3","x",0,1);
    TGaxis *A3 = new TGaxis(0,-0.12,0.5,-0.12,"f3",506);
    A3->SetTitle("d_{z} (cm)");
    A3->SetTitleOffset(0.80);
    A3->SetTitleSize(0.03);
    A3 ->SetTitleFont(132);
    A3->SetLabelSize(0.03);
    A3->SetLabelOffset(0.00);
    A3->SetLabelFont(132);
    A3->Draw();
    
    TF1 *f4=new TF1("f4","x",0,40);
    TGaxis *A4 = new TGaxis(0,-0.25,0.5,-0.25,"f4",506);
    A4->SetTitle("3DIPsig");
    A4->SetTitleOffset(0.80);
    A4->SetTitleSize(0.03);
    A4 ->SetTitleFont(132);
    A4->SetLabelSize(0.03);
    A4->SetLabelOffset(0.00);
    A4->SetLabelFont(132);
    A4->Draw();
    
    
    leg_n_ossf_pt_eff ->Draw();
    c_pt_eff->Print("PLOT_PT/c_pt_eff_max.root");
    delete  c_pt_eff ;
    
    
    TCanvas *c_pt_eff_2d= new TCanvas("c_pt_eff_2d","c_pt_eff_2d");
    c_pt_eff_2d->Divide(2,1);
    c_pt_eff_2d->cd(1);
    gPad->SetBottomMargin(0.2);
    
    g_pt_s_b_6_ossf->Draw("ap");
    g_pt_s_b_1_ossf->Draw("psame");
    g_pt_s_b_2_ossf->Draw("psame");
    g_pt_s_b_5_ossf->Draw("psame");
    g_pt_s_b_4_ossf->Draw("psame");
    g_pt_s_b_7_ossf->Draw("psame");
    g_pt_s_b_8_ossf->Draw("psame");
    g_pt_s_b_9_ossf->Draw("psame");
    g_pt_s_b_3_ossf->Draw("psame");
    leg_ossf_pt_eff_2d  = new TLegend(0.1,0.77,0.56,0.98);
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_1_ossf,"m_{N} = 1 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_2_ossf,"m_{N} = 2 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_3_ossf,"m_{N} = 3 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_4_ossf,"m_{N} = 4 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_5_ossf,"m_{N} = 5 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_6_ossf,"m_{N} = 6 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_7_ossf,"m_{N} = 7 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_8_ossf,"m_{N} = 8 GeV","p");
    leg_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_9_ossf,"m_{N} = 9 GeV","p");
    leg_ossf_pt_eff_2d ->Draw();
    
    TF1 *f5=new TF1("f5","x",0,1);
    TGaxis *A5 = new TGaxis(0,-0.12,0.5,-0.12,"f5",506);
    A5->SetTitle("d_{z} (cm)");
    A5->SetTitleOffset(0.80);
    A5->SetTitleSize(0.03);
    A5 ->SetTitleFont(132);
    A5->SetLabelSize(0.03);
    A5->SetLabelOffset(0.00);
    A5->SetLabelFont(132);
    A5->Draw();
    
    TF1 *f6=new TF1("f6","x",0,40);
    TGaxis *A6 = new TGaxis(0,-0.25,0.5,-0.25,"f6",506);
    A6->SetTitle("3DIPsig");
    A6->SetTitleOffset(0.80);
    A6->SetTitleSize(0.03);
    A6 ->SetTitleFont(132);
    A6->SetLabelSize(0.03);
    A6->SetLabelOffset(0.00);
    A6->SetLabelFont(132);
    A6->Draw();
    
    
    
    c_pt_eff_2d->cd(2);
    gPad->SetBottomMargin(0.2);
    
    g_pt_s_b_6_n_ossf->Draw("ap");
    g_pt_s_b_1_n_ossf->Draw("psame");
    g_pt_s_b_2_n_ossf->Draw("psame");
    g_pt_s_b_5_n_ossf->Draw("psame");
    g_pt_s_b_4_n_ossf->Draw("psame");
    g_pt_s_b_7_n_ossf->Draw("psame");
    g_pt_s_b_8_n_ossf->Draw("psame");
    g_pt_s_b_9_n_ossf->Draw("psame");
    g_pt_s_b_3_n_ossf->Draw("psame");
    leg_n_ossf_pt_eff_2d  = new TLegend(0.1,0.77,0.56,0.98);
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_1_n_ossf,"m_{N} = 1 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_2_n_ossf,"m_{N} = 2 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_3_n_ossf,"m_{N} = 3 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_4_n_ossf,"m_{N} = 4 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_5_n_ossf,"m_{N} = 5 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_6_n_ossf,"m_{N} = 6 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_7_n_ossf,"m_{N} = 7 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_8_n_ossf,"m_{N} = 8 GeV","p");
    leg_n_ossf_pt_eff_2d ->AddEntry(g_pt_s_b_9_n_ossf,"m_{N} = 9 GeV","p");
    leg_n_ossf_pt_eff_2d ->Draw();
    
    TF1 *f7=new TF1("f7","x",0,1);
    TGaxis *A7 = new TGaxis(0,-0.12,0.5,-0.12,"f7",506);
    A7->SetTitle("d_{z} (cm)");
    A7->SetTitleOffset(0.80);
    A7->SetTitleSize(0.03);
    A7 ->SetTitleFont(132);
    A7->SetLabelSize(0.03);
    A7->SetLabelOffset(0.00);
    A7->SetLabelFont(132);
    A7->Draw();
    
    TF1 *f8=new TF1("f8","x",0,40);
    TGaxis *A8 = new TGaxis(0,-0.25,0.5,-0.25,"f8",506);
    A8->SetTitle("3DIPsig");
    A8->SetTitleOffset(0.80);
    A8->SetTitleSize(0.03);
    A8 ->SetTitleFont(132);
    A8->SetLabelSize(0.03);
    A8->SetLabelOffset(0.00);
    A8->SetLabelFont(132);
    A8->Draw();
    
    c_pt_eff_2d->cd(2);
    c_pt_eff_2d->Print("PLOT_PT/c_pt_eff_2d_max.root");
    delete  c_pt_eff_2d ;
    
    
}// end analisi
