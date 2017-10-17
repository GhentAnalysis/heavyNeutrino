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
#include <Analysis_mc.h>
#include <THStack.h>
#include <TPaveText.h>
#include <THStack.h>

void majo (string inputRootFile= "file_mva_gen1.root");


// ********************************************************************
void majo(string inputRootFile){
    Double_t pigreco= TMath::ACos(-1);
    
    cout<<"in analisi"<<endl;
    cout<<"---------------------------"<<endl;
    

//==========================================================================================
    Analysis_mc all("/Users/Martina/Desktop/CMS/file_bck/zg.root");
    all.analisi(1);


}
