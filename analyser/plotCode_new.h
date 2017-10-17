#ifndef PlotScript
#define PlotScript
 
//import ROOT classes
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH2D.h"
 
const double xPad = 0.05;
const Color_t colors[8] = {kCyan +2, kYellow +1,  kBlue, kGreen +2, kMagenta, kOrange +1, kMagenta -1, kCyan};
const Color_t sigCols[11] = {kRed, kBlue, kGreen, kYellow, kCyan, kMagenta, kMagenta, kRed-3, kBlue-3, kRed, kYellow + 4 };
//Set histogram colors and lines
void histcol(TH1D *, const Color_t);
//Return histogram divided by other histogram (both are normalized)
TH1D *HistDiv(TH1D *, TH1D *, const bool abs = false);
//Set Histogram labelsizes
void HistLabelSizes(TH1D *h, const double xlabel = 0.045, const double xtitle = 0.05, const double ylabel = 0.045, const double ytitle = 0.045);
void HistLabelSizes(TH2D *h, const double xlabel = 0.045, const double xtitle = 0.05, const double ylabel = 0.045, const double ytitle = 0.045);
//Set Stack colors
void StackCol(TH1D *h, const Color_t);
//Order histograms in terms of yields (biggest first)
void yieldOrder(TH1D**&, unsigned*, const unsigned);
void plotDataVSMC(TH1D* data, TH1D** bkg, const TString* names, const unsigned nHist, const TString& file, const bool ylog = false, const unsigned widthopt = 0, const bool plotsig = false, TH1D** signal = nullptr, const TString* signames =nullptr, const unsigned nSig = 0, const bool signorm=false);
#endif
