#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include "TAxis.h"
#include "TColor.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGeoManager.h"
#include "TGeoGlobalMagField.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraphPolar.h"
#include "TGrid.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLinearFitter.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TRegexp.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TUUID.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TVirtualPad.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TFitter.h"
#include "TMinuit.h"

using namespace std;

class PlotUtils
{
  public:
  double PrintStat(const TString tag, TH1 *hh, const double val0, const double val1, const double oldsel);
  void FillHist(TH1 * hh,  double xx, const double yy);
  void ProcessHist(TList *lout, const bool kMC);
  void DrawHist(TList *lout, const double plotscale, TList * overlayList, const TString outdir);
  THStack * ConvertToStack(const TH2D * hh, const bool kMC);
  TH1D * GetStackedSum(THStack *stk);
  void ScaleStack(THStack *stk, const double scale);
  TH2D * NormalHist(const TH2D *hraw, const Double_t thres, const Bool_t kmax);
  THStack * NormalizeStack(THStack * hstk);
  void getProfileFit(TH2D * h2d);
  TH1D * GetCDF(const TH2D *hraw, const TString hname);
  void DrawOverlay(TH1D *holay);
  void IniColorCB();
  void SetColor();
  int GetColor(const int col);
  int * GetColorArray(const int minsize = -999);
  void PadSetup(TPad *currentPad, const Double_t currentLeft=0.12, const Double_t currentTop=0.09, const Double_t currentRight=0.13, const Double_t currentBottom=0.14);
  void gStyleSetup();
  TLegend * DrawLegend(const vector<TString> &entries, const vector<TString>& htype, const TString tag, const int *tmpcol=0x0, const int * tmpmkr=0x0, const int ncol = 1);
  void getSliceXDrawY(TH2D * h2d);
  void xSlicedEnergyCorrection(TH2D * h2d);

  static Double_t CauchyDens(Double_t *x, Double_t *par)
  {
    Double_t pi   = TMath::Pi();
    Double_t mean = par[0];
    Double_t fwhm = par[1];
    Double_t height = par[2];

    Double_t arg = x[0]-mean;
    Double_t top = fwhm;
    Double_t bot = pi*(arg*arg+top*top);

    Double_t func = height*(top/bot);
    return func;
  }

  static Double_t ShowerEenergyFCN(Double_t *x, Double_t *par)
  {
    return par[3] + (par[0] - par[3])/(1 + pow((x[0]/par[2]),par[1]));
  }


  private:
};
#endif
