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

#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGaxis.h"

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

#include "TLegend.h"
#include "TLegendEntry.h"

#include <iostream>
#include <random>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TH2D.h"


using namespace std;
//double plotScale = 0.561778;

void DrawHistOutput(TH1D *h_truth, TH1D *h_reco, TH1D *h_data, TH1D *h_unfold, TH1D *h_unfold_data, const bool &rebin);
TH1D * GetIncidentHist(TH1D * InitialHist, TH1D * InteractingHist);
void DrawTotalXSPlots(TList *lout, TH1D *hbeamint_truth, TH1D *hinc_truth, TH1D *hint_truth, TH1D *hbeamint_reco, TH1D *hinc_reco, TH1D *hint_reco, TH1D *hbeamint_data, TH1D *hinc_data, TH1D *hint_data, TGraph* g_cex, double & totalCEXXS_800_truth, double & totalCEXXS_800_reco, double & totalCEXXS_800_data, double & totalCEXXS_800_truth_error, double & totalCEXXS_800_reco_error, double & totalCEXXS_800_data_error, bool doTruth = true, bool doReco = true, bool doData = true);
void DrawDiffXSPlots(TList *lout, TH1D *hpi0KE_truth, TH1D *hpi0KE_reco, TH1D *hpi0KE_data, TGraph* g_775, double Int_800_truth, double Int_800err_truth, double Int_800_reco, double Int_800err_reco, double Int_800_data, double Int_800err_data, double totalCEXXS_800_truth, double totalCEXXS_800_reco, double totalCEXXS_800_data, double totalCEXXS_800_truth_error, double totalCEXXS_800_reco_error, double totalCEXXS_800_data_error, TString name, bool doTruth = true, bool doReco = true, bool doData = true);
TH1D * TotalCEXXSCal(TList *lout, TH1D * hbeamint, TH1D * hh, TH1D * InteractingHist, const bool & Eslice, const bool & widerBin, const bool & newMethod, const TString tag);
TH1D * DiffCEXXSCal(TList *lout, TH1D * DiffCEXInteractingHist, const double &diffInt,  const double &diffInterror, const double &scale_bin, const double &scale_totalXS, const double &scale_totalXS_error, const TString tag);
void SetTitleFormat(TH1 * hh);
void DrawDataMCRatio(TH1D * hratio, bool xsec, const TString tag);
void WeightBeamData(TH1D * hbeam_data, const bool IsRange, const TString name);
void WeightPi0Data(TH1D * hpi0KE_data, const bool IsRange, const bool IsKF, const TString name);
double CalIniWeight(const int & binx, const bool IsRange);
double CalBeamIntWeight(const int & binx, const bool IsRange);
double CalCEXIntWeight(const int & binx, const bool IsRange);
double CalCEXPi0KEWeight(const int & binx, const bool IsRange, const bool IsKF);
double CalCEXPi0CosThetaWeight(const int & binx, const bool IsRange, const bool IsKF);
double CalCEXPi0ThetaWeight(const int & binx, const bool IsRange, const bool IsKF);

TH1D * RebinThetaHist(TH1D * hist);

void DrawSystError(TH1D * hh);

const bool IsRange = true;
const bool IsKF = true;

int main(int argc, char * argv[])
{
  gStyle->SetOptStat(0);

  //TString finName = "input/outana_trkLenhigh.root";
  //TString finName = "input/outana_April_last_DataTree_final.root";
  //TString finName = "input/outana_MCXSnorminal.root";
  //TString finName = "input/outana_withThresholdnew.root";
  TString finName = "input/outana.root";
  if(IsRange && !IsKF) finName = "input/outana_noKF.root";

  //if(IsRange && IsKF) finName = "input/outana.root";

  TFile *file = TFile::Open(finName);

  if(!file->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }


  // ========== Initial Histogram ========== //
  RooUnfoldResponse *response_SliceID_Ini = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Ini");
  TH1D *hini_truth = (TH1D*) file->Get("mc/i000hTruthBeamInitialHist"); // Truth
  TH1D *hini = (TH1D*) file->Get("mc/i016hRecoInitialHist"); // Fake data
  TH1D *hini_data = (TH1D*) file->Get("data/i019hRecoBeamInitialHistData"); // Data

  // Data bck sub
  //WeightBeamData(hini_data, IsRange, "Ini");

  RooUnfoldBayes unfold_Ini (response_SliceID_Ini, hini, 4); // Unfold fake data
  RooUnfoldBayes unfold_Ini_Data (response_SliceID_Ini, hini_data, 4); // Unfold data

  TH1D *hini_uf = (TH1D*)unfold_Ini.Hreco();
  TH1D *hini_uf_data = (TH1D*)unfold_Ini_Data.Hreco();

  // Draw Hist using tmp hist to rebin
  TH1D *hini_truth_tmp = (TH1D*)hini_truth->Clone();
  TH1D *hini_tmp = (TH1D*)hini->Clone();
  TH1D *hini_data_tmp = (TH1D*)hini_data->Clone();

  TH1D *hini_uf_tmp = (TH1D*)hini_uf->Clone();
  TH1D *hini_uf_data_tmp = (TH1D*)hini_uf_data->Clone();

  DrawHistOutput(hini_truth_tmp, hini_tmp, hini_data_tmp, hini_uf_tmp, hini_uf_data_tmp, true); // rebin == true

  // ========== BeamInt Histogram ========== //
  RooUnfoldResponse *response_SliceID_BeamInt = (RooUnfoldResponse*)file->Get("mc/response_SliceID_BeamInt");
  TH1D *hbeamint_truth = (TH1D*) file->Get("mc/i002hTruthBeamInteractingHist");
  TH1D *hbeamint = (TH1D*) file->Get("mc/i016hRecoBeamInteractingHist");
  TH1D *hbeamint_data = (TH1D*) file->Get("data/i020hRecoBeamInteractingHistData");

  // Data bck sub
  //WeightBeamData(hbeamint_data, IsRange, "BeamInt");

  RooUnfoldBayes unfold_BeamInt (response_SliceID_BeamInt, hbeamint, 4);
  RooUnfoldBayes unfold_BeamInt_Data (response_SliceID_BeamInt, hbeamint_data, 4);

  TH1D *hbeamint_uf = (TH1D*)unfold_BeamInt.Hreco();
  TH1D *hbeamint_uf_data = (TH1D*)unfold_BeamInt_Data.Hreco();

  // Draw Hist
  TH1D *hbeamint_truth_tmp = (TH1D*)hbeamint_truth->Clone();
  TH1D *hbeamint_tmp = (TH1D*)hbeamint->Clone();
  TH1D *hbeamint_data_tmp = (TH1D*)hbeamint_data->Clone();

  TH1D *hbeamint_uf_tmp = (TH1D*)hbeamint_uf->Clone();
  TH1D *hbeamint_uf_data_tmp = (TH1D*)hbeamint_uf_data->Clone();

  DrawHistOutput(hbeamint_truth_tmp, hbeamint_tmp, hbeamint_data_tmp, hbeamint_uf_tmp, hbeamint_uf_data_tmp, true); // rebin == true

  // ========== Int Histogram ========== //
  RooUnfoldResponse *response_SliceID_Int = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Int");
  TH1D *hint_truth = (TH1D*) file->Get("mc/i004hTruthCEXInteractingHist");
  TH1D *hint = (TH1D*) file->Get("mc/i018hRecoInteractingHist");
  TH1D *hint_data = (TH1D*) file->Get("data/i018hRecoInteractingHistData");

  // Data bck sub
  WeightBeamData(hint_data, IsRange, "CEXInt");

  RooUnfoldBayes unfold_Int (response_SliceID_Int, hint, 4);
  RooUnfoldBayes unfold_Int_Data (response_SliceID_Int, hint_data, 4);

  TH1D *hint_uf = (TH1D*)unfold_Int.Hreco();
  TH1D *hint_uf_data = (TH1D*)unfold_Int_Data.Hreco();
  hint_uf->Scale(0.5); // Need to scale down the unfolded fake data since used half truth MC
  //hint_uf_data->Scale(0.5);


  DrawHistOutput(hint_truth, hint, hint_data, hint_uf, hint_uf_data, false); // rebin == false

  // ======= GEANT4 Input XS ======== //
  const TString xcesfinName = "../Analysis/input/LAr_PiPlus_exclusive_cross_section.root";
  TFile *f_CrossSection = TFile::Open(xcesfinName);
  if(!f_CrossSection->IsOpen()){
    cout << "xsce file not open" << endl;
    exit(1);
  }

  TGraph *g_cex;

  f_CrossSection->GetObject("cex_KE",g_cex);


  //const TString DiffxcesfinName = "../Analysis/input/cross_section_out_withNewVar.root";
  const TString DiffxcesfinName = "input/cross_section_out.root";
  TFile *f_DiffCrossSection = TFile::Open(DiffxcesfinName);
  if(!f_DiffCrossSection->IsOpen()) {
    cout << "Diffxsce file not open" << endl;
    exit(1);
  }

  TH1D * g_cex_775MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dKEpi01775_MeV");
  const Int_t x0_diff775 = g_cex_775MeV->GetXaxis()->GetFirst();
  const Int_t x1_diff775 = g_cex_775MeV->GetXaxis()->GetLast();
  double x_775[x1_diff775], y_775[x1_diff775];
  for(Int_t ix=x0_diff775; ix<=x1_diff775; ix++){
    //cout << "Bin 775: "  << ix << " " << g_cex_775MeV->GetBinContent(ix) << endl;
    x_775[ix] = g_cex_775MeV->GetBinCenter(ix);
    y_775[ix] = g_cex_775MeV->GetBinContent(ix);
  }

  TGraph *g_775 = new TGraph(40,x_775,y_775);

  // ======== Calculate incident histogram ======== //
  // Truth
  TH1D * hinc_truth = GetIncidentHist(hini_truth,hbeamint_truth);
  // Unfolded reco
  TH1D * hinc_uf = GetIncidentHist(hini_uf,hbeamint_uf);
  // Data
  TH1D * hinc_uf_data = GetIncidentHist(hini_uf_data,hbeamint_uf_data);
  // Need to rebin the histogram after calculation for XS measurements
  hinc_truth->Rebin(50); //hint_truth->Rebin(50);
  hinc_uf->Rebin(50); //hint_uf->Rebin(50);
  hinc_uf_data->Rebin(50); //hint_uf_data->Rebin(50);

    // Declare output list
  TList * datalout = 0x0;
  datalout = new TList;

  // ======== Draw Total CEX XS plots ======== //
  double totalCEXXS_800_truth = 0, totalCEXXS_800_reco = 0, totalCEXXS_800_data = 0;
  double totalCEXXS_800_truth_error = 0, totalCEXXS_800_reco_error = 0, totalCEXXS_800_data_error = 0;

  // Overlay
  DrawTotalXSPlots(datalout,hbeamint_truth,hinc_truth,hint_truth,hbeamint_uf,hinc_uf,hint_uf,hbeamint_uf_data,hinc_uf_data,hint_uf_data,g_cex,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data,totalCEXXS_800_truth_error,totalCEXXS_800_reco_error,totalCEXXS_800_data_error,true,true,true);
  // Truth
  DrawTotalXSPlots(datalout,hbeamint_truth,hinc_truth,hint_truth,hbeamint_uf,hinc_uf,hint_uf,hbeamint_uf_data,hinc_uf_data,hint_uf_data,g_cex,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data,totalCEXXS_800_truth_error,totalCEXXS_800_reco_error,totalCEXXS_800_data_error,true,false,false);
  // Reco
  DrawTotalXSPlots(datalout,hbeamint_truth,hinc_truth,hint_truth,hbeamint_uf,hinc_uf,hint_uf,hbeamint_uf_data,hinc_uf_data,hint_uf_data,g_cex,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data,totalCEXXS_800_truth_error,totalCEXXS_800_reco_error,totalCEXXS_800_data_error,false,true,false);
  // Data
  DrawTotalXSPlots(datalout,hbeamint_truth,hinc_truth,hint_truth,hbeamint_uf,hinc_uf,hint_uf,hbeamint_uf_data,hinc_uf_data,hint_uf_data,g_cex,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data,totalCEXXS_800_truth_error,totalCEXXS_800_reco_error,totalCEXXS_800_data_error,false,false,true);

  // ======== Pi0 KE Histogram (KE PI+ = 650to800 MeV)======= //
  RooUnfoldResponse *response_SliceID_Pi0KE = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Pi0KE");
  TH1D *hpi0KE_truth = (TH1D*) file->Get("mc/i006hTruthDiffCEXInteractingHist_800MeV");
  //if(IsRange) hpi0KE_truth = (TH1D*) file->Get("mc/i005hTruthDiffCEXInteractingHist_700MeV");
  if(IsRange) hpi0KE_truth = (TH1D*) file->Get("mc/i006hTruthDiffCEXInteractingHist_650to800MeV");
  TH1D *hpi0KE = (TH1D*) file->Get("mc/i030hRecoPi0KEHist");
  TH1D *hpi0KE_data = (TH1D*) file->Get("data/i030hRecoPi0KEHistData");

  // Data bck sub
  WeightPi0Data(hpi0KE_data, IsRange, IsKF, "KE");

  RooUnfoldBayes unfold_Pi0KE (response_SliceID_Pi0KE, hpi0KE, 4);
  RooUnfoldBayes unfold_Pi0KE_Data (response_SliceID_Pi0KE, hpi0KE_data, 4);

  TH1D *hpi0KE_uf = (TH1D*)unfold_Pi0KE.Hreco();
  TH1D *hpi0KE_uf_data = (TH1D*)unfold_Pi0KE_Data.Hreco();
  hpi0KE_uf->Scale(0.5);
  //hpi0KE_uf_data->Scale(0.5);
  //hpi0KE_truth->Rebin(2); hpi0KE->Rebin(2); hpi0KE_uf->Rebin(2); hpi0KE_uf_data->Rebin(2);
  TH1D *hpi0KE_truth_tmp = (TH1D*)hpi0KE_truth->Clone();
  TH1D *hpi0KE_tmp = (TH1D*)hpi0KE->Clone();
  TH1D *hpi0KE_data_tmp = (TH1D*)hpi0KE_data->Clone();

  TH1D *hpi0KE_uf_tmp = (TH1D*)hpi0KE_uf->Clone();
  TH1D *hpi0KE_uf_data_tmp = (TH1D*)hpi0KE_uf_data->Clone();
  hpi0KE_truth_tmp->Rebin(2); hpi0KE_tmp->Rebin(2); hpi0KE_data_tmp->Rebin(2); hpi0KE_uf_tmp->Rebin(2); hpi0KE_uf_data_tmp->Rebin(2);

  DrawHistOutput(hpi0KE_truth_tmp, hpi0KE_tmp, hpi0KE_data_tmp, hpi0KE_uf_tmp, hpi0KE_uf_data_tmp, false); // rebin == false

  double Int_800_truth = hint_truth->GetBinContent(16);
  if(IsRange) Int_800_truth = hint_truth->GetBinContent(14) + hint_truth->GetBinContent(15) + hint_truth->GetBinContent(16);
  //double Int_800err_truth = hint_truth->GetBinError(16);
  double Int_800err_truth = sqrt(pow(hint_truth->GetBinError(14),2)+pow(hint_truth->GetBinError(15),2)+pow(hint_truth->GetBinError(16),2));

  double Int_800_reco = hint_uf->GetBinContent(16);
  if(IsRange) Int_800_reco = hint_uf->GetBinContent(14) + hint_uf->GetBinContent(15) + hint_uf->GetBinContent(16);
  //double Int_800err_reco = hint_uf->GetBinError(16);
  double Int_800err_reco = sqrt(pow(hint_uf->GetBinError(14),2)+pow(hint_uf->GetBinError(15),2)+pow(hint_uf->GetBinError(16),2));

  double Int_800_data = hint_uf_data->GetBinContent(16);
  if(IsRange) Int_800_data = hint_uf_data->GetBinContent(14) + hint_uf_data->GetBinContent(15) + hint_uf_data->GetBinContent(16);
  //double Int_800err_data = hint_uf_data->GetBinError(16);
  double Int_800err_data = sqrt(pow(hint_uf_data->GetBinError(14),2)+pow(hint_uf_data->GetBinError(15),2)+pow(hint_uf_data->GetBinError(16),2));

  if(!IsRange) DrawDiffXSPlots(datalout,hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_775,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "775");

  // ======= Plot Pi0KE ======== //

  const TString finNameXS = "input/cross_section_out.root";
  TFile *fileXS = TFile::Open(finNameXS);

  TH1D * g_cex_650to800MeV = (TH1D*)fileXS->Get("inel_cex_1dKEpi0650to800_MeV");

  const Int_t x0_diff650to800 = g_cex_650to800MeV->GetXaxis()->GetFirst();
  const Int_t x1_diff650to800 = g_cex_650to800MeV->GetXaxis()->GetLast();
  double x_650to800[x1_diff650to800], y_650to800[x1_diff650to800];
  for(Int_t ix=x0_diff650to800; ix<=x1_diff650to800; ix++){
    x_650to800[ix] = g_cex_650to800MeV->GetBinCenter(ix);
    y_650to800[ix] = g_cex_650to800MeV->GetBinContent(ix);
    //cout << "******* x_650to800 center = "<< g_cex_650to800MeV->GetBinCenter(ix) << " y_650to800 content " << g_cex_650to800MeV->GetBinContent(ix) << endl;
  }

  TGraph *g_650to800 = new TGraph(40,x_650to800,y_650to800);

  // overlay
  if(IsRange) DrawDiffXSPlots(datalout,hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_650to800,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800");
  // truth
  if(IsRange) DrawDiffXSPlots(datalout,hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_650to800,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800",true,false,false);
  // reco
  if(IsRange) DrawDiffXSPlots(datalout,hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_650to800,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800",false,true,false);
  // data
  if(IsRange) DrawDiffXSPlots(datalout,hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_650to800,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800",false,false,true);


  // ======== Pi0 CosTheta Histogram (KE PI+ = 650to800 MeV)======= //
  RooUnfoldResponse *response_SliceID_Pi0CosTheta = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Pi0CosTheta");
  TH1D *hpi0CosTheta_truth = (TH1D*) file->Get("mc/i666hTruthDiffCEXInteractingHistCosTheta_650to800MeV");
  TH1D *hpi0CosTheta = (TH1D*) file->Get("mc/i032hRecoPi0CosThetaHist");
  TH1D *hpi0CosTheta_data = (TH1D*) file->Get("data/i032hRecoPi0CosThetaHistData");

  // Data bck sub
  WeightPi0Data(hpi0CosTheta_data, IsRange, IsKF, "CosTheta");

  RooUnfoldBayes unfold_Pi0CosTheta (response_SliceID_Pi0CosTheta, hpi0CosTheta, 4);
  RooUnfoldBayes unfold_Pi0CosTheta_Data (response_SliceID_Pi0CosTheta, hpi0CosTheta_data, 4);

  TH1D *hpi0CosTheta_uf = (TH1D*)unfold_Pi0CosTheta.Hreco();
  TH1D *hpi0CosTheta_uf_data = (TH1D*)unfold_Pi0CosTheta_Data.Hreco();
  hpi0CosTheta_uf->Scale(0.5);
  //hpi0CosTheta_uf_data->Scale(0.5);

  DrawHistOutput(hpi0CosTheta_truth, hpi0CosTheta, hpi0CosTheta_data, hpi0CosTheta_uf, hpi0CosTheta_uf_data, false); // rebin == false

  TH1D * g_cex_650to800MeV_CosTheta = (TH1D*)fileXS->Get("inel_cex_1dcosThetapi0650to800_MeV");
  const Int_t x0_diff650to800_CosTheta = g_cex_650to800MeV_CosTheta->GetXaxis()->GetFirst();
  const Int_t x1_diff650to800_CosTheta = g_cex_650to800MeV_CosTheta->GetXaxis()->GetLast();
  double x_650to800_CosTheta[x1_diff650to800_CosTheta], y_650to800_CosTheta[x1_diff650to800_CosTheta];
  for(Int_t ix=x0_diff650to800_CosTheta; ix<=x1_diff650to800_CosTheta; ix++){
    //cout << "Bin 650to800: "  << ix << " " << g_cex_650to800MeV_CosTheta->GetBinContent(ix) << endl;
    x_650to800_CosTheta[ix] = g_cex_650to800MeV_CosTheta->GetBinCenter(ix);
    y_650to800_CosTheta[ix] = g_cex_650to800MeV_CosTheta->GetBinContent(ix);
  }

  TGraph *g_650to800_CosTheta = new TGraph(40,x_650to800_CosTheta,y_650to800_CosTheta);
  g_650to800_CosTheta->RemovePoint(0);
  g_650to800_CosTheta->RemovePoint(1);

  // overlay
  if(IsRange) DrawDiffXSPlots(datalout,hpi0CosTheta_truth,hpi0CosTheta_uf,hpi0CosTheta_uf_data,g_650to800_CosTheta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800CosTheta");
  // truth
  if(IsRange) DrawDiffXSPlots(datalout,hpi0CosTheta_truth,hpi0CosTheta_uf,hpi0CosTheta_uf_data,g_650to800_CosTheta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800CosTheta",true,false,false);
  // reco
  if(IsRange) DrawDiffXSPlots(datalout,hpi0CosTheta_truth,hpi0CosTheta_uf,hpi0CosTheta_uf_data,g_650to800_CosTheta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800CosTheta",false,true,false);
  // data
  if(IsRange) DrawDiffXSPlots(datalout,hpi0CosTheta_truth,hpi0CosTheta_uf,hpi0CosTheta_uf_data,g_650to800_CosTheta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800CosTheta",false,false,true);


  // ======== Pi0 Theta Histogram (KE PI+ = 650to800 MeV)======= //
  RooUnfoldResponse *response_SliceID_Pi0Theta = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Pi0Theta");
  TH1D *hpi0Theta_truth = (TH1D*) file->Get("mc/i066hTruthDiffCEXInteractingHistTheta_650to800MeV");
  TH1D *hpi0Theta = (TH1D*) file->Get("mc/i034hRecoPi0ThetaHist");
  TH1D *hpi0Theta_data = (TH1D*) file->Get("data/i034hRecoPi0ThetaHistData");

  // Data bck sub
  WeightPi0Data(hpi0Theta_data, IsRange, IsKF, "Theta");

  RooUnfoldBayes unfold_Pi0Theta (response_SliceID_Pi0Theta, hpi0Theta, 4);
  RooUnfoldBayes unfold_Pi0Theta_Data (response_SliceID_Pi0Theta, hpi0Theta_data, 4);

  TH1D *hpi0Theta_uf = (TH1D*)unfold_Pi0Theta.Hreco();
  TH1D *hpi0Theta_uf_data = (TH1D*)unfold_Pi0Theta_Data.Hreco();
  hpi0Theta_uf->Scale(0.5);
  //hpi0Theta_uf_data->Scale(0.5);

  DrawHistOutput(hpi0Theta_truth, hpi0Theta, hpi0Theta_data, hpi0Theta_uf, hpi0Theta_uf_data, false); // rebin == false

  TH1D * g_cex_650to800MeV_Theta = (TH1D*)fileXS->Get("inel_cex_1dThetapi0650to800_MeV");
  const Int_t x0_diff650to800_Theta = g_cex_650to800MeV_Theta->GetXaxis()->GetFirst();
  const Int_t x1_diff650to800_Theta = g_cex_650to800MeV_Theta->GetXaxis()->GetLast();
  double x_650to800_Theta[x1_diff650to800_Theta], y_650to800_Theta[x1_diff650to800_Theta];
  for(Int_t ix=x0_diff650to800_Theta; ix<=x1_diff650to800_Theta; ix++){
    //cout << "Bin 650to800: "  << ix << " " << g_cex_650to800MeV_Theta->GetBinContent(ix) << endl;
    x_650to800_Theta[ix] = g_cex_650to800MeV_Theta->GetBinCenter(ix);
    y_650to800_Theta[ix] = g_cex_650to800MeV_Theta->GetBinContent(ix);
  }

  TGraph *g_650to800_Theta = new TGraph(40,x_650to800_Theta,y_650to800_Theta);

  // overlay
  if(IsRange) DrawDiffXSPlots(datalout,hpi0Theta_truth,hpi0Theta_uf,hpi0Theta_uf_data,g_650to800_Theta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800Theta");
  // truth
  if(IsRange) DrawDiffXSPlots(datalout,hpi0Theta_truth,hpi0Theta_uf,hpi0Theta_uf_data,g_650to800_Theta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800Theta",true,false,false);
  // reco
  if(IsRange) DrawDiffXSPlots(datalout,hpi0Theta_truth,hpi0Theta_uf,hpi0Theta_uf_data,g_650to800_Theta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800Theta",false,true,false);
  // data
  if(IsRange) DrawDiffXSPlots(datalout,hpi0Theta_truth,hpi0Theta_uf,hpi0Theta_uf_data,g_650to800_Theta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, totalCEXXS_800_truth_error, totalCEXXS_800_reco_error, totalCEXXS_800_data_error, "650to800Theta",false,false,true);


  // Declare output root file
  TFile * fout = new TFile("output/outXShists.root","recreate");
  // Create data subdirectory
  TDirectory * ld_data = gDirectory->mkdir("data");
  ld_data->cd();
  // Write list to the root file
  datalout->Write();
  gDirectory->cd("../");
  // Save the info
  fout->Save();
  fout->Close();


  // Close files
  f_CrossSection->Close();
  f_DiffCrossSection->Close();
  file->Close();
  fileXS->Close();

}

void DrawHistOutput(TH1D *h_truth, TH1D *h_reco, TH1D *h_data, TH1D *h_unfold, TH1D *h_unfold_data, const bool &rebin){
  // Rebin the histogram
  if(rebin){
    h_truth->Rebin(50); h_reco->Rebin(50); h_data->Rebin(50); h_unfold->Rebin(50); h_unfold_data->Rebin(50);
  }
  TString tag = h_truth->GetName();

  TLatex tt;
  tt.SetNDC();

  TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 1200, 800);

  SetTitleFormat(h_truth);

  h_truth->SetLineColor(kRed);
  h_truth->SetLineWidth(2);
  h_truth->SetMaximum(h_truth->GetMaximum()*1.5);
  if(tag.Contains("i006")) h_truth->SetMaximum(1500);
  h_truth->SetMinimum(0.0);
  h_truth->GetYaxis()->SetTitle("Candidates");
  h_truth->GetXaxis()->SetTitle("T_{#pi^{+}} (MeV)");
  if(tag.Contains("i006")) h_truth->GetXaxis()->SetTitle("T_{#pi^{0}} (MeV)");
  if(tag.Contains("i066")) h_truth->GetXaxis()->SetTitle("#theta_{#pi^{0}} (deg.)");
  if(tag.Contains("i666")) h_truth->GetXaxis()->SetTitle("cos#theta_{#pi^{0}}");

  h_truth->SetFillStyle(1001);
  h_truth->SetLineWidth(3);
  h_truth->SetFillColorAlpha(kRed,0.3);

  h_truth->Draw("hist");
/*
  h_reco->SetMarkerStyle(22);
  h_reco->SetMarkerSize(1.5);
  h_reco->SetMarkerColor(kBlue);
  h_reco->SetLineColor(kBlue);
  h_reco->SetLineWidth(1);
  h_reco->Draw("e1 sames");
*/
  h_data->SetFillStyle(1001);
  h_data->SetLineWidth(3);
  h_data->SetLineStyle(7);
  h_data->SetFillColorAlpha(kBlue,0.3);
  h_data->Draw("e1 hist sames");


  cout << "tag: " << tag << endl;
  if(tag.Contains("i000") || tag.Contains("i002")) cout << "eff: " << h_reco->Integral(0,10000)/h_truth->Integral(0,10000) << endl;
  else cout << "eff: " << h_reco->Integral(0,10000)/2/h_truth->Integral(0,10000) << endl;

/*
  h_unfold->SetMarkerStyle(23);
  h_unfold->SetMarkerSize(1.5);
  h_unfold->SetMarkerColor(kGreen+3);
  h_unfold->SetLineColor(kGreen+3);
  h_unfold->SetLineWidth(1);
  h_unfold->Draw("e1 sames");
*/

  h_unfold_data->SetMarkerStyle(8);
  h_unfold_data->SetMarkerSize(1);
  h_unfold_data->SetMarkerColor(kBlack);
  h_unfold_data->SetLineColor(kBlack);
  h_unfold_data->SetLineWidth(1);
  h_unfold_data->Draw("e1 sames");

  TLegend * legend = new TLegend(0.15, 0.62, 0.48, 0.88);
  if(tag.Contains("650to800") && !tag.Contains("CosTheta")) legend = new TLegend(0.55, 0.62, 0.88, 0.88);
  legend->AddEntry(h_truth, "MC Truth", "lf");
  legend->AddEntry(h_data, "Data - BckSub.", "lf");
  //legend->AddEntry(h_unfold, "MC Reco - Unfolded", "lp");
  legend->AddEntry(h_unfold_data, "Data - Unfolded", "lp");
  legend->SetBorderSize(0);
  legend->Draw("same");

  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  if(tag.Contains("i006")||tag.Contains("i066")||tag.Contains("i666")) tt.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
  else tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
  if(tag.Contains("i006")||tag.Contains("i066")) {
    tt.DrawLatex(0.155,0.825,"#bf{#it{Preliminary}}");
    tt.DrawLatex(0.155,0.775,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");
  }
  else {
    tt.DrawLatex(0.705,0.825,"#bf{#it{Preliminary}}");
    tt.DrawLatex(0.635,0.775,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");
  }
  //c1->Print("output/hint_uf.pdf");
  c1->Print("output/"+tag+".pdf");

}

void DrawTotalXSPlots(TList *lout, TH1D *hbeamint_truth, TH1D *hinc_truth, TH1D *hint_truth, TH1D *hbeamint_reco, TH1D *hinc_reco, TH1D *hint_reco, TH1D *hbeamint_data, TH1D *hinc_data, TH1D *hint_data, TGraph* g_cex, double & totalCEXXS_800_truth, double & totalCEXXS_800_reco, double & totalCEXXS_800_data, double & totalCEXXS_800_truth_error, double & totalCEXXS_800_reco_error, double & totalCEXXS_800_data_error, bool doTruth, bool doReco, bool doData){

  TLatex tt;
  tt.SetNDC();
  TString label = "overlay"; TString name = "Total";
  if(doTruth && !doReco && !doData) label = "truth";
  if(!doTruth && doReco && !doData) label = "reco";
  if(!doTruth && !doReco && doData) label = "data";


  TCanvas * c1 = new TCanvas(Form("c1_%s_%s",name.Data(),label.Data()), Form("c1_%s_%s",name.Data(),label.Data()), 1200, 800);

  // Truth
  TH1D * hcex_truth = TotalCEXXSCal(lout,hbeamint_truth,hinc_truth,hint_truth,true,false,false,"truth");
  // Reco
  TH1D * hcex_reco = TotalCEXXSCal(lout,hbeamint_reco,hinc_reco,hint_reco,true,false,false,"reco");
  // Data
  TH1D * hcex_data = TotalCEXXSCal(lout,hbeamint_data,hinc_data,hint_data,true,false,false,"data");

  totalCEXXS_800_truth = (hcex_truth->GetBinContent(14) + hcex_truth->GetBinContent(15) + hcex_truth->GetBinContent(16))/3.0;
  totalCEXXS_800_reco = (hcex_reco->GetBinContent(14) + hcex_reco->GetBinContent(15) + hcex_reco->GetBinContent(16))/3.0;
  totalCEXXS_800_data = (hcex_data->GetBinContent(14) + hcex_data->GetBinContent(15) + hcex_data->GetBinContent(16))/3.0;

  totalCEXXS_800_truth_error = sqrt(pow(hcex_truth->GetBinError(14),2)+pow(hcex_truth->GetBinError(15),2)+pow(hcex_truth->GetBinError(16),2))/3.0;
  totalCEXXS_800_reco_error = sqrt(pow(hcex_reco->GetBinError(14),2)+pow(hcex_reco->GetBinError(15),2)+pow(hcex_reco->GetBinError(16),2))/3.0;
  totalCEXXS_800_data_error = sqrt(pow(hcex_data->GetBinError(14),2)+pow(hcex_data->GetBinError(15),2)+pow(hcex_data->GetBinError(16),2))/3.0;

  SetTitleFormat(hcex_truth);
  SetTitleFormat(hcex_reco);
  SetTitleFormat(hcex_data);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 0.98);
  // Set 0.01 will hide axis label
  pad1->SetBottomMargin(0.01);
  //  pad1->SetGridx();
  //  pad1->SetGridy();
  pad1->Draw();
  pad1->cd();

  hcex_truth->SetMarkerStyle(22);
  hcex_truth->SetMarkerSize(1.5);
  hcex_truth->SetMarkerColor(kBlue);
  hcex_truth->SetLineColor(kBlue);
  hcex_truth->SetLineWidth(1);
  //hcex_truth->SetMaximum(0.4);
  hcex_truth->SetMinimum(0.);

  hcex_truth->GetYaxis()->SetTitle("#sigma_{CEX} (mb)");
  if(doTruth) hcex_truth->Draw("E1 sames");

  hcex_reco->SetMarkerStyle(23);
  hcex_reco->SetMarkerSize(1.5);
  hcex_reco->SetMarkerColor(kGreen+3);
  hcex_reco->SetLineColor(kGreen+3);
  hcex_reco->SetLineWidth(1);
  //hcex_reco->SetMaximum(0.4);
  hcex_reco->SetMinimum(0.);

  hcex_reco->GetYaxis()->SetTitle("#sigma_{CEX} (mb)");
  if(doReco) hcex_reco->Draw("E1 sames");


  hcex_data->SetMarkerStyle(8);
  hcex_data->SetMarkerSize(1);
  hcex_data->SetMarkerColor(kBlack);
  hcex_data->SetLineColor(kBlack);
  hcex_data->SetLineWidth(1);
  //hcex_data->SetMaximum(0.4);
  hcex_data->SetMinimum(0.);
  hcex_data->GetYaxis()->SetTitle("#sigma_{CEX} (mb)");

  if(doData) DrawSystError(hcex_data); // it works but not done

  if(doData) hcex_data->Draw("E1 sames");
  if(doData) {
    hcex_reco->SetLineWidth(2);
    hcex_reco->SetLineColor(kGray+2);
    hcex_reco->SetMarkerSize(0);
    hcex_reco->SetFillColorAlpha(kGray+1, 0.4);
    hcex_reco->Draw("e2 sames");
    TH1D *hcex_reco_overlay = (TH1D*)hcex_reco->Clone("hcex_reco");
    hcex_reco_overlay->SetFillStyle(0);
    hcex_reco_overlay->SetLineWidth(2);
    hcex_reco_overlay->Draw("HIST sames");
  }

  g_cex->SetLineColor(kRed);
  g_cex->SetLineWidth(2);
  g_cex->Draw("sames C");

  TLegend * legend = new TLegend(0.15, 0.48, 0.48, 0.88);
  legend->AddEntry(g_cex, "Geant4 v4.10.6", "l");
  legend->AddEntry(hcex_truth, "MC Truth", "lep");
  legend->AddEntry(hcex_reco, "MC Reco", "lep");
  legend->AddEntry(hcex_data, "Data", "lep");
  //TString lheader("1GeV Pion Data");
  //legend->SetHeader(lheader);
  legend->SetBorderSize(0);
  legend->Draw("same");


  c1->Update();
  c1->cd();

  // ======= pad2 =======//
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.2);

  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.42);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  // ratio with data
  TH1D *hratio = (TH1D*)hcex_data->Clone("hcex_data");
  if(doTruth && !doReco && !doData) hratio = (TH1D*)hcex_truth->Clone("hcex_truth");
  if(!doTruth && doReco && !doData) hratio = (TH1D*)hcex_reco->Clone("hcex_reco");

  const Int_t x0_x = hcex_data->GetXaxis()->GetFirst();
  const Int_t x1_x = hcex_data->GetXaxis()->GetLast();

  hratio->Scale(0);

  for(Int_t ix=x0_x; ix<=x1_x; ix++){
    double bin_center = hcex_data->GetBinCenter(ix);

    double MC = g_cex->Eval(bin_center);

    double data = hcex_data->GetBinContent(ix);
    double edata = hcex_data->GetBinError(ix);

    if(doTruth && !doReco && !doData) {data = hcex_truth->GetBinContent(ix); edata = hcex_truth->GetBinError(ix);}
    if(!doTruth && doReco && !doData) {data = hcex_reco->GetBinContent(ix); edata = hcex_reco->GetBinError(ix);}

    if(MC != 0){
      double ratio = data/MC;
      hratio->SetBinContent(ix,ratio);
      double error = sqrt(ratio*ratio*(pow(edata/data,2)));
      hratio->SetBinError(ix,error);
    }
  }
  hratio->SetTitle(" ");

  DrawDataMCRatio(hratio,true,"");
  c1->cd();

  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
  tt.DrawLatex(0.705,0.825,"#bf{#it{Preliminary}}");
  //tt.DrawLatex(0.705,0.775,"#bf{#it{(Stats. Error Only)}}");
  tt.DrawLatex(0.685,0.775,"#bf{#it{(Stats. + Syst. Error)}}");

  if(doTruth && !doReco && !doData) c1->Print("output/truth_TotalCEXXS.pdf");
  if(!doTruth && doReco && !doData) c1->Print("output/reco_TotalCEXXS.pdf");

  if(!doTruth && !doReco && doData) c1->Print("output/data_TotalCEXXS.pdf");

  if(doTruth && doReco && doData) c1->Print("output/overlay_TotalCEXXS.pdf");


}


void DrawDiffXSPlots(TList *lout, TH1D *hpi0KE_truth, TH1D *hpi0KE_reco, TH1D *hpi0KE_data, TGraph* g_775, double Int_800_truth, double Int_800err_truth, double Int_800_reco, double Int_800err_reco, double Int_800_data, double Int_800err_data, double totalCEXXS_800_truth, double totalCEXXS_800_reco, double totalCEXXS_800_data, double totalCEXXS_800_truth_error, double totalCEXXS_800_reco_error, double totalCEXXS_800_data_error, TString name, bool doTruth, bool doReco, bool doData){

  TLatex tt;
  tt.SetNDC();

  TString label = "overlay";
  if(doTruth && !doReco && !doData) label = "truth";
  if(!doTruth && doReco && !doData) label = "reco";
  if(!doTruth && !doReco && doData) label = "data";

  TCanvas * c2 = new TCanvas(Form("c2_%s_%s",name.Data(),label.Data()), Form("c2_%s_%s",name.Data(),label.Data()), 1200, 800);

  double scale_bin = 100.0;
  if(name.Contains("CosTheta")) scale_bin = 0.2;
  if(name.Contains("Theta") && !name.Contains("Cos")) scale_bin = 18.0;

  double scale_totalXS_truth = totalCEXXS_800_truth;
  double scale_totalXS_reco = totalCEXXS_800_reco;
  double scale_totalXS_data = totalCEXXS_800_data;

  double scale_totalXS_truth_error = totalCEXXS_800_truth_error;
  double scale_totalXS_reco_error = totalCEXXS_800_reco_error;
  double scale_totalXS_data_error = totalCEXXS_800_data_error;

  TH1D * hdiffcex_truth = DiffCEXXSCal(lout,hpi0KE_truth,Int_800_truth,Int_800err_truth, scale_bin, scale_totalXS_truth, scale_totalXS_truth_error, "truth");
  TH1D * hdiffcex_reco = DiffCEXXSCal(lout,hpi0KE_reco,Int_800_reco,Int_800err_reco, scale_bin, scale_totalXS_reco, scale_totalXS_reco_error, "reco");
  TH1D * hdiffcex_data = DiffCEXXSCal(lout,hpi0KE_data,Int_800_data,Int_800err_data, scale_bin, scale_totalXS_data, scale_totalXS_data_error, "data");
  if(!name.Contains("Theta")){
    hdiffcex_truth->Rebin(2);
    hdiffcex_reco->Rebin(2);
    hdiffcex_data->Rebin(2);
  }

  hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/dT_{#pi^{0}} (mb/MeV)");
  if(name.Contains("CosTheta")) hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi^{0}} (mb)");
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/d#theta_{#pi^{0}} (mb/deg.)");

  SetTitleFormat(hdiffcex_truth);
  SetTitleFormat(hdiffcex_reco);
  SetTitleFormat(hdiffcex_data);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 0.98);
  // Set 0.01 will hide axis label
  pad1->SetBottomMargin(0.01);
  //  pad1->SetGridx();
  //  pad1->SetGridy();
  pad1->Draw();
  pad1->cd();

  hdiffcex_truth->SetMarkerStyle(22);
  hdiffcex_truth->SetMarkerSize(1.5);
  hdiffcex_truth->SetMarkerColor(kBlue);
  hdiffcex_truth->SetLineColor(kBlue);
  hdiffcex_truth->SetLineWidth(1);
  hdiffcex_truth->SetMaximum(0.4);
  if(name.Contains("CosTheta")) hdiffcex_truth->SetMaximum(250);
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_truth->SetMaximum(1.8);

  hdiffcex_truth->SetMinimum(0.);

  if(doTruth) hdiffcex_truth->Draw("E1 sames");


  hdiffcex_reco->GetYaxis()->SetTitle("d#sigma/dT_{#pi^{0}} (mb/MeV)");
  if(name.Contains("CosTheta")) hdiffcex_reco->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi^{0}} (mb)");
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_reco->GetYaxis()->SetTitle("d#sigma/d#theta_{#pi^{0}} (mb/deg.)");


  hdiffcex_reco->SetMarkerStyle(23);
  hdiffcex_reco->SetMarkerSize(1.5);
  hdiffcex_reco->SetMarkerColor(kGreen+3);
  hdiffcex_reco->SetLineColor(kGreen+3);
  hdiffcex_reco->SetLineWidth(1);
  hdiffcex_reco->SetMaximum(0.4);
  hdiffcex_reco->SetMinimum(0.);
  if(name.Contains("CosTheta")) hdiffcex_reco->SetMaximum(250);
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_reco->SetMaximum(1.8);

  if(doReco) hdiffcex_reco->Draw("E1 sames");

  hdiffcex_data->GetYaxis()->SetTitle("d#sigma/dT_{#pi^{0}} (mb/MeV)");
  if(name.Contains("CosTheta")) hdiffcex_data->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi^{0}} (mb)");
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_data->GetYaxis()->SetTitle("d#sigma/d#theta_{#pi^{0}} (mb/deg.)");



  hdiffcex_data->SetMarkerStyle(8);
  hdiffcex_data->SetMarkerSize(1);
  hdiffcex_data->SetMarkerColor(kBlack);
  hdiffcex_data->SetLineColor(kBlack);
  hdiffcex_data->SetLineWidth(1);
  hdiffcex_data->SetMaximum(0.4);
  hdiffcex_data->SetMinimum(0.);
  if(name.Contains("CosTheta")) hdiffcex_data->SetMaximum(250);
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_data->SetMaximum(1.8);


  if(doData) hdiffcex_data->Draw("E1 sames");
  if(doData) {
    hdiffcex_reco->SetLineWidth(2);
    hdiffcex_reco->SetLineColor(kGray+2);
    hdiffcex_reco->SetMarkerSize(0);
    hdiffcex_reco->SetFillColorAlpha(kGray+1, 0.4);
    hdiffcex_reco->Draw("e2 sames");
    TH1D *hdiffcex_reco_overlay = (TH1D*)hdiffcex_reco->Clone("hdiffcex_reco");
    hdiffcex_reco_overlay->SetFillStyle(0);
    hdiffcex_reco_overlay->SetLineWidth(2);
    hdiffcex_reco_overlay->Draw("HIST sames");
  }

  g_775->SetLineColor(kRed);
  g_775->SetLineWidth(2);
  g_775->Draw("sames C");

  TLegend * legend = new TLegend(0.55, 0.48, 0.88, 0.88);
  if(name.Contains("CosTheta")) legend = new TLegend(0.15, 0.48, 0.48, 0.88);

  legend->AddEntry(g_775, "Geant4 v4.10.6", "l");
  legend->AddEntry(hdiffcex_truth, "MC Truth", "lep");
  legend->AddEntry(hdiffcex_reco, "MC Reco", "lep");
  legend->AddEntry(hdiffcex_data, "Data", "lep");
  //TString lheader("1GeV Pion Data");
  //legend->SetHeader(lheader);
  legend->SetBorderSize(0);
  legend->Draw("same");


  c2->Update();
  c2->cd();

  // ======= pad2 =======//
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.2);

  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.42);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  // ratio with data
  TH1D *hratio = (TH1D*)hdiffcex_data->Clone("hdiffcex_data");
  if(doTruth && !doReco && !doData) hratio = (TH1D*)hdiffcex_truth->Clone("hdiffcex_truth");
  if(!doTruth && doReco && !doData) hratio = (TH1D*)hdiffcex_reco->Clone("hdiffcex_reco");


  const Int_t x0_x = hdiffcex_data->GetXaxis()->GetFirst();
  const Int_t x1_x = hdiffcex_data->GetXaxis()->GetLast();

  hratio->Scale(0);

  for(Int_t ix=x0_x; ix<=x1_x; ix++){
    double bin_center = hdiffcex_data->GetBinCenter(ix);
    double MC = g_775->Eval(bin_center);

    double data = hdiffcex_data->GetBinContent(ix);
    double edata = hdiffcex_data->GetBinError(ix);

    if(doTruth && !doReco && !doData) {data = hdiffcex_truth->GetBinContent(ix); edata = hdiffcex_truth->GetBinError(ix);}
    if(!doTruth && doReco && !doData) {data = hdiffcex_reco->GetBinContent(ix); edata = hdiffcex_reco->GetBinError(ix);}

    if(MC != 0){
      double ratio = data/MC;
      hratio->SetBinContent(ix,ratio);
      double error = sqrt(ratio*ratio*(pow(edata/data,2)));
      hratio->SetBinError(ix,error);
    }
  }
  hratio->SetTitle(" ");
  DrawDataMCRatio(hratio,false,name);
  c2->cd();

  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  tt.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
  if(!name.Contains("CosTheta")) tt.DrawLatex(0.22,0.825,"#bf{#it{Preliminary}}");
  if(name.Contains("CosTheta")) tt.DrawLatex(0.605,0.825,"#bf{#it{Preliminary}}");

  if(!name.Contains("CosTheta")) tt.DrawLatex(0.22,0.775,"#bf{#it{(Stats. Error Only)}}");
  if(name.Contains("CosTheta")) tt.DrawLatex(0.605,0.775,"#bf{#it{(Stats. Error Only)}}");

  //c2->Print("output/plot_Diff"+name+"CEXXS.pdf");

  if(doTruth && !doReco && !doData) c2->Print("output/truth_Diff"+name+"CEXXS.pdf");
  if(!doTruth && doReco && !doData) c2->Print("output/reco_Diff"+name+"CEXXS.pdf");

  if(!doTruth && !doReco && doData) c2->Print("output/data_Diff"+name+"CEXXS.pdf");

  if(doTruth && doReco && doData) c2->Print("output/overlay_Diff"+name+"CEXXS.pdf");

}


void SetTitleFormat(TH1 * hh){
  hh->SetTitle(" ");
  hh->GetYaxis()->CenterTitle();
  hh->GetYaxis()->SetTitleFont(22);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetYaxis()->SetTitleOffset(0.92);
  hh->GetXaxis()->CenterTitle();
  hh->GetXaxis()->SetTitleFont(22);
  hh->GetXaxis()->SetTitleSize(0.05);
  hh->GetXaxis()->SetTitleOffset(0.9);
  //hh->GetXaxis()->SetLabelOffset(999);
  //hh->GetXaxis()->SetLabelSize(0);
}



TH1D * TotalCEXXSCal(TList *lout, TH1D * hbeamint, TH1D * hh, TH1D * InteractingHist, const bool & Eslice, const bool & widerBin, const bool & newMethod, const TString tag)
{
  TH1D *xsec = (TH1D*)hh->Clone();
  // Scale down it to zero for later calculation
  xsec->Scale(0);

  double avogadro_constant = 6.02214076e23;  // 1 / mol
  double argon_molar_mass = 39.95;           // g / mol
  double liquid_argon_density = 1.39;        // g / cm^3
  double fiducial_thickness = 0.479;//slice_width;//0.479;         // cm wire spacing
  if(Eslice) fiducial_thickness = 1.0;
  if(Eslice && newMethod) fiducial_thickness = 1.0*50.0;
  if(Eslice && widerBin) fiducial_thickness = 1.0*50.0;
  //double fiducial_thickness = 20;         // cm wire spacing
  double sigma_factor = argon_molar_mass / (avogadro_constant * liquid_argon_density * fiducial_thickness); // pre-factor
  double meandEdx[] = {5.06921, 2.78604, 2.35879, 2.20657, 2.14137, 2.11321, 2.10322, 2.10293, 2.10805, 2.11628, 2.12627, 2.13724, 2.14869, 2.16033, 2.17195,
                       2.18344, 2.19472, 2.20573, 2.21646, 2.2269};
  // Define start and end bins
  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();

  // Loop over each bin and combine the incident and interacting histograms to get xsec
  for(Int_t ix=x0+9; ix<x1; ix++){
    double correction = 1.0;//InteractingHist->GetBinContent(ix)/hbeamint->GetBinContent(ix);
    cout << "correction: " << correction << endl;
    // Read off the entry for each bin of the two histograms
    double incE = hh->GetBinContent(ix);
    double intE = InteractingHist->GetBinContent(ix);
    //double intE = hbeamint->GetBinContent(ix);
    // Calculate the ratio in the full log form
    double ratio = 0;
    if(incE > intE) ratio = log(incE/(incE-intE));//*meandEdx[ix-1]*sigma_factor*1e27;

    if(Eslice) ratio = correction * ratio * meandEdx[ix-1]*sigma_factor*1e27;

    if(incE != 0) cout << "ix: "<< ix << "incE: " << incE << "intE: " << intE << " ratio: " << ratio << "meandEdx[ix]: " << meandEdx[ix-1] << endl;

    // Simple ratio form
    //double ratio = intE/incE;

    // If the incE entry is not zero set the bin content
    if(incE != 0) xsec->SetBinContent(ix,ratio);
    //if( Eslice && ratio > 150) xsec->SetBinContent(ix,9999);
    //if( Eslice && ratio < 50) xsec->SetBinContent(ix,9999);

    // Error propagation method 1
    //double error = sqrt(intE+pow(intE,2)/incE)/incE;
    // Error propagation method 2
    double einc = hh->GetBinError(ix);
    double eint = InteractingHist->GetBinError(ix);
    double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)));

    if( Eslice && error > 100) xsec->SetBinContent(ix,9999);

    // If the ratio is not zero set the error
    if(ratio != 0 ) xsec->SetBinError(ix,error);

  }
  // The xsec histogram entry is now set
  if(!Eslice) xsec->Scale(sigma_factor*1e27);
  xsec->SetMaximum(300);

  xsec->SetMinimum(0);
   /*for(Int_t ix=x0; ix<=x1; ix++){
    cout << "ix: "<< ix << "bin content: " << xsec->GetBinContent(ix) << endl;
   }*/
  if(tag.Contains("data")){
    xsec->SetName("data_TotalCEXxsec");
    lout->Add(xsec);
  }
  if(tag.Contains("reco")){
    xsec->SetName("reco_TotalCEXxsec");
    lout->Add(xsec);
  }
  if(tag.Contains("truth")){
    xsec->SetName("truth_TotalCEXxsec");
    lout->Add(xsec);
  }
  return xsec;
}

TH1D * DiffCEXXSCal(TList *lout, TH1D * DiffCEXInteractingHist, const double &diffInt,  const double &diffInterror, const double &scale_bin, const double &scale_totalXS, const double &scale_totalXS_error, const TString tag)
{
  TH1D *DiffCEXxsec = (TH1D*)DiffCEXInteractingHist->Clone();
  DiffCEXxsec->Scale(0);

  const Int_t x0_diff = DiffCEXInteractingHist->GetXaxis()->GetFirst();
  const Int_t x1_diff = DiffCEXInteractingHist->GetXaxis()->GetLast();

  for(Int_t ix=x0_diff; ix<=x1_diff; ix++){
    double incE = diffInt;
    double intE = DiffCEXInteractingHist->GetBinContent(ix);
    if(intE != 0 ){
      //double ratio = log(incE/(incE-intE));
      double ratio = intE/incE*scale_totalXS;
      DiffCEXxsec->SetBinContent(ix,ratio);
      //double error = sqrt(intE+pow(intE,2)/incE)/incE;
      double einc = diffInterror;
      double eint = DiffCEXInteractingHist->GetBinError(ix);
      double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)+pow(scale_totalXS_error/scale_totalXS,2)));
      DiffCEXxsec->SetBinError(ix,error);
      //cout << "error: " << error << endl;
      //cout << "ix: " << ix << endl;
      //cout << "ratio: " << ratio << endl;
    }
  }
  // Normalise to total cross section at that slice
  DiffCEXxsec->Scale(1/scale_bin);
  //DiffCEXxsec->Scale(scale_totalXS);
  if(tag.Contains("data")){
    TString name = DiffCEXxsec->GetName();
    if(name.Contains("Pi0KE")) DiffCEXxsec->SetName("data_DiffCEXxsec_Pi0KE");
    if(name.Contains("Pi0CosTheta")) DiffCEXxsec->SetName("data_DiffCEXxsec_Pi0CosTheta");
    if(name.Contains("Pi0Theta")) DiffCEXxsec->SetName("data_DiffCEXxsec_Pi0Theta");
    lout->Add(DiffCEXxsec);
  }
  if(tag.Contains("reco")){
    TString name = DiffCEXxsec->GetName();
    if(name.Contains("Pi0KE")) DiffCEXxsec->SetName("reco_DiffCEXxsec_Pi0KE");
    if(name.Contains("Pi0CosTheta")) DiffCEXxsec->SetName("reco_DiffCEXxsec_Pi0CosTheta");
    if(name.Contains("Pi0Theta")) DiffCEXxsec->SetName("reco_DiffCEXxsec_Pi0Theta");
    lout->Add(DiffCEXxsec);
  }
  if(tag.Contains("truth")){
    TString name = DiffCEXxsec->GetName();
    if(name.Contains("Pi0KE")) DiffCEXxsec->SetName("truth_DiffCEXxsec_Pi0KE");
    if(name.Contains("Pi0CosTheta")) DiffCEXxsec->SetName("truth_DiffCEXxsec_Pi0CosTheta");
    if(name.Contains("Pi0Theta")) DiffCEXxsec->SetName("truth_DiffCEXxsec_Pi0Theta");
    lout->Add(DiffCEXxsec);
  }
  return DiffCEXxsec;
}


TH1D * GetIncidentHist(TH1D * InitialHist, TH1D * InteractingHist)
{
  TH1D *IncidentHist = (TH1D*)InitialHist->Clone();
  IncidentHist->Scale(0);

  //double IniSum = 0, IntSum = 0;
  const Int_t x0 = InteractingHist->GetXaxis()->GetFirst();
  const Int_t x1 = InteractingHist->GetXaxis()->GetLast();
/*  for(Int_t ix=x1; ix>=x0; ix--){
    IniSum += InitialHist->GetBinContent(ix+1);
    //double error = InitialHist->GetBinError(ix);
    IntSum += InteractingHist->GetBinContent(ix+1);
    IncidentHist->SetBinContent(ix, IniSum-IntSum);
  }
*/

  for(Int_t ix=x1; ix>=x0; ix--){
    double iniN = 0;
    double intN = 0;
    for(Int_t jx=x1; jx>=ix; jx--){
      iniN += InitialHist->GetBinContent(jx);
    }
    for(Int_t jx=x1; jx>=ix+1; jx--){
      intN += InteractingHist->GetBinContent(jx);
    }
    /*for(Int_t jx=ix; jx>=x0; jx--){
      intN += InteractingHist->GetBinContent(jx);
    }
    for(Int_t jx=ix-1; jx>=x0; jx--){
      iniN += InitialHist->GetBinContent(jx);
    }*/
    //cout << "iniN: " << iniN << endl;
    //cout << "intN: " << intN << endl;
    //cout << "iniN - intN: " << iniN - intN << endl;

    //cout << "intN: " << intN << endl;
    //cout << "iniN: " << iniN << endl;
    //cout << "intN - iniN: " << intN - iniN << endl;

    //double entry = intN - iniN;
    double entry = iniN - intN;
    if (entry < 0) entry = 0;
    IncidentHist->SetBinContent(ix,entry);
    IncidentHist->SetBinError(ix,InitialHist->GetBinError(ix));

  }

  return IncidentHist;
}


void DrawDataMCRatio(TH1D * hratio, bool xsec, const TString tag){

  hratio->GetYaxis()->SetTitle("Data/MC");
  if(xsec) hratio->GetXaxis()->SetTitle("T_{#pi^{+}} (MeV)");
  if(!xsec && tag.Contains("CosTheta")) hratio->GetXaxis()->SetTitle("cos#theta_{#pi^{0}}");
  if(!xsec && tag.Contains("Theta") && !tag.Contains("Cos")) hratio->GetXaxis()->SetTitle("#theta_{#pi^{0}} (deg.)");

  hratio->GetXaxis()->SetLabelSize(0.15);
  hratio->GetXaxis()->SetTitleSize(0.15);
  hratio->GetXaxis()->SetTitleOffset(1.);
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetYaxis()->SetTitleOffset(.3);
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->GetYaxis()->SetRangeUser(0,2.1);
  hratio->GetYaxis()->SetRangeUser(0.5,1.6);
  hratio->GetXaxis()->CenterTitle();
  hratio->GetYaxis()->CenterTitle();

  hratio->GetXaxis()->SetTitleFont(22);
  hratio->GetYaxis()->SetTitleFont(22);
  hratio->GetXaxis()->SetTitleOffset(0.7);
  hratio->GetXaxis()->SetTitleSize(0.25);

  hratio->Draw();
}


// Ini Hist
double CalIniWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  /*if(IsRange){
    double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.801123, 0.859177, 0.86919, 0.880718, 0.886977, 0.889779, 0.887032};
    weight *= IntWeight[binx-1];
  }

  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517,
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }*/
  weight = 1.;
  return weight;
}

double CalBeamIntWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  /*if(IsRange){
    double IntWeight[] = {0.439768, 1.0, 0.0, 0.240104, 0.325475, 0.547874, 0.525402, 0.569011, 0.663258, 0.755661, 0.820706, 0.858978,
                     0.886342, 0.901809, 0.920025, 0.931599, 0.945724, 0.945179, 0.943406, 0.938913};
    weight *= IntWeight[binx-1];
  }

  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517,
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }*/
  weight = 1.;
  return weight;
}

double CalCEXIntWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  /*if(IsRange){
    //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.68125, 0.811044, 0.829454, 0.84448, 0.600005,
    //                0.731605, 0.683442, 0.650626, 0.651113, 0.668391, 0.633431, 0.65286, 1.0};
    // Feb
    //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.681394, 0.811332, 0.768364, 0.771537, 0.483201,
    //                0.597788, 0.496897, 0.699494, 0.705544, 0.718769, 0.687874, 0.693459, 1.0};
    // April
    double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.69794, 0.818919, 0.775483, 0.800256, 0.49334,
                    0.653571, 0.504115, 0.740037, 0.727999, 0.725943, 0.717317, 0.742088, 1.0};

    weight *= IntWeight[binx-1];
  }

  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517,
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }*/
  weight = 1.;
  return weight;
}


double CalCEXPi0KEWeight(const int & binx, const bool IsRange, const bool IsKF){

  double weight = 1.;
  /*if(IsRange && IsKF){
    //==  small bin beam wt
    //double IntWeight[] = {0.0, 0.409994, 0.357135, 0.382372, 0.465516, 0.597296, 0.672419, 0.564701, 0.61694, 0.571107, 0.791561,
    //                  0.733454, 0.945928, 0.941008, 0.785822, 1.0, 1,0, 0.0, 0.0, 0.0};
    //double IntWeight[] = {0.0, 0.539643, 0.378196, 0.398281, 0.504129, 0.658985, 0.715681, 0.638402, 0.66754, 0.70573, 0.791561,
    //                  0.733454, 0.945928, 0.941008, 0.785822, 1.0, 1,0, 0.0, 0.0, 0.0};
    // no wt
    //double IntWeight[] = {0.0, 0.627507, 0.429589, 0.45521, 0.562079, 0.702244, 0.776978, 0.695155, 0.6979, 0.763299, 0.812858,
    //                  0.780085, 0.950007, 0.946627, 0.863655, 1.0, 1,0, 0.0, 0.0, 0.0};
    // beam wt
    //double IntWeight[] = {0.0, 0.585378, 0.422816, 0.443579, 0.550453, 0.699467, 0.75336, 0.681768, 0.709005, 0.744256, 0.82169,
    //                  0.769536, 0.955012, 0.950875, 0.816588, 1.0, 1,0, 0.0, 0.0, 0.0};

    double IntWeight[] = {0.0, 0.586469, 0.388012, 0.412948, 0.519354, 0.665045, 0.740563, 0.651374, 0.654317, 0.72544, 0.780647,
                      0.74401, 0.93965, 0.935617, 0.838449, 1.0, 1,0, 0.0, 0.0, 0.0};


    weight *= IntWeight[binx-1];
  }
  else if(IsRange && !IsKF){

    // no weight
    double IntWeight[] = {0.0, 0.0, 0.559718, 0.528286, 0.422081, 0.487519, 0.666412, 0.544518, 0.627939, 0.604775, 0.780329,
                    0.621141, 0.934835, 0.934835, 0.836084, 1.0, 1.0, 0.0, 0.0, 0.0};
    // beam weight
    //double IntWeight[] = {0.0, 0.0, 0.616256, 0.51759, 0.375281, 0.435319, 0.50785, 0.531448, 0.597912, 0.536377, 0.797653,
    //                0.642434, 0.949061, 0.955738, 0.803652, 1.0, 1.0, 0.0, 0.0, 0.0};

    weight *= IntWeight[binx-1];

  }
  else{
    double IntWeight[] = {0.0, 0.369517, 0.409675, 0.415022, 0.373534, 0.348398, 0.459325, 0.422841, 0.481166, 0.426144, 0.680265,
                    0.633115, 0.933773, 0.909191,0.747696, 0.9244, 0.864349, 0.0, 1.0, 0.402144};
    weight *= IntWeight[binx-1];
  }*/
  weight = 1.;
  return weight;
}
// CosTheta
double CalCEXPi0CosThetaWeight(const int & binx, const bool IsRange, const bool IsKF){

  double weight = 1.;
  /*if(IsRange && IsKF){
    // beam wt
    //double IntWeight[] = {0.367657, 0.577443, 0.446119, 0.738655, 0.459862, 0.467202, 0.518055, 0.435467, 0.522927, 0.702523};
    //double IntWeight[] = {0.367657, 0.577443, 0.596655, 0.738655, 0.459862, 0.467202, 0.547708, 0.513704, 0.586615, 0.732868};
    // no wt
    //double IntWeight[] = {0.402573, 0.57405, 0.627507, 0.787512, 0.528985, 0.519493, 0.639224, 0.597047, 0.627748, 0.773835};
    // beam wt
    //double IntWeight[] = {0.411856, 0.622053, 0.640499, 0.772938, 0.506272, 0.513948, 0.593623, 0.560292, 0.631538, 0.768512};

    // pi0 wt 0410
    double IntWeight[] = {0.361952, 0.531519, 0.586469, 0.757284, 0.485984, 0.475453, 0.597117, 0.553845, 0.58418, 0.739109};
    // no weight
    //double IntWeight[] = {0.5, 0.666667, 0.714286, 0.846154, 0.625, 0.586207, 0.685714, 0.655738, 0.651786, 0.767717};

    weight *= IntWeight[binx-1];
  }
  else if(IsRange && !IsKF){
    // No KF only beam weight
    //double IntWeight[] = {0.509613, 0.709514, 0.524426, 0.804804, 0.586394, 0.536927, 0.596884, 0.549961, 0.625892, 0.720649};
    // No KF both weight
    double IntWeight[] = {0.517665, 0.583722, 0.599486, 0.868833, 0.679566, 0.503973, 0.573794, 0.566587, 0.621706, 0.706317};

    weight *= IntWeight[binx-1];

  }
  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517,
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }*/
  weight = 1.;
  return weight;
}
// Theta
double CalCEXPi0ThetaWeight(const int & binx, const bool IsRange, const bool IsKF){

  double weight = 1.;
  /*if(IsRange && IsKF){
    // beam wt
    //double IntWeight[] = {0.763616, 0.6884, 0.51305, 0.494907, 0.424725, 0.535193, 0.567364, 0.577443, 0.367657, 0.0};
    //double IntWeight[] = {0.783151, 0.72413, 0.578315, 0.551875, 0.44398, 0.535193, 0.632291, 0.577443, 0.367657, 0.0};
    // no wt
    //double IntWeight[] = {0.818232, 0.767555, 0.619857, 0.635357, 0.514164, 0.611244, 0.669043, 0.57405, 0.402573, 0.0};
    // beam wt
    //double IntWeight[] = {0.813937, 0.760516, 0.623557, 0.597683, 0.490585, 0.581027, 0.674375, 0.622053, 0.411856, 0.0};

    // pi0 wt 0410
    double IntWeight[] = {0.787775, 0.732404, 0.576011, 0.593378, 0.469896, 0.569643, 0.629882, 0.531519, 0.361952, 0.0};
    // no weight
    //double IntWeight[] = {0.8, 0.765363, 0.644628, 0.689189, 0.574468, 0.7, 0.75, 0.666667, 0.5, 0.0};

    weight *= IntWeight[binx-1];
  }
  else if(IsRange && !IsKF){
    // No KF only beam weight
    double IntWeight[] = {0.720196, 0.714559, 0.549956, 0.647322, 0.560627, 0.743115, 0.657209, 0.765056, 0.0};
    // No KF both weight
    //double IntWeight[] = {0.753768, 0.69482, 0.570662, 0.584574, 0.553642, 0.769505, 0.523964, 0.7708, 0.0};

    weight *= IntWeight[binx-1];

  }
  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517,
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }*/
  weight = 1.;
  return weight;
}

void WeightBeamData(TH1D * hbeam_data, const bool IsRange, const TString name){

  const Int_t x0 = hbeam_data->GetXaxis()->GetFirst();
  const Int_t x1 = hbeam_data->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double ientry = hbeam_data->GetBinContent(ix);
    double weight = -999;
    //double weight = 1;
    if(name.Contains("Ini")) weight = CalIniWeight(ix,IsRange);
    if(name.Contains("BeamInt")) weight = CalBeamIntWeight(ix,IsRange);
    if(name.Contains("CEXInt")) weight = CalCEXIntWeight(ix,IsRange);
    hbeam_data->SetBinContent(ix, ientry*weight);
  }

}

void WeightPi0Data(TH1D * hpi0_data, const bool IsRange, const bool IsKF,const TString name){

  const Int_t x0 = hpi0_data->GetXaxis()->GetFirst();
  const Int_t x1 = hpi0_data->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double ientry = hpi0_data->GetBinContent(ix);
    double weight = -999;
    //double weight = 1;
    if(name.Contains("KE")) weight = CalCEXPi0KEWeight(ix,IsRange,IsKF);
    if(name.Contains("CosTheta")) weight = CalCEXPi0CosThetaWeight(ix,IsRange,IsKF);
    if(name.Contains("Theta") && !name.Contains("Cos")) weight = CalCEXPi0ThetaWeight(ix,IsRange,IsKF);
    hpi0_data->SetBinContent(ix, ientry*weight);
  }

}

void DrawSystError(TH1D * hh){

  TH1D *hsyst = (TH1D*)hh->Clone("hsyst");
  const Int_t x0 = hsyst->GetXaxis()->GetFirst();
  const Int_t x1 = hsyst->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double syst[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 35.574,25.271,16.661,17.593,15.142,13.621,14.135,14.484,21.067,18.678};
    double stats = hsyst->GetBinError(ix);
    double systs = syst[ix-1];
    //hsyst->SetBinError(ix, sqrt(stats*stats + systs*systs));
    hsyst->SetBinError(ix, sqrt(stats*stats + systs*systs));

  }
  //hsyst->SetMarkerColor(kRed);
  //hsyst->SetLineColor(kRed);
  //hsyst->GetXaxis()->SetRangeUser(450,950);
  hsyst->Draw("E1 sames");
  //gStyle->SetErrorX(0);
}


TH1D * RebinThetaHist(TH1D * hist){

  Int_t nbins = hist->GetNbinsX();
  Double_t xbins[nbins+1];
  hist->GetXaxis()->GetLowEdge(xbins);
  xbins[nbins] = hist->GetXaxis()->GetBinUpEdge(nbins);

  Double_t newbins[nbins+1];
  for (Int_t i=0; i<=nbins; i++) {
      if (i < nbins-3) {
          newbins[i] = xbins[i]; // leave first nbins-4 bins unchanged
      } else {
          newbins[i] = xbins[nbins]; // set remaining bins to the same upper edge as last bin
      }
  }

  TH1D *rebinHist = new TH1D(hist->GetName(), hist->GetTitle(), nbins, xbins);
  for (Int_t i=1; i<=nbins; i++) {
      rebinHist->SetBinContent(i, hist->GetBinContent(i));
      rebinHist->SetBinError(i, hist->GetBinError(i));
  }
  rebinHist = (TH1D*)rebinHist->Rebin(nbins, "", newbins);

  for (Int_t i=1; i<=nbins-3; i++) {
      if(i == nbins-3) rebinHist->SetBinContent(i, rebinHist->GetBinContent(i)/4);
      if(i == nbins-3) rebinHist->SetBinError(i, rebinHist->GetBinError(i)/2);
  }

  return rebinHist;
}
