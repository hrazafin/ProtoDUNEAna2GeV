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

void DrawHistOutput(TH1D *h_truth, TH1D *h_reco, TH1D *h_unfold, TH1D *h_unfold_data, const bool &rebin);
TH1D * GetIncidentHist(TH1D * InitialHist, TH1D * InteractingHist);
void DrawTotalXSPlots(TH1D *hinc_truth, TH1D *hint_truth, TH1D *hinc_reco, TH1D *hint_reco, TH1D *hinc_data, TH1D *hint_data, TGraph* g_cex, double & totalCEXXS_800_truth, double & totalCEXXS_800_reco, double & totalCEXXS_800_data);
void DrawDiffXSPlots(TH1D *hpi0KE_truth, TH1D *hpi0KE_reco, TH1D *hpi0KE_data, TGraph* g_775, double Int_800_truth, double Int_800err_truth, double Int_800_reco, double Int_800err_reco, double Int_800_data, double Int_800err_data, double totalCEXXS_800_truth, double totalCEXXS_800_reco, double totalCEXXS_800_data, TString name);
TH1D * TotalCEXXSCal(TH1D * hh, TH1D * InteractingHist, const bool & Eslice, const bool & widerBin, const bool & newMethod);
TH1D * DiffCEXXSCal(TH1D * DiffCEXInteractingHist, const double &diffInt,  const double &diffInterror, const double &scale_bin, const double &scale_totalXS);
void SetTitleFormat(TH1 * hh);
void DrawDataMCRatio(TH1D * hratio, bool xsec, const TString tag);
void WeightPi0Data(TH1D * hpi0KE_data, const bool IsRange, const TString name);
double CalCEXPi0KEWeight(const int & binx, const bool IsRange);
double CalCEXPi0CosThetaWeight(const int & binx, const bool IsRange);
double CalCEXPi0ThetaWeight(const int & binx, const bool IsRange);

void DrawSystError(TH1D * hh);


const bool IsRange = true;

int main(int argc, char * argv[])
{ 
  gStyle->SetOptStat(0);

  //const TString finName = "input/outana.root";
  //const TString finName = "input/outana_dataNobckSub.root";
  TString finName = "input/outana_dataNobckSub.root";
  //if(IsRange) finName = "input/outana_wholeRange.root";
  if(IsRange) finName = "input/outana_wholeRange_0123.root";

  TFile *file = TFile::Open(finName);

  if(!file->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // ========== Initial Histogram ========== //
  RooUnfoldResponse *response_SliceID_Ini = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Ini");
  TH1D *hini_truth = (TH1D*) file->Get("mc/i000hTruthBeamInitialHist"); // Truth
  TH1D *hini = (TH1D*) file->Get("mc/i016hRecoInitialHist"); // Fake data
  TH1D *hini_data = (TH1D*) file->Get("mc/i019hRecoBeamInitialHistData"); // Data

  RooUnfoldBayes unfold_Ini (response_SliceID_Ini, hini, 4); // Unfold fake data
  RooUnfoldBayes unfold_Ini_Data (response_SliceID_Ini, hini_data, 4); // Unfold data

  TH1D *hini_uf = (TH1D*)unfold_Ini.Hreco();
  TH1D *hini_uf_data = (TH1D*)unfold_Ini_Data.Hreco();

  // Draw Hist using tmp hist to rebin
  TH1D *hini_truth_tmp = (TH1D*)hini_truth->Clone();
  TH1D *hini_tmp = (TH1D*)hini->Clone();
  TH1D *hini_uf_tmp = (TH1D*)hini_uf->Clone();
  TH1D *hini_uf_data_tmp = (TH1D*)hini_uf_data->Clone();

  DrawHistOutput(hini_truth_tmp, hini_tmp, hini_uf_tmp, hini_uf_data_tmp, true); // rebin == true

  // ========== BeamInt Histogram ========== //
  RooUnfoldResponse *response_SliceID_BeamInt = (RooUnfoldResponse*)file->Get("mc/response_SliceID_BeamInt");
  TH1D *hbeamint_truth = (TH1D*) file->Get("mc/i002hTruthBeamInteractingHist");
  TH1D *hbeamint = (TH1D*) file->Get("mc/i016hRecoBeamInteractingHist");
  TH1D *hbeamint_data = (TH1D*) file->Get("mc/i020hRecoBeamInteractingHistData");

  RooUnfoldBayes unfold_BeamInt (response_SliceID_BeamInt, hbeamint, 4);
  RooUnfoldBayes unfold_BeamInt_Data (response_SliceID_BeamInt, hbeamint_data, 4);

  TH1D *hbeamint_uf = (TH1D*)unfold_BeamInt.Hreco();
  TH1D *hbeamint_uf_data = (TH1D*)unfold_BeamInt_Data.Hreco();

  // Draw Hist
  TH1D *hbeamint_truth_tmp = (TH1D*)hbeamint_truth->Clone();
  TH1D *hbeamint_tmp = (TH1D*)hbeamint->Clone();
  TH1D *hbeamint_uf_tmp = (TH1D*)hbeamint_uf->Clone();
  TH1D *hbeamint_uf_data_tmp = (TH1D*)hbeamint_uf_data->Clone();

  DrawHistOutput(hbeamint_truth_tmp, hbeamint_tmp, hbeamint_uf_tmp, hbeamint_uf_data_tmp, true); // rebin == true

  // ========== Int Histogram ========== //
  RooUnfoldResponse *response_SliceID_Int = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Int");
  TH1D *hint_truth = (TH1D*) file->Get("mc/i004hTruthCEXInteractingHist");
  TH1D *hint = (TH1D*) file->Get("mc/i018hRecoInteractingHist");
  TH1D *hint_data = (TH1D*) file->Get("mc/i018hRecoInteractingHistData");

  RooUnfoldBayes unfold_Int (response_SliceID_Int, hint, 4);
  RooUnfoldBayes unfold_Int_Data (response_SliceID_Int, hint_data, 4);

  TH1D *hint_uf = (TH1D*)unfold_Int.Hreco();
  TH1D *hint_uf_data = (TH1D*)unfold_Int_Data.Hreco();
  hint_uf->Scale(0.5); // Need to scale down the unfolded fake data since used half truth MC
  //hint_uf_data->Scale(0.5);

  DrawHistOutput(hint_truth, hint, hint_uf, hint_uf_data, false); // rebin == false

  // ======= GEANT4 Input XS ======== //
  const TString xcesfinName = "../Analysis/input/LAr_PiPlus_exclusive_cross_section.root";
  TFile *f_CrossSection = TFile::Open(xcesfinName);
  if(!f_CrossSection->IsOpen()){
    cout << "xsce file not open" << endl;
    exit(1);
  }

  TGraph *g_cex;

  f_CrossSection->GetObject("cex_KE",g_cex);
  

  const TString DiffxcesfinName = "../Analysis/input/cross_section_out_withNewVar.root";
  TFile *f_DiffCrossSection = TFile::Open(DiffxcesfinName);
  if(!f_DiffCrossSection->IsOpen()) {
    cout << "Diffxsce file not open" << endl;
    exit(1);
  }

  TH1D * g_cex_775MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dKEpi0775_MeV");
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
  hinc_truth->Rebin(50); hint_truth->Rebin(50);
  hinc_uf->Rebin(50); hint_uf->Rebin(50);
  hinc_uf_data->Rebin(50); hint_uf_data->Rebin(50);

  // ======== Draw Total CEX XS plots ======== //
  double totalCEXXS_800_truth = 0, totalCEXXS_800_reco = 0, totalCEXXS_800_data = 0;
  DrawTotalXSPlots(hinc_truth,hint_truth,hinc_uf,hint_uf,hinc_uf_data,hint_uf_data,g_cex,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data);

  // ======== Pi0 KE Histogram (KE PI+ = 650to800 MeV)======= //
  RooUnfoldResponse *response_SliceID_Pi0KE = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Pi0KE");
  TH1D *hpi0KE_truth = (TH1D*) file->Get("mc/i006hTruthDiffCEXInteractingHist_800MeV");
  if(IsRange) hpi0KE_truth = (TH1D*) file->Get("mc/i006hTruthDiffCEXInteractingHist_650to800MeV");
  TH1D *hpi0KE = (TH1D*) file->Get("mc/i030hRecoPi0KEHist");
  TH1D *hpi0KE_data = (TH1D*) file->Get("mc/i030hRecoPi0KEHistData");

  // Data weight
  WeightPi0Data(hpi0KE_data, IsRange, "KE");

  RooUnfoldBayes unfold_Pi0KE (response_SliceID_Pi0KE, hpi0KE, 4);
  RooUnfoldBayes unfold_Pi0KE_Data (response_SliceID_Pi0KE, hpi0KE_data, 4);

  TH1D *hpi0KE_uf = (TH1D*)unfold_Pi0KE.Hreco();
  TH1D *hpi0KE_uf_data = (TH1D*)unfold_Pi0KE_Data.Hreco();
  hpi0KE_uf->Scale(0.5);
  //hpi0KE_uf_data->Scale(0.5);

  DrawHistOutput(hpi0KE_truth, hpi0KE, hpi0KE_uf, hpi0KE_uf_data, false); // rebin == false

  double Int_800_truth = hint_truth->GetBinContent(16);
  if(IsRange) Int_800_truth = hint_truth->GetBinContent(14) + hint_truth->GetBinContent(15) + hint_truth->GetBinContent(16);
  double Int_800err_truth = hint_truth->GetBinError(16);

  double Int_800_reco = hint_uf->GetBinContent(16);
  if(IsRange) Int_800_reco = hint_uf->GetBinContent(14) + hint_uf->GetBinContent(15) + hint_uf->GetBinContent(16);
  double Int_800err_reco = hint_uf->GetBinError(16);

  double Int_800_data = hint_uf_data->GetBinContent(16);
  if(IsRange) Int_800_data = hint_uf_data->GetBinContent(14) + hint_uf_data->GetBinContent(15) + hint_uf_data->GetBinContent(16);
  double Int_800err_data = hint_uf_data->GetBinError(16);

  if(!IsRange) DrawDiffXSPlots(hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_775,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, "775");

  // ======= Plot Pi0KE ======== //
  const TString finNameNew = "input/outana_600to800.root";
  TFile *fileNew = TFile::Open(finNameNew);

  const TString finNameXS = "input/cross_section_out.root";
  TFile *fileXS = TFile::Open(finNameXS);

  TH1D * g_cex_650to800MeV = (TH1D*)fileXS->Get("inel_cex_1dKEpi0650to800_MeV");
  const Int_t x0_diff650to800 = g_cex_650to800MeV->GetXaxis()->GetFirst();
  const Int_t x1_diff650to800 = g_cex_650to800MeV->GetXaxis()->GetLast();
  double x_650to800[x1_diff650to800], y_650to800[x1_diff650to800];
  for(Int_t ix=x0_diff650to800; ix<=x1_diff650to800; ix++){
    //cout << "Bin 650to800: "  << ix << " " << g_cex_650to800MeV->GetBinContent(ix) << endl;
    x_650to800[ix] = g_cex_650to800MeV->GetBinCenter(ix);
    y_650to800[ix] = g_cex_650to800MeV->GetBinContent(ix);
  }

  TGraph *g_650to800 = new TGraph(40,x_650to800,y_650to800);

  if(IsRange) DrawDiffXSPlots(hpi0KE_truth,hpi0KE_uf,hpi0KE_uf_data,g_650to800,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, "650to800");


  // ======== Pi0 CosTheta Histogram (KE PI+ = 650to800 MeV)======= //
  RooUnfoldResponse *response_SliceID_Pi0CosTheta = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Pi0CosTheta");
  TH1D *hpi0CosTheta_truth = (TH1D*) file->Get("mc/i666hTruthDiffCEXInteractingHistCosTheta_650to800MeV");
  TH1D *hpi0CosTheta = (TH1D*) file->Get("mc/i032hRecoPi0CosThetaHist");
  TH1D *hpi0CosTheta_data = (TH1D*) file->Get("mc/i032hRecoPi0CosThetaHistData");

  // Data weight
  WeightPi0Data(hpi0CosTheta_data, IsRange, "CosTheta");

  RooUnfoldBayes unfold_Pi0CosTheta (response_SliceID_Pi0CosTheta, hpi0CosTheta, 4);
  RooUnfoldBayes unfold_Pi0CosTheta_Data (response_SliceID_Pi0CosTheta, hpi0CosTheta_data, 4);

  TH1D *hpi0CosTheta_uf = (TH1D*)unfold_Pi0CosTheta.Hreco();
  TH1D *hpi0CosTheta_uf_data = (TH1D*)unfold_Pi0CosTheta_Data.Hreco();
  hpi0CosTheta_uf->Scale(0.5);
  //hpi0CosTheta_uf_data->Scale(0.5);

  DrawHistOutput(hpi0CosTheta_truth, hpi0CosTheta, hpi0CosTheta_uf, hpi0CosTheta_uf_data, false); // rebin == false

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

  if(IsRange) DrawDiffXSPlots(hpi0CosTheta_truth,hpi0CosTheta_uf,hpi0CosTheta_uf_data,g_650to800_CosTheta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, "650to800CosTheta");


  // ======== Pi0 Theta Histogram (KE PI+ = 650to800 MeV)======= //
  RooUnfoldResponse *response_SliceID_Pi0Theta = (RooUnfoldResponse*)file->Get("mc/response_SliceID_Pi0Theta");
  TH1D *hpi0Theta_truth = (TH1D*) file->Get("mc/i066hTruthDiffCEXInteractingHistTheta_650to800MeV");
  TH1D *hpi0Theta = (TH1D*) file->Get("mc/i034hRecoPi0ThetaHist");
  TH1D *hpi0Theta_data = (TH1D*) file->Get("mc/i034hRecoPi0ThetaHistData");

  // Data weight
  WeightPi0Data(hpi0Theta_data, IsRange, "Theta");

  RooUnfoldBayes unfold_Pi0Theta (response_SliceID_Pi0Theta, hpi0Theta, 4);
  RooUnfoldBayes unfold_Pi0Theta_Data (response_SliceID_Pi0Theta, hpi0Theta_data, 4);

  TH1D *hpi0Theta_uf = (TH1D*)unfold_Pi0Theta.Hreco();
  TH1D *hpi0Theta_uf_data = (TH1D*)unfold_Pi0Theta_Data.Hreco();
  hpi0Theta_uf->Scale(0.5);
  //hpi0Theta_uf_data->Scale(0.5);

  DrawHistOutput(hpi0Theta_truth, hpi0Theta, hpi0Theta_uf, hpi0Theta_uf_data, false); // rebin == false

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

  if(IsRange) DrawDiffXSPlots(hpi0Theta_truth,hpi0Theta_uf,hpi0Theta_uf_data,g_650to800_Theta,Int_800_truth,Int_800err_truth, Int_800_reco, Int_800err_reco, Int_800_data, Int_800err_data, totalCEXXS_800_truth, totalCEXXS_800_reco, totalCEXXS_800_data, "650to800Theta");

  

  // ======= Seperate Truth Study ======== //
  TH1D *hpi0KETotal_truth = (TH1D*) fileNew->Get("mc/i006hTruthDiffCEXInteractingHist_650to800MeV");

  double scale_bin = 50.0;
  double scale_totalXS_truth = totalCEXXS_800_truth;

  double Int_total_truth = hint_truth->GetBinContent(14) + hint_truth->GetBinContent(15) + hint_truth->GetBinContent(16);
  double Int_totalerr_truth = hint_truth->GetBinError(16);

  TH1D * hdiffcex_truth = DiffCEXXSCal(hpi0KETotal_truth,Int_total_truth,Int_totalerr_truth, scale_bin, scale_totalXS_truth);
  
  TCanvas * ctmp = new TCanvas("ctmp", "", 1200, 800);

  hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/dT_{#pi^{0}} (mb/MeV)");
  SetTitleFormat(hdiffcex_truth);

  hdiffcex_truth->SetMarkerStyle(8);
  hdiffcex_truth->SetMarkerSize(1);
  hdiffcex_truth->SetMarkerColor(kBlue);
  hdiffcex_truth->SetLineColor(kBlue);
  hdiffcex_truth->SetLineWidth(1);
  hdiffcex_truth->SetMaximum(0.6);
  hdiffcex_truth->SetMinimum(0.);

  hdiffcex_truth->Draw("E1 sames");

  g_650to800->SetLineColor(kRed);
  g_650to800->Draw("sames C");

  ctmp->Print("output/truth_Diff600to800CEXXS.png");


  
  // Close files
  f_CrossSection->Close();
  f_DiffCrossSection->Close();
  file->Close();
  fileNew->Close();
  fileXS->Close();

}

void DrawHistOutput(TH1D *h_truth, TH1D *h_reco, TH1D *h_unfold, TH1D *h_unfold_data, const bool &rebin){
  // Rebin the histogram
  if(rebin){
    h_truth->Rebin(50); h_reco->Rebin(50); h_unfold->Rebin(50); h_unfold_data->Rebin(50); 
  }
  TString tag = h_truth->GetName();

  TLatex tt;
  tt.SetNDC();

  TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 1200, 800);
  
  SetTitleFormat(h_truth);

  h_truth->SetLineColor(kRed);
  h_truth->SetMaximum(h_truth->GetMaximum()*1.5);
  h_truth->SetMinimum(0.0);
  h_truth->GetYaxis()->SetTitle("Candidates");
  h_truth->GetXaxis()->SetTitle("T_{#pi^{+}} (MeV)");
  if(tag.Contains("i006")) h_truth->GetXaxis()->SetTitle("T_{#pi^{0}} (MeV)");
  if(tag.Contains("i066")) h_truth->GetXaxis()->SetTitle("#theta_{#pi^{0}} (deg.)");
  if(tag.Contains("i666")) h_truth->GetXaxis()->SetTitle("cos#theta_{#pi^{0}}");

  h_truth->Draw("hist");

  h_reco->SetMarkerStyle(8);
  h_reco->SetMarkerSize(1);
  h_reco->SetMarkerColor(kBlue);
  h_reco->SetLineColor(kBlue);
  h_reco->SetLineWidth(1);
  h_reco->Draw("e1 sames");

  h_unfold->SetMarkerStyle(8);
  h_unfold->SetMarkerSize(1);
  h_unfold->SetMarkerColor(kGreen+3);
  h_unfold->SetLineColor(kGreen+3);
  h_unfold->SetLineWidth(1);
  h_unfold->Draw("e1 sames");

  h_unfold_data->SetMarkerStyle(8);
  h_unfold_data->SetMarkerSize(1);
  h_unfold_data->SetMarkerColor(kBlack);
  h_unfold_data->SetLineColor(kBlack);
  h_unfold_data->SetLineWidth(1);
  h_unfold_data->Draw("e1 sames");

  TLegend * legend = new TLegend(0.15, 0.48, 0.48, 0.88);
  if(tag.Contains("650to800") && !tag.Contains("CosTheta")) legend = new TLegend(0.55, 0.48, 0.88, 0.88);
  legend->AddEntry(h_truth, "MC Truth", "l");
  legend->AddEntry(h_reco, "MC Reco. - BckSub", "lep");
  legend->AddEntry(h_unfold, "MC Reco - Unfolded", "lep");
  legend->AddEntry(h_unfold_data, "Data - Unfolded", "lep");
  legend->SetBorderSize(0);
  legend->Draw("same");

  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  if(tag.Contains("i006")||tag.Contains("i066")) tt.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
  else tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
  if(tag.Contains("i006")||tag.Contains("i066")) tt.DrawLatex(0.22,0.825,"#bf{#it{Preliminary}}");
  else tt.DrawLatex(0.725,0.825,"#bf{#it{Preliminary}}");

  //c1->Print("output/hint_uf.png");
  c1->Print("output/"+tag+".png");

}

void DrawTotalXSPlots(TH1D *hinc_truth, TH1D *hint_truth, TH1D *hinc_reco, TH1D *hint_reco, TH1D *hinc_data, TH1D *hint_data, TGraph* g_cex, double & totalCEXXS_800_truth, double & totalCEXXS_800_reco, double & totalCEXXS_800_data){

  TLatex tt;
  tt.SetNDC();

  TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
  
  // Truth
  TH1D * hcex_truth = TotalCEXXSCal(hinc_truth,hint_truth,true,false,false);
  // Reco
  TH1D * hcex_reco = TotalCEXXSCal(hinc_reco,hint_reco,true,false,false);
  // Data
  TH1D * hcex_data = TotalCEXXSCal(hinc_data,hint_data,true,false,false);

  totalCEXXS_800_truth = hcex_truth->GetBinContent(16);
  totalCEXXS_800_reco = hcex_reco->GetBinContent(16);
  totalCEXXS_800_data = hcex_data->GetBinContent(16);


  SetTitleFormat(hcex_truth);
  SetTitleFormat(hcex_reco);

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

  hcex_truth->Draw("E1 sames");

  hcex_reco->SetMarkerStyle(23);
  hcex_reco->SetMarkerSize(1.5);
  hcex_reco->SetMarkerColor(kGreen+3);
  hcex_reco->SetLineColor(kGreen+3);
  hcex_reco->SetLineWidth(1);
  //hcex_reco->SetMaximum(0.4);
  hcex_reco->SetMinimum(0.);

  hcex_reco->Draw("E1 sames");


  hcex_data->SetMarkerStyle(8);
  hcex_data->SetMarkerSize(1);
  hcex_data->SetMarkerColor(kBlack);
  hcex_data->SetLineColor(kBlack);
  hcex_data->SetLineWidth(1);
  //hcex_data->SetMaximum(0.4);
  hcex_data->SetMinimum(0.);
  DrawSystError(hcex_data); // it works but not done

  hcex_data->Draw("E1 sames");

  g_cex->SetLineColor(kRed);
  g_cex->Draw("sames C");

  TLegend * legend = new TLegend(0.15, 0.48, 0.48, 0.88);
  legend->AddEntry(g_cex, "Geant4 Prediction", "l");
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
  
  const Int_t x0_x = hcex_data->GetXaxis()->GetFirst();
  const Int_t x1_x = hcex_data->GetXaxis()->GetLast();

  hratio->Scale(0);

  for(Int_t ix=x0_x; ix<=x1_x; ix++){
    double bin_center = hcex_data->GetBinCenter(ix);
    double MC = g_cex->Eval(bin_center);
    
    double data = hcex_data->GetBinContent(ix);
    double edata = hcex_data->GetBinError(ix);
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
  tt.DrawLatex(0.705,0.775,"#bf{#it{(Stats. + Syst. Error)}}");

  c1->Print("output/truth_TotalCEXXS.png");

  
}


void DrawDiffXSPlots(TH1D *hpi0KE_truth, TH1D *hpi0KE_reco, TH1D *hpi0KE_data, TGraph* g_775, double Int_800_truth, double Int_800err_truth, double Int_800_reco, double Int_800err_reco, double Int_800_data, double Int_800err_data, double totalCEXXS_800_truth, double totalCEXXS_800_reco, double totalCEXXS_800_data, TString name){

  TLatex tt;
  tt.SetNDC();

  TCanvas * c2 = new TCanvas("c2", "", 1200, 800);
  
  double scale_bin = 50.0;
  if(name.Contains("CosTheta")) scale_bin = 0.2;
  if(name.Contains("Theta") && !name.Contains("Cos")) scale_bin = 20.0;

  double scale_totalXS_truth = totalCEXXS_800_truth;
  double scale_totalXS_reco = totalCEXXS_800_reco;
  double scale_totalXS_data = totalCEXXS_800_data;

  TH1D * hdiffcex_truth = DiffCEXXSCal(hpi0KE_truth,Int_800_truth,Int_800err_truth, scale_bin, scale_totalXS_truth);
  TH1D * hdiffcex_reco = DiffCEXXSCal(hpi0KE_reco,Int_800_reco,Int_800err_reco, scale_bin, scale_totalXS_reco);
  TH1D * hdiffcex_data = DiffCEXXSCal(hpi0KE_data,Int_800_data,Int_800err_data, scale_bin, scale_totalXS_data);

  hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/dT_{#pi^{0}} (mb/MeV)");
  if(name.Contains("CosTheta")) hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi^{0}} (mb)");
  if(name.Contains("Theta") && !name.Contains("Cos")) hdiffcex_truth->GetYaxis()->SetTitle("d#sigma/d#theta_{#pi^{0}} (mb/deg.)");
  SetTitleFormat(hdiffcex_truth);
  SetTitleFormat(hdiffcex_reco);

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

  hdiffcex_truth->Draw("E1 sames");

  hdiffcex_reco->SetMarkerStyle(23);
  hdiffcex_reco->SetMarkerSize(1.5);
  hdiffcex_reco->SetMarkerColor(kGreen+3);
  hdiffcex_reco->SetLineColor(kGreen+3);
  hdiffcex_reco->SetLineWidth(1);
  hdiffcex_reco->SetMaximum(0.4);
  hdiffcex_reco->SetMinimum(0.);

  hdiffcex_reco->Draw("E1 sames");

  hdiffcex_data->SetMarkerStyle(8);
  hdiffcex_data->SetMarkerSize(1);
  hdiffcex_data->SetMarkerColor(kBlack);
  hdiffcex_data->SetLineColor(kBlack);
  hdiffcex_data->SetLineWidth(1);
  hdiffcex_data->SetMaximum(0.4);
  hdiffcex_data->SetMinimum(0.);

  hdiffcex_data->Draw("E1 sames");

  g_775->SetLineColor(kRed);
  g_775->Draw("sames C");

  TLegend * legend = new TLegend(0.55, 0.48, 0.88, 0.88);
  if(name.Contains("CosTheta")) legend = new TLegend(0.15, 0.48, 0.48, 0.88);

  legend->AddEntry(g_775, "Geant4 Prediction", "l");
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
  
  const Int_t x0_x = hdiffcex_data->GetXaxis()->GetFirst();
  const Int_t x1_x = hdiffcex_data->GetXaxis()->GetLast();

  hratio->Scale(0);

  for(Int_t ix=x0_x; ix<=x1_x; ix++){
    double bin_center = hdiffcex_data->GetBinCenter(ix);
    double MC = g_775->Eval(bin_center);
    
    double data = hdiffcex_data->GetBinContent(ix);
    double edata = hdiffcex_data->GetBinError(ix);
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

  c2->Print("output/plot_Diff"+name+"CEXXS.png");

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



TH1D * TotalCEXXSCal(TH1D * hh, TH1D * InteractingHist, const bool & Eslice, const bool & widerBin, const bool & newMethod)
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
  for(Int_t ix=x0; ix<=x1; ix++){

    // Read off the entry for each bin of the two histograms
    double incE = hh->GetBinContent(ix);
    double intE = InteractingHist->GetBinContent(ix);
    // Calculate the ratio in the full log form
    double ratio = 0;
    if(incE > intE) ratio = log(incE/(incE-intE));//*meandEdx[ix-1]*sigma_factor*1e27;

    if(Eslice) ratio = ratio * meandEdx[ix-1]*sigma_factor*1e27;

    //if(incE != 0) cout << "ix: "<< ix << "incE: " << incE << "intE: " << intE << " ratio: " << ratio << "meandEdx[ix]: " << meandEdx[ix-1] << endl;

    // Simple ratio form 
    //double ratio = intE/incE;

    // If the incE entry is not zero set the bin content
    if(incE != 0) xsec->SetBinContent(ix,ratio);
    if( Eslice && ratio > 150) xsec->SetBinContent(ix,9999);

    // Error propagation method 1
    //double error = sqrt(intE+pow(intE,2)/incE)/incE;
    // Error propagation method 2
    double einc = hh->GetBinError(ix);
    double eint = InteractingHist->GetBinError(ix);
    double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)));
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
  return xsec;
}

TH1D * DiffCEXXSCal(TH1D * DiffCEXInteractingHist, const double &diffInt,  const double &diffInterror, const double &scale_bin, const double &scale_totalXS)
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
      double ratio = intE/incE;
      DiffCEXxsec->SetBinContent(ix,ratio);
      //double error = sqrt(intE+pow(intE,2)/incE)/incE;
      double einc = diffInterror;
      double eint = DiffCEXInteractingHist->GetBinError(ix);
      double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)));
      DiffCEXxsec->SetBinError(ix,error);
      //cout << "error: " << error << endl;
      //cout << "ratio: " << ratio << endl;

    }
  }
  // Normalise to total cross section at that slice
  DiffCEXxsec->Scale(1/scale_bin);
  DiffCEXxsec->Scale(scale_totalXS);
  
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

double CalCEXPi0KEWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  if(IsRange){
    double IntWeight[] = {0.0, 0.0, 0.652846, 0.596701, 0.521758, 0.579765, 0.689989, 0.544885, 0.607432, 0.591426, 0.837472, 
                      0.624503, 0.944532, 0.940689, 0.778484, 1.0, 1,0, 0.0, 0.0, 0.0};
    //double IntWeight[] = {0.0, 0.0, 0.652846, 0.596701, 0.521758, 0.579765, 0.689989, 0.544885, 0.607432, 0.591426, 0.837472, 
    //                  0.624503, 0.944532, 1.340689, 1.778484, 1.3, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}
// CosTheta
double CalCEXPi0CosThetaWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  if(IsRange){
    double IntWeight[] = {0.509613, 0.709514, 0.665169, 1, 0.633065, 0.581289, 0.671715, 0.554388, 0.592632, 0.727953};
    
    weight *= IntWeight[binx-1];
  }
  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}
// Theta
double CalCEXPi0ThetaWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  if(IsRange){
    double IntWeight[] = {0.712032, 0.72855, 0.540438, 0.680656, 0.574381, 0.827391, 0.657209, 0.765056, 0.0};
    
    weight *= IntWeight[binx-1];
  }
  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}

void WeightPi0Data(TH1D * hpi0_data, const bool IsRange, const TString name){

  const Int_t x0 = hpi0_data->GetXaxis()->GetFirst();
  const Int_t x1 = hpi0_data->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double ientry = hpi0_data->GetBinContent(ix);
    double weight = -999;
    if(name.Contains("KE")) weight = CalCEXPi0KEWeight(ix,IsRange);
    if(name.Contains("CosTheta")) weight = CalCEXPi0CosThetaWeight(ix,IsRange);
    if(name.Contains("Theta") && !name.Contains("Cos")) weight = CalCEXPi0ThetaWeight(ix,IsRange);
    hpi0_data->SetBinContent(ix, ientry*weight);
  }

}

void DrawSystError(TH1D * hh){

  TH1D *hsyst = (TH1D*)hh->Clone("hsyst");
  const Int_t x0 = hsyst->GetXaxis()->GetFirst();
  const Int_t x1 = hsyst->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double syst[] = {0, 0, 0, 0, 0, 0, 0, 0, 4.69, 4.69, 2.37, 1.51, 0.67, 0.34, 0.28, 0.45, 0.44, 4.88};
    double stats = hsyst->GetBinError(ix);
    double systs = syst[ix-1];
    //hsyst->SetBinError(ix, sqrt(stats*stats + systs*systs));
    hsyst->SetBinError(ix, stats + systs);

  }
  //hsyst->SetMarkerColor(kRed);
  //hsyst->SetLineColor(kRed);
  hsyst->Draw("E1 sames");
  //gStyle->SetErrorX(0);
}

