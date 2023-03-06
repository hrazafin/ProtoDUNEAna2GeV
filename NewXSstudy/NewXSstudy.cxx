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
TH1D * GetNewIncidentHist(TH1D * InitialHist, TH1D * InteractingHist);
void DrawTotalXSPlots(TH1D *hinc_truth, TH1D *hint_truth, TH1D *hinc_reco, TH1D *hint_reco, TH1D *hinc_data, TH1D *hint_data, TGraph* g_cex, double & totalCEXXS_800_truth, double & totalCEXXS_800_reco, double & totalCEXXS_800_data, const bool kTotal);
void DrawDiffXSPlots(TH1D *hpi0KE_truth, TH1D *hpi0KE_reco, TH1D *hpi0KE_data, TGraph* g_775, double Int_800_truth, double Int_800err_truth, double Int_800_reco, double Int_800err_reco, double Int_800_data, double Int_800err_data, double totalCEXXS_800_truth, double totalCEXXS_800_reco, double totalCEXXS_800_data, TString name);
TH1D * TotalCEXXSCal(TH1D * hh, TH1D * InteractingHist, const bool & Eslice, const bool & widerBin, const bool & newMethod, const bool kTotal);
TH1D * DiffCEXXSCal(TH1D * DiffCEXInteractingHist, const double &diffInt,  const double &diffInterror, const double &scale_bin, const double &scale_totalXS);
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
void DrawSystError(TH1D * hh);
double Density_Correction(double beta, double gamma);
double dEdx_Bethe_Bloch(double KE, double mass);


double mass_muon = 105.658; // [MeV]
double mass_pion = 139.57; // [MeV]
double mass_proton = 938.272; // [MeV]

const double K = 0.307075; // [MeV cm2 / mol]
const double I = 188.0e-6; // [MeV], mean excitation energy
const double Me = 0.511; // [Mev], mass of electron
// == Parameters for the density correction
const double density_C = 5.2146;
const double density_y0 = 0.2;
const double density_y1 = 3.0;
const double density_a = 0.19559;
const double density_k = 3.0;

const double LAr_density = 1.39; // [g/cm3]
const double N_A = 6.02; // [10^23 / mol]
const double M_Ar = 39.948; // [g / mol]
const double xsec_unit = M_Ar / (LAr_density * N_A); // [10^4 mb cm]

const bool IsRange = true;
const bool IsKF = false;

int main(int argc, char * argv[])
{ 
  gStyle->SetOptStat(0);

  //const TString finName = "input/outana.root";
  const TString finName = "input/ori_nokefflim.root";

  TFile *file = TFile::Open(finName);

  if(!file->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

 
  // ======= GEANT4 Input XS ======== //
  const TString xcesfinName = "../Analysis/input/LAr_PiPlus_exclusive_cross_section.root";
  TFile *f_CrossSection = TFile::Open(xcesfinName);
  if(!f_CrossSection->IsOpen()){
    cout << "xsce file not open" << endl;
    exit(1);
  }

  TGraph *g_cex;
  TGraph *g_totalinel;

  f_CrossSection->GetObject("cex_KE",g_cex);
  f_CrossSection->GetObject("total_inel_KE",g_totalinel);

  // ============= Sungbin’s new method ============= //
  TH1D *hNewini_truth = (TH1D*) file->Get("mc/i000hTruthBeamInitialHist");
  TH1D *hNewbeamint_truth = (TH1D*) file->Get("mc/i002hTruthBeamInteractingHist");
  TH1D *hNewint_truth = (TH1D*) file->Get("mc/i002hTruthInteractingHist");
  TH1D *hNewCEXint_truth = (TH1D*) file->Get("mc/i004hTruthCEXInteractingHist");
  //TH1D *hNewinc_truth_dir = (TH1D*) file->Get("mc/i000hNewTruthBeamIncidentHist");


  TH1D *hbeamint_truth_tmp = (TH1D*)hNewini_truth->Clone();
  

  // ======== Calculate incident histogram ======== //
  // Truth
  TH1D * hNewinc_truth = GetIncidentHist(hNewini_truth,hNewbeamint_truth);

  hNewinc_truth->Rebin(50);

  hNewini_truth->Rebin(50); hNewbeamint_truth->Rebin(50);

  DrawHistOutput(hNewini_truth, hbeamint_truth_tmp, hbeamint_truth_tmp, hbeamint_truth_tmp, false); // rebin == false
  DrawHistOutput(hNewbeamint_truth, hbeamint_truth_tmp, hbeamint_truth_tmp, hbeamint_truth_tmp, false); // rebin == false
  DrawHistOutput(hNewint_truth, hbeamint_truth_tmp, hbeamint_truth_tmp, hbeamint_truth_tmp, false); // rebin == false
  DrawHistOutput(hNewCEXint_truth, hbeamint_truth_tmp, hbeamint_truth_tmp, hbeamint_truth_tmp, false); // rebin == false
  //DrawHistOutput(hNewinc_truth_dir, hbeamint_truth_tmp, hbeamint_truth_tmp, hbeamint_truth_tmp, false); // rebin == false

  DrawHistOutput(hNewinc_truth, hbeamint_truth_tmp, hbeamint_truth_tmp, hbeamint_truth_tmp, false); // rebin == false

  
  double totalCEXXS_800_truth = 0, totalCEXXS_800_reco = 0, totalCEXXS_800_data = 0;
  
  DrawTotalXSPlots(hNewinc_truth,hNewCEXint_truth,hbeamint_truth_tmp,hbeamint_truth_tmp,hbeamint_truth_tmp,hbeamint_truth_tmp,g_cex,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data,false);
  
  DrawTotalXSPlots(hNewinc_truth,hNewint_truth,hbeamint_truth_tmp,hbeamint_truth_tmp,hbeamint_truth_tmp,hbeamint_truth_tmp,g_totalinel,totalCEXXS_800_truth,totalCEXXS_800_reco,totalCEXXS_800_data,true);

  //vector<double> binning_100MeV = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000.};
  //cout << "bin size: " << binning_100MeV.size() << endl;
  //cout << "bin at 0: " << binning_100MeV.at(0) << endl;
  //cout << "bin at 19: " << binning_100MeV.at(19) << endl;
  //cout << "bin at 20: " << binning_100MeV.at(20) << endl;

  

  // Close files
  f_CrossSection->Close();
  file->Close();

}

void DrawHistOutput(TH1D *h_truth, TH1D *h_reco, TH1D *h_unfold, TH1D *h_unfold_data, const bool &rebin){
  // Rebin the histogram
  if(rebin){
    //h_truth->Rebin(50); h_reco->Rebin(50); h_unfold->Rebin(50); h_unfold_data->Rebin(50);
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
  //h_reco->Draw("e1 sames");

  h_unfold->SetMarkerStyle(8);
  h_unfold->SetMarkerSize(1);
  h_unfold->SetMarkerColor(kGreen+3);
  h_unfold->SetLineColor(kGreen+3);
  h_unfold->SetLineWidth(1);
  //h_unfold->Draw("e1 sames");

  h_unfold_data->SetMarkerStyle(8);
  h_unfold_data->SetMarkerSize(1);
  h_unfold_data->SetMarkerColor(kBlack);
  h_unfold_data->SetLineColor(kBlack);
  h_unfold_data->SetLineWidth(1);
  //h_unfold_data->Draw("e1 sames");

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

void DrawTotalXSPlots(TH1D *hinc_truth, TH1D *hint_truth, TH1D *hinc_reco, TH1D *hint_reco, TH1D *hinc_data, TH1D *hint_data, TGraph* g_cex, double & totalCEXXS_800_truth, double & totalCEXXS_800_reco, double & totalCEXXS_800_data, const bool kTotal){

  TLatex tt;
  tt.SetNDC();

  TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
  
  // Truth
  TH1D * hcex_truth = TotalCEXXSCal(hinc_truth,hint_truth,true,false,false,kTotal);
  // Reco
  TH1D * hcex_reco = TotalCEXXSCal(hinc_reco,hint_reco,true,false,false,kTotal);
  // Data
  TH1D * hcex_data = TotalCEXXSCal(hinc_data,hint_data,true,false,false,kTotal);

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
  TH1D *hratio = (TH1D*)hcex_truth->Clone("hcex_truth");
  
  const Int_t x0_x = hcex_truth->GetXaxis()->GetFirst();
  const Int_t x1_x = hcex_truth->GetXaxis()->GetLast();

  hratio->Scale(0);

  for(Int_t ix=x0_x; ix<=x1_x; ix++){
    double bin_center = hcex_truth->GetBinCenter(ix);
    double MC = g_cex->Eval(bin_center);
    
    double data = hcex_truth->GetBinContent(ix);
    double edata = hcex_truth->GetBinError(ix);
    if(MC != 0){
      double ratio = data/MC;
      hratio->SetBinContent(ix,ratio);
      cout << "\nratio: " << ratio << endl;
      cout << "\ndiff: " << MC - data << endl;

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

  if(!kTotal) c1->Print("output/truth_TotalCEXXS.png");
  else c1->Print("output/truth_TotalInelXS.png");

  
}


void DrawDiffXSPlots(TH1D *hpi0KE_truth, TH1D *hpi0KE_reco, TH1D *hpi0KE_data, TGraph* g_775, double Int_800_truth, double Int_800err_truth, double Int_800_reco, double Int_800err_reco, double Int_800_data, double Int_800err_data, double totalCEXXS_800_truth, double totalCEXXS_800_reco, double totalCEXXS_800_data, TString name){

  TLatex tt;
  tt.SetNDC();

  TCanvas * c2 = new TCanvas("c2", "", 1200, 800);
  
  double scale_bin = 100.0;
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



TH1D * TotalCEXXSCal(TH1D * hh, TH1D * InteractingHist, const bool & Eslice, const bool & widerBin, const bool & newMethod, const bool kTotal)
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
  //double meandEdx[] = {5.06921, 2.78604, 2.35879, 2.20657, 2.14137, 2.11321, 2.10322, 2.10293, 2.10805, 2.11628, 2.12627, 2.13724, 2.14869, 2.16033, 2.17195,
   //                    2.18344, 2.19472, 2.20573, 2.21646, 2.2269};
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
    vector<double> KE_binning = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1050, 1100, 1150, 1200};
    double this_KE = (KE_binning[ix] + KE_binning[ix - 1]) / 2.0;
    double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
    //cout << "this KE: " << this_KE << "this_dEdx: " << this_dEdx << endl;
    //if(Eslice) ratio = ratio * meandEdx[ix-1]*sigma_factor*1e27;
    if(Eslice) ratio = ratio * this_dEdx*sigma_factor*1e27;

    //if(incE != 0 && kTotal) cout << "ix: "<< ix << "incE: " << incE << "intE: " << intE << " ratio: " << ratio << "meandEdx[ix]: " << meandEdx[ix-1] << endl;

    // Simple ratio form 
    //double ratio = intE/incE;

    // If the incE entry is not zero set the bin content
    if(incE != 0) xsec->SetBinContent(ix,ratio);
    if( Eslice && !kTotal && ratio > 150) xsec->SetBinContent(ix,9999);
    if( Eslice && kTotal && ratio > 750) xsec->SetBinContent(ix,9999);

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
  if(!kTotal)xsec->SetMaximum(300);
  if(kTotal)xsec->SetMaximum(2200);
   
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
  IncidentHist->SetName("i001hTruthBeamIncidentHist");
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
TH1D * GetNewIncidentHist(TH1D * InitialHist, TH1D * InteractingHist){

  TH1D *IncidentHist = (TH1D*)InitialHist->Clone();
  const Int_t N_KE_bins = InteractingHist->GetXaxis()->GetLast();
  
  for(int i = 1; i < N_KE_bins; i++){
    double sum_N_end = 0.;
    double sum_N_init = 0.;
    for(int j = 1; j < i + 1; j++){
      cout << "[sum_N_end] (i, j) : (" << i << ", " << j << ")" << endl;
      sum_N_end += InteractingHist->GetBinContent(j);
    }
    for(int j = 1; j < i + 0; j++){
      cout << "[sum_N_init] (i, j) : (" << i << ", " << j << ")" << endl;
      sum_N_init += InitialHist->GetBinContent(j);
    }
    double this_content = sum_N_end - sum_N_init;
    double this_err = sqrt(this_content);
    IncidentHist->SetBinContent(i, this_content);
    IncidentHist->SetBinError(i, this_err);
    cout << "[Draw_Eslice_distribution] " << ", " << i << " : sum_N_end : " << sum_N_end << ", sum_N_init : " << sum_N_init << ", this_content : " << this_content << endl;
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
  if(IsRange){
    double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     0.0, 0.801123, 0.859177, 0.86919, 0.880718, 0.886977, 0.889779, 0.887032};
    weight *= IntWeight[binx-1];
  }

  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}

double CalBeamIntWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  if(IsRange){
    double IntWeight[] = {0.439768, 1.0, 0.0, 0.240104, 0.325475, 0.547874, 0.525402, 0.569011, 0.663258, 0.755661, 0.820706, 0.858978, 
                     0.886342, 0.901809, 0.920025, 0.931599, 0.945724, 0.945179, 0.943406, 0.938913};
    weight *= IntWeight[binx-1];
  }

  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}

double CalCEXIntWeight(const int & binx, const bool IsRange){

  double weight = 1.;
  if(IsRange){
    double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.68125, 0.811044, 0.829454, 0.84448, 0.600005, 
                    0.731605, 0.683442, 0.650626, 0.651113, 0.668391, 0.633431, 0.65286, 1.0};
    
    weight *= IntWeight[binx-1];
  }

  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}


double CalCEXPi0KEWeight(const int & binx, const bool IsRange, const bool IsKF){

  double weight = 1.;
  if(IsRange && IsKF){
    // ======================= Old weight =================// 
    // With KF
    //double IntWeight[] = {0.0, 0.331871, 0.49273, 0.585248, 0.5774, 0.733449, 0.627068, 0.610235, 0.560636, 0.599749, 0.854405, 
    //                  0.643868, 0.943692, 0.948054, 0.767483, 1.0, 1,0, 0.0, 0.0, 0.0};
    // No KF another weight
    //double IntWeight[] = {0.0, 0.0, 0.652846, 0.596701, 0.521758, 0.579765, 0.689989, 0.544885, 0.607432, 0.591426, 0.837472, 
    //                  0.624503, 0.944532, 0.940689, 0.778484, 1.0, 1,0, 0.0, 0.0, 0.0};

    //With KF both weight
    //double IntWeight[] = {0.0, 0.333956, 0.521258, 0.557553, 0.501853, 0.74563, 0.648465, 0.569635, 0.556742, 0.593001, 0.868087, 
    //                  0.693125, 0.954237, 0.967019, 0.798809, 1.0, 1,0, 0.0, 0.0, 0.0};
    // test
    double IntWeight[] = {0.333956, 0.542595, 0.617736, 0.6163, 0.57292, 0.804125, 0.95878, 0.886411, 1.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  else if(IsRange && !IsKF){
    // No KF only beam weight
    //double IntWeight[] = {0.0, 0.0, 0.652846, 0.596701, 0.521758, 0.579765, 0.689989, 0.544885, 0.607432, 0.591426, 0.837472, 
    //                  0.624503, 0.944532, 0.940689, 0.778484, 1.0, 1,0, 0.0, 0.0, 0.0};
    // No KF both weight
    //double IntWeight[] = {0.0, 0.0, 0.71829, 0.61267, 0.480149, 0.505302, 0.710389, 0.532536, 0.606493, 0.573254, 0.869121, 
    //                  0.668477, 0.954928, 0.962236, 0.808697, 1.0, 1,0, 0.0, 0.0, 0.0};
    // test
    double IntWeight[] = {0.0, 0.648345, 0.492083, 0.631761, 0.590167, 0.79236, 0.957278, 0.886136, 1.0, 0.0};

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
double CalCEXPi0CosThetaWeight(const int & binx, const bool IsRange, const bool IsKF){

  double weight = 1.;
  if(IsRange && IsKF){
    // ======================= Old weight =================// 
    // With KF
    //double IntWeight[] = {0.509613, 0.709514, 0.665169, 1, 0.633065, 0.581289, 0.671715, 0.554388, 0.592632, 0.727953};
    // NoKF another weight
    //double IntWeight[] = {0.517665, 0.583722, 0.599486, 0.868833, 0.679566, 0.503973, 0.573794, 0.566587, 0.621706, 0.706317};

    // With KF both weight
    double IntWeight[] = {0.517665, 0.583722, 0.599486, 0.868833, 0.679566, 0.503973, 0.573794, 0.566587, 0.621706, 0.706317};
    
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
  }
  return weight;
}
// Theta
double CalCEXPi0ThetaWeight(const int & binx, const bool IsRange, const bool IsKF){

  double weight = 1.;
  if(IsRange && IsKF){
    // ======================= Old weight =================// 
    // With KF
    //double IntWeight[] = {0.712032, 0.72855, 0.540438, 0.680656, 0.574381, 0.827391, 0.657209, 0.765056, 0.0};
    // NoKF another weight
    //double IntWeight[] = {0.753768, 0.69482, 0.570662, 0.584574, 0.553642, 0.769505, 0.523964, 0.7708, 0.0};

    // With KF both weight
    double IntWeight[] = {0.753768, 0.69482, 0.570662, 0.584574, 0.553642, 0.769505, 0.523964, 0.7708, 0.0};
    
    weight *= IntWeight[binx-1];
  }
  else if(IsRange && !IsKF){
    // No KF only beam weight
    //double IntWeight[] = {0.715657, 0.718822, 0.571186, 0.627761, 0.519831, 0.695564, 0.657209, 0.765056, 0.0};
    // No KF both weight
    double IntWeight[] = {0.753768, 0.69482, 0.570662, 0.584574, 0.553642, 0.769505, 0.523964, 0.7708, 0.0};
    
    weight *= IntWeight[binx-1];

  }
  else{
    double IntWeight[] = {0.0, 0.0, 0.636412, 0.570768, 0.348301, 0.633125, 0.747292, 0.461547, 0.711167, 0.603851, 0.660517, 
                    0.747348, 0.859409, 1, 0.742225, 1.0, 1,0, 0.0, 0.0, 0.0};
    weight *= IntWeight[binx-1];
  }
  return weight;
}

void WeightBeamData(TH1D * hbeam_data, const bool IsRange, const TString name){

  const Int_t x0 = hbeam_data->GetXaxis()->GetFirst();
  const Int_t x1 = hbeam_data->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double ientry = hbeam_data->GetBinContent(ix);
    double weight = -999;
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

double Density_Correction(double beta, double gamma){
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

double dEdx_Bethe_Bloch(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = Density_Correction(beta, gamma);

  // == dE/dx with the density correction
  double f = LAr_density * K * (18.0 / M_Ar) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm]

  return this_dEdx;
}




