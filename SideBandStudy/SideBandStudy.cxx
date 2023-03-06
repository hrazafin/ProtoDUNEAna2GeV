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

#include "TemplateFitter.h"
#include "TemplateFitter.cxx"


using namespace std;
double plotScale = 0.499375;

double GetChi2(TH1D * hdata, TH1D * hmc);
void DrawOutput(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h2_low, TH1D *h2_high, TH1D *h2_mid, const double &lowScale, const double &highScale, TString tag);
void GetPi0EnergyScale(const TString tag, double & low_scale, double & high_scale);
void SidebandDefinition();

int main(int argc, char * argv[])
{ 
  double low_low = 1.0, low_high = 1.0;
  GetPi0EnergyScale("lowPi0E",low_low,low_high);
  double high_low = 1.0, high_high = 1.0;
  GetPi0EnergyScale("highPi0E",high_low,high_high);

  SidebandDefinition();
  
}

void GetPi0EnergyScale(const TString tag, double & low_scale, double & high_scale){
  
  TH1D * h0 = new TH1D(Form("h0_%s",tag.Data()),";Pi0 Kinetic Energy (MeV); Candidates", 20, 0, 0.5); 
  TH1D * h1 = new TH1D(Form("h1_%s",tag.Data()),";Pi0 Kinetic Energy (MeV); Candidates", 20, 0, 0.5); 
  TH1D * h2 = new TH1D(Form("h2_%s",tag.Data()),";Pi0 Kinetic Energy (MeV); Candidates", 20, 0, 0.5);

  // Scale components
  TH1D * h2_low = new TH1D(Form("h2_low_%s",tag.Data()),";Pi0 Kinetic Energy (MeV); Candidates", 20, 0, 0.5); 
  TH1D * h2_high = new TH1D(Form("h2_high_%s",tag.Data()),";Pi0 Kinetic Energy (MeV); Candidates", 20, 0, 0.5); 
  TH1D * h2_mid = new TH1D(Form("h2_mid_%s",tag.Data()),";Pi0 Kinetic Energy (MeV); Candidates", 20, 0, 0.5); 

  //const TString finName = "input/outana_noKF.root";
  const TString finName = "input/outana_withKF_new.root";
    
  TFile *file = TFile::Open(finName);

  if(!file->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }
  TH2D * h2d_mc = (TH2D*) file->Get("mc/testb022hRecPi0Mass_COMPOSE_EVT;1");
  TH2D * h2d_data = (TH2D*) file->Get("data/testb022hRecPi0Mass_COMPOSE_EVT;1");

  if(tag.Contains("highPi0E")){
    h2d_mc = (TH2D*) file->Get("mc/test1b022hRecPi0Mass_COMPOSE_EVT;1");
    h2d_data = (TH2D*) file->Get("data/test1b022hRecPi0Mass_COMPOSE_EVT;1");
  }

  //TH2D * h2d_mc_high = (TH2D*) file->Get("mc/test1b022hRecPi0Mass_COMPOSE;1");
  //TH2D * h2d_data_high = (TH2D*) file->Get("data/test1b022hRecPi0Mass_COMPOSE;1");

  int nx = h2d_mc->GetNbinsX();
  int ny = h2d_data->GetNbinsY();

  for(int ix=0; ix<=nx+1; ix++){
    // Data
    double ientry_data = h2d_data->GetBinContent(ix, ny);
    double ientry_dataerr = h2d_data->GetBinError(ix, ny);
    // MC component
    double ientry_mc_sig = h2d_mc->GetBinContent(ix, 1);
    double ientry_mc_abs = h2d_mc->GetBinContent(ix, 2);
    double ientry_mc_0pi0 = h2d_mc->GetBinContent(ix, 3);
    double ientry_mc_1pi0 = h2d_mc->GetBinContent(ix, 4);
    double ientry_mc_mulpi0 = h2d_mc->GetBinContent(ix, 5);
    double ientry_mc_beambck = h2d_mc->GetBinContent(ix, 6);
    // MC component error
    double ientry_mc_sigerr = h2d_mc->GetBinError(ix, 1);
    double ientry_mc_abserr = h2d_mc->GetBinError(ix, 2);
    double ientry_mc_0pi0err = h2d_mc->GetBinError(ix, 3);
    double ientry_mc_1pi0err = h2d_mc->GetBinError(ix, 4);
    double ientry_mc_mulpi0err = h2d_mc->GetBinError(ix, 5);
    double ientry_mc_beambckerr = h2d_mc->GetBinError(ix, 6);

    h0->SetBinContent(ix, ientry_data);
    h0->SetBinError(ix, ientry_dataerr);

    double ifixed = ientry_mc_sig; 
    double ifixederr = ientry_mc_sigerr;//sqrt(pow(ientry_mc_sigerr,2)+pow(ientry_mc_abserr,2)+pow(ientry_mc_0pi0err,2)+pow(ientry_mc_beambckerr,2));
    h1->SetBinContent(ix, ifixed);
    h1->SetBinError(ix, ifixederr);

    double ivaried = ientry_mc_abs + ientry_mc_0pi0 + ientry_mc_beambck + ientry_mc_1pi0 + ientry_mc_mulpi0;
    double ivariederr = sqrt(pow(ientry_mc_abserr,2)+pow(ientry_mc_0pi0err,2)+pow(ientry_mc_beambckerr,2)+pow(ientry_mc_1pi0err,2)+pow(ientry_mc_mulpi0err,2));

    h2->SetBinContent(ix, ivaried);
    h2->SetBinError(ix, ivariederr);
  }
  
  file->Close();

  TemplateFitter fitter;   
  // Scale MC to data
  h1->Scale(plotScale); h2->Scale(plotScale);

  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(2, 11);
  fitter.Fit();
  double par_low = fitter.GetPar();
  double parerr_low = fitter.GetParError();
  cout << "===== Low Region ==== " << endl;
  cout << par_low << " " << parerr_low << endl;
  low_scale = par_low;
/*
  fitter.SetFitRange(11, 19);
  fitter.Fit();
  double par_high = fitter.GetPar();
  double parerr_high = fitter.GetParError();
  cout << "===== High Region ==== " << endl;
  cout << par_high << " " << parerr_high << endl;
  high_scale = par_high;
*/
  high_scale = par_low;

  for(int ix=0; ix<=nx+1; ix++){
    h2_low->SetBinContent(ix,0);
    h2_high->SetBinContent(ix,0);
    h2_mid->SetBinContent(ix,0);
    double ientry = h2->GetBinContent(ix);
    if(ix <= 2) h2_low->SetBinContent(ix, ientry);
    else if(ix >= 11) h2_high->SetBinContent(ix, ientry);
    else h2_mid->SetBinContent(ix, ientry);
  }
  // Save the optimised scaling factor and error
  //Par.push_back(par_low); Par.push_back(par_high); 
  //Parerr.push_back(parerr_low); Parerr.push_back(parerr_high); 
  DrawOutput(h0, h1, h2, h2_low, h2_high, h2_mid, par_low, par_low, tag);
  
}


double GetChi2(TH1D * hdata, TH1D * hmc){
    int nx = hdata->GetNbinsX();
    double chisq = 0;
    for(int ix=1; ix<=nx; ix++){
      double x = hdata->GetBinContent(ix);
      double y = hmc->GetBinContent(ix);
        
      double ex = hdata->GetBinError(ix);
      double ey = hmc->GetBinError(ix);

      if (!ex) ex = 1;
      chisq += pow(x-y,2)/(pow(ex,2)+pow(ey,2));
    }
    return chisq;
}

void DrawOutput(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h2_low, TH1D *h2_high, TH1D *h2_mid, const double &lowScale, const double &highScale, TString tag)
{
    TLatex tt;
    tt.SetNDC();

    TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 1200, 800);
    gStyle->SetOptStat(0);

    h0->SetMarkerStyle(25);
    h0->SetMarkerSize(1);
    h0->SetMarkerColor(kBlack);
    h0->SetLineColor(kBlack);

    h0->Draw("e1");
    TH1D * hmc = (TH1D*)h1->Clone();
    hmc->Add(h2);
    hmc->SetLineColor(kRed);

    auto lg = new TLegend(0.65,0.7,0.85,0.88);
    lg->SetHeader(Form("#chi^{2}: %f",GetChi2(h0,hmc)));
    lg->AddEntry(h0,"data","lp");
    lg->AddEntry(hmc,"MC","l");
    lg->SetBorderSize(0);
    lg->Draw("sames");
    hmc->Draw("hists sames");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    c1->Print("output/ori_"+tag+".png");

    TCanvas * c2 = new TCanvas(Form("c2_%s",tag.Data()), "", 1200, 800);

    h0->SetMarkerStyle(8);
    h0->SetMarkerSize(1);
    h0->SetMarkerColor(kBlack);
    h0->SetLineColor(kBlack);
    h0->SetMaximum(h0->GetMaximum()*1.5);

    h0->Draw("e1");
    TH1D * hmc_fit = (TH1D*)h1->Clone();
    hmc_fit->Add(h2_low,lowScale);
    hmc_fit->Add(h2_high,highScale);
    hmc_fit->Add(h2_mid,1.0);

    hmc_fit->SetLineColor(kRed);
    hmc_fit->SetLineWidth(3);
    
    hmc_fit->Draw("hists sames");

    hmc->SetLineColor(kGreen+3);
    hmc->SetLineStyle(9);
    hmc->SetLineWidth(3);
    hmc->Draw("hists sames");

    auto lgfit = new TLegend(0.65,0.60,0.85,0.88);
    //lgfit->SetHeader(Form("#chi^{2}: %f",GetChi2(h0,hmc_fit)));
    lgfit->SetHeader(Form("#chi^{2}_{def}: %.2f, #chi^{2}_{fit}: %.2f", GetChi2(h0,hmc), GetChi2(h0,hmc_fit)));
    lgfit->AddEntry(h0,"data","lp");
    lgfit->AddEntry(hmc,"MC Default","l");
    lgfit->AddEntry(hmc_fit,"MC Fit","l");
    lgfit->SetBorderSize(0);
    lgfit->Draw("sames");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    c2->Print("output/fit_"+tag+".png");

}

void SidebandDefinition(){
  const TString finName = "input/outana_sideband.root";
    
  TFile *file = TFile::Open(finName);

  if(!file->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  TH2D * sideband = (TH2D*) file->Get("mc/b099hRecPi0KineticEnergyVSPi0Mass;1");

  TLatex tt;
  tt.SetNDC();

  TCanvas * c1 = new TCanvas("c1_sideband", "", 1200, 800);
  

  sideband->SetFillColor(kBlue-3);
  
  sideband->GetYaxis()->SetTitle("#pi^{0} Mass (GeV/c^{2})");
  sideband->GetXaxis()->SetTitle("Kinetic Energy (GeV)");

  sideband->Draw("BOX");

  TLine *line = new TLine(0.3,0,0.3,0.05);

  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->SetLineWidth(3);

  line->Draw("sames");

  TLine *line3 = new TLine(0.3,0.25,0.3,0.5);

  line3->SetLineColor(kRed);
  line3->SetLineStyle(2);
  line3->SetLineWidth(3);

  line3->Draw("sames");

  TLine *line1 = new TLine(0,0.05,1,0.05);

  line1->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line1->SetLineWidth(3);

  line1->Draw("sames");

  TLine *line2 = new TLine(0,0.25,1,0.25);

  line2->SetLineColor(kRed);
  line2->SetLineStyle(2);
  line2->SetLineWidth(3);

  line2->Draw("sames");

  tt.DrawLatex(0.755,0.255,"#color[3]{#bf{#it{Signal}}}");

  tt.DrawLatex(0.705,0.785,"#color[2]{#bf{#it{#alpha_{2} = 0.96}}}");
  tt.DrawLatex(0.705,0.625,"#color[6]{#bf{#it{Sideband 2}}}");
  tt.DrawLatex(0.705,0.125,"#color[6]{#bf{#it{Sideband 2}}}");

  tt.DrawLatex(0.135,0.785,"#color[2]{#bf{#it{#alpha_{1} = 1.79}}}");

  tt.DrawLatex(0.135,0.625,"#color[6]{#bf{#it{Sideband 1}}}");
  tt.DrawLatex(0.135,0.125,"#color[6]{#bf{#it{Sideband 1}}}");


  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

  c1->Print("output/sideband.png");

}