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

#include "TEllipse.h"

using namespace std;

void DrawOutput(TH2D* h2d_beam, TH2D* h2d_scraper, const double & x_mean, const double & y_mean, const double & radius, const double & multi);
TF1 * FitGausFunc(TH1D * hh, const double &min, const double &max);
void Draw1D(TH1D* hh, TF1* f1);
void DoPlots(TH2D* h2d_beam, TH2D* h2d_scraper,TH1D* h1d_beamx, TH1D* h1d_beamy);
void DrawBeamScraperDef(TH1D* hh, bool kCut = false);
void DrawUpStreamELossAfterCuts(TH2D * h2d_InstE);
void SetProtoDUNELabel();
void SetBeamInstELabel(const TString name);
void SetTitleFormat(TH1 * hh);
void SetTitleFormat(TH2 * hh);
TGraphErrors * ConvertTH1DtoGraph(TH1D * hist);

int main(int argc, char * argv[])
{
    
    const TString finName = "input/outana.root";
    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    // ========== Beam Inst P Before Smearing ========== //
    TH2D * h2d_InstP_mc_raw = (TH2D*) file->Get("mc/i057hRecPiPlusInstMomentumNoSmearing_STK;1");
    TH1D * h1d_InstP_data_raw = (TH1D*) file->Get("data/i057hRecPiPlusInstMomentumNoSmearing_STK_sum;1");
    TH1D * h2d_InstP_mc_projX_raw = h2d_InstP_mc_raw->ProjectionX();
    TF1 *f1_raw = FitGausFunc(h2d_InstP_mc_projX_raw,850.0,1150.0);
    double mc_sigma_raw = f1_raw->GetParameter("Sigma");
    double mc_mean_raw = f1_raw->GetParameter("Mean");
    TF1 *f2_raw = FitGausFunc(h1d_InstP_data_raw,850.0,1150.0);
    double data_sigma_raw = f2_raw->GetParameter("Sigma");
    double data_mean_raw = f2_raw->GetParameter("Mean");
    cout << "mu_raw: " << (data_mean_raw - mc_mean_raw)/1000.0 << " sigma_raw: " << sqrt(pow(data_sigma_raw/1000.0,2)-pow(mc_sigma_raw/1000.0,2)) << endl;
    Draw1D(h2d_InstP_mc_projX_raw,f1_raw);
    Draw1D(h1d_InstP_data_raw,f2_raw);

    // ========== Beam Inst P After Smearing ========== //
    TH2D * h2d_InstP_mc = (TH2D*) file->Get("mc/i051hRecPiPlusInstMomentum_STK;1");
    TH1D * h1d_InstP_data = (TH1D*) file->Get("data/i051hRecPiPlusInstMomentum_STK_sum;1");
    TH1D * h2d_InstP_mc_projX = h2d_InstP_mc->ProjectionX();
    TF1 *f1 = FitGausFunc(h2d_InstP_mc_projX,850.0,1150.0);
    double mc_sigma = f1->GetParameter("Sigma");
    double mc_mean = f1->GetParameter("Mean");
    TF1 *f2 = FitGausFunc(h1d_InstP_data,850.0,1150.0);
    double data_sigma = f2->GetParameter("Sigma");
    double data_mean = f2->GetParameter("Mean");
    cout << "mu: " << (data_mean - mc_mean)/1000.0 << " sigma: " << sqrt(pow(data_sigma/1000.0,2)-pow(mc_sigma/1000.0,2)) << endl;
    Draw1D(h2d_InstP_mc_projX,f1);
    Draw1D(h1d_InstP_data,f2);

    // ========== Define Beam Scraper Events ========== //
    TH1D * h1d_upEloss700MeV = (TH1D*) file->Get("mc/j003hUpStreamELoss700MeV;1");
    TH1D * h1d_upEloss800MeV = (TH1D*) file->Get("mc/j003hUpStreamELoss800MeV;1");
    TH1D * h1d_upEloss900MeV = (TH1D*) file->Get("mc/j003hUpStreamELoss900MeV;1");
    TH1D * h1d_upEloss1000MeV = (TH1D*) file->Get("mc/j003hUpStreamELoss1000MeV;1");

    DrawBeamScraperDef(h1d_upEloss700MeV);
    DrawBeamScraperDef(h1d_upEloss800MeV);
    DrawBeamScraperDef(h1d_upEloss900MeV);
    DrawBeamScraperDef(h1d_upEloss1000MeV);


    // ========== Beam Scraper Cut ========== //

    TH2D * h2d_beam700MeV = (TH2D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam700MeV_RES;1");
    TH2D * h2d_scraper700MeV = (TH2D*) file->Get("mc/j005hBeamInstXVSBeamInstYScraper700MeV_RES;1");

    TH1D * h1d_beamx700MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam700MeV_RES_projX;1");
    TH1D * h1d_beamy700MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam700MeV_RES_projY;1");

    DoPlots(h2d_beam700MeV,h2d_scraper700MeV,h1d_beamx700MeV,h1d_beamy700MeV);

    TH2D * h2d_beam800MeV = (TH2D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam800MeV_RES;1");
    TH2D * h2d_scraper800MeV = (TH2D*) file->Get("mc/j005hBeamInstXVSBeamInstYScraper800MeV_RES;1");

    TH1D * h1d_beamx800MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam800MeV_RES_projX;1");
    TH1D * h1d_beamy800MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam800MeV_RES_projY;1");

    DoPlots(h2d_beam800MeV,h2d_scraper800MeV,h1d_beamx800MeV,h1d_beamy800MeV);

    TH2D * h2d_beam900MeV = (TH2D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam900MeV_RES;1");
    TH2D * h2d_scraper900MeV = (TH2D*) file->Get("mc/j005hBeamInstXVSBeamInstYScraper900MeV_RES;1");

    TH1D * h1d_beamx900MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam900MeV_RES_projX;1");
    TH1D * h1d_beamy900MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam900MeV_RES_projY;1");

    DoPlots(h2d_beam900MeV,h2d_scraper900MeV,h1d_beamx900MeV,h1d_beamy900MeV);

    TH2D * h2d_beam1000MeV = (TH2D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam1000MeV_RES;1");
    TH2D * h2d_scraper1000MeV = (TH2D*) file->Get("mc/j005hBeamInstXVSBeamInstYScraper1000MeV_RES;1");

    TH1D * h1d_beamx1000MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam1000MeV_RES_projX;1");
    TH1D * h1d_beamy1000MeV = (TH1D*) file->Get("mc/j004hBeamInstXVSBeamInstYBeam1000MeV_RES_projY;1");

    DoPlots(h2d_beam1000MeV,h2d_scraper1000MeV,h1d_beamx1000MeV,h1d_beamy1000MeV);



    // ========== After Beam Scraper Cut (Check Performance) ========== //
    TH1D * h1d_upEloss700MeV_afterCuts = (TH1D*) file->Get("mc/j007hUpStreamELoss700MeVAfterScraperCuts;1");
    TH1D * h1d_upEloss800MeV_afterCuts = (TH1D*) file->Get("mc/j007hUpStreamELoss800MeVAfterScraperCuts;1");
    TH1D * h1d_upEloss900MeV_afterCuts = (TH1D*) file->Get("mc/j007hUpStreamELoss900MeVAfterScraperCuts;1");
    TH1D * h1d_upEloss1000MeV_afterCuts = (TH1D*) file->Get("mc/j007hUpStreamELoss1000MeVAfterScraperCuts;1");
    
    DrawBeamScraperDef(h1d_upEloss700MeV_afterCuts,true);
    DrawBeamScraperDef(h1d_upEloss800MeV_afterCuts,true);
    DrawBeamScraperDef(h1d_upEloss900MeV_afterCuts,true);
    DrawBeamScraperDef(h1d_upEloss1000MeV_afterCuts,true);

    // ========== Upstream Energy Loss (after smearing and beam reweight) ========== //

    TH2D * h2d_InstE = (TH2D*) file->Get("mc/j008hUpStreamELossAfterSmearingAndWeight;1");
    
    DrawUpStreamELossAfterCuts(h2d_InstE);
    
    file->Close();

}

void DoPlots(TH2D* h2d_beam, TH2D* h2d_scraper,TH1D* h1d_beamx, TH1D* h1d_beamy){

/*
    // Get my own circle but now used Sungbin's cut below
    TF1 *f1 = FitGausFunc(h1d_beamx,-50.0,0.0);
    double x_sigma = f1->GetParameter("Sigma");
    double x_mean = f1->GetParameter("Mean");

    TF1 *f2 = FitGausFunc(h1d_beamy,400.0,450.0);
    double y_sigma = f2->GetParameter("Sigma");
    double y_mean = f2->GetParameter("Mean");

    double radius = sqrt(pow(x_sigma,2) + pow(y_sigma,2)); 

    cout << "x_mean: " << x_mean << endl;
    cout << "y_mean: " << y_mean << endl;
    cout << "radius: " << radius << endl;
*/
    // Draw projection plots
    //Draw1D(h1d_beamx,f1);
    //Draw1D(h1d_beamy,f2);

    // Get my own circle but now used Sungbin's cut below
    //DrawOutput(h2d_beam, h2d_scraper, x_mean, y_mean, radius, 1.4);
    // Sungbin's circle for MC
    DrawOutput(h2d_beam, h2d_scraper, -29.6, 422, 4.8, 1.4);
}

TF1 * FitGausFunc(TH1D * hh, const double &min, const double &max){
  TString tag = hh->GetName();
  TF1 *f1 = new TF1(Form("f1%s",tag.Data()),"gaus",min,max);    
  if (hh->GetEntries() != 0) {
    f1->SetParameters(hh->GetMaximum(), hh->GetMean(), hh->GetRMS() );
    hh->Fit(Form("f1%s",tag.Data()),"q");
  }
  return f1;
}

void Draw1D(TH1D* hh, TF1* f1){
  TString tag = hh->GetName();
  TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 800, 800);
  const Double_t currentLeft=0.12, currentTop=0.09, currentRight=0.13, currentBottom=0.14;
  c1->SetTicks(1,1);
  c1->SetLeftMargin(currentLeft);
  c1->SetTopMargin(currentTop);
  c1->SetRightMargin(currentRight);
  c1->SetBottomMargin(currentBottom);
  c1->SetFillColor(0);

  hh->Draw("hists");
  f1->Draw("sames C");
  
  c1->Print("output/"+tag+".pdf");

}
void DrawOutput(TH2D* h2d_beam, TH2D* h2d_scraper, const double & x_mean, const double & y_mean, const double & radius, const double & multi){

    TLatex tt;
    tt.SetNDC();

    TString tagbeam = h2d_beam->GetName();

    TCanvas * c1 = new TCanvas(Form("c1_%s",tagbeam.Data()), "", 800, 800);
    const Double_t currentLeft=0.12, currentTop=0.09, currentRight=0.13, currentBottom=0.14;
    c1->SetTicks(1,1);
    c1->SetLeftMargin(currentLeft);
    c1->SetTopMargin(currentTop);
    c1->SetRightMargin(currentRight);
    c1->SetBottomMargin(currentBottom);
    c1->SetFillColor(0);

    gStyle->SetOptStat(0);
    
    //h2d_beam -> SetMarkerColor(kGreen);
    //h2d_scraper -> SetMarkerColor(kRed);
    //h2d_beam -> SetLineColor(kGreen);
    //h2d_scraper -> SetLineColor(kRed);
    //h2d_beam -> SetFillColor(kWhite);
    //h2d_scraper -> SetFillColor(kWhite);

    h2d_beam->SetMarkerStyle(7);
    h2d_beam->SetFillColor(kGreen);
    h2d_beam->SetMarkerColor(kGreen);
    h2d_beam->SetLineColor(kGreen);
    h2d_beam->GetXaxis()->SetTitle("X_{inst} (cm)");
    h2d_beam->GetYaxis()->SetTitle("Y_{inst} (cm)");
    SetTitleFormat(h2d_beam);
    h2d_scraper->SetMarkerStyle(7);
    h2d_scraper->SetMarkerColor(kRed);
    h2d_scraper->SetLineColor(kRed);

    h2d_scraper->SetFillColor(kRed);

    h2d_beam->Draw("");
    h2d_scraper->Draw("same");
    TString tit = "";
    if(tagbeam.Contains("700")) tit = "KE_{Beam  Inst.}^{reco.} (700-800 MeV)";
    else if(tagbeam.Contains("800")) tit = "KE_{Beam  Inst.}^{reco.} (800-900 MeV)";
    else if(tagbeam.Contains("900")) tit = "KE_{Beam  Inst.}^{reco.} (900-1000 MeV)";
    else if(tagbeam.Contains("1000")) tit = "KE_{Beam  Inst.}^{reco.} (1000-1100 MeV)";

/*
    h2d_beam -> SetMarkerColor(kGreen);
    h2d_scraper -> SetMarkerColor(kRed);
    h2d_beam -> SetLineColor(kGreen);
    h2d_scraper -> SetLineColor(kRed);
    h2d_beam -> SetFillColor(kWhite);
    h2d_scraper -> SetFillColor(kWhite);

    h2d_beam->Draw("box");
    h2d_scraper->Draw("boxsame");
*/
    TEllipse *el1 = new TEllipse(x_mean,y_mean,radius*multi,radius*multi);
    el1->SetFillStyle(0);
    el1->SetLineColor(kBlue);
    el1->SetLineWidth(3);
    el1->SetLineStyle(kDashed);
    el1->Draw("same");

    auto lg = new TLegend(0.45,0.6,0.85,0.88);
    lg->SetHeader(tit);
    lg->AddEntry(h2d_beam,"Non-Scraper","l");
    lg->AddEntry(h2d_scraper,"Scraper","l");
    //lg->AddEntry(el1,"Sungbin's cut","l");
    lg->AddEntry(el1,"Beam Scraper cut","l");
    lg->SetBorderSize(0);
    lg->SetTextFont(22);
    lg->Draw("sames");


    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    //tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
    c1->Print("output/"+tagbeam+".pdf");
}   

void DrawBeamScraperDef(TH1D* hh, bool kCut){
  
    TString tag = hh->GetName();
    TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 1200, 800);
    // Show the fitted pars
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    // Remove the borderSize size
    gStyle->SetStatBorderSize(-1);
    SetTitleFormat(hh);
    hh->Draw();
    
    TF1 *f1 = FitGausFunc(hh,-300,300);

    gPad->Update();
    TPaveStats *st = (TPaveStats*) gPad->GetPrimitive("stats");
    gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

    st->SetY1NDC(0.65); 
    st->SetY2NDC(0.85);
    st->SetX1NDC(0.65); 
    st->SetX2NDC(0.88);
    st->SetTextColor(kBlack);

    f1->Draw("sames C");

    double sigma = f1->GetParameter("Sigma");
    double mean = f1->GetParameter("Mean");
    // 3 sigma value
    double scraperValue = mean+sigma*3.0;
    //cout << "scraperValue: " << scraperValue << endl;
    //cout << "sigma: " << sigma << " mean: " << mean << endl;

    // Draw the line
    TLine *line = new TLine(scraperValue,0,scraperValue,f1->GetMaximum());
    line->SetLineColor(kBlue);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(3);
    line->Draw("sames");

    auto lg = new TLegend(0.15,0.6,0.4,0.88);
    lg->AddEntry(f1,"Gaus Fit","l");
    lg->AddEntry(line,"Scraper Definition","l");
    lg->SetBorderSize(0);
    //lg->SetTextFont(22);
    lg->Draw("sames");

    TLatex latex(.18,.52,Form("#bf{#mu + 3#sigma = %1.1f [MeV]}",scraperValue)); 
    latex.SetTextSize(0.04);
    latex.SetNDC(kTRUE);
    if(!kCut)latex.Draw();

    SetProtoDUNELabel();
    SetBeamInstELabel(tag);
    
    c1->Print("output/"+tag+".pdf");
}

void DrawUpStreamELossAfterCuts(TH2D * h2d_InstE){

  const int ny = h2d_InstE->GetNbinsY();
  const int nx = h2d_InstE->GetNbinsX();
  const double xmin = h2d_InstE->GetXaxis()->GetBinLowEdge(1);
  const double xmax = h2d_InstE->GetXaxis()->GetBinUpEdge(nx);

  //TH1D * hFitElossSW = new TH1D("hFitElossSW",";KE_{Beam  Inst.}^{reco.} (MeV);#mu(KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true}) (MeV)", 16, 400, 1200); 
  TH1D * hFitElossSW = new TH1D("hFitElossSW",";Beam Inst. KE (MeV);Upstream Energy Loss (MeV)", 16, 400, 1200); 
  TH1D * hFitElossSW_statsE = new TH1D("hFitElossSW_statsE",";Beam Inst. KE (MeV);UpStream Energy Loss (MeV)", 16, 400, 1200); 

  SetTitleFormat(hFitElossSW);

  for(int iy = 1; iy <= ny-1; iy++){

    TH1D * htmp = new TH1D(Form("tmp%d",iy), "", nx, xmin, xmax);
    for(int ix=0; ix<=nx+1; ix++){
      const double ientry = h2d_InstE->GetBinContent(ix, iy);
      htmp->SetBinContent(ix, ientry);
    }
    TF1 *f1 = FitGausFunc(htmp,-300,300);
    hFitElossSW->SetBinContent(iy+6,f1->GetParameter("Mean"));
    hFitElossSW_statsE->SetBinContent(iy+6,f1->GetParameter("Mean"));
    // Sigma band
    hFitElossSW->SetBinError(iy+6,f1->GetParameter("Sigma"));
    // Stats Error
    hFitElossSW_statsE->SetBinError(iy+6,f1->GetParError(1));

  }

 
  TCanvas * ctmp1 = new TCanvas("ctmp1", "", 1200, 800);
  TF1 *f_pol2 = new TF1("f_pol2", "pol2", 700, 1050);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  hFitElossSW->SetMaximum(120);
  hFitElossSW->SetMinimum(-55);

  hFitElossSW->SetMarkerColor(kBlack);
  hFitElossSW->SetMarkerStyle(kFullCircle);
  hFitElossSW->SetMaximum(120);
  hFitElossSW->SetStats(0);
  hFitElossSW->Draw("AXIS");

  auto graph = ConvertTH1DtoGraph(hFitElossSW);
  graph->SetLineWidth(0);
  graph->SetLineColor(kOrange);
  graph->SetMarkerSize(0);
  graph->SetFillColorAlpha(kOrange, 0.35);
  graph->Draw("3");

  hFitElossSW->Fit("f_pol2","0");
/*
  gPad->Update();
  TPaveStats *st = (TPaveStats*) gPad->GetPrimitive("stats");
  gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

  st->SetY1NDC(0.68); 
  st->SetY2NDC(0.88);
  st->SetX1NDC(0.65); 
  st->SetX2NDC(0.88);
  st->SetTextColor(kBlue);
*/
  hFitElossSW_statsE->SetMarkerColor(kBlack);
  hFitElossSW_statsE->SetMarkerStyle(kFullCircle);
  hFitElossSW_statsE->SetMaximum(120);
  hFitElossSW_statsE->Draw("e1 sames");

  f_pol2->Draw("sames C");

  TLine *line = new TLine(400,0,1200,0);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1);
  line->Draw("sames");


  //TLine *line1 = new TLine(400,11.966,1200,11.966);
  TLine *line1 = new TLine(400,13.583,1200,13.583);

  line1->SetLineColor(kBlue);
  line1->SetLineStyle(kDashed);
  line1->SetLineWidth(2);
  line1->Draw("sames");
  
  //TLine *line2 = new TLine(400,2.652,1200,2.652);
  TLine *line2 = new TLine(400,3.830,1200,3.830);

  line2->SetLineColor(kGreen+3);
  line2->SetLineStyle(kDashed);
  line2->SetLineWidth(2);
  //line2->Draw("sames");

  auto lg = new TLegend(0.15,0.68,0.55,0.88);
  lg->AddEntry(hFitElossSW,"#Delta KE_{upstream}","lp");
  lg->AddEntry(graph,"1 #sigma band","f");
  lg->AddEntry(f_pol2,"Fitted Poly2","l");
  lg->AddEntry(line1,"Mean Energy Loss (13.58 MeV)","l");
  //lg->AddEntry(line1,"Constant Eloss (13.583 MeV) - Before Smearing","l");
  //lg->AddEntry(line2,"Constant Eloss (3.702 MeV) - After Smearing","l");

  //lg->AddEntry(line1,"Constant Eloss (11.966 MeV) - Before Smearing","l");
  //lg->AddEntry(line2,"Constant Eloss (2.652 MeV) - After Smearing","l");

  lg->SetBorderSize(0);
  lg->Draw("sames");

  TLatex tt;
  tt.SetNDC();
  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");

  tt.DrawLatex(0.665,0.845,"#color[4]{p_{0} = 171.8 #pm 13.1}");
  tt.DrawLatex(0.665,0.785,"#color[4]{p_{1} = -0.575 #pm 0.030}");
  tt.DrawLatex(0.665,0.725,"#color[4]{p_{2} = 4.3e-4 #pm 1.71e-5}");


  ctmp1->Print("output/ElossSW.pdf");


  // Print const E loss before and after smearing
  TH1D * hLoss = h2d_InstE->ProjectionX();
 
  TF1 *f11 = new TF1("f11","gaus",-100,100);    
  cout << "Const. Energy loss from histogram (no fit): " << hLoss->GetMean() << endl;
  hLoss->Fit("f11","0");
  cout << "Const. Energy loss (after fit and smearing): " << f11->GetParameter("Mean") << endl;


}

void SetProtoDUNELabel(){
  TLatex tt;
  tt.SetNDC();
  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  //tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
}

void SetBeamInstELabel(const TString name){
  TLatex tt;
  tt.SetNDC();
  tt.SetTextSize(0.035);
  if(name.Contains("700")) tt.DrawLatex(0.600,0.925,"700 MeV < E_{inst}^{reco.} < 800 MeV");
  if(name.Contains("800")) tt.DrawLatex(0.600,0.925,"800 MeV < E_{inst}^{reco.} < 900 MeV");
  if(name.Contains("900")) tt.DrawLatex(0.600,0.925,"900 MeV < E_{inst}^{reco.} < 1000 MeV");
  if(name.Contains("1000")) tt.DrawLatex(0.600,0.925,"1000 MeV < E_{inst}^{reco.} < 1100 MeV");
}

void SetTitleFormat(TH1 * hh){
  hh->SetTitle(" ");
  hh->GetYaxis()->CenterTitle();
  // Bold text
  hh->GetYaxis()->SetTitleFont(22);
  hh->GetYaxis()->SetTitleSize(0.04);
  hh->GetYaxis()->SetTitleOffset(0.99);
  hh->GetXaxis()->CenterTitle();
  hh->GetXaxis()->SetTitleFont(22);
  hh->GetXaxis()->SetTitleSize(0.04);
  hh->GetXaxis()->SetTitleOffset(0.99);
  //hh->GetXaxis()->SetLabelOffset(999);
  //hh->GetXaxis()->SetLabelSize(0);
}

void SetTitleFormat(TH2 * hh){
  hh->SetTitle(" ");
  hh->GetYaxis()->CenterTitle();
  // Bold text
  hh->GetYaxis()->SetTitleFont(22);
  hh->GetYaxis()->SetTitleSize(0.04);
  hh->GetYaxis()->SetTitleOffset(1.35);
  hh->GetXaxis()->CenterTitle();
  hh->GetXaxis()->SetTitleFont(22);
  hh->GetXaxis()->SetTitleSize(0.04);
  hh->GetXaxis()->SetTitleOffset(0.99);
  //hh->GetXaxis()->SetLabelOffset(999);
  //hh->GetXaxis()->SetLabelSize(0);
}

TGraphErrors * ConvertTH1DtoGraph(TH1D * hist){

    int n = hist->GetNbinsX(); 

    double x[n], y[n], ex[n], ey[n]; // define arrays to store the x and y values, as well as the error bars

    for (int i=1; i<=n; i++) {
        if(i==7) x[i-1] = hist->GetBinLowEdge(i);
        else if (i==13) x[i-1] = hist->GetBinLowEdge(i)+50.0;
        else x[i-1] = hist->GetBinCenter(i); // get the bin center for x
        y[i-1] = hist->GetBinContent(i); // get the bin content for y
        ex[i-1] = hist->GetBinWidth(i)*0.0; // set the error bar for x to half the bin width
        ey[i-1] = hist->GetBinError(i);//
    }

    TGraphErrors *graph = new TGraphErrors(n,x,y,ex,ey);

    Int_t npoints = graph->GetN();
    for (Int_t i=0; i<npoints; i++) {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        if (y == 0) {
            graph->RemovePoint(i);
            npoints--;
            i--;
        }
    }
    

    return graph;
}
