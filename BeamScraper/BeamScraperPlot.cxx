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
void DrawEffResponsePlots(TH2D * h2d, TH1D * h1d, TString tag, bool rebin, double overall);
double GetOverallEff(TH1D * h1d); // This function is incorrect


int main(int argc, char * argv[])
{
    
    //const TString finName = "input/outana.root"; // Scaper plots
    //const TString finName = "input/out_smearing_new.root";
    
    //const TString finName = "input/out_InstEstudy.root";
    //const TString finName = "input/out_UpstreamFit.root"; // Eloss

    //const TString finName = "input/out_new1115.root"; // Eff and response
    //const TString finName = "input/out_new1116.root"; // Eff and response
    //const TString finName = "input/out_new1117.root"; // Eff and response
    //const TString finName = "input/out_new1118.root"; // Eff and response
    const TString finName = "input/outana_wholeRange_0123.root"; // Eff and response

    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }
/* == change to newer file
    // ========== Beam Inst P Smearing ========== //

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
*/
/*
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


    // ========== Upstream Energy Loss ========== //

    TH1D * h1d_InstP_700 = (TH1D*) file->Get("mc/j007hUpStreamELoss700MeVAfterScraperCuts;1");
    TH1D * h1d_InstP_800 = (TH1D*) file->Get("mc/j007hUpStreamELoss800MeVAfterScraperCuts;1");
    TH1D * h1d_InstP_900 = (TH1D*) file->Get("mc/j007hUpStreamELoss900MeVAfterScraperCuts;1");
    TH1D * h1d_InstP_1000 = (TH1D*) file->Get("mc/j007hUpStreamELoss1000MeVAfterScraperCuts;1");

    TH1D * hFitEloss = new TH1D("hFitEloss","", 8, 400, 1200); 
    hFitEloss->SetBinContent(3,h1d_InstP_700->GetMean());
    hFitEloss->SetBinError(3,h1d_InstP_700->GetRMS());

    hFitEloss->SetBinContent(4,h1d_InstP_800->GetMean());
    hFitEloss->SetBinError(4,h1d_InstP_800->GetRMS());

    hFitEloss->SetBinContent(5,h1d_InstP_900->GetMean());
    hFitEloss->SetBinError(5,h1d_InstP_900->GetRMS());

    hFitEloss->SetBinContent(6,h1d_InstP_1000->GetMean());
    hFitEloss->SetBinError(6,h1d_InstP_1000->GetRMS());


    TCanvas * ctmp = new TCanvas("ctmp", "", 1200, 800);
    hFitEloss->Fit("pol2");
    hFitEloss->Draw("e1");
    ctmp->Print("output/Eloss.pdf");

*/
    TH2D * h2d_InstE = (TH2D*) file->Get("mc/j008hUpStreamELossAfterSmearingAndWeight;1");

    const int ny = h2d_InstE->GetNbinsY();
    const int nx = h2d_InstE->GetNbinsX();
    const double xmin = h2d_InstE->GetXaxis()->GetBinLowEdge(1);
    const double xmax = h2d_InstE->GetXaxis()->GetBinUpEdge(nx);

    TH1D * hFitElossSW = new TH1D("hFitElossSW",";KE_{Beam  Inst.}^{reco.} (MeV);#mu(KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true}) (MeV)", 16, 400, 1200); 


    for(int iy = 1; iy <= ny-1; iy++){

      TH1D * htmp = new TH1D(Form("tmp%d",iy), "", nx, xmin, xmax);
      for(int ix=0; ix<=nx+1; ix++){
        const double ientry = h2d_InstE->GetBinContent(ix, iy);
        htmp->SetBinContent(ix, ientry);
      }
      TF1 *f1 = FitGausFunc(htmp,-300,300);
      hFitElossSW->SetBinContent(iy+6,f1->GetParameter("Mean"));
      // Sigma band
      //hFitElossSW->SetBinError(iy+6,f1->GetParameter("Sigma"));
      // Stats Error
      hFitElossSW->SetBinError(iy+6,1/sqrt(htmp->Integral()));

    }

    TLatex tt;
    tt.SetNDC();
    gStyle->SetOptStat(0);

    TCanvas * ctmp1 = new TCanvas("ctmp1", "", 1200, 800);
    TF1 *f_pol2 = new TF1("f_pol2", "pol2", 650, 1100);
    hFitElossSW->Fit("f_pol2","0");
    hFitElossSW->SetMarkerColor(kBlack);
    hFitElossSW->SetMarkerStyle(kFullCircle);
    hFitElossSW->Draw("e1");
    f_pol2->Draw("sames C");

    TLine *line = new TLine(400,0,1200,0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->Draw("sames");


    //TLine *line1 = new TLine(400,12.74,1200,12.74);
    TLine *line1 = new TLine(400,11.966,1200,11.966);

    line1->SetLineColor(kBlue);
    line1->SetLineStyle(kDashed);
    line1->SetLineWidth(2);
    line1->Draw("sames");

    TLine *line2 = new TLine(400,2.652,1200,2.652);

    line2->SetLineColor(kGreen+3);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);
    line2->Draw("sames");

    auto lg = new TLegend(0.15,0.55,0.62,0.88);
    lg->AddEntry(hFitElossSW,"#Delta KE_{upstream}","lp");
    lg->AddEntry(f_pol2,"Fitted Poly2","l");
    lg->AddEntry(line1,"Constant Eloss (11.966 MeV) - Before Smearing","l");
    lg->AddEntry(line2,"Constant Eloss (2.652 MeV) - After Smearing","l");

    lg->SetBorderSize(0);
    lg->Draw("sames");


    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    ctmp1->Print("output/ElossSW.pdf");


    TH1D * hLoss = h2d_InstE->ProjectionX();
    TCanvas * ctmp2 = new TCanvas("ctmp2", "", 1200, 800);
    TF1 *f11 = new TF1("f11","gaus",-100,100);    

    cout << "hLoss mean: " << hLoss->GetMean() << endl;
    hLoss->Fit("f11","0");
    hLoss->Draw("hist");
    f11->Draw("sames");
    cout << "f11 mean: " << f11->GetParameter("Mean") << endl;

    ctmp2->Print("output/hLoss.pdf");

    // ========== Efficiency and Response Matrix ========== //

    TH2D * h2d_beamint = (TH2D*) file->Get("mc/i902response_BeamInt;1");
    TH1D * h1d_beamint = (TH1D*) file->Get("mc/i913eff_BeamInt;1");
    TString name_beamint = h2d_beamint->GetName();
    DrawEffResponsePlots(h2d_beamint, h1d_beamint, name_beamint, true, 0.499024);

    TH2D * h2d_int = (TH2D*) file->Get("mc/i901response_Int;1");
    TH1D * h1d_int = (TH1D*) file->Get("mc/i910eff_Int;1");
    TString name_int = h2d_int->GetName();
    DrawEffResponsePlots(h2d_int, h1d_int, name_int, false, 0.0482581);

    TH2D * h2d_ini = (TH2D*) file->Get("mc/i904response_Ini;1");
    TH1D * h1d_ini = (TH1D*) file->Get("mc/i912eff_Ini;1");
    TString name_ini = h2d_ini->GetName();
    DrawEffResponsePlots(h2d_ini, h1d_ini, name_ini, false, 0.499024);

    // Pi0 related
    TH2D * h2d_Pi0KE = (TH2D*) file->Get("mc/i905response_pi0KE;1");
    TH1D * h1d_Pi0KE = (TH1D*) file->Get("mc/i914eff_Pi0KE;1");
    TString name_Pi0KE = h2d_Pi0KE->GetName();
    DrawEffResponsePlots(h2d_Pi0KE, h1d_Pi0KE, name_Pi0KE, false, 0.059414);

    TH2D * h2d_Pi0CosTheta = (TH2D*) file->Get("mc/i906response_pi0CosTheta;1");
    TH1D * h1d_Pi0CosTheta = (TH1D*) file->Get("mc/i915eff_Pi0CosTheta;1");
    TString name_Pi0CosTheta = h2d_Pi0CosTheta->GetName();
    DrawEffResponsePlots(h2d_Pi0CosTheta, h1d_Pi0CosTheta, name_Pi0CosTheta, false, 0.059414);

    TH2D * h2d_Pi0Theta = (TH2D*) file->Get("mc/i907response_pi0Theta;1");
    TH1D * h1d_Pi0Theta = (TH1D*) file->Get("mc/i916eff_Pi0Theta;1");
    TString name_Pi0Theta = h2d_Pi0Theta->GetName();
    DrawEffResponsePlots(h2d_Pi0Theta, h1d_Pi0Theta, name_Pi0Theta, false, 0.059414);


    file->Close();
}
double GetOverallEff(TH1D * h1d){
  
  const int nx = h1d->GetNbinsX();
  double eff = 0;
  double nn = 0;
  for(int ix=0; ix<=nx+1; ix++){
    eff += h1d->GetBinContent(ix);
    if(h1d->GetBinContent(ix)!=0) nn++;
  }
  return eff/nn;
}

void DrawEffResponsePlots(TH2D * h2d, TH1D * h1d, TString tag, bool rebin, double overall){
    gStyle->SetOptStat(0);

    TLatex tt;
    tt.SetNDC();

    TCanvas * c11 = new TCanvas(Form("c11%s",tag.Data()), "", 1200, 800);
    if(rebin) h2d->Rebin2D(50,50);
    h2d->SetTitle(" ");
    h2d->Draw("box");
    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    c11->Print("output/"+tag+"_response.pdf");

    TCanvas * c12 = new TCanvas(Form("c12%s",tag.Data()), "", 1200, 800);
    h1d->SetTitle("");
    h1d->GetXaxis()->SetTitle("True Pion Energy (MeV)");
    if(tag.Contains("pi0Theta")) h1d->GetXaxis()->SetTitle("True Pion Theta (deg.)");
    if(tag.Contains("pi0CosTheta")) h1d->GetXaxis()->SetTitle("True Pion Cos Theta");
    h1d->GetYaxis()->SetTitle("Efficiency");

    h1d->SetMarkerColor(kBlack);
    h1d->SetLineColor(kBlack);
    h1d->SetMarkerStyle(kFullSquare);
    // This is incorrect -- see main Analysis
    //cout << "GetOverallEff: "  << h1d->GetName() << " "<< GetOverallEff(h1d) << endl;
    h1d->Draw("e1");

    TLine *line = new TLine(0,overall,1000,overall);
    if(tag.Contains("pi0Theta")) line = new TLine(0,overall,180,overall);
    if(tag.Contains("pi0CosTheta")) line = new TLine(-1,overall,1,overall);

    line->SetLineColor(kBlue);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line->Draw("sames");

    auto lg = new TLegend(0.15,0.58,0.4,0.68);
    if(tag.Contains("pi0Theta")) lg = new TLegend(0.6,0.58,0.85,0.68);
    //if(tag.Contains("Pi0CosTheta")) new TLegend(0.15,0.58,0.4,0.68);

    lg->AddEntry(line,"Overall Efficiency","l");
    lg->SetBorderSize(0);
    lg->Draw("sames");



    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    c12->Print("output/"+tag+"_eff.pdf");


}

void DoPlots(TH2D* h2d_beam, TH2D* h2d_scraper,TH1D* h1d_beamx, TH1D* h1d_beamy){

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

    Draw1D(h1d_beamx,f1);
    Draw1D(h1d_beamy,f2);

    //DrawOutput(h2d_beam, h2d_scraper, x_mean, y_mean, radius, 1.4);
    DrawOutput(h2d_beam, h2d_scraper, -29.6, 422, 4.8, 1.4);
}

TF1 * FitGausFunc(TH1D * hh, const double &min, const double &max){
  TString tag = hh->GetName();
  TF1 *f1 = new TF1(Form("f1%s",tag.Data()),"gaus",min,max);    
  if (hh->GetEntries() != 0) {
    f1->SetParameters(hh->GetMaximum(), hh->GetMean(), hh->GetRMS() );
    hh->Fit(Form("f1%s",tag.Data()));
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
    lg->AddEntry(el1,"Sungbin's cut","l");
    lg->SetBorderSize(0);
    lg->SetTextFont(22);
    lg->Draw("sames");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    c1->Print("output/"+tagbeam+".pdf");
}   