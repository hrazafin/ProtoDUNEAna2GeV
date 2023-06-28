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

void DrawEffResponsePlots(TH2D * h2d, TH1D * h1d, TString tag, bool rebin, double overall);
double GetOverallEff(TH1D * h1d); // This function is incorrect


int main(int argc, char * argv[])
{
    const TString finName = "input/outana.root"; // Eff and response

    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    // ========== Efficiency and Response Matrix ========== //

    TH2D * h2d_ini = (TH2D*) file->Get("mc/i904response_Ini;1");
    TH1D * h1d_ini = (TH1D*) file->Get("mc/i912eff_Ini;1");
    TString name_ini = h2d_ini->GetName();
    //DrawEffResponsePlots(h2d_ini, h1d_ini, name_ini, false, 0.499024);
    // 0410
    DrawEffResponsePlots(h2d_ini, h1d_ini, name_ini, false, 0.433599);

    TH2D * h2d_beamint = (TH2D*) file->Get("mc/i902response_BeamInt;1");
    TH1D * h1d_beamint = (TH1D*) file->Get("mc/i913eff_BeamInt;1");
    TString name_beamint = h2d_beamint->GetName();
    //DrawEffResponsePlots(h2d_beamint, h1d_beamint, name_beamint, true, 0.499024);
    // 0410
    DrawEffResponsePlots(h2d_beamint, h1d_beamint, name_beamint, true, 0.433599);

    TH2D * h2d_int = (TH2D*) file->Get("mc/i901response_Int;1");
    TH1D * h1d_int = (TH1D*) file->Get("mc/i910eff_Int;1");
    TString name_int = h2d_int->GetName();
    //DrawEffResponsePlots(h2d_int, h1d_int, name_int, false, 0.0482581);
    // 0410
    DrawEffResponsePlots(h2d_int, h1d_int, name_int, false, 0.0463485);


    // Pi0 related
    TH2D * h2d_Pi0KE = (TH2D*) file->Get("mc/i905response_pi0KE;1");
    TH1D * h1d_Pi0KE = (TH1D*) file->Get("mc/i914eff_Pi0KE;1");
    TString name_Pi0KE = h2d_Pi0KE->GetName();
    //DrawEffResponsePlots(h2d_Pi0KE, h1d_Pi0KE, name_Pi0KE, false, 0.059414);
    // 0410
    DrawEffResponsePlots(h2d_Pi0KE, h1d_Pi0KE, name_Pi0KE, false, 0.0534021);

    TH2D * h2d_Pi0CosTheta = (TH2D*) file->Get("mc/i906response_pi0CosTheta;1");
    TH1D * h1d_Pi0CosTheta = (TH1D*) file->Get("mc/i915eff_Pi0CosTheta;1");
    TString name_Pi0CosTheta = h2d_Pi0CosTheta->GetName();
    //DrawEffResponsePlots(h2d_Pi0CosTheta, h1d_Pi0CosTheta, name_Pi0CosTheta, false, 0.059414);
    // 0410
    DrawEffResponsePlots(h2d_Pi0CosTheta, h1d_Pi0CosTheta, name_Pi0CosTheta, false, 0.0534021);

    TH2D * h2d_Pi0Theta = (TH2D*) file->Get("mc/i907response_pi0Theta;1");
    TH1D * h1d_Pi0Theta = (TH1D*) file->Get("mc/i916eff_Pi0Theta;1");
    TString name_Pi0Theta = h2d_Pi0Theta->GetName();
    //DrawEffResponsePlots(h2d_Pi0Theta, h1d_Pi0Theta, name_Pi0Theta, false, 0.059414);
    // 0410
    DrawEffResponsePlots(h2d_Pi0Theta, h1d_Pi0Theta, name_Pi0Theta, false, 0.0534021);

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
    //gStyle->SetPalette(kBlackBody);//only 56 available

    TLatex tt;
    tt.SetNDC();

    TCanvas * c11 = new TCanvas(Form("c11%s",tag.Data()), "", 1200, 800);
    if(rebin) h2d->Rebin2D(50,50);
    h2d->SetTitle(" ");
    h2d->GetYaxis()->CenterTitle();
    h2d->GetYaxis()->SetTitleFont(22);
    h2d->GetYaxis()->SetTitleSize(0.05);
    h2d->GetYaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->CenterTitle();
    h2d->GetXaxis()->SetTitleFont(22);
    h2d->GetXaxis()->SetTitleSize(0.05);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    if(tag.Contains("pi0KE")){
      h2d->GetYaxis()->SetTitle("Reco T_{#pi^{0}} (MeV) ");
      h2d->GetXaxis()->SetTitle("True T_{#pi^{0}} (MeV) ");
    }
    else if(tag.Contains("pi0Theta")){
      h2d->GetYaxis()->SetTitle("Reco #theta_{#pi^{0}} (deg.)");
      h2d->GetXaxis()->SetTitle("True #theta_{#pi^{0}} (deg.)");
    }
    else if(tag.Contains("pi0CosTheta")) {
      h2d->GetYaxis()->SetTitle("Reco cos#theta_{#pi^{0}} ");
      h2d->GetXaxis()->SetTitle("True cos#theta_{#pi^{0}} ");
    }
    else{
      h2d->GetYaxis()->SetTitle("Reco T_{#pi^{+}} (MeV) ");
      h2d->GetXaxis()->SetTitle("True T_{#pi^{+}} (MeV) ");
    }
    h2d->Draw("box");
    //h2d->Draw("colz");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    if(!tag.Contains("pi0")) tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
    else tt.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
    tt.DrawLatex(0.155,0.825,"#bf{#it{Preliminary}}");
    tt.DrawLatex(0.155,0.775,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");

    c11->Print("output/"+tag+"_response.pdf");

    TCanvas * c12 = new TCanvas(Form("c12%s",tag.Data()), "", 1200, 800);
    h1d->SetTitle("");
    h1d->GetXaxis()->SetTitle("True T_{#pi^{+}} (MeV)");
    if(tag.Contains("pi0KE")) h1d->GetXaxis()->SetTitle("True T_{#pi^{0}} (MeV)");
    if(tag.Contains("pi0Theta")) h1d->GetXaxis()->SetTitle("True #theta_{#pi^{0}} (deg.)");
    if(tag.Contains("pi0CosTheta")) h1d->GetXaxis()->SetTitle("True cos#theta_{#pi^{0}} ");
    h1d->GetYaxis()->SetTitle("Efficiency");
    
    if(tag.Contains("Ini") || tag.Contains("BeamInt")) h1d->SetMaximum(h1d->GetMaximum()*1.3);

    h1d->SetMarkerColor(kBlack);
    h1d->SetLineColor(kBlack);
    h1d->SetMarkerStyle(kFullSquare);
    h1d->GetYaxis()->CenterTitle();
    h1d->GetYaxis()->SetTitleFont(22);
    h1d->GetYaxis()->SetTitleSize(0.05);
    h1d->GetYaxis()->SetTitleOffset(0.9);
    h1d->GetXaxis()->CenterTitle();
    h1d->GetXaxis()->SetTitleFont(22);
    h1d->GetXaxis()->SetTitleSize(0.05);
    h1d->GetXaxis()->SetTitleOffset(0.9);

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

    auto lg = new TLegend(0.15,0.67,0.4,0.77);
    if(tag.Contains("pi0Theta")) lg = new TLegend(0.6,0.67,0.85,0.77);
    //if(tag.Contains("Pi0CosTheta")) new TLegend(0.15,0.58,0.4,0.68);

    lg->AddEntry(line,"Overall Efficiency","l");
    lg->SetBorderSize(0);
    lg->Draw("sames");



    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    if(!tag.Contains("pi0")) tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
    else tt.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
    if(!tag.Contains("pi0Theta")){
      tt.DrawLatex(0.155,0.825,"#bf{#it{Preliminary}}");
      tt.DrawLatex(0.155,0.775,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");
    }
    else{
      tt.DrawLatex(0.675,0.825,"#bf{#it{Preliminary}}");
      tt.DrawLatex(0.605,0.775,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");
    }

    c12->Print("output/"+tag+"_eff.pdf");

}
