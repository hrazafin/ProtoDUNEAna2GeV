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
void Get1pi0BckScale(vector<double> &Par, vector<double> &Parerr);
void Get0pi0BckScale(vector<double> &Par, vector<double> &Parerr, const double &low1pi0Scale, const double &high1pi0Scale);
void GetSignalResults(const double &low1pi0Scale, const double &high1pi0Scale, const double &low0pi0Scale, const double &high0pi0Scale, TString tag, const int &break_s1, const int &break_s2);
void DrawOutput(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h2_low, TH1D *h2_high, const double &lowScale, const double &highScale, TString tag);
void DrawOutputSignal(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h0pi0, TH1D *h1pi0_s1, TH1D *h1pi0_s2, TString tag);

void GetSimulFitScale(vector<double> &Par, vector<double> &Parerr, int &bs1, int &bs2);

int main(int argc, char * argv[])
{
    // 1pi0 Scaling factors
    vector<double> bck1pi0Scale, bck1pi0ScaleError;
    Get1pi0BckScale(bck1pi0Scale, bck1pi0ScaleError);

    // 0pi0 Scaling factors
    vector<double> bck0pi0Scale, bck0pi0ScaleError;
    Get0pi0BckScale(bck0pi0Scale, bck0pi0ScaleError, bck1pi0Scale[0], bck1pi0Scale[1]);

    GetSignalResults(bck1pi0Scale[0], bck1pi0Scale[1], bck0pi0Scale[0], bck0pi0Scale[1], "signal", 13, 13);

    //======= Simul Fit ======//
    vector<double> simulFitScale, simulFitScaleError;
    int bs1 = 0, bs2 = 0;
    GetSimulFitScale(simulFitScale, simulFitScaleError, bs1, bs2);

    GetSignalResults(simulFitScale[0], simulFitScale[1], simulFitScale[2], simulFitScale[3], "Sim_signal", bs1, bs2);

/*
    cout << "bck1pi0Scale: " << bck1pi0Scale[0] << " " << bck1pi0Scale[1] << endl;
    cout << "bck1pi0ScaleError: " << bck1pi0ScaleError[0] << " " << bck1pi0ScaleError[1] << endl;

    cout << "bck0pi0Scale: " << bck0pi0Scale[0] << " " << bck0pi0Scale[1] << endl;
    cout << "bck0pi0ScaleError: " << bck0pi0ScaleError[0] << " " << bck0pi0ScaleError[1] << endl;

    cout << "simulFitScale: " << simulFitScale[0] << " " << simulFitScale[1] << " " << simulFitScale[2] << " " << simulFitScale[3] << endl;
    cout << "simulFitScaleError: " << simulFitScaleError[0] << " " << simulFitScaleError[1] << " " << simulFitScaleError[2] << " " << simulFitScaleError[3] << endl;
*/  

    // ==== Scale factors in summary ====== //
    const double SFBin[] = {0, 700, 1000};

    TH1D * hSacleFactor_1pi0 = new TH1D("hSacleFactor_1pi0",";Interacting Energy (MeV); Background Scaling Factors", sizeof(SFBin)/sizeof(double)-1, SFBin); 
    //hSacleFactor_1pi0->SetBinContent(1,2.13924); hSacleFactor_1pi0->SetBinError(1,0.573771);
    //hSacleFactor_1pi0->SetBinContent(2,0.688961); hSacleFactor_1pi0->SetBinError(2,0.161465);
    hSacleFactor_1pi0->SetBinContent(1,2.54864); hSacleFactor_1pi0->SetBinError(1,0.745619);
    hSacleFactor_1pi0->SetBinContent(2,0.622455); hSacleFactor_1pi0->SetBinError(2,0.153662);
    TH1D * hSacleFactor_0pi0 = new TH1D("hSacleFactor_0pi0",";Interacting Energy (MeV); Background Scaling Factors", sizeof(SFBin)/sizeof(double)-1, SFBin); 
    //hSacleFactor_0pi0->SetBinContent(1,1.01219); hSacleFactor_0pi0->SetBinError(1,0.0304626);
    //hSacleFactor_0pi0->SetBinContent(2,0.966249); hSacleFactor_0pi0->SetBinError(2,0.0179907);
    hSacleFactor_0pi0->SetBinContent(1,1.0101); hSacleFactor_0pi0->SetBinError(1,0.0316139);
    hSacleFactor_0pi0->SetBinContent(2,0.968285); hSacleFactor_0pi0->SetBinError(2,0.0177684);
    TLatex tt;
    tt.SetNDC();

    TCanvas * ctmp1 = new TCanvas("ctmp1", "", 1200, 800);
    
    hSacleFactor_1pi0->SetMarkerColor(kRed);
    hSacleFactor_1pi0->SetLineColor(kRed);
    hSacleFactor_1pi0->SetMarkerStyle(kFullCircle);
    hSacleFactor_1pi0->SetMaximum(3.0);
    hSacleFactor_1pi0->SetMinimum(0.0);

    //hSacleFactor_1pi0->Draw("ex0p");
    hSacleFactor_1pi0->Draw("e1");

    hSacleFactor_0pi0->SetMarkerColor(kGreen+3);
    hSacleFactor_0pi0->SetLineColor(kGreen+3);
    hSacleFactor_0pi0->SetMarkerStyle(kFullCircle);
    //hSacleFactor_0pi0->Draw("same ex0p");
    hSacleFactor_0pi0->Draw("same e1");

    TLine *line1 = new TLine(0,1.0,1000,1.0);
    line1->SetLineColor(kBlue);
    line1->SetLineStyle(kDashed);
    line1->SetLineWidth(1);
    line1->Draw("sames");

    TLine *line2 = new TLine(700,0,700,3);
    line2->SetLineColor(kBlack);
    line2->SetLineWidth(1);
    //line2->Draw("sames");

    auto lg = new TLegend(0.65,0.60,0.85,0.88);
    lg->AddEntry(hSacleFactor_1pi0,"sideband 1 (#pi^{0} #geq 1)","p");
    lg->AddEntry(hSacleFactor_0pi0,"sideband 2 (#pi^{0} = 0)","p");
    lg->SetBorderSize(0);
    lg->Draw("sames");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    
    ctmp1->Print("output/BckFactor.png");


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

void Get1pi0BckScale(vector<double> &Par, vector<double> &Parerr)
{
    // Clear vector first
    Par.clear(); Parerr.clear(); 

    // Declare data and MC histograms (placeholder)
    TH1D * h0 = new TH1D("h01pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h1 = new TH1D("h11pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2 = new TH1D("h21pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000);
    // Scale components
    TH1D * h2_low = new TH1D("h2_low1pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2_high = new TH1D("h2_high1pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 

    //const TString finName = "input/outana_1pi0New1117.root";
    const TString finName = "input/outana_Bck1pi0_benchmark.root";
    //const TString finName = "input/outana_Bck1pi0_benchmark_12MeV.root";
    
    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    TH2D * h2d_mc = (TH2D*) file->Get("mc/i076PiPlusInteractingEnergyEvt_COMPOSE;1");
    TH2D * h2d_data = (TH2D*) file->Get("data/i076PiPlusInteractingEnergyEvt_COMPOSE;1");

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

      double ifixed = ientry_mc_sig + ientry_mc_abs + ientry_mc_0pi0 + ientry_mc_beambck;
      double ifixederr = sqrt(pow(ientry_mc_sigerr,2)+pow(ientry_mc_abserr,2)+pow(ientry_mc_0pi0err,2)+pow(ientry_mc_beambckerr,2));
      h1->SetBinContent(ix, ifixed);
      h1->SetBinError(ix, ifixederr);

      double ivaried = ientry_mc_1pi0 + ientry_mc_mulpi0;
      double ivariederr = sqrt(pow(ientry_mc_1pi0err,2)+pow(ientry_mc_mulpi0err,2));

      h2->SetBinContent(ix, ivaried);
      h2->SetBinError(ix, ivariederr);
    }
    
    file->Close();

    TemplateFitter fitter;   
    // Scale MC to data
    h1->Scale(plotScale); h2->Scale(plotScale);

    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(1, 13);
    fitter.Fit();
    double par_low = fitter.GetPar();
    double parerr_low = fitter.GetParError();
    cout << "===== 1pi0 Low Region ==== " << endl;
    cout << par_low << " " << parerr_low << endl;

    fitter.SetFitRange(14, 20);
    fitter.Fit();
    double par_high = fitter.GetPar();
    double parerr_high = fitter.GetParError();
    cout << "===== 1pi0 High Region ==== " << endl;
    cout << par_high << " " << parerr_high << endl;

    for(int ix=0; ix<=nx+1; ix++){
      h2_low->SetBinContent(ix,0);
      h2_high->SetBinContent(ix,0);
      double ientry = h2->GetBinContent(ix);
      if(ix <= 13) h2_low->SetBinContent(ix, ientry);
      if(ix > 13) h2_high->SetBinContent(ix, ientry);
    }
    // Save the optimised scaling factor and error
    Par.push_back(par_low); Par.push_back(par_high); 
    Parerr.push_back(parerr_low); Parerr.push_back(parerr_high); 

    DrawOutput(h0, h1, h2, h2_low, h2_high, par_low, par_high, "1pi0");

}

void Get0pi0BckScale(vector<double> &Par, vector<double> &Parerr, const double &low1pi0Scale, const double &high1pi0Scale)
{
    // Clear vector first
    Par.clear(); Parerr.clear(); 

    // Declare data and MC histograms (placeholder)
    TH1D * h0 = new TH1D("h00pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h1 = new TH1D("h10pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2 = new TH1D("h20pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000);
    // Scale components
    TH1D * h2_low = new TH1D("h2_low0pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2_high = new TH1D("h2_high0pi0",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 

    //const TString finName = "input/outana_0pi0New1117.root";
    const TString finName = "input/outana_Bck0pi0_benchmark.root";
    //const TString finName = "input/outana_Bck0pi0_benchmark_12MeV.root";
    
    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    TH2D * h2d_mc = (TH2D*) file->Get("mc/i076PiPlusInteractingEnergyEvt_COMPOSE;1");
    TH2D * h2d_data = (TH2D*) file->Get("data/i076PiPlusInteractingEnergyEvt_COMPOSE;1");

    int nx = h2d_mc->GetNbinsX();
    int ny = h2d_data->GetNbinsY();

    for(int ix=0; ix<=nx+1; ix++){
      double weight_1pi0sample = low1pi0Scale;
      if(ix>13) weight_1pi0sample = high1pi0Scale;

      // Data
      double ientry_data = h2d_data->GetBinContent(ix, ny);
      double ientry_dataerr = h2d_data->GetBinError(ix, ny);
      // MC component
      double ientry_mc_sig = h2d_mc->GetBinContent(ix, 1);
      double ientry_mc_abs = h2d_mc->GetBinContent(ix, 2);
      double ientry_mc_0pi0 = h2d_mc->GetBinContent(ix, 3);
      double ientry_mc_1pi0 = h2d_mc->GetBinContent(ix, 4) * weight_1pi0sample;
      double ientry_mc_mulpi0 = h2d_mc->GetBinContent(ix, 5) * weight_1pi0sample;
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

      double ifixed = ientry_mc_sig + ientry_mc_abs + ientry_mc_1pi0 + ientry_mc_mulpi0 + ientry_mc_beambck;
      double ifixederr = sqrt(pow(ientry_mc_sigerr,2)+pow(ientry_mc_abserr,2)+pow(ientry_mc_1pi0err,2)+pow(ientry_mc_mulpi0err,2)+pow(ientry_mc_beambckerr,2));

      h1->SetBinContent(ix, ifixed);
      h1->SetBinError(ix, ifixederr);

      double ivaried = ientry_mc_0pi0; 
      double ivariederr = ientry_mc_0pi0err;

      h2->SetBinContent(ix, ivaried);
      h2->SetBinError(ix, ivariederr);
    }
    
    file->Close();

    TemplateFitter fitter;   
    // Scale MC to data
    h1->Scale(plotScale); h2->Scale(plotScale);

    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(1, 13);
    fitter.Fit();
    double par_low = fitter.GetPar();
    double parerr_low = fitter.GetParError();
    cout << "===== 0pi0 Low Region ==== " << endl;
    cout << par_low << " " << parerr_low << endl;

    fitter.SetFitRange(14, 20);
    fitter.Fit();
    double par_high = fitter.GetPar();
    double parerr_high = fitter.GetParError();
    cout << "===== 0pi0 High Region ==== " << endl;
    cout << par_high << " " << parerr_high << endl;

    for(int ix=0; ix<=nx+1; ix++){
      h2_low->SetBinContent(ix,0);
      //h2_mid->SetBinContent(ix,0);
      h2_high->SetBinContent(ix,0);
      double ientry = h2->GetBinContent(ix);
      if(ix <= 13) h2_low->SetBinContent(ix, ientry);
      //if(ix > 13 && ix <= 17) h2_mid->SetBinContent(ix, ientry);
      if(ix > 13) h2_high->SetBinContent(ix, ientry);
    }
    // Save the optimised scaling factor and error
    Par.push_back(par_low); Par.push_back(par_high); 
    Parerr.push_back(parerr_low); Parerr.push_back(parerr_high); 

    DrawOutput(h0, h1, h2, h2_low, h2_high, par_low, par_high, "0pi0");

}

void GetSignalResults(const double &low1pi0Scale, const double &high1pi0Scale, const double &low0pi0Scale, const double &high0pi0Scale, TString tag, const int &break_s1, const int &break_s2)
{
    // Declare data and MC histograms (placeholder)
    TH1D * h0 = new TH1D(Form("h0_%s",tag.Data()),";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h1 = new TH1D(Form("h1_%s",tag.Data()),";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2 = new TH1D(Form("h2_%s",tag.Data()),";Interacting Energy (MeV); Candidates", 20, 0, 1000);

    TH1D * h0pi0 = new TH1D(Form("h0pi0_%s",tag.Data()),";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h1pi0_s1 = new TH1D(Form("h1pi0_s1_%s",tag.Data()),";Interacting Energy (MeV); Candidates", 20, 0, 1000);
    TH1D * h1pi0_s2 = new TH1D(Form("h1pi0_s2_%s",tag.Data()),";Interacting Energy (MeV); Candidates", 20, 0, 1000);
    
    //const TString finName = "input/outana_sigNew1117.root";
    const TString finName = "input/outana_Sig_benchmark.root";
    //const TString finName = "input/outana_Sig_benchmark_12MeV.root";
    
    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    TH2D * h2d_mc = (TH2D*) file->Get("mc/i076PiPlusInteractingEnergyEvt_COMPOSE;1");
    TH2D * h2d_data = (TH2D*) file->Get("data/i076PiPlusInteractingEnergyEvt_COMPOSE;1");

    int nx = h2d_mc->GetNbinsX();
    int ny = h2d_data->GetNbinsY();

    for(int ix=0; ix<=nx+1; ix++){
      double weight_1pi0sample = low1pi0Scale;
      if(ix>break_s1) weight_1pi0sample = high1pi0Scale;

      double weight_0pi0sample = low0pi0Scale;
      if(ix>break_s2) weight_0pi0sample = high0pi0Scale;

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
      // Ori MC
      double ifixed = ientry_mc_sig + ientry_mc_abs + ientry_mc_0pi0  + ientry_mc_1pi0 + ientry_mc_mulpi0 + ientry_mc_beambck;
      double ifixederr = sqrt(pow(ientry_mc_sigerr,2)+pow(ientry_mc_0pi0err,2)+pow(ientry_mc_abserr,2)+pow(ientry_mc_1pi0err,2)+pow(ientry_mc_mulpi0err,2)+pow(ientry_mc_beambckerr,2));
      h1->SetBinContent(ix, ifixed);
      h1->SetBinError(ix, ifixederr);

      // Fit MC
      double ivaried = ientry_mc_sig + ientry_mc_abs + ientry_mc_0pi0 * weight_0pi0sample + ientry_mc_1pi0 * weight_1pi0sample + ientry_mc_mulpi0 * weight_1pi0sample + ientry_mc_beambck;; 
      double ivariederr = ifixederr;
      h2->SetBinContent(ix, ivaried);
      h2->SetBinError(ix, ivariederr);

      cout << "ix: " << ix << " ratio: " << ientry_mc_sig/ivaried << endl;
      //cout << "ivaried: " << ivaried << " ratio: " << ivaried/ivaried << endl;
      //cout << "0pi0: " << ientry_mc_0pi0 * weight_0pi0sample << " ratio: " << (ientry_mc_0pi0 * weight_0pi0sample)/ivaried << endl;
      //cout << "1pi0: " << ientry_mc_1pi0 * weight_1pi0sample + ientry_mc_mulpi0 * weight_1pi0sample  << " ratio: " << (ientry_mc_1pi0 * weight_1pi0sample + ientry_mc_mulpi0 * weight_1pi0sample)/ivaried << endl;
      //cout << "beambck: " << ientry_mc_beambck << " ratio: " << ientry_mc_beambck/ivaried << endl;
      //cout << endl;

      double ih0pi0 = ientry_mc_0pi0 * weight_0pi0sample;
      double ih0pi0err = ientry_mc_0pi0err;
      h0pi0->SetBinContent(ix, ih0pi0);
      
      h0pi0->SetBinError(ix, ih0pi0err);

      double ih1pi0_s1 = ientry_mc_1pi0 * weight_1pi0sample;
      double ih1pi0err_s1 = ientry_mc_1pi0err;
      h1pi0_s1->SetBinContent(ix, ih1pi0_s1);
      h1pi0_s1->SetBinError(ix, ih1pi0err_s1);

      double ih1pi0_s2 = ientry_mc_mulpi0 * weight_1pi0sample;
      double ih1pi0err_s2 = ientry_mc_mulpi0err;
      h1pi0_s2->SetBinContent(ix, ih1pi0_s2);
      h1pi0_s2->SetBinError(ix, ih1pi0err_s2);

    }
    
    file->Close();

    // Scale MC to data
    h1->Scale(plotScale); h2->Scale(plotScale);
    h0pi0->Scale(plotScale); h1pi0_s1->Scale(plotScale); h1pi0_s2->Scale(plotScale);


    DrawOutputSignal(h0, h1, h2, h0pi0, h1pi0_s1, h1pi0_s2, tag.Data());

}

void DrawOutput(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h2_low, TH1D *h2_high, const double &lowScale, const double &highScale, TString tag)
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

    auto lg = new TLegend(0.15,0.7,0.45,0.88);
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

    hmc_fit->SetLineColor(kRed);
    hmc_fit->SetLineWidth(3);
    
    hmc_fit->Draw("hists sames");

    hmc->SetLineColor(kGreen+3);
    hmc->SetLineStyle(9);
    hmc->SetLineWidth(3);
    hmc->Draw("hists sames");

    auto lgfit = new TLegend(0.15,0.60,0.45,0.88);
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

void DrawOutputSignal(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h0pi0, TH1D *h1pi0_s1, TH1D *h1pi0_s2, TString tag)
{
    TLatex tt;
    tt.SetNDC();

    TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 1200, 800);
    gStyle->SetOptStat(0);
    
    h0->SetMarkerStyle(8);
    h0->SetMarkerSize(1);
    h0->SetMarkerColor(kBlack);
    h0->SetLineColor(kBlack);
    h0->SetMaximum(h0->GetMaximum()*1.5);

    h0->Draw("e1");
    TH1D * hmc = (TH1D*)h1->Clone();
    hmc->SetLineColor(kGreen+3);
    hmc->SetLineStyle(9);
    hmc->SetLineWidth(3);

    hmc->Draw("hists sames");

    TH1D * hmc_fit = (TH1D*)h2->Clone();
    hmc_fit->SetLineColor(kRed);
    hmc_fit->SetLineWidth(3);
    hmc_fit->Draw("hists sames");

    h0pi0->SetLineColor(kMagenta-3);
    h0pi0->SetLineWidth(3);
    //h0pi0->Draw("hists sames");

    h1pi0_s1->SetLineColor(kOrange);
    h1pi0_s1->SetLineWidth(3);
    //h1pi0_s1->Draw("hists sames");

    h1pi0_s2->SetLineColor(kBlue);
    h1pi0_s2->SetLineWidth(3);
    //h1pi0_s2->Draw("hists sames");

    auto lg = new TLegend(0.15,0.60,0.45,0.88);
    lg->SetHeader(Form("#chi^{2}_{def}: %.2f, #chi^{2}_{fit}: %.2f", GetChi2(h0,hmc), GetChi2(h0,hmc_fit)));
    lg->AddEntry(h0,"data","lp");
    lg->AddEntry(hmc,"MC Default","l");
    lg->AddEntry(hmc_fit,"MC Fit","l");
    //lg->AddEntry(h0pi0,"MC Fit (0pi0)","l");
    //lg->AddEntry(h1pi0_s1,"MC Fit (1pi0)","l");
    //lg->AddEntry(h1pi0_s2,"MC Fit (mul. pi0)","l");
    lg->SetBorderSize(0);

    lg->Draw("sames");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    c1->Print("output/ori_"+tag+".png");

    TCanvas * c2 = new TCanvas(Form("c2_%s",tag.Data()), "", 1200, 800);

    h0->SetMarkerStyle(25);
    h0->SetMarkerSize(1);
    h0->SetMarkerColor(kBlack);
    h0->SetLineColor(kBlack);
    h0->Draw("e1");
    auto lgfit = new TLegend(0.15,0.7,0.45,0.88);
    lgfit->SetHeader(Form("#chi^{2}: %.2f",GetChi2(h0,hmc_fit)));
    lgfit->AddEntry(h0,"data","lp");
    lgfit->AddEntry(hmc_fit,"MC","l");
    lgfit->SetBorderSize(0);

    lgfit->Draw("sames");

    hmc_fit->Draw("hists sames");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");

    c2->Print("output/fit_"+tag+".png");

}

void GetSimulFitScale(vector<double> &Par, vector<double> &Parerr, int &bs1, int &bs2)
{
    // Clear vector first
    Par.clear(); Parerr.clear(); 

    // Declare data and MC histograms (placeholder)
    TH1D * h0_sample1 = new TH1D("h0_s1",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h1_sample1 = new TH1D("h1_s1",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2_sample1 = new TH1D("h2_s1",";Interacting Energy (MeV); Candidates", 20, 0, 1000);

    TH1D * h0_sample2 = new TH1D("h0_s2",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h1_sample2 = new TH1D("h1_s2",";Interacting Energy (MeV); Candidates", 20, 0, 1000); 
    TH1D * h2_sample2 = new TH1D("h2_s2",";Interacting Energy (MeV); Candidates", 20, 0, 1000);
  
    //const TString finName_1pi0 = "input/outana_1pi0New1117.root";  // Same with benchmark
    const TString finName_1pi0 = "input/outana_Bck1pi0_benchmark.root";  
    //const TString finName_1pi0 = "input/outana_Bck1pi0_benchmark_12MeV.root"; // Hadron Ana Meeting Plots

    TFile *file_1pi0 = TFile::Open(finName_1pi0);   

    if(!file_1pi0->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    TH2D * h2d_mc_1pi0 = (TH2D*) file_1pi0->Get("mc/i076PiPlusInteractingEnergyEvt_COMPOSE;1");
    TH2D * h2d_data_1pi0 = (TH2D*) file_1pi0->Get("data/i076PiPlusInteractingEnergyEvt_COMPOSE;1");

    int nx = h2d_mc_1pi0->GetNbinsX();
    int ny = h2d_data_1pi0->GetNbinsY();

    for(int ix=0; ix<=nx+1; ix++){

      // Data
      double ientry_data = h2d_data_1pi0->GetBinContent(ix, ny);
      double ientry_dataerr = h2d_data_1pi0->GetBinError(ix, ny);
      // MC component
      double ientry_mc_sig = h2d_mc_1pi0->GetBinContent(ix, 1);
      double ientry_mc_abs = h2d_mc_1pi0->GetBinContent(ix, 2);
      double ientry_mc_0pi0 = h2d_mc_1pi0->GetBinContent(ix, 3);
      double ientry_mc_1pi0 = h2d_mc_1pi0->GetBinContent(ix, 4);
      double ientry_mc_mulpi0 = h2d_mc_1pi0->GetBinContent(ix, 5);
      double ientry_mc_beambck = h2d_mc_1pi0->GetBinContent(ix, 6);
      // MC component error
      double ientry_mc_sigerr = h2d_mc_1pi0->GetBinError(ix, 1);
      double ientry_mc_abserr = h2d_mc_1pi0->GetBinError(ix, 2);
      double ientry_mc_0pi0err = h2d_mc_1pi0->GetBinError(ix, 3);
      double ientry_mc_1pi0err = h2d_mc_1pi0->GetBinError(ix, 4);
      double ientry_mc_mulpi0err = h2d_mc_1pi0->GetBinError(ix, 5);
      double ientry_mc_beambckerr = h2d_mc_1pi0->GetBinError(ix, 6);

      h0_sample1->SetBinContent(ix, ientry_data);
      h0_sample1->SetBinError(ix, ientry_dataerr);

      double ifixed = ientry_mc_sig + ientry_mc_abs + ientry_mc_0pi0+ ientry_mc_beambck;
      double ifixederr = sqrt(pow(ientry_mc_sigerr,2)+pow(ientry_mc_abserr,2)+pow(ientry_mc_0pi0err,2)+pow(ientry_mc_beambckerr,2));

      h1_sample1->SetBinContent(ix, ifixed);
      h1_sample1->SetBinError(ix, ifixederr);

      double ivaried = ientry_mc_1pi0 + ientry_mc_mulpi0; 
      double ivariederr = sqrt(pow(ientry_mc_1pi0err,2)+pow(ientry_mc_mulpi0err,2));

      h2_sample1->SetBinContent(ix, ivaried);
      h2_sample1->SetBinError(ix, ivariederr);
    }
    
    file_1pi0->Close();

    //const TString finName_0pi0 = "input/outana_0pi0New1117.root";
    const TString finName_0pi0 = "input/outana_Bck0pi0_benchmark.root";
    //const TString finName_0pi0 = "input/outana_Bck0pi0_benchmark_12MeV.root";
    
    TFile *file_0pi0 = TFile::Open(finName_0pi0);

    if(!file_0pi0->IsOpen()){
      cout << "file file not open" << endl;
      exit(1);
    }

    TH2D * h2d_mc_0pi0 = (TH2D*) file_0pi0->Get("mc/i076PiPlusInteractingEnergyEvt_COMPOSE;1");
    TH2D * h2d_data_0pi0 = (TH2D*) file_0pi0->Get("data/i076PiPlusInteractingEnergyEvt_COMPOSE;1");

    for(int ix=0; ix<=nx+1; ix++){

      // Data
      double ientry_data = h2d_data_0pi0->GetBinContent(ix, ny);
      double ientry_dataerr = h2d_data_0pi0->GetBinError(ix, ny);
      // MC component
      double ientry_mc_sig = h2d_mc_0pi0->GetBinContent(ix, 1);
      double ientry_mc_abs = h2d_mc_0pi0->GetBinContent(ix, 2);
      double ientry_mc_0pi0 = h2d_mc_0pi0->GetBinContent(ix, 3);
      double ientry_mc_1pi0 = h2d_mc_0pi0->GetBinContent(ix, 4);
      double ientry_mc_mulpi0 = h2d_mc_0pi0->GetBinContent(ix, 5);
      double ientry_mc_beambck = h2d_mc_0pi0->GetBinContent(ix, 6);
      // MC component error
      double ientry_mc_sigerr = h2d_mc_0pi0->GetBinError(ix, 1);
      double ientry_mc_abserr = h2d_mc_0pi0->GetBinError(ix, 2);
      double ientry_mc_0pi0err = h2d_mc_0pi0->GetBinError(ix, 3);
      double ientry_mc_1pi0err = h2d_mc_0pi0->GetBinError(ix, 4);
      double ientry_mc_mulpi0err = h2d_mc_0pi0->GetBinError(ix, 5);
      double ientry_mc_beambckerr = h2d_mc_0pi0->GetBinError(ix, 6);

      h0_sample2->SetBinContent(ix, ientry_data);
      h0_sample2->SetBinError(ix, ientry_dataerr);

      double ifixed = ientry_mc_sig + ientry_mc_abs + ientry_mc_1pi0 + ientry_mc_mulpi0 + ientry_mc_beambck;
      double ifixederr = sqrt(pow(ientry_mc_sigerr,2)+pow(ientry_mc_abserr,2)+pow(ientry_mc_1pi0err,2)+pow(ientry_mc_mulpi0err,2)+pow(ientry_mc_beambckerr,2));

      h1_sample2->SetBinContent(ix, ifixed);
      h1_sample2->SetBinError(ix, ifixederr);

      double ivaried = ientry_mc_0pi0; 
      double ivariederr = ientry_mc_0pi0err;

      h2_sample2->SetBinContent(ix, ivaried);
      h2_sample2->SetBinError(ix, ivariederr);
    }
    
    file_0pi0->Close();

    
    TemplateFitter fitter(4,2);   
    // Scale MC to data
    h1_sample1->Scale(plotScale); h2_sample1->Scale(plotScale);
    h1_sample2->Scale(plotScale); h2_sample2->Scale(plotScale);

    // Some tests to see which break is the best
    // Fix sample 2 = 13
    //bs1 = 13; bs2 = 13; // chisq = 11.95;
    //bs1 = 12; bs2 = 13; // chisq = 16.37;
    //bs1 = 11; bs2 = 13; // chisq = 21.87;

    //bs1 = 14; bs2 = 13; // chisq = 15.94;
    //bs1 = 15; bs2 = 13; // chisq = 19.80;

    // Fix sample 1 = 13
    //bs1 = 13; bs2 = 12; // chisq = 12.04;
    //bs1 = 13; bs2 = 11; // chisq = 12.04;
    //bs1 = 13; bs2 = 10; // chisq = 11.82;
    //bs1 = 13; bs2 = 9; // chisq = 12.25;
    //bs1 = 13; bs2 = 8; // chisq = 12.39;
    //bs1 = 13; bs2 = 7; // chisq = 12.03;

    //bs1 = 13; bs2 = 14; // chisq = 12.04;
    //bs1 = 13; bs2 = 15; // chisq = 12.10;
    //bs1 = 13; bs2 = 16; // chisq = 12.11;

    bs1 = 14; bs2 = 13; // chisq = 11.98;
    
    

    fitter.SetSimulHistograms(h0_sample1, h1_sample1, h2_sample1, h0_sample2, h1_sample2, h2_sample2, bs1, bs2);


    fitter.SimulFit();

    vector<double> ParVect = fitter.GetSimulPar();
    vector<double> ParVectError = fitter.GetSimulParError();
    
    cout << "===== Simul Fit Pars ==== " << endl;
    for(unsigned int ii = 0; ii < ParVect.size(); ii++){
      cout << ParVect[ii] << " " << ParVectError[ii] << endl;
    }

    Par = ParVect;
    Parerr = ParVectError;
      
}