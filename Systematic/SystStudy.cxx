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

using namespace std;
TMatrixD GetCovMatrix(TH1D* h1d_cv, TH1D* h1d_upper, TH1D* h1d_lower);
void SetTitleFormat(TH1 * hh);

int main(int argc, char * argv[])
{ 
  // Open file
  const TString finName_CV = "input/outana_CV.root";
  
  TFile *file_CV = TFile::Open(finName_CV);

  if(!file_CV->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // Open file
  const TString finName_Upper = "input/outana_IniE_upper25.root";//outana_IniE_upperLimit.root";
  //const TString finName_Upper = "input/outana_IniE_upperLimit.root";
  
  TFile *file_Upper = TFile::Open(finName_Upper);

  if(!file_Upper->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // Open file
  const TString finName_Lower = "input/outana_IniE_lower25.root";//outana_IniE_lowerLimit.root";
  //const TString finName_Lower = "input/outana_IniE_lowerLimit.root";
  
  TFile *file_Lower = TFile::Open(finName_Lower);

  if(!file_Lower->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // Get Histograms
  TH1D * h1d_cv = (TH1D*) file_CV->Get("mc/i027hUnFoldedBeamIncidentHist_xsec;1");
  TH1D * h1d_upper = (TH1D*) file_Upper->Get("mc/i027hUnFoldedBeamIncidentHist_xsec;1");
  TH1D * h1d_lower = (TH1D*) file_Lower->Get("mc/i027hUnFoldedBeamIncidentHist_xsec;1");
  
  TMatrixD CovMatrix = GetCovMatrix(h1d_cv, h1d_upper, h1d_lower);

  cout << "CVM: " << endl;
  CovMatrix.Print();


  TLatex tt;
  tt.SetNDC();

  TCanvas * ctmp = new TCanvas("ctmp", "", 1200, 800);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlackBody);//only 56 available

  TH2D * hCVM = new TH2D("hCVM", "", 9, 10, 19, 9, 10, 19);
  hCVM->GetYaxis()->SetTitle("T_{#pi^{+}} Bin");
  hCVM->GetXaxis()->SetTitle("T_{#pi^{+}} Bin");
  SetTitleFormat(hCVM);


  // Correlation
  //hCVM->SetMaximum(+1.); 
  //hCVM->SetMinimum(-1.);
  //CVM
  hCVM->SetMaximum(+3.); 
  hCVM->SetMinimum(-3.);


  for(Int_t ix=10; ix<=19; ix++){
    for(Int_t jx=10; jx<=19; jx++){
      double content = CovMatrix[ix-1][jx-1]; 
      hCVM->SetBinContent(ix-9,jx-9,content);  
      cout << "19 ix: " << ix << "jx: " << jx  << "CVM: " <<  CovMatrix[ix-1][jx-1] << endl; 
    }
  }
  hCVM->Draw("colz");
  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");

  ctmp->Print("output/CVM.pdf");

  // Example
  TCanvas * cdemo = new TCanvas("cdemo", "", 1200, 800);
  TH1D * hCV = new TH1D("hCV", "", 3, 0, 3);
  TH1D * hU1 = new TH1D("hU1", "", 3, 0, 3);
  TH1D * hU2 = new TH1D("hU2", "", 3, 0, 3);
  
  hCV->GetYaxis()->SetTitle("Candidates");
  hCV->GetXaxis()->SetTitle("Variable of Interest");
  hCV->SetMaximum(5.0); 
  hCV->SetMinimum(0.0);
  
  hCV->SetBinContent(1,1); hCV->SetBinContent(2,3); hCV->SetBinContent(3,2);
  hCV->SetLineColor(kBlue);
  hCV->SetLineWidth(2);

  hU1->SetBinContent(1,1.5); hU1->SetBinContent(2,2.5); hU1->SetBinContent(3,1.5);
  hU1->SetLineColor(kGreen+3);
  hU1->SetLineStyle(kDashed);
  hU1->SetLineWidth(2);

  hU2->SetBinContent(1,0.5); hU2->SetBinContent(2,4); hU2->SetBinContent(3,2.5);
  hU2->SetLineColor(kRed-3);
  hU2->SetLineStyle(10);
  hU2->SetLineWidth(2);

  SetTitleFormat(hCV);
  hCV->Draw("HIST");
  hU1->Draw("HIST SAMES");
  hU2->Draw("HIST SAMES");


  TLegend * legend = new TLegend(0.68, 0.68, 0.89, 0.88);
  legend->AddEntry(hCV, "Central Value", "l");
  legend->AddEntry(hU1, "Universe 1", "l");
  legend->AddEntry(hU2, "Universe 2", "l");
  legend->SetBorderSize(0);
  legend->Draw("same");

  cdemo->Print("output/Demo.pdf");


  // Close files
  file_CV->Close();
  file_Upper->Close();
  file_Lower->Close();
  
}

// Get the covariance matirx and extract the error bar
TMatrixD GetCovMatrix(TH1D* h1d_cv, TH1D* h1d_upper, TH1D* h1d_lower){
  
  // Get the bin number
  const Int_t x0 = h1d_cv->GetXaxis()->GetFirst();
  const Int_t x1 = h1d_cv->GetXaxis()->GetLast();
  // Declare the CVM
  TMatrixD CovMatrix(x1,x1);

  vector<double> xCV_vect, xUp_vect, xLow_vect, xbar_vect;
  vector<vector<double>> Universe;

  // Loop over all bins and get bin content
  for(Int_t ix=x0; ix<=x1; ix++){
    double bin_cv = 0.0, bin_upper = 0.0, bin_lower = 0.0, xbar = 0.0;
    if(h1d_cv->GetBinContent(ix)!= 0) bin_cv = h1d_cv->GetBinContent(ix);
    if(h1d_upper->GetBinContent(ix)!= 0) bin_upper = h1d_upper->GetBinContent(ix);
    if(h1d_lower->GetBinContent(ix)!= 0) bin_lower = h1d_lower->GetBinContent(ix);

    xbar = (bin_cv + bin_upper + bin_lower)/3.0;

    cout << "ix: " << ix << endl;
    cout << "bin_cv: " << bin_cv << endl;
    cout << "bin_upper: " << bin_upper << endl;
    cout << "bin_lower: " << bin_lower << endl;

    xCV_vect.push_back(bin_cv);
    xUp_vect.push_back(bin_upper);
    xLow_vect.push_back(bin_lower);
    xbar_vect.push_back(xbar);
  }

  Universe.push_back(xCV_vect);
  Universe.push_back(xUp_vect);
  Universe.push_back(xLow_vect);

  for(Int_t ix=x0; ix<=x1; ix++){
    for(Int_t jx=x0; jx<=x1; jx++){
      double sigmax = 0.0, sigmay = 0.0;
      double CVM_ele = 0.0;
      
      /*CovMatrix[ix-1][jx-1] = ((xCV_vect[ix]-xbar_vect[ix])*(xCV_vect[jx]-xbar_vect[jx])+(xUp_vect[ix]-xbar_vect[ix])*(xUp_vect[jx]-xbar_vect[jx])+(xLow_vect[ix]-xbar_vect[ix])*(xLow_vect[jx]-xbar_vect[jx]))/3.0;
      sigmax = ((xCV_vect[ix]-xbar_vect[ix])*(xCV_vect[ix]-xbar_vect[ix])+(xUp_vect[ix]-xbar_vect[ix])*(xUp_vect[ix]-xbar_vect[ix])+(xLow_vect[ix]-xbar_vect[ix])*(xLow_vect[ix]-xbar_vect[ix]))/3.0;
      sigmay = ((xCV_vect[jx]-xbar_vect[jx])*(xCV_vect[jx]-xbar_vect[jx])+(xUp_vect[jx]-xbar_vect[jx])*(xUp_vect[jx]-xbar_vect[jx])+(xLow_vect[jx]-xbar_vect[jx])*(xLow_vect[jx]-xbar_vect[jx]))/3.0;
      CovMatrix[ix-1][jx-1] = CovMatrix[ix-1][jx-1]/(sqrt(sigmax)*sqrt(sigmay));
      */
      int ent = 0;
      for(auto uni : Universe){
        ent++;
        double size = Universe.size();

        cout << "ix: " << ix << " jx: " << jx << endl;
        cout << "ent: " << ent << endl;
        cout << "uni[ix]: " << uni[ix] << "xbar_vect[ix]: " << xbar_vect[ix] << endl;
        cout << "uni[jx]: " << uni[jx] << "xbar_vect[jx]: " << xbar_vect[jx] << endl;

        CVM_ele += (uni[ix]-xbar_vect[ix])*(uni[jx]-xbar_vect[jx])/size;
        sigmax += (uni[ix]-xbar_vect[ix])*(uni[ix]-xbar_vect[ix])/size;
        sigmay += (uni[jx]-xbar_vect[jx])*(uni[jx]-xbar_vect[jx])/size;
      }
      // Correlation matrix
      //CovMatrix[ix-1][jx-1] = CVM_ele/(sqrt(sigmax)*sqrt(sigmay));
      //cout << "ix: " << ix << " jx: " << jx << endl;
      //cout << "sigmax: " << sigmax << " sigmay: " << sigmay << endl;
      // CVM
      CovMatrix[ix-1][jx-1] = CVM_ele;
    }
  }


  // Clear the vectors
  xCV_vect.clear();
  xUp_vect.clear();
  xLow_vect.clear();
  Universe.clear();
  xbar_vect.clear();

  return CovMatrix;

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