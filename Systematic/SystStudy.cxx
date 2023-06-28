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
void DrawErrorSum(TH1D * hCV_Total,TH1D * hMCXS,TH1D * hInstP,TH1D * hEloss,TH1D * htrkLen,TH1D * hstats,TH1D * htotalerr,THStack *hstckErrorPercent,double sys_err_MCXS[],double sys_err_InstP[],double sys_err_Eloss[],double sys_err_trkLen[],double total_err[], const TString tag);

int main(int argc, char * argv[])
{ 
  // Open file
  const TString finName_CV = "input/outXShists_CV.root";
  
  TFile *file_CV = TFile::Open(finName_CV);

  if(!file_CV->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // Open file
  const TString finName_Upper = "input/outXShists_trkLenhigh.root";
  
  TFile *file_Upper = TFile::Open(finName_Upper);

  if(!file_Upper->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // Open file
  const TString finName_Lower = "input/outXShists_trkLenlow.root";
  
  TFile *file_Lower = TFile::Open(finName_Lower);

  if(!file_Lower->IsOpen()){
    cout << "file file not open" << endl;
    exit(1);
  }

  // Get Histograms
  TH1D * h1d_cv = (TH1D*) file_CV->Get("data/data_TotalCEXxsec;1");
  TH1D * h1d_upper = (TH1D*) file_Upper->Get("data/data_TotalCEXxsec;1");
  TH1D * h1d_lower = (TH1D*) file_Lower->Get("data/data_TotalCEXxsec;1");

  //TH1D * h1d_cv = (TH1D*) file_CV->Get("data/data_DiffCEXxsec_Pi0KE;1");
  //TH1D * h1d_upper = (TH1D*) file_Upper->Get("data/data_DiffCEXxsec_Pi0KE;1");
  //TH1D * h1d_lower = (TH1D*) file_Lower->Get("data/data_DiffCEXxsec_Pi0KE;1");
  
  //TH1D * h1d_cv = (TH1D*) file_CV->Get("data/data_DiffCEXxsec_Pi0Theta;1");
  //TH1D * h1d_upper = (TH1D*) file_Upper->Get("data/data_DiffCEXxsec_Pi0Theta;1");
  //TH1D * h1d_lower = (TH1D*) file_Lower->Get("data/data_DiffCEXxsec_Pi0Theta;1");
   
  TMatrixD CovMatrix = GetCovMatrix(h1d_cv, h1d_upper, h1d_lower);

  cout << "CVM: " << endl;
  CovMatrix.Print();


  TLatex tt;
  tt.SetNDC();

  TCanvas * ctmp = new TCanvas("ctmp", "", 1200, 800);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlackBody);//only 56 available

  TH2D * hCVM = new TH2D("hCVM", "", 20, 1, 21, 20, 1, 21);
  //TH2D * hCVM = new TH2D("hCVM", "", 10, 1, 11, 10, 1, 11);

  hCVM->GetYaxis()->SetTitle("T_{#pi^{+}} Bin");
  hCVM->GetXaxis()->SetTitle("T_{#pi^{+}} Bin");
  SetTitleFormat(hCVM);


  // Correlation
  //hCVM->SetMaximum(+1.); 
  //hCVM->SetMinimum(-1.);
  //CVM
  hCVM->SetMaximum(+100.); 
  hCVM->SetMinimum(-100.);


  /*for(Int_t ix=9; ix<=19; ix++){
    for(Int_t jx=9; jx<=19; jx++){
      double content = CovMatrix[ix-1][jx-1]; 
      hCVM->SetBinContent(ix-8,jx-8,sqrt(content)); 
      cout << "content: " << content << endl;
      cout << "19 ix: " << ix << "jx: " << jx  << "CVM: " <<  CovMatrix[ix-1][jx-1] << endl; 
    }
  }*/

  for(Int_t ix=1; ix<=20; ix++){
    for(Int_t jx=1; jx<=20; jx++){
      double content = CovMatrix[ix-1][jx-1]; 
      hCVM->SetBinContent(ix,jx,sqrt(content)); 
      cout << "content: " << content << endl;
      cout << "19 ix: " << ix << "jx: " << jx  << "CVM: " <<  CovMatrix[ix-1][jx-1] << endl; 
    }
  }

  hCVM->Draw("colz");

  for (int i = 1; i <= hCVM->GetNbinsX(); i++) {
    for (int j = 1; j <= hCVM->GetNbinsY(); j++) {
      if (i == j) {
        float x = hCVM->GetXaxis()->GetBinCenter(i);
        float y = hCVM->GetYaxis()->GetBinCenter(j);
        float val = hCVM->GetBinContent(i, j);
        TLatex *tex = new TLatex(x, y, Form("%.2f", val)); // create a TLatex object with the number of entries
        cout << "i: " << i << "i syst: " << val << endl;
        tex->SetTextAlign(22); // set the text alignment to center
        tex->SetTextSize(0.02); // set the text size
        tex->Draw(); // draw the TLatex object on top of the histogram
      }
    }
  }

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

  tt.DrawLatex(0.755,0.925,"Toy Example");

  cdemo->Print("output/Demo.pdf");


  TCanvas * cSyst_Total = new TCanvas("cSyst_Total", "", 1200, 800);

  TH1D * hCV_Total = (TH1D*) file_CV->Get("data/data_TotalCEXxsec;1");
  TH1D * hMCXS = new TH1D("hMCXS", "", 20, 0, 1000);
  TH1D * hInstP = new TH1D("hInstP", "", 20, 0, 1000);
  TH1D * hEloss = new TH1D("hEloss", "", 20, 0, 1000);
  TH1D * htrkLen = new TH1D("htrkLen", "", 20, 0, 1000);
  TH1D * hstats = new TH1D("hstats", "", 20, 0, 1000);
  TH1D * htotalerr = new TH1D("htotalerr", "", 20, 0, 1000);
  THStack * hstckErrorPercent = new THStack("hstckErrorPercent","");

  // Input errors
  double sys_err_MCXS[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 6.72699,8.80628,10.5353,9.98448,10.75,11.5913,13.7413,12.561,15.4875,15.8728,0};
  double sys_err_InstP[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 34.7672,23.1971,12.3552,13.5464,7.34627,3.70384,3.04849,6.46126,14.0317,5.02415,0};
  double sys_err_Eloss[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1.94672,1.76384,1.33447,1.32445,2.53843,1.24545,0.72337,0.880737,0.634467,6.1957,0};
  double sys_err_trkLen[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 2.77684,4.45628,3.48969,4.95647,7.30045,5.99141,1.08109,3.0805,2.58046,5.76954,0};
  double total_err[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 39.1958,29.6274,24.4663,24.1578,21.1504,18.3059,18.6177,18.5062,28.3507,39.0797,0};

  DrawErrorSum(hCV_Total,hMCXS,hInstP,hEloss,htrkLen,hstats,htotalerr,hstckErrorPercent,sys_err_MCXS,sys_err_InstP,sys_err_Eloss,sys_err_trkLen,total_err,"Total");
  
  cSyst_Total->Print("output/1GeV_Pion_TotalCEX_Error.pdf");

  // ============= Pi0 KE ==============//

  TCanvas * cSyst_Pi0KE = new TCanvas("cSyst_Pi0KE", "", 1200, 800);

  TH1D * hCV_Pi0KE = (TH1D*) file_CV->Get("data/data_DiffCEXxsec_Pi0KE;1");
  TH1D * hMCXS_Pi0KE = new TH1D("hMCXS_Pi0KE", "", 8, 0, 800);
  TH1D * hInstP_Pi0KE = new TH1D("hInstP_Pi0KE", "", 8, 0, 800);
  TH1D * hEloss_Pi0KE = new TH1D("hEloss_Pi0KE", "", 8, 0, 800);
  TH1D * htrkLen_Pi0KE = new TH1D("htrkLen_Pi0KE", "", 8, 0, 800);
  TH1D * hstats_Pi0KE = new TH1D("hstats_Pi0KE", "", 8, 0, 800);
  TH1D * htotalerr_Pi0KE = new TH1D("htotalerr_Pi0KE", "", 8, 0, 800);
  THStack * hstckErrorPercent_Pi0KE = new THStack("hstckErrorPercent_Pi0KE","");

  // Input errors
  double sys_err_MCXS_Pi0KE[] = {0.0387427,0.0193442,0.0255634,0.0172559,0.0117158,0.0096578,0.00831253,0.00176299,0};

  double sys_err_InstP_Pi0KE[] = {0.0279626,0.0141663,0.0191867,0.0132152,0.00918648,0.00747897,0.0067952,0.00136478,0};
  double sys_err_Eloss_Pi0KE[] = {0.00172252,0.000872654,0.00118192,0.000814065,0.000565894,0.000460711,0.00041859,8.40718e-05,0};

  double sys_err_trkLen_Pi0KE[] = {0.00675079,0.00374966,0.00388837,0.0045393,0.00195442,0.00162326,0.000257695,0.000122726,0};

  double total_err_Pi0KE[] = {0.176816,0.0443948,0.0417384,0.0329094,0.0230236,0.0202084,0.0167687,0.00508369,0};

  DrawErrorSum(hCV_Pi0KE,hMCXS_Pi0KE,hInstP_Pi0KE,hEloss_Pi0KE,htrkLen_Pi0KE,hstats_Pi0KE,htotalerr_Pi0KE,hstckErrorPercent_Pi0KE,sys_err_MCXS_Pi0KE,sys_err_InstP_Pi0KE,sys_err_Eloss_Pi0KE,sys_err_trkLen_Pi0KE,total_err_Pi0KE,"Pi0KE");
  
  cSyst_Pi0KE->Print("output/Pi0KE_DifferentialCEX_Error.pdf");

  // ============= Pi0 Theta ==============//

  TCanvas * cSyst_Pi0Theta = new TCanvas("cSyst_Pi0Theta", "", 1200, 800);

  TH1D * hCV_Pi0Theta = (TH1D*) file_CV->Get("data/data_DiffCEXxsec_Pi0Theta;1");
  TH1D * hMCXS_Pi0Theta = new TH1D("hMCXS_Pi0Theta", "", 10, 0, 180);
  TH1D * hInstP_Pi0Theta = new TH1D("hInstP_Pi0Theta", "", 10, 0, 180);
  TH1D * hEloss_Pi0Theta = new TH1D("hEloss_Pi0Theta", "", 10, 0, 180);
  TH1D * htrkLen_Pi0Theta = new TH1D("htrkLen_Pi0Theta", "", 10, 0, 180);
  TH1D * hstats_Pi0Theta = new TH1D("hstats_Pi0Theta", "", 10, 0, 180);
  TH1D * htotalerr_Pi0Theta = new TH1D("htotalerr_Pi0Theta", "", 10, 0, 180);
  THStack * hstckErrorPercent_Pi0Theta = new THStack("hstckErrorPercent_Pi0Theta","");

  // Input errors
  double sys_err_MCXS_Pi0Theta[] = {0.0647985,0.118325,0.115649,0.105839,0.0643276,0.0685551,0.0433671,0.0188147,0.049131,0};

  double sys_err_InstP_Pi0Theta[] = {0.0515921,0.0927121,0.0889976,0.0791777,0.0479888,0.0497884,0.0313087,0.0137916,0.0353815,0};

  double sys_err_Eloss_Pi0Theta[] = {0.00317811,0.00571114,0.00548232,0.00487741,0.00295615,0.00306701,0.00192864,0.000849572,0.00217953,0};

  double sys_err_trkLen_Pi0Theta[] = {0.00948149,0.0132473,0.0338935,0.0210826,0.00631149,0.021305,0.0273017,0.00295538,0.0020769,0};

  double total_err_Pi0Theta[] = {0.142981,0.237759,0.247947,0.253884,0.173279,0.238311,0.217488,0.127014,0.355222,0};

  DrawErrorSum(hCV_Pi0Theta,hMCXS_Pi0Theta,hInstP_Pi0Theta,hEloss_Pi0Theta,htrkLen_Pi0Theta,hstats_Pi0Theta,htotalerr_Pi0Theta,hstckErrorPercent_Pi0Theta,sys_err_MCXS_Pi0Theta,sys_err_InstP_Pi0Theta,sys_err_Eloss_Pi0Theta,sys_err_trkLen_Pi0Theta,total_err_Pi0Theta,"Pi0Theta");
  
  cSyst_Pi0Theta->Print("output/Pi0Theta_DifferentialCEX_Error.pdf");

  // ============= Pi0 CosTheta ==============//

  TCanvas * cSyst_Pi0CosTheta = new TCanvas("cSyst_Pi0CosTheta", "", 1200, 800);

  TH1D * hCV_Pi0CosTheta = (TH1D*) file_CV->Get("data/data_DiffCEXxsec_Pi0CosTheta;1");
  TH1D * hMCXS_Pi0CosTheta = new TH1D("hMCXS_Pi0CosTheta", "", 10, -1, 1);
  TH1D * hInstP_Pi0CosTheta = new TH1D("hInstP_Pi0CosTheta", "", 10, -1, 1);
  TH1D * hEloss_Pi0CosTheta = new TH1D("hEloss_Pi0CosTheta", "", 10, -1, 1);
  TH1D * htrkLen_Pi0CosTheta = new TH1D("htrkLen_Pi0CosTheta", "", 10, -1, 1);
  TH1D * hstats_Pi0CosTheta = new TH1D("hstats_Pi0CosTheta", "", 10, -1, 1);
  TH1D * htotalerr_Pi0CosTheta = new TH1D("htotalerr_Pi0CosTheta", "", 10, -1, 1);
  THStack * hstckErrorPercent_Pi0CosTheta = new THStack("hstckErrorPercent_Pi0CosTheta","");

  // Input errors
  double sys_err_MCXS_Pi0CosTheta[] = {2.74602,1.4075,1.74537,3.04156,4.73819,3.03186,4.73649,7.77741,9.38828,17.2057};

  double sys_err_InstP_Pi0CosTheta[] = {1.97616,1.03686,1.25381,2.21927,3.43312,2.2611,3.52764,5.82954,7.25274,13.5235};

  double sys_err_Eloss_Pi0CosTheta[] = {0.121733,0.0638712,0.0772357,0.136709,0.211483,0.139286,0.217305,0.359104,0.446774,0.833061};

  double sys_err_trkLen_Pi0CosTheta[] = {0.116001,0.590888,1.46974,0.20461,1.67914,0.845945,2.38934,1.90939,2.1647,1.539};

  double total_err_Pi0CosTheta[] = {17.3039,9.59741,12.0657,9.16669,15.8287,9.6101,12.9418,18.3434,20.3443,30.9255};
  DrawErrorSum(hCV_Pi0CosTheta,hMCXS_Pi0CosTheta,hInstP_Pi0CosTheta,hEloss_Pi0CosTheta,htrkLen_Pi0CosTheta,hstats_Pi0CosTheta,htotalerr_Pi0CosTheta,hstckErrorPercent_Pi0CosTheta,sys_err_MCXS_Pi0CosTheta,sys_err_InstP_Pi0CosTheta,sys_err_Eloss_Pi0CosTheta,sys_err_trkLen_Pi0CosTheta,total_err_Pi0CosTheta,"Pi0CosTheta");
  
  cSyst_Pi0CosTheta->Print("output/Pi0CosTheta_DifferentialCEX_Error.pdf");



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
  cout << "x1 total: " << x1 << endl;
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
    if(bin_cv == 9999){
      bin_cv = 0; bin_upper = 0; bin_lower = 0;
    }
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
        //if(ix > 9 && jx > 9) {
        cout << "ix: " << ix-1 << " jx: " << jx-1 << endl;
        cout << "ent: " << ent << endl;
        cout << "uni[ix]: " << uni[ix-1] << "xbar_vect[ix]: " << xbar_vect[ix-1] << endl;
        cout << "uni[jx]: " << uni[jx-1] << "xbar_vect[jx]: " << xbar_vect[jx-1] << endl;
        cout << "(uni[ix]-xbar_vect[ix]): " << (uni[ix-1]-xbar_vect[ix-1]) << endl;
        cout << "(uni[jx]-xbar_vect[jx]): " << (uni[jx-1]-xbar_vect[jx-1]) << endl;
        //}
        CVM_ele += (uni[ix-1]-xbar_vect[ix-1])*(uni[jx-1]-xbar_vect[jx-1])/size;
        sigmax += (uni[ix-1]-xbar_vect[ix-1])*(uni[ix-1]-xbar_vect[ix-1])/size;
        sigmay += (uni[jx-1]-xbar_vect[jx-1])*(uni[jx-1]-xbar_vect[jx-1])/size;
      }
      if(CVM_ele > 999) CVM_ele = 0.0;
      // Correlation matrix
      //CovMatrix[ix-1][jx-1] = CVM_ele/(sqrt(sigmax)*sqrt(sigmay));
      //cout << "ix: " << ix << " jx: " << jx << endl;
      cout << "sigmax: " << sigmax << " sigmay: " << sigmay << endl;
      if(ix == jx) cout << "error" << sqrt(sigmax) << endl; 
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

void DrawErrorSum(TH1D * hCV_Total,TH1D * hMCXS,TH1D * hInstP,TH1D * hEloss,TH1D * htrkLen,TH1D * hstats,TH1D * htotalerr,THStack *hstckErrorPercent,double sys_err_MCXS[],double sys_err_InstP[],double sys_err_Eloss[],double sys_err_trkLen[],double total_err[], const TString tag)
{
  const Int_t x0 = hCV_Total->GetXaxis()->GetFirst();
  const Int_t x1 = hCV_Total->GetXaxis()->GetLast();
  for(Int_t ix=x0; ix<=x1; ix++){
    double value = hCV_Total->GetBinContent(ix);
    double statsError = hCV_Total->GetBinError(ix);
    if(value <= 0) value = 9999999;
    hMCXS->SetBinContent(ix,sys_err_MCXS[ix-1]/value);
    hInstP->SetBinContent(ix,sys_err_InstP[ix-1]/value);
    hEloss->SetBinContent(ix,sys_err_Eloss[ix-1]/value);
    htrkLen->SetBinContent(ix,sys_err_trkLen[ix-1]/value);
    hstats->SetBinContent(ix,statsError/value);
    htotalerr->SetBinContent(ix,total_err[ix-1]/value);
    if(ix ==20) hstats->SetBinContent(ix,0.0);
  }
  htotalerr->SetLineWidth(2);
  htotalerr->SetLineColor(kBlack);
  hstckErrorPercent->Add(htotalerr);
  hstats->SetLineWidth(2);
  hstats->SetLineColor(kBlack);
  hstats->SetLineStyle(kDashed);
  hstckErrorPercent->Add(hstats);
  hMCXS->SetLineWidth(2);
  hMCXS->SetLineColor(kRed);
  hstckErrorPercent->Add(hMCXS);
  hInstP->SetLineWidth(2);
  hInstP->SetLineColor(kBlue);
  hstckErrorPercent->Add(hInstP);
  hEloss->SetLineWidth(2);
  hEloss->SetLineColor(kGreen);
  hstckErrorPercent->Add(hEloss);
  htrkLen->SetLineWidth(2);
  htrkLen->SetLineColor(kOrange);
  hstckErrorPercent->Add(htrkLen);

  hstckErrorPercent->Draw();
  hstckErrorPercent->GetYaxis()->CenterTitle();
  hstckErrorPercent->GetYaxis()->SetTitleFont(22);
  hstckErrorPercent->GetYaxis()->SetTitleSize(0.05);
  hstckErrorPercent->GetYaxis()->SetTitleOffset(0.92);
  hstckErrorPercent->GetXaxis()->CenterTitle();
  hstckErrorPercent->GetXaxis()->SetTitleFont(22);
  hstckErrorPercent->GetXaxis()->SetTitleSize(0.05);
  hstckErrorPercent->GetXaxis()->SetTitleOffset(0.9);
  

  hstckErrorPercent->SetMaximum(1.0);
  if(tag.Contains("Pi0Theta")) hstckErrorPercent->SetMaximum(1.5);
  if(tag.Contains("Pi0CosTheta")) hstckErrorPercent->SetMaximum(1.5);

  hstckErrorPercent->GetYaxis()->SetTitle("Fractional Uncertainty");
  hstckErrorPercent->GetXaxis()->SetTitle("T_{#pi^{+}} (MeV)");
  if(tag.Contains("Pi0KE")) hstckErrorPercent->GetXaxis()->SetTitle("T_{#pi^{0}} (MeV)");
  if(tag.Contains("Pi0Theta")) hstckErrorPercent->GetXaxis()->SetTitle("#theta_{#pi^{0}} (deg.)");
  if(tag.Contains("Pi0CosTheta")) hstckErrorPercent->GetXaxis()->SetTitle("cos#theta_{#pi^{0}}");
  
  hstckErrorPercent->GetXaxis()->SetLimits(450, 950);
  if(tag.Contains("Pi0KE")) hstckErrorPercent->GetXaxis()->SetLimits(0, 800);
  if(tag.Contains("Pi0Theta")) hstckErrorPercent->GetXaxis()->SetLimits(0, 160);
  if(tag.Contains("Pi0CosTheta")) hstckErrorPercent->GetXaxis()->SetLimits(-1, 1);

  hstckErrorPercent->Draw("no stack");

  TLegend * lg = new TLegend(0.45, 0.48, 0.88, 0.88);
  if(tag.Contains("Pi0Theta")) lg = new TLegend(0.15, 0.48, 0.58, 0.88);

  lg->SetBorderSize(0);
  lg->AddEntry(htotalerr, "Total Uncertainty", "l");
  lg->AddEntry(hstats, "Statistical", "l");
  lg->AddEntry(hMCXS, "Cross-Section Model", "l");
  lg->AddEntry(hInstP, "Beam Momentum", "l");
  lg->AddEntry(hEloss, "Upstream Energy Loss", "l");
  lg->AddEntry(htrkLen, "Track Length", "l");
  lg->SetNColumns(2);
  lg->Draw("sames");
  TLatex tt1;
  tt1.SetNDC();
  tt1.SetTextSize(0.035);
  tt1.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  if(tag.Contains("Total")) tt1.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
  else tt1.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
  if(!tag.Contains("Pi0Theta")){
    tt1.DrawLatex(0.145,0.835,"#bf{#it{Preliminary}}");
    tt1.DrawLatex(0.145,0.785,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");
  }
  else{
    tt1.DrawLatex(0.635,0.835,"#bf{#it{Preliminary}}");
    tt1.DrawLatex(0.635,0.785,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");
  }
}
