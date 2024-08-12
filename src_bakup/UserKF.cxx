#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH1D.h"

#include "TLegend.h"
#include "TLegendEntry.h"

#include <iostream>
#include <random>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "../include/VariadicTable.h"

using namespace std;


class UserKF{

public:
  static void SetPar(const double par[], const int npar){
    if(npar!=fNpar-1){ //npar inside TMinuit is only for free par
      printf("npar mismatch! %d %d\n", npar, fNpar); exit(1);
    }
    fLambda = par[0];
    fOptX   = par[1];
    fOptY   = par[2];
    fOptZ   = par[3];
  }

  static double Constraint(){
    //return fOptX + fOptY + fOptZ - 26;
    return 2*fOptX*fOptY*(1-cos(fOptZ)) - 0.134977*0.134977;
  }
 
  static double FullLagrangian(){
    return CoreLikelihood() + fLambda * Constraint();
  }
  
  static void IniCoreMIN(TMinuit * mnt, const double inputl, const int flag){

    //cout << "inputl: " << inputl << endl;
    //cout << "fOptX: " << fOptX << endl;
    //cout << "fOptY: " << fOptY << endl;
    //cout << "fOptZ: " << fOptZ << endl;

    mnt->DefineParameter(0, "lambda", inputl, 1e-2, -1e6, 1e6);//overwrite lambda from input
    mnt->FixParameter(0);  // Fix Lambda

    //cout << "flag: " << flag << endl;

    // Use user previous fit result as initial values
    mnt->DefineParameter(1, "x", fOptX, 1e-2, -1e6, 1e6);
    mnt->DefineParameter(2, "y", fOptY, 1e-2, -1e6, 1e6);
    mnt->DefineParameter(3, "z", fOptZ, 1e-2, -1e6, 1e6);
    //mnt->FixParameter(1);
    //cout << "fOptX: " << fOptX << endl;
    //cout << "fOptY: " << fOptY << endl;
    //cout << "fOptZ: " << fOptZ << endl;
    //cout << "Mass Opt: " << sqrt(2*fOptX*fOptY*(1-cos(fOptZ))) << endl;
    //cout << endl;
  }

  static bool IsConstraintGood(){
    // This is a loose constraint, can be set to 1e-4
    const double eps = 1e-2;
    return (TMath::Abs(Constraint())<eps);
  }

  static void Print(const TString tag){
    // Print the values
    cout << "\n============================ " << tag << " ============================" << endl;
    VariadicTable<double, double, double, double, double, double, double, double, double, double> vt({"Lambda", "fX", "fY", "fZ", "fOptX", "fOptY", "fOptZ", "CoreLikelihood", "Constraint", "Lagrangian"},10);
    vt.addRow(fLambda,fX,fY,fZ,fOptX,fOptY,fOptZ,CoreLikelihood(),TMath::Abs(Constraint()),FullLagrangian());
    vt.print(std::cout);
  }
  static double GetfOptX(){return fOptX;}
  
  static double GetfOptY(){return fOptY;}

  static double GetfOptZ(){return fOptZ;}

  static double GetfX(){return fX;}
  
  static double GetfY(){return fY;}

  static double GetfZ(){return fZ;}
  
  static double GetfLambda(){return fLambda;}
  
  static double GetfNpar(){return fNpar;}
  
  static void SetVars(const double iniVarX, const double iniVarY, const double iniVarZ){
    fX = iniVarX;
    fY = iniVarY;
    fZ = iniVarZ;
  }
  static void SetOptVars(const double iniVarX, const double iniVarY, const double iniVarZ){
    fOptX = iniVarX;
    fOptY = iniVarY;
    fOptZ = iniVarZ;
  }
  static void SetCVM(const vector<double> V){
    int size = V.size();
    int dim = sqrt(size);
    for(int i=0;i<size;i++) {
      int x=floor(i/dim);
      int y=i%dim;
      fCovMatrix[x][y]=V[i];
    }
  }
private:
  static double fLambda;
  static double fOptX;
  static double fOptY;
  static double fOptZ;
  static const int fNpar;

  static double fX;
  static double fY;
  static double fZ;

  static TMatrixD fCovMatrix;

  static double CoreLikelihood(){

    // Get the Covariance matrix (CVM)
    TMatrixD CovMatrixtmp = fCovMatrix;
    // Note function Invert() will change the original matrix, need to use a copy instead
    // Get the inverse of CVM
    TMatrixD CovMatrixInverse = CovMatrixtmp.Invert();
    // Define the dimension of parameter matrix 
    int nparameters = fNpar-1;
    // Define the column vector in matrix form 
    TMatrixD Diff(nparameters,1);
    // Calculate the difference between optimal value and measured value
    Diff[0][0]=fOptX-fX;
    Diff[1][0]=fOptY-fY;
    Diff[2][0]=fOptZ-fZ;
    // Define the row vector in matrix form
    TMatrixD DiffT(1,nparameters);
    DiffT.Transpose(Diff);
    // Calculate the value
    TMatrixD Chi2 = (DiffT*CovMatrixInverse*Diff);
    
    return Chi2[0][0];
  }
};

double UserKF::fLambda = -999;
double UserKF::fOptX = -999;
double UserKF::fOptY = -999;
double UserKF::fOptZ = -999;
const int UserKF::fNpar = 4;

double UserKF::fX = -999;
double UserKF::fY = -999;
double UserKF::fZ = -999;

const int npars = UserKF::GetfNpar() - 1; 
TMatrixD UserKF::fCovMatrix(npars,npars);


vector<double> IniE(const double &bias, const double &sigma1, const double &sigma2, const double &sigma3){

  std::default_random_engine engine; 
  engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
  
  // Gaussian mean followed by stdiv
  std::normal_distribution<double> nd1(3*bias, sigma1); 
  std::normal_distribution<double> nd2(8*bias, sigma2); 
  std::normal_distribution<double> nd3(15*bias, sigma3); 
  // Generate the intial E's value
  double iniX = nd1(engine);
  double iniY = nd2(engine);
  double iniZ = nd3(engine);
  vector<double> vet;
  vet.push_back(iniX);
  vet.push_back(iniY);
  vet.push_back(iniZ);
  return vet;

}

vector<double> GetCVM(const double &sigma1, const double &sigma2, const double &sigma3){
  // Defines the covariance matrix and the variables
  double CVM[9];
  vector<double> CVM_Vet;
  // Set the covariance matrix elements 
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      CVM[3*i+j]= 10E-5; //off diagonal elements have no co-correlation here but a small number here does no effect the algorithm
      if(i==j && (i==0)) {CVM[3*i+j]=TMath::Power(sigma1,2.);}
      if(i==j && (i==1)) {CVM[3*i+j]=TMath::Power(sigma2,2.);}
      if(i==j && (i==2)) {CVM[3*i+j]=TMath::Power(sigma3,2.);}
    }
  }
  for(int i=0;i<9;i++){
    CVM_Vet.push_back(CVM[i]);
  }
  return CVM_Vet;
}

void CoreFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  UserKF::SetPar(par, npars);
  
  value = UserKF::FullLagrangian();
}

void LambdaFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  // Parameters and results are passed in this order
  // Par/UserKF -> CoreMIN -> UserKF -> value
  
  //UserKF::Print("CoreMIN before fit");

  // Second Minimization
  TMinuit * CoreMIN = new TMinuit(3);
  CoreMIN->SetPrintLevel(-1);
  
  CoreMIN->SetFCN(CoreFCN);

  UserKF::IniCoreMIN(CoreMIN, par[0],flag);

  int flagL = CoreMIN->Command("MIGRAD");
  //UserKF::Print("CoreMIN after fit");
  
  int irun = 1;
  const int maxnrun = 2;//no need to try many times, fail always if lambda is bad
  while(flagL!=0){
    printf("CoreMIN bad fit! %d ---- run once more [%d]\n", flagL, irun++); 
    flagL = CoreMIN->Command("MIGRAD");
    //UserKF::Print("CoreMIN after fit");
    if(irun>=maxnrun){
      break;
    }
  }
  
  if(flagL==0){
    value = TMath::Abs(UserKF::Constraint());
  }
  else{
    printf("CoreMIN giving up now... %d\n", flagL);
    value = 1E50;//must be large enough wrt possible FullLagrangian when fit fail
  }

  delete CoreMIN;
}

void SetIniValues(const vector<double> iniVar, const vector<double> CVM){
  int VarsSize = iniVar.size();
  int CVMDim = sqrt(CVM.size());
  if(VarsSize != CVMDim) {
    cout <<  "Variable size is not the same as CVM dimension! " << VarsSize << " " << CVMDim << endl;
    exit(1);
  }
  UserKF::SetVars(iniVar[0],iniVar[1],iniVar[2]);
  UserKF::SetCVM(CVM);
}


bool DoubleMin(const double iniLambda, const double lmin, const double lmax, vector<double> iniVar, vector<double> CVM)
{
  //UserKF::Print("DoubleMin before fit");
  // Set the value of measured quantities and CVM (optimised values are not set yet)
  SetIniValues(iniVar, CVM);

  TMinuit * LambdaMIN = new TMinuit(1);
  LambdaMIN->SetPrintLevel(-1);

  LambdaMIN->SetFCN(LambdaFCN);

  LambdaMIN->DefineParameter(0, "lambda", iniLambda, 1e-2, lmin, lmax);
  
  int flag = LambdaMIN->Command("MIGRAD");
  //UserKF::Print("LambdaMIN after fit");
  
  int irun = 1;
  const int maxnrun=20; //if time permits, the larger the better  
  std::srand(std::time(nullptr)); // Seed the time

  while(flag != 0 || ! UserKF::IsConstraintGood()){
    printf("LambdaMIN bad fit! %d %e ------- run once more! [%d]\n", flag, UserKF::Constraint(), irun++);
    // Random choose lambda between min and max
    //srand((unsigned)time(NULL)); // Seed the time
    
    //cout << "rand(): " << rand() << endl;
    //double finalNum = (lmax - lmin) * ((double)rand()/(double)RAND_MAX) + lmin;
    //cout << "finalNum: " << finalNum << endl;
    //LambdaMIN->DefineParameter(0, "lambda", finalNum, 1e-2, lmin, lmax);
    UserKF::SetOptVars(UserKF::GetfX(),UserKF::GetfY(),UserKF::GetfZ());

    flag = LambdaMIN->Command("MIGRAD");
    //UserKF::Print("LambdaMIN after fit");

    if(irun>=maxnrun){
      break;
    }
  }

  delete LambdaMIN;

  if(flag==0 && UserKF::IsConstraintGood()){
    printf("DoubleMin finishes: it works for %f %f %f\n", UserKF::GetfLambda(), lmin, lmax);
    UserKF::Print("DoubleMin Outputs (successful)");
    return true;
  }
  else{
    printf("DoubleMin finishes: giving up now... %d %f for  %f %f %f\n", flag, UserKF::Constraint(), UserKF::GetfLambda(), lmin, lmax);
    UserKF::Print("DoubleMin Outputs (failed)");
    return false;
  }
}

vector<double> DoKF(const double &LdShowerEnergyRaw, const double &SlShowerEnergyRaw, const double &OpenAngle, vector<double> CVM, bool &GoodFit)
{
  // Declare variables holder
  vector<double> IniVars;
  vector<double> FittedVars;
  IniVars.push_back(LdShowerEnergyRaw);
  IniVars.push_back(SlShowerEnergyRaw);
  IniVars.push_back(OpenAngle); // in radians
  UserKF::SetOptVars(LdShowerEnergyRaw,SlShowerEnergyRaw,OpenAngle);

  GoodFit = DoubleMin(5, -100, 100, IniVars, CVM);
  // Check if good fit
  if(GoodFit){
    // Save the fitted variables
    FittedVars.push_back(UserKF::GetfOptX());
    FittedVars.push_back(UserKF::GetfOptY());
    FittedVars.push_back(UserKF::GetfOptZ());
  }
  else {
    FittedVars.push_back(UserKF::GetfX());
    FittedVars.push_back(UserKF::GetfY());
    FittedVars.push_back(UserKF::GetfZ());
  }

  // Clean the vector
  IniVars.clear();

  return FittedVars; 
}

/*
// Only for toy example
int main()
{
  TH1F *hBefore = new TH1F("hBefore","Energy - Before and After Fitting;Energy (MeV);Candidates",80,0,20);
  TH1F *hAfter = new TH1F("hAfter","Energy - After Fitting;Energy (MeV);Candidates",80,0,20);

  TH1F *hBefore_Sum = new TH1F("hBefore_Sum","Energy Sum - Before Fitting",80,22,42);
  TH1F *hAfter_Sum = new TH1F("hAfter_Sum","Energy Sum - After Fitting",80,22,42);

  double FailEvt = 0;

  int ientries = 1000;

  for(int i = 0; i < ientries; i++){
    // Bias, Sigmas
    double Bias = 1;
    double Sigma1 = 0.4;
    double Sigma2 = 0.8;
    double Sigma3 = 1.2;

    vector<double> vet = IniE(Bias,Sigma1,Sigma2,Sigma3);
    vector<double> CVMvet = GetCVM(Sigma1,Sigma2,Sigma3);

    UserKF::SetOptVars(vet[0],vet[1],vet[2]);
    bool GoodFit = DoubleMin(5, -10, 100, vet, CVMvet);

    if(GoodFit){
      // Fill histograms for good fit
      hBefore->Fill(vet[0]);
      hBefore->Fill(vet[1]);
      hBefore->Fill(vet[2]);
      hBefore_Sum->Fill(vet[0]+vet[1]+vet[2]);
      double finX = UserKF::GetfOptX();
      double finY = UserKF::GetfOptY();
      double finZ = UserKF::GetfOptZ();

      hAfter->Fill(finX);
      hAfter->Fill(finY);
      hAfter->Fill(finZ);
      hAfter_Sum->Fill(finX+finY+finZ);
    }
    else {
      cout << "This event is not well fitted, discard!" << endl;
      FailEvt++;
    }
  }


  cout << "FailEvt: " << FailEvt << endl;

  TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
  auto legend = new TLegend(0.6,0.6,0.88,0.88);
  hBefore->SetMaximum(3330);
  hBefore->SetStats(0);
  hBefore->SetFillStyle(4050);
  hBefore->SetFillColor(24);
  hBefore->SetLineColor(24);
  hBefore->Draw("hist");
  hBefore->SetTitle(" ");
  hBefore->GetYaxis()->CenterTitle();
  hBefore->GetYaxis()->SetTitleFont(22);
  hBefore->GetYaxis()->SetTitleSize(0.05);
  hBefore->GetYaxis()->SetTitleOffset(0.9);
  hBefore->GetXaxis()->CenterTitle();
  hBefore->GetXaxis()->SetTitleFont(22);
  hBefore->GetXaxis()->SetTitleSize(0.05);
  hBefore->GetXaxis()->SetTitleOffset(0.9);

  hAfter->SetFillStyle(3001);
  hAfter->SetFillColor(46);
  hAfter->SetLineColor(46);
  hAfter->Draw("SAMES hist");
  hAfter->SetStats(0);
  TF1 *f1 = new TF1(Form("f1%d",1),"gaus",1,5);
  f1->SetParameters(200, 4, 0.5);
  f1->SetLineColor(kBlue);
  TF1 *f2 = new TF1(Form("f1%d",2),"gaus",5,12);
  f2->SetParameters(200, 4, 0.5);
  f2->SetLineColor(kBlue);
  TF1 *f3 = new TF1(Form("f1%d",3),"gaus",12,19);
  f3->SetParameters(200, 4, 0.5);
  f3->SetLineColor(kBlue);

  hAfter->Fit(Form("f1%d",1),"","",1,5);
  f1->Draw("same");
  
  hAfter->Fit(Form("f1%d",2),"","",5,12);
  f2->Draw("same");

  hAfter->Fit(Form("f1%d",3),"","",12,19);
  f3->Draw("same");

  Double_t par1[3];
  Double_t par2[3];
  Double_t par3[3];
  
  // writes the fit results into the par array
  f1->GetParameters(par1);
  f2->GetParameters(par2);
  f3->GetParameters(par3);
  cout << "f1 mean: " << par1[1] << " sigma: " << par1[2] << endl;
  cout << "f2 mean: " << par2[1] << " sigma: " << par2[2] << endl;
  cout << "f3 mean: " << par3[1] << " sigma: " << par3[2] << endl;

  legend->AddEntry(hBefore,"Before Fitting","f");
  legend->AddEntry(hAfter,"After Fitting","f");
  legend->AddEntry(f1,"After Fitting Gaus Fit","l");
  legend->Draw("same");

  c1->Print("hEnergyFitting.pdf");

  TCanvas * c2 = new TCanvas("c2", "", 1200, 800);
  auto legend_Sum = new TLegend(0.5,0.7,0.68,0.88);
  hBefore_Sum->SetMaximum(3330);
  hBefore_Sum->SetStats(0);
  hBefore_Sum->SetFillStyle(4050);
  hBefore_Sum->SetFillColor(24);
  hBefore_Sum->SetLineColor(24);
  hBefore_Sum->Draw("hist");
  hAfter_Sum->SetFillStyle(3001);
  hAfter_Sum->SetFillColor(46);
  hAfter_Sum->SetLineColor(46);
  hAfter_Sum->Draw("SAMES hist");
  legend_Sum->AddEntry(hBefore_Sum,"Before Fitting","f");
  legend_Sum->AddEntry(hAfter_Sum,"After Fitting","f");
  legend_Sum->Draw("same");
  c2->Print("hSumFitting.pdf");

  return 0;
}
*/

