#include "util.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>

using namespace std;

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype, bool fill_underflow, bool fill_overflow){
  //hist[0]->Fill(value);
  double weight = 1.;
  if (partype == 3) // kMuon
    weight = 1.;//59;
  if (partype>=0 && partype < pi::nIntTypes+1){
    if (value < hist[partype]->GetXaxis()->GetXmin()){ // underflow values
      if (fill_underflow) hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmin(), weight);
    }
    else if (value < hist[partype]->GetXaxis()->GetXmax()){
      hist[partype]->Fill(value, weight);
    }
    else{ // overflow values
      if (fill_overflow) hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmax()-0.000001, weight);
    }
  }
}

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype){
  double weight = 1.;
  if (partype == 3) // kMuon
    weight = 1.;//59;
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < pi::nIntTypes+1){
    hist[partype]->Fill(value1, value2, weight);
  }
}

void FillHist1D(TH1D *hist, const double &value, const double &wei){
  if (value<hist->GetXaxis()->GetXmin()){
    hist->Fill(hist->GetXaxis()->GetXmin(), wei);
  }
  else if (value<hist->GetXaxis()->GetXmax()){
    hist->Fill(value, wei);
  }
  else{
    hist->Fill(hist->GetXaxis()->GetXmax()-0.000001, wei);
  }
}

double GetPionKE(double length){
  if (length>0){
    // https://arxiv.org/abs/1306.1712
    return (8/(1-0.37))*pow(length,1-0.37);
  }
  else{
    cout<<"Invalid length "<<length<<", returning -1"<<endl;
    return -1;
  }
}

double GetTheta(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return v.Theta()*180/TMath::Pi();
  }
}

double GetPhi(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return v.Phi()*180/TMath::Pi();
  }
}

double GetThetaxz(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return atan2(v.X(), v.Z())*180/TMath::Pi();
  }
}

double GetThetayz(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return atan2(v.Y(), v.Z())*180/TMath::Pi();
  }
}
