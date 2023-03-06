#include "TemplateFitter.h"
#include "TH1D.h"
#include <iostream>
#include <vector>

TemplateFitter::TemplateFitter(){
  gMinuit = new TMinuit(1);
  gMinuit->SetFCN(fcn);
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
}

TemplateFitter::TemplateFitter(const int &nvars, const int &nsams){
  nparameters = nvars;
  nsamples = nsams;
  gMinuit = new TMinuit(nparameters);
  gMinuit->SetFCN(Simulfcn);
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
}

void TemplateFitter::SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2){
  h0 = hist0;
  h1 = hist1;
  h2 = hist2;

  i0 = 1;
  i1 = h0->GetNbinsX();

}

void TemplateFitter::SetSimulHistograms(TH1D *hist0_s1, TH1D *hist1_s1, TH1D *hist2_s1, TH1D *hist0_s2, TH1D *hist1_s2, TH1D *hist2_s2, const int &break_s1, const int &break_s2){
  h0_sample1 = hist0_s1;
  h1_sample1 = hist1_s1;
  h2_sample1 = hist2_s1;

  h0_sample2 = hist0_s2;
  h1_sample2 = hist1_s2;
  h2_sample2 = hist2_s2;

  iend_sample1 = h0_sample1->GetNbinsX();
  iend_sample2 = h0_sample2->GetNbinsX();

  ibreak_sample1 = break_s1;
  ibreak_sample2 = break_s2;

}


void TemplateFitter::SetFitRange(int imin, int imax){
  if (imin>=0) i0 = imin;
  if (imax>=0) i1 = imax;
}

void TemplateFitter::fcn(int &npar, double *gin, double &f, double *par, int iflag){
  
  double chisq = 0;

  for (int i = 0; i<=i0; ++i){
    double x = h0->GetBinContent(i);
    double y = h1->GetBinContent(i);
    double z = h2->GetBinContent(i);
    
    double ex = h0->GetBinError(i);
    double ey = h1->GetBinError(i);
    double ez = h2->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[0],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[0],2));
  }

  for (int i = i1; i<=19; ++i){
    double x = h0->GetBinContent(i);
    double y = h1->GetBinContent(i);
    double z = h2->GetBinContent(i);
    
    double ex = h0->GetBinError(i);
    double ey = h1->GetBinError(i);
    double ez = h2->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[0],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[0],2));
  }

  f = chisq;

}

void TemplateFitter::Simulfcn(int &npar, double *gin, double &f, double *par, int iflag){
  
  double chisq = 0;
  
  //Low sample1
  for (int i = 1; i<=ibreak_sample1; ++i){
    double x = h0_sample1->GetBinContent(i);
    double y = h1_sample1->GetBinContent(i);
    double z = h2_sample1->GetBinContent(i);
    
    double ex = h0_sample1->GetBinError(i);
    double ey = h1_sample1->GetBinError(i);
    double ez = h2_sample1->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[0],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[0],2));

  }
  
  //High sample1
  for (int i = ibreak_sample1+1; i<=iend_sample1; ++i){
    double x = h0_sample1->GetBinContent(i);
    double y = h1_sample1->GetBinContent(i);
    double z = h2_sample1->GetBinContent(i);
    
    double ex = h0_sample1->GetBinError(i);
    double ey = h1_sample1->GetBinError(i);
    double ez = h2_sample1->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[1],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[1],2));
  }

  //Low sample2
  for (int i = 1; i<=ibreak_sample2; ++i){
    double x = h0_sample2->GetBinContent(i);
    double y = h1_sample2->GetBinContent(i);
    double z = h2_sample2->GetBinContent(i);
    
    double ex = h0_sample2->GetBinError(i);
    double ey = h1_sample2->GetBinError(i);
    double ez = h2_sample2->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[2],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[2],2));
  }
  //High sample2
  for (int i = ibreak_sample2+1; i<=iend_sample2; ++i){
    double x = h0_sample2->GetBinContent(i);
    double y = h1_sample2->GetBinContent(i);
    double z = h2_sample2->GetBinContent(i);
    
    double ex = h0_sample2->GetBinError(i);
    double ey = h1_sample2->GetBinError(i);
    double ez = h2_sample2->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[3],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[3],2));
  }
  
  f = chisq; 

}

void TemplateFitter::Fit(){

  double arglist[10];
  int ierflg = 0;

  gMinuit->mncler();
  double vstart = 1;
  double step = 0.01;
  gMinuit->mnparm(0,"corr_fact",vstart,step,0,10,ierflg);

  fitsuccess = false;

  if (h0->Integral(i0,i1) && h2->Integral(i0,i1)){
    arglist[0] = 500;
    arglist[1] = 1;
    gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
    
    double par, epar;
    gMinuit->GetParameter(0,par,epar);

    TString test =  gMinuit->fCstatu.Data(); 
    if (test.EqualTo("CONVERGED ")){
      std::cout<<"Best fit = "<<par<<" error = "<<epar<<std::endl;
      fitsuccess = true;
    }
  }
  else{
    std::cout<<"No fit was done because data and/or template are empty."<<std::endl;
  }
}

void TemplateFitter::SimulFit(){

  double arglist[10];
  int ierflg = 0;

  gMinuit->mncler();
  double vstart = 1;
  double step = 0.01;
  gMinuit->mnparm(0,"corr_fact_1",vstart,step,0,10,ierflg);
  gMinuit->mnparm(1,"corr_fact_2",vstart,step,0,10,ierflg);
  gMinuit->mnparm(2,"corr_fact_3",vstart,step,0,10,ierflg);
  gMinuit->mnparm(3,"corr_fact_4",vstart,step,0,10,ierflg);

  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);


  fitsuccess = false;

  if (h0_sample1->Integral(0,iend_sample1) && h2_sample2->Integral(0,iend_sample2)){
    arglist[0] = 500;
    arglist[1] = 1;
    gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
    
    double par[nparameters], epar[nparameters];
    gMinuit->GetParameter(0,par[0],epar[0]);
    gMinuit->GetParameter(1,par[1],epar[1]);
    gMinuit->GetParameter(2,par[2],epar[2]);
    gMinuit->GetParameter(3,par[3],epar[3]);

    TString test =  gMinuit->fCstatu.Data(); 
    if (test.EqualTo("CONVERGED ")){
      std::cout<<"Best fit 1 = "<<par[0]<<" error = "<<epar[0]<<std::endl;
      std::cout<<"Best fit 2 = "<<par[1]<<" error = "<<epar[1]<<std::endl;
      std::cout<<"Best fit 3 = "<<par[2]<<" error = "<<epar[2]<<std::endl;
      std::cout<<"Best fit 4 = "<<par[3]<<" error = "<<epar[3]<<std::endl;
      fitsuccess = true;
    }
  }
  else{
    std::cout<<"No fit was done because data and/or template are empty."<<std::endl;
  }
}

double TemplateFitter::GetPar(){
  if (fitsuccess){
    double par, epar;
    gMinuit->GetParameter(0,par,epar);
    return par;
  }
  else{
    return 1;
  }
}

vector<double> TemplateFitter::GetSimulPar(){
  vector<double> parVect; 
  if (fitsuccess){
    double par[4], epar[4];
    gMinuit->GetParameter(0,par[0],epar[0]); parVect.push_back(par[0]);
    gMinuit->GetParameter(1,par[1],epar[1]); parVect.push_back(par[1]);
    gMinuit->GetParameter(2,par[2],epar[2]); parVect.push_back(par[2]);
    gMinuit->GetParameter(3,par[3],epar[3]); parVect.push_back(par[3]);
    return parVect;
  }
  else{
    parVect.push_back(1.0);
    parVect.push_back(1.0);
    parVect.push_back(1.0);
    parVect.push_back(1.0); 
    return parVect;
  }
}

double TemplateFitter::GetParError(){
  if (fitsuccess){
    double par, epar;
    gMinuit->GetParameter(0,par,epar);
    return epar;
  }
  else{
    return 2;
  }
}

vector<double> TemplateFitter::GetSimulParError(){
  vector<double> parVectError; 
  if (fitsuccess){
    double par[4], epar[4];
    gMinuit->GetParameter(0,par[0],epar[0]); parVectError.push_back(epar[0]);
    gMinuit->GetParameter(1,par[1],epar[1]); parVectError.push_back(epar[1]);
    gMinuit->GetParameter(2,par[2],epar[2]); parVectError.push_back(epar[2]);
    gMinuit->GetParameter(3,par[3],epar[3]); parVectError.push_back(epar[3]);
    return parVectError;
  }
  else{
    parVectError.push_back(2.0);
    parVectError.push_back(2.0);
    parVectError.push_back(2.0);
    parVectError.push_back(2.0); 
    return parVectError;
  }
}

bool TemplateFitter::GetFitStatus(){
  return fitsuccess;
}
