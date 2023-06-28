#ifndef TEMPLATEFITTER_H
#define TEMPLATEFITTER_H

#include "TMinuit.h"

using namespace std;

class TH1D;

class TemplateFitter {
  
 public:
  
  TemplateFitter();
  TemplateFitter(const int &nvars, const int &nsams);
  
  void SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2);
  void SetSimulHistograms(TH1D *hist0_s1, TH1D *hist1_s1, TH1D *hist2_s1, TH1D *hist0_s2, TH1D *hist1_s2, TH1D *hist2_s2, const int &break_s1, const int &break_s2);

  static void fcn(int &npar, double *gin, double &f, double *par, int iflag);
  static void Simulfcn(int &npar, double *gin, double &f, double *par, int iflag);
  
  void SetFitRange(int imin, int imax);

  void Fit();
  void SimulFit();

  double GetPar();
  vector<double> GetSimulPar();
  double GetParError();
  vector<double> GetSimulParError();

  bool GetFitStatus();

  static TH1D *h0, *h1, *h2;
  static int i0, i1;

  // Simul fit histograms and parameters
  static TH1D *h0_sample1, *h1_sample1, *h2_sample1;
  static TH1D *h0_sample2, *h1_sample2, *h2_sample2;
  static int ibreak_sample1, ibreak_sample2;
  static int iend_sample1, iend_sample2;
  static int nparameters, nsamples; 

  TMinuit *gMinuit;

 private:

  bool fitsuccess;
  
};

TH1D* TemplateFitter::h0;
TH1D* TemplateFitter::h1;
TH1D* TemplateFitter::h2;

int TemplateFitter::i0;
int TemplateFitter::i1;

// Simul fit histograms and parameters
TH1D* TemplateFitter::h0_sample1;
TH1D* TemplateFitter::h1_sample1;
TH1D* TemplateFitter::h2_sample1;
TH1D* TemplateFitter::h0_sample2;
TH1D* TemplateFitter::h1_sample2;
TH1D* TemplateFitter::h2_sample2;

int TemplateFitter::ibreak_sample1;
int TemplateFitter::ibreak_sample2;
int TemplateFitter::iend_sample1;
int TemplateFitter::iend_sample2;

int TemplateFitter::nparameters;
int TemplateFitter::nsamples;

#endif
