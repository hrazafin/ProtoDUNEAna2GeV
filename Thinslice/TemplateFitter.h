#ifndef TEMPLATEFITTER_H
#define TEMPLATEFITTER_H

#include "TMinuit.h"

class TH1D;

class TemplateFitter {
  
 public:
  
  TemplateFitter();
  
  void SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2);
  static void fcn(int &npar, double *gin, double &f, double *par, int iflag);
  void SetFitRange(int imin, int imax);

  void Fit();
  
  double GetPar();
  double GetParError();
  bool GetFitStatus();

  static TH1D *h0, *h1, *h2;
  static int i0, i1;

  TMinuit *gMinuit;

 private:

  bool fitsuccess;
  
};

TH1D* TemplateFitter::h0;
TH1D* TemplateFitter::h1;
TH1D* TemplateFitter::h2;

int TemplateFitter::i0;
int TemplateFitter::i1;

#endif
