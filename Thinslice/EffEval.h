#ifndef EFFEVAL_H
#define EFFEVAL_H

#include "TFile.h"
#include "SliceParams.h"
#include "HadAna.h"

class anavar;
class HadAna;
class TEfficiency;
class TH1D;
class TH2D;

class EffEval {

 public:

  EffEval();

  HadAna hadana;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(const anavar & evt);
  void SaveHistograms();

  void ProcessEvent(const anavar & evt);

  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;
  TEfficiency *h_Eff_Ppi;
  TEfficiency *h_Eff_Ppi_michel;
  TEfficiency *h_Eff_thetapi;
  TEfficiency *h_Eff_thetapi_michel;
  TEfficiency *h_Eff_phipi;
  TEfficiency *h_Eff_phipi_michel;
  TH1D *h_true_Ppi_all;
  TH1D *h_true_Ppi_sel;
  TH2D *h_reco_true_Ppi_sel;
  TH1D *h_res_Ppi_sel;
  TH2D *h_res_Ppi_michelscore;
  TH1D *h_true_Ppi_michel;
  TH2D *h_reco_true_Ppi_michel;
  TH1D *h_res_Ppi_michel;

  TH1D *h_true_thetapi_all;
  TH1D *h_true_thetapi_sel;
  TH2D *h_reco_true_thetapi_sel;
  TH1D *h_res_thetapi_sel;
  TH1D *h_true_thetapi_michel;
  TH2D *h_reco_true_thetapi_michel;
  TH1D *h_res_thetapi_michel;

  TH1D *h_true_phipi_all;
  TH1D *h_true_phipi_sel;
  TH2D *h_reco_true_phipi_sel;
  TH1D *h_res_phipi_sel;
  TH1D *h_true_phipi_michel;
  TH2D *h_reco_true_phipi_michel;
  TH1D *h_res_phipi_michel;

  ///////////////////////////////

  TEfficiency *h_Eff_Pp;
  TEfficiency *h_Eff_thetap;
  TEfficiency *h_Eff_phip;
  TH1D *h_true_Pp_all;
  TH1D *h_true_Pp_sel;
  TH2D *h_reco_true_Pp_sel;
  TH1D *h_res_Pp_sel;

  TH1D *h_true_thetap_all;
  TH1D *h_true_thetap_sel;
  TH2D *h_reco_true_thetap_sel;
  TH1D *h_res_thetap_sel;

  TH1D *h_true_phip_all;
  TH1D *h_true_phip_sel;
  TH2D *h_reco_true_phip_sel;
  TH1D *h_res_phip_sel;

};

#endif
