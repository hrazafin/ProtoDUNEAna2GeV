#ifndef PROTONINEL_H
#define PROTONINEL_H

#include "TFile.h"
#include "SliceParams.h"
#include "HadAna.h"
#include "Unfold.h"

class anavar;
class HadAna;

class ProtonInel {

 public:

  ProtonInel();

  HadAna hadana;

  int reco_sliceID;
  int true_sliceID;

  bool isTestSample;

  TH1D *reco_incE[p::nthinslices];
  TH1D *true_incE[p::nthinslices];
  TH1D *reco_AngCorr;
  TH1D *true_AngCorr;

  TH1D *h_truesliceid_pion_all;
  TH1D *h_truesliceid_pion_uf;
  TH1D *h_truesliceid_pion_cuts;
  TH1D *h_truesliceid_pioninelastic_all;
  TH1D *h_truesliceid_pioninelastic_uf;
  TH1D *h_truesliceid_pioninelastic_cuts;
  TH1D *h_recosliceid_allevts_cuts;
  TH1D *h_recosliceid_pion_cuts;
  TH1D *h_recosliceid_pioninelastic_cuts;

  double true_interactions[p::nthinslices];
  double true_incidents[p::nthinslices];

  TH1D *htrue_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_true_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ[p::nCuts][p::nIntTypes+1];

  TH1D *htrue_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_true_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];

  TH1D *htrue_sliceID[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_sliceID[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_true_sliceID[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_sliceID[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_true_vs_true_sliceID[p::nCuts][p::nIntTypes+1];

  TH1D *hmediandEdx[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_score[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_scoreMu[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_scorePi[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_score2Mu[p::nCuts][p::nIntTypes+1];
  TH1D *htrackscore[p::nCuts][p::nIntTypes+1];
  TH1D *hemscore[p::nCuts][p::nIntTypes+1];
  TH1D *hdEdx_5cm[p::nCuts][p::nIntTypes+1];

  TH1D *hmediandEdxSlice[p::nthinslices][p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_scoreSlice[p::nthinslices][p::nCuts][p::nIntTypes+1];

  TH1D *hdeltax[p::nCuts][p::nIntTypes+1];
  TH1D *hdeltay[p::nCuts][p::nIntTypes+1];
  TH1D *hdeltaz[p::nCuts][p::nIntTypes+1];
  TH1D *hcostheta[p::nCuts][p::nIntTypes+1];

  TH1D *htrklen[p::nCuts][p::nIntTypes+1];

  TH1D *hreco_beam_startX_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_startY_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_startZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_dcosX_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_dcosY_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_dcosZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_angleX_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_angleY_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_angleZ_SCE[p::nCuts][p::nIntTypes+1];

  TH2D *hreco_beam_startXY_SCE[p::nCuts][p::nIntTypes+1];

  std::string fOutputFileName;
  TFile *outputFile;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(int cut, const anavar & evt);
  void SaveHistograms();

  void ProcessEvent(const anavar & evt, Unfold & uf);
  void CalcXS(const Unfold & uf);

  void Run(anavar & evt, Unfold & uf);

};

#endif
