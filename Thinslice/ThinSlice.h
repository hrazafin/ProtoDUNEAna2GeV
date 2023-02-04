#ifndef THINSLICE_H
#define THINSLICE_H

#include "TFile.h"
#include "SliceParams.h"
#include "HadAna.h"
#include "Unfold.h"

class anavar;
class HadAna;

class ThinSlice {

 public:

  ThinSlice();

  HadAna hadana;

  int reco_sliceID;
  int true_sliceID;

  bool isTestSample;

  bool selectCosmics;

  TH1D *reco_incE[pi::nthinslices];
  TH1D *true_incE[pi::nthinslices];
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

  double true_interactions[pi::nthinslices];
  double true_incidents[pi::nthinslices];

  TH1D *hreco_beam_type[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_reconstructable_beam_event[pi::nCuts][pi::nIntTypes+1];
  
  TH1D *htrue_beam_endZ[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_endZ[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_true_beam_endZ[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_vs_true_beam_endZ[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ[pi::nCuts][pi::nIntTypes+1];

  // after SCE correction
  TH1D *htrue_beam_endZ_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_endZ_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_true_beam_endZ_SCE[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_vs_true_beam_endZ_SCE[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ_SCE[pi::nCuts][pi::nIntTypes+1];

  TH1D *htrue_sliceID[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_sliceID[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_true_sliceID[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_vs_true_sliceID[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_true_vs_true_sliceID[pi::nCuts][pi::nIntTypes+1];

  TH1D *hmediandEdx[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdaughter_michel_score[pi::nCuts][pi::nIntTypes+1];
  TH1D *henergy_calorimetry_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdEdx_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdaughter_michel_scoreMu[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdaughter_michel_scorePi[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdaughter_michel_score2Mu[pi::nCuts][pi::nIntTypes+1];
  TH1D *htrackscore[pi::nCuts][pi::nIntTypes+1];
  TH1D *hemscore[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdEdx_5cm[pi::nCuts][pi::nIntTypes+1];

  TH1D *hmediandEdx_bkg[pi::nCuts][pi::nIntTypes+1];
  TH1D *hChi2_proton_bkg[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdaughter_michel_score_bkg[pi::nCuts][pi::nIntTypes+1];
  TH1D *hcostheta_bkg[pi::nCuts][pi::nIntTypes+1];
  TH1D *hmediandEdxSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hChi2_protonSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hdaughter_michel_scoreSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hcosthetaSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];

  TH1D *hdeltax[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdeltay[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdeltaz[pi::nCuts][pi::nIntTypes+1];
  TH1D *hcostheta[pi::nCuts][pi::nIntTypes+1];

  TH1D *hreco_beam_true_byE_matched[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_trklen[pi::nCuts][pi::nIntTypes+1];
  TH1D *htrue_trklen[pi::nCuts][pi::nIntTypes+1];
  TH1D *hdiff_trklen[pi::nCuts][pi::nIntTypes+1];
  TH2D *hreco_vs_true_trklen[pi::nCuts][pi::nIntTypes+1];
  TH1D *hbeam_score[pi::nCuts][pi::nIntTypes+1];
  TH2D *beam_score_vs_hreco_trklen[pi::nCuts][pi::nIntTypes+1];

  TH1D *hreco_beam_startX_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_startY_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_startZ_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_dcosX_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_dcosY_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_dcosZ_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_angleX_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_angleY_SCE[pi::nCuts][pi::nIntTypes+1];
  TH1D *hreco_beam_angleZ_SCE[pi::nCuts][pi::nIntTypes+1];

  TH2D *hreco_beam_startXY_SCE[pi::nCuts][pi::nIntTypes+1];
  
  TH1D *htrklen_csda_proton[pi::nCuts][pi::nIntTypes+1];
  TH1D *hChi2_proton[pi::nCuts][pi::nIntTypes+1];

  std::string fOutputFileName;
  TFile *outputFile;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(int cut, const anavar & evt);
  void FillSliceHist(const anavar & evt, int constraint_type, int cut=6);
  void SaveHistograms();

  void ProcessEvent(const anavar & evt, Unfold & uf);
  void CalcXS(const Unfold & uf);

  void Run(anavar & evt, Unfold & uf, Long64_t nentries);

  void SetSelectCosmics(bool sc);

};

#endif
