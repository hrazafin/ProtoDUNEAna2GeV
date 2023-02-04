#include "ThinSlice.h"
#include "HadAna.h"
#include "anavar.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "util.h"

#include <iostream>

ThinSlice::ThinSlice(){
  hadana.InitPi();
  selectCosmics = false;
}

void ThinSlice::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");
  
  for (int i = 0; i<pi::nthinslices; ++i){ // energy distribution in each thin slice
    reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %.1f < z < %.1f (cm)",i*pi::thinslicewidth, (i+1)*pi::thinslicewidth), pi::nbinse, 0, 1200.);
    true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %.1f < z < %.1f (cm)",i*pi::thinslicewidth, (i+1)*pi::thinslicewidth), pi::nbinse, 0, 1200.);
    reco_incE[i]->Sumw2();
    true_incE[i]->Sumw2();
  }

  reco_AngCorr = new TH1D("reco_AngCorr","Reco angle correction", 100, 0, 1.);
  true_AngCorr = new TH1D("true_AngCorr","true angle correction", 100, 0, 1.);
  reco_AngCorr->Sumw2();
  true_AngCorr->Sumw2();

  h_truesliceid_pion_all = new TH1D("h_truesliceid_pion_all","h_truesliceid_pion_all;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_truesliceid_pion_cuts = new TH1D("h_truesliceid_pion_cuts","h_truesliceid_pion_cuts;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_truesliceid_pioninelastic_all = new TH1D("h_truesliceid_pioninelastic_all","h_truesliceid_pioninelastic_all;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_truesliceid_pioninelastic_cuts = new TH1D("h_truesliceid_pioninelastic_cuts","h_truesliceid_pioninelastic_cuts;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recosliceid_pion_cuts = new TH1D("h_recosliceid_pion_cuts","h_recosliceid_pion_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recosliceid_pioninelastic_cuts = new TH1D("h_recosliceid_pioninelastic_cuts","h_recosliceid_pioninelastic_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);

  h_truesliceid_pion_all->Sumw2();
  h_truesliceid_pion_cuts->Sumw2();
  h_truesliceid_pioninelastic_all->Sumw2();
  h_truesliceid_pioninelastic_cuts->Sumw2();
  h_recosliceid_allevts_cuts->Sumw2();
  h_recosliceid_pion_cuts->Sumw2();
  h_recosliceid_pioninelastic_cuts->Sumw2();

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      hreco_beam_type[i][j] = new TH1D(Form("hreco_beam_type_%d_%d",i,j),Form("hreco_beam_type, %s, %s;hreco_beam_type", pi::cutName[i], pi::intTypeName[j]), 21, -1, 20);
      hreco_beam_type[i][j]->Sumw2();
      hreco_reconstructable_beam_event[i][j] = new TH1D(Form("hreco_reconstructable_beam_event_%d_%d",i,j),Form("hreco_reconstructable_beam_event, %s, %s;hreco_reconstructable_beam_event", pi::cutName[i], pi::intTypeName[j]), 21, -1, 20);
      hreco_reconstructable_beam_event[i][j]->Sumw2();
      
      htrue_beam_endZ[i][j] = new TH1D(Form("htrue_beam_endZ_%d_%d",i,j),Form("true_beam_endZ, %s, %s;true_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      htrue_beam_endZ[i][j]->Sumw2();
      hreco_beam_endZ[i][j] = new TH1D(Form("hreco_beam_endZ_%d_%d",i,j),Form("reco_beam_endZ, %s, %s;reco_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      hreco_beam_endZ[i][j]->Sumw2();
      hreco_true_beam_endZ[i][j] = new TH1D(Form("hreco_true_beam_endZ_%d_%d",i,j), Form("reco_true_beam_endZ, %s, %s;reco_beam_endZ - true_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ[i][j]->Sumw2();
      hreco_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 70, -100, 600);
      hreco_true_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco - true_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 100, -100, 100);

      htrue_beam_endZ_SCE[i][j] = new TH1D(Form("htrue_beam_endZ_SCE_%d_%d",i,j),Form("true_beam_endZ_SCE, %s, %s;true_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      htrue_beam_endZ_SCE[i][j]->Sumw2();
      hreco_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_beam_endZ_SCE_%d_%d",i,j),Form("reco_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      hreco_beam_endZ_SCE[i][j]->Sumw2();
      hreco_true_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_true_beam_endZ_SCE_%d_%d",i,j), Form("reco_true_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE - true_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ_SCE[i][j]->Sumw2();
      hreco_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 70, -100, 600);
      hreco_true_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco - true_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 100, -100, 100);

      htrue_sliceID[i][j] = new TH1D(Form("htrue_sliceID_%d_%d",i,j),Form("true_sliceID, %s, %s;true_sliceID (cm)", pi::cutName[i], pi::intTypeName[j]), 26, -1, 25);
      htrue_sliceID[i][j]->Sumw2();
      hreco_sliceID[i][j] = new TH1D(Form("hreco_sliceID_%d_%d",i,j),Form("reco_sliceID, %s, %s;reco_sliceID", pi::cutName[i], pi::intTypeName[j]), 26, -1, 25);
      hreco_sliceID[i][j]->Sumw2();
      hreco_true_sliceID[i][j] = new TH1D(Form("hreco_true_sliceID_%d_%d",i,j), Form("reco_true_sliceID, %s, %s;reco_sliceID - true_sliceID", pi::cutName[i], pi::intTypeName[j]), 20, -10, 10);
      hreco_true_sliceID[i][j]->Sumw2();
      hreco_vs_true_sliceID[i][j]= new TH2D(Form("hreco_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID", pi::cutName[i], pi::intTypeName[j]), 26, -1, 25, 26, -1, 25);
      hreco_true_vs_true_sliceID[i][j]= new TH2D(Form("hreco_true_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID - true_sliceID", pi::cutName[i], pi::intTypeName[j]), 26, -1, 25, 20, -10, 10);

      hmediandEdx[i][j] = new TH1D(Form("hmediandEdx_%d_%d",i,j), Form("mediandEdx, %s, %s;Median dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 5);
      hmediandEdx[i][j]->Sumw2();

      hdaughter_michel_score[i][j] = new TH1D(Form("hdaughter_michel_score_%d_%d",i,j), Form("daughter_michel_score, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 110, -0.1, 1);
      hdaughter_michel_score[i][j]->Sumw2();
      
      henergy_calorimetry_SCE[i][j] = new TH1D(Form("henergy_calorimetry_SCE_%d_%d",i,j), Form("Energy_calorimetry_SCE_corrected, %s, %s;Energy (MeV)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 500);
      henergy_calorimetry_SCE[i][j]->Sumw2();
      hdEdx_SCE[i][j] = new TH1D(Form("hdEdx_SCE_%d_%d",i,j), Form("dEdx_SCE_corrected, %s, %s; dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 1.7, 4.2);
      hdEdx_SCE[i][j]->Sumw2();

      hdaughter_michel_scoreMu[i][j] = new TH1D(Form("hdaughter_michel_scoreMu_%d_%d",i,j), Form("daughter_michel_scoreMu, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 10, 0, 1);
      hdaughter_michel_scoreMu[i][j]->Sumw2();

      hdaughter_michel_score2Mu[i][j] = new TH1D(Form("hdaughter_michel_score2Mu_%d_%d",i,j), Form("daughter_michel_score2Mu, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 10, 0, 1);
      hdaughter_michel_score2Mu[i][j]->Sumw2();

      hdaughter_michel_scorePi[i][j] = new TH1D(Form("hdaughter_michel_scorePi_%d_%d",i,j), Form("daughter_michel_scorePi, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 10, 0, 1);
      hdaughter_michel_scorePi[i][j]->Sumw2();

      hmediandEdx_bkg[i][j] = new TH1D(Form("hmediandEdx_bkg_%d_%d",i,j), Form("mediandEdx_bkg, %s, %s;Median dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 8);
      hmediandEdx_bkg[i][j]->Sumw2();
      hChi2_proton_bkg[i][j] = new TH1D(Form("hChi2_proton_bkg_%d_%d",i,j), Form("Chi2_proton_bkg, %s, %s;Chi2/Ndof", pi::cutName[i], pi::intTypeName[j]), 100, 0, 100);
      hChi2_proton_bkg[i][j]->Sumw2();
      hdaughter_michel_score_bkg[i][j] = new TH1D(Form("hdaughter_michel_score_bkg_%d_%d",i,j), Form("daughter_michel_score_bkg, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 100, 0, 1);
      hdaughter_michel_score_bkg[i][j]->Sumw2();
      hcostheta_bkg[i][j] = new TH1D(Form("hcostheta_bkg_%d_%d",i,j), Form("costheta_bkg, %s, %s;cos#theta", pi::cutName[i], pi::intTypeName[j]), 75, 0.85, 1);
      hcostheta_bkg[i][j]->Sumw2();
      for (int k = 0; k<pi::nthinslices; ++k){
        hmediandEdxSlice[k][i][j] = new TH1D(Form("hmediandEdxSlice_%d_%d_%d",k,i,j), Form("mediandEdx, %s, %s, sliceID = %d;Median dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j], k), 14, 1, 8);
        hmediandEdxSlice[k][i][j]->Sumw2();
        
        hChi2_protonSlice[k][i][j] = new TH1D(Form("hChi2_protonSlice_%d_%d_%d",k,i,j), Form("Chi2_proton, %s, %s, sliceID = %d;Chi2_proton/Ndf", pi::cutName[i], pi::intTypeName[j], k), 10, 0, 100);
        hChi2_protonSlice[k][i][j]->Sumw2();

        hdaughter_michel_scoreSlice[k][i][j] = new TH1D(Form("hdaughter_michel_scoreSlice_%d_%d_%d",k,i,j), Form("daughter_michel_score, %s, %s, sliceID = %d;Michel score", pi::cutName[i], pi::intTypeName[j], k), 10, 0, 1);
        hdaughter_michel_scoreSlice[k][i][j]->Sumw2();
        
        hcosthetaSlice[k][i][j] = new TH1D(Form("hcosthetaSlice_%d_%d_%d",k,i,j), Form("costheta, %s, %s, sliceID = %d;Cos(theta)", pi::cutName[i], pi::intTypeName[j], k), 15, 0.85, 1);
        hcosthetaSlice[k][i][j]->Sumw2();
      }        

      htrackscore[i][j] = new TH1D(Form("htrackscore_%d_%d",i,j), Form("trackscore, %s, %s;Track score", pi::cutName[i], pi::intTypeName[j]), 110, -0.1, 1);
      htrackscore[i][j]->Sumw2();

      hemscore[i][j] = new TH1D(Form("hemscore_%d_%d",i,j), Form("emscore, %s, %s;Em score", pi::cutName[i], pi::intTypeName[j]), 50, 0, 1);
      hemscore[i][j]->Sumw2();

      hdEdx_5cm[i][j] = new TH1D(Form("hdEdx_5cm_%d_%d",i,j), Form("dEdx_5cm, %s, %s;dE/dx_5cm (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 5);
      hdEdx_5cm[i][j]->Sumw2();

      hdeltax[i][j] = new TH1D(Form("hdeltax_%d_%d",i,j), Form("deltax, %s, %s;#Deltax/#sigma_{x}", pi::cutName[i], pi::intTypeName[j]), 100, -10, 10);
      hdeltax[i][j]->Sumw2();

      hdeltay[i][j] = new TH1D(Form("hdeltay_%d_%d",i,j), Form("deltay, %s, %s;#Deltay/#sigma_{y}", pi::cutName[i], pi::intTypeName[j]), 100, -10, 10);
      hdeltay[i][j]->Sumw2();

      hdeltaz[i][j] = new TH1D(Form("hdeltaz_%d_%d",i,j), Form("deltaz, %s, %s;#Deltaz/#sigma_{z}", pi::cutName[i], pi::intTypeName[j]), 100, -10, 10);
      hdeltaz[i][j]->Sumw2();

      hcostheta[i][j] = new TH1D(Form("hcostheta_%d_%d",i,j), Form("costheta, %s, %s;cos#theta", pi::cutName[i], pi::intTypeName[j]), 100, 0.9, 1);
      hcostheta[i][j]->Sumw2();

      hreco_beam_true_byE_matched[i][j] = new TH1D(Form("hreco_beam_true_byE_matched_%d_%d",i,j), Form("reco_beam_true_byE_matched, %s, %s;Truth matched", pi::cutName[i], pi::intTypeName[j]), 2, 0, 2);
      hreco_beam_true_byE_matched[i][j]->Sumw2();
      //const double xbins[25] = {-10.,0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,600}; // a user-defined binning
      //hreco_trklen[i][j] = new TH1D(Form("hreco_trklen_%d_%d",i,j), Form("reco_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 24, xbins);
      hreco_trklen[i][j] = new TH1D(Form("hreco_trklen_%d_%d",i,j), Form("reco_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 61, -10, 600);
      hreco_trklen[i][j]->Sumw2();
      htrue_trklen[i][j] = new TH1D(Form("htrue_trklen_%d_%d",i,j), Form("true_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 61, -10, 600);
      htrue_trklen[i][j]->Sumw2();
      hdiff_trklen[i][j] = new TH1D(Form("hdiff_trklen_%d_%d",i,j), Form("diff_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 60, -600, 600);
      hdiff_trklen[i][j]->Sumw2();
      hreco_vs_true_trklen[i][j]= new TH2D(Form("hreco_vs_true_trklen_%d_%d",i,j), Form("%s, %s;true_trklen (cm);reco_trklen (cm)", pi::cutName[i], pi::intTypeName[j]), 61, -10, 600, 61, -10, 600);
      hbeam_score[i][j] = new TH1D(Form("hbeam_score_%d_%d",i,j), Form("Beam_score, %s, %s;Beam score", pi::cutName[i], pi::intTypeName[j]), 140, -0.23, 0.12);
      hbeam_score[i][j]->Sumw2();
      beam_score_vs_hreco_trklen[i][j]= new TH2D(Form("beam_score_vs_hreco_trklen_%d_%d",i,j), Form("%s, %s;reco trklen (cm);beam_score", pi::cutName[i], pi::intTypeName[j]), 50, -100, 400, 140, -0.23, 0.12);

      hreco_beam_startX_SCE[i][j] = new TH1D(Form("hreco_beam_startX_SCE_%d_%d",i,j), Form("reco_beam_startX_SCE, %s, %s; reco_beam_startX_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -80, 20);
      hreco_beam_startX_SCE[i][j]->Sumw2();

      hreco_beam_startY_SCE[i][j] = new TH1D(Form("hreco_beam_startY_SCE_%d_%d",i,j), Form("reco_beam_startY_SCE, %s, %s; reco_beam_startY_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, 350, 500);
      hreco_beam_startY_SCE[i][j]->Sumw2();

      hreco_beam_startZ_SCE[i][j] = new TH1D(Form("hreco_beam_startZ_SCE_%d_%d",i,j), Form("reco_beam_startZ_SCE, %s, %s; reco_beam_startZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -5, 10);
      hreco_beam_startZ_SCE[i][j]->Sumw2();

      hreco_beam_dcosX_SCE[i][j] = new TH1D(Form("hreco_beam_dcosX_SCE_%d_%d",i,j), Form("hreco_beam_dcosX_SCE, %s, %s; reco_beam_dcosX_SCE", pi::cutName[i], pi::intTypeName[j]), 100, -1, 1);
      hreco_beam_dcosX_SCE[i][j]->Sumw2();

      hreco_beam_dcosY_SCE[i][j] = new TH1D(Form("hreco_beam_dcosY_SCE_%d_%d",i,j), Form("hreco_beam_dcosY_SCE, %s, %s; reco_beam_dcosY_SCE", pi::cutName[i], pi::intTypeName[j]), 100, -1, 1);
      hreco_beam_dcosY_SCE[i][j]->Sumw2();

      hreco_beam_dcosZ_SCE[i][j] = new TH1D(Form("hreco_beam_dcosZ_SCE_%d_%d",i,j), Form("hreco_beam_dcosZ_SCE, %s, %s; reco_beam_dcosZ_SCE", pi::cutName[i], pi::intTypeName[j]), 100, -1, 1);
      hreco_beam_dcosZ_SCE[i][j]->Sumw2();

      hreco_beam_angleX_SCE[i][j] = new TH1D(Form("hreco_beam_angleX_SCE_%d_%d",i,j), Form("hreco_beam_angleX_SCE, %s, %s; #theta_{x} (deg)", pi::cutName[i], pi::intTypeName[j]), 180, 0, 180);
      hreco_beam_angleX_SCE[i][j]->Sumw2();

      hreco_beam_angleY_SCE[i][j] = new TH1D(Form("hreco_beam_angleY_SCE_%d_%d",i,j), Form("hreco_beam_angleY_SCE, %s, %s; #theta_{y} (deg)", pi::cutName[i], pi::intTypeName[j]), 180, 0, 180);
      hreco_beam_angleY_SCE[i][j]->Sumw2();

      hreco_beam_angleZ_SCE[i][j] = new TH1D(Form("hreco_beam_angleZ_SCE_%d_%d",i,j), Form("hreco_beam_angleZ_SCE, %s, %s; #theta_{z} (deg)", pi::cutName[i], pi::intTypeName[j]), 90, 0, 90);
      hreco_beam_angleZ_SCE[i][j]->Sumw2();

      hreco_beam_startXY_SCE[i][j] = new TH2D(Form("hreco_beam_startXY_SCE_%d_%d",i,j), Form("reco_beam_startXY_SCE, %s, %s;reco_beam_startX_SCE (cm);reco_beam_startY_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 1000, -360, 360, 1000, 0, 700);
      
      htrklen_csda_proton[i][j] = new TH1D(Form("htrklen_csda_proton_%d_%d",i,j), Form("trklen_csda_proton, %s, %s;Track length / CSDA", pi::cutName[i], pi::intTypeName[j]), 61, -0.1, 6);
      htrklen_csda_proton[i][j]->Sumw2();
      hChi2_proton[i][j] = new TH1D(Form("hChi2_proton_%d_%d",i,j), Form("Chi2_proton, %s, %s;Chi2/Ndof", pi::cutName[i], pi::intTypeName[j]), 101, -1, 100);
      hChi2_proton[i][j]->Sumw2();
    }
  }

   for (int i = 0; i<pi::nthinslices; ++i){
     true_interactions[i] = 0;
     true_incidents[i] = 0;
   }
   
   //response_SliceID_Pion = new RooUnfoldResponse(pi::nthinslices+2, -1, pi::nthinslices+1, "response_SliceID_Pion");
   //response_SliceID_PionInEl = new RooUnfoldResponse(pi::nthinslices+2, -1, pi::nthinslices+1, "response_SliceID_PionInEl");

}

void ThinSlice::ProcessEvent(const anavar & evt, Unfold & uf){
  hadana.ProcessEvent(evt);
  reco_sliceID = -1;
  true_sliceID = -1;

  isTestSample = (hadana.pitype == pi::kData); // fake data
  //if (evt.MC && evt.event%2 == 0) isTestSample = false;

  if (hadana.fAllTrackCheck) {//using all-track reconstruction
    if (evt.MC){
      //true_sliceID = int(evt.true_beam_endZ/pi::thinslicewidth);
      true_sliceID = int(hadana.true_trklen/pi::thinslicewidth);
      if (true_sliceID < 0) true_sliceID = -1;
      if (evt.true_beam_endZ < 0) true_sliceID = -1;
      if (true_sliceID >= pi::nthinslices) true_sliceID = pi::nthinslices;
      if (evt.true_beam_PDG == 211){
        for (int i = 0; i<=true_sliceID; ++i){
          if (i<pi::nthinslices) ++true_incidents[i]; // count incident events
        }
      }
      if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
        if (true_sliceID < pi::nthinslices && true_sliceID>=0){
          ++true_interactions[true_sliceID]; // count interaction events
        }
        // Reco info
        if (!(evt.reco_beam_calo_wire_allTrack->empty()) && evt.reco_beam_true_byE_matched){ // truth matched. so it must be the true track?
          std::vector<std::vector<double>> vincE(pi::nthinslices);
          for (size_t i = 0; i<evt.reco_beam_calo_wire_allTrack->size(); ++i){
            //int this_sliceID = int((*evt.reco_beam_calo_Z)[i]/pi::thinslicewidth);
            int this_sliceID = int((hadana.reco_trklen_accum)[i]/pi::thinslicewidth);
            if (this_sliceID>=pi::nthinslices) continue;
            if (this_sliceID<0) continue;
            double this_incE = (*evt.reco_beam_incidentEnergies_allTrack)[i];
            vincE[this_sliceID].push_back(this_incE);
          }
          for (size_t i = 0; i<vincE.size(); ++i){
            if (!vincE[i].empty()){
              double sum_incE = 0;
              for (size_t j = 0; j<vincE[i].size(); ++j){
                sum_incE += vincE[i][j];
              }
              reco_incE[i]->Fill(sum_incE/vincE[i].size());
            }
          }
          TVector3 pt0(evt.reco_beam_calo_startX_allTrack,
                       evt.reco_beam_calo_startY_allTrack,
                       evt.reco_beam_calo_startZ_allTrack);
          TVector3 pt1(evt.reco_beam_calo_endX_allTrack,
                       evt.reco_beam_calo_endY_allTrack,
                       evt.reco_beam_calo_endZ_allTrack);
          TVector3 dir = pt1 - pt0; // direction of the track is determine by the start/end point
          dir = dir.Unit();
          reco_AngCorr->Fill(dir.Z()); // projection to Z of the direction of the track
        }

        // True info
        if (!(evt.true_beam_traj_Z->empty())){
          std::vector<std::vector<double>> vincE(pi::nthinslices);
          for (size_t i = 0; i<evt.true_beam_traj_Z->size()-1; ++i){//last point always has KE = 0
            //int this_sliceID = int((*evt.true_beam_traj_Z)[i]/pi::thinslicewidth);
            int this_sliceID = int((hadana.true_trklen_accum)[i]/pi::thinslicewidth);
            double this_incE = (*evt.true_beam_traj_KE)[i];
            if (this_sliceID>=pi::nthinslices) continue;
            if (this_sliceID<0) continue;
            vincE[this_sliceID].push_back(this_incE);
          }
          for (size_t i = 0; i<vincE.size(); ++i){
            if (!vincE[i].empty()){
              double sum_incE = 0;
              for (size_t j = 0; j<vincE[i].size(); ++j){
                sum_incE += vincE[i][j];
              }
              true_incE[i]->Fill(sum_incE/vincE[i].size());
            }
          }
          TVector3 pt0(evt.true_beam_startX,
                       evt.true_beam_startY,
                       evt.true_beam_startZ);
          TVector3 pt1(evt.true_beam_endX,
                       evt.true_beam_endY,
                       evt.true_beam_endZ);
          TVector3 dir = pt1 - pt0;
          dir = dir.Unit();
          true_AngCorr->Fill(dir.Z());
        }
      }
    }

    if (!evt.reco_beam_calo_wire_allTrack->empty()){
      //reco_sliceID = int(evt.reco_beam_calo_endZ_allTrack/pi::thinslicewidth);
      reco_sliceID = int(hadana.reco_trklen/pi::thinslicewidth);
      if (reco_sliceID < 0) reco_sliceID = -1;
      if (evt.reco_beam_calo_endZ_allTrack < 0) reco_sliceID = -1;
      if (reco_sliceID >= pi::nthinslices) reco_sliceID = pi::nthinslices;
    }
  }
  else {//not using all track reconstruction
    if (evt.MC){
      //true_sliceID = int(evt.true_beam_endZ/pi::thinslicewidth);
      true_sliceID = int(hadana.true_trklen/pi::thinslicewidth);
      if (true_sliceID < 0) true_sliceID = -1;
      if (evt.true_beam_endZ < 0) true_sliceID = -1;
      if (true_sliceID >= pi::nthinslices) true_sliceID = pi::nthinslices;
      if (evt.true_beam_PDG == 211){
        for (int i = 0; i<=true_sliceID; ++i){
          if (i<pi::nthinslices) ++true_incidents[i]; // count incident events
        }
      }
      if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
        if (true_sliceID < pi::nthinslices && true_sliceID>=0){
          ++true_interactions[true_sliceID]; // count interaction events
        }
        // Reco info
        if (!(evt.reco_beam_calo_wire->empty()) && evt.reco_beam_true_byE_matched){ // truth matched. so it must be the true track?
          std::vector<std::vector<double>> vincE(pi::nthinslices);
          for (size_t i = 0; i<evt.reco_beam_calo_wire->size(); ++i){
            //int this_sliceID = int((*evt.reco_beam_calo_Z)[i]/pi::thinslicewidth);
            int this_sliceID = int((hadana.reco_trklen_accum)[i]/pi::thinslicewidth);
            if (this_sliceID>=pi::nthinslices) continue;
            if (this_sliceID<0) continue;
            double this_incE = (*evt.reco_beam_incidentEnergies)[i];
            vincE[this_sliceID].push_back(this_incE);
          }
          for (size_t i = 0; i<vincE.size(); ++i){
            if (!vincE[i].empty()){
              double sum_incE = 0;
              for (size_t j = 0; j<vincE[i].size(); ++j){
                sum_incE += vincE[i][j];
              }
              reco_incE[i]->Fill(sum_incE/vincE[i].size());
            }
          }
          TVector3 pt0(evt.reco_beam_calo_startX,
                       evt.reco_beam_calo_startY,
                       evt.reco_beam_calo_startZ);
          TVector3 pt1(evt.reco_beam_calo_endX,
                       evt.reco_beam_calo_endY,
                       evt.reco_beam_calo_endZ);
          TVector3 dir = pt1 - pt0; // direction of the track is determine by the start/end point
          dir = dir.Unit();
          reco_AngCorr->Fill(dir.Z()); // projection to Z of the direction of the track
        }

        // True info
        if (!(evt.true_beam_traj_Z->empty())){
          std::vector<std::vector<double>> vincE(pi::nthinslices);
          for (size_t i = 0; i<evt.true_beam_traj_Z->size()-1; ++i){//last point always has KE = 0
            //int this_sliceID = int((*evt.true_beam_traj_Z)[i]/pi::thinslicewidth);
            int this_sliceID = int((hadana.true_trklen_accum)[i]/pi::thinslicewidth);
            double this_incE = (*evt.true_beam_traj_KE)[i];
            if (this_sliceID>=pi::nthinslices) continue;
            if (this_sliceID<0) continue;
            vincE[this_sliceID].push_back(this_incE);
          }
          for (size_t i = 0; i<vincE.size(); ++i){
            if (!vincE[i].empty()){
              double sum_incE = 0;
              for (size_t j = 0; j<vincE[i].size(); ++j){
                sum_incE += vincE[i][j];
              }
              true_incE[i]->Fill(sum_incE/vincE[i].size());
            }
          }
          TVector3 pt0(evt.true_beam_startX,
                       evt.true_beam_startY,
                       evt.true_beam_startZ);
          TVector3 pt1(evt.true_beam_endX,
                       evt.true_beam_endY,
                       evt.true_beam_endZ);
          TVector3 dir = pt1 - pt0;
          dir = dir.Unit();
          true_AngCorr->Fill(dir.Z());
        }
      }
    }

    if (!evt.reco_beam_calo_wire->empty()){
      //reco_sliceID = int(evt.reco_beam_calo_endZ/pi::thinslicewidth);
      reco_sliceID = int(hadana.reco_trklen/pi::thinslicewidth);
      if (reco_sliceID < 0) reco_sliceID = -1;
      if (evt.reco_beam_calo_endZ < 0) reco_sliceID = -1;
      if (reco_sliceID >= pi::nthinslices) reco_sliceID = pi::nthinslices;
    }
  }

  if (evt.MC){
    if (evt.true_beam_PDG == 211){ // true pion beam incident event
      if (isTestSample){ // fake data
        h_truesliceid_pion_all->Fill(true_sliceID);
      }
      else{
        uf.eff_den_Inc->Fill(true_sliceID);
      }
      if (hadana.PassPiCuts(evt) && evt.reco_beam_true_byE_matched){ // the beam pion passed full selections (reco_beam_true_byE_matched is used to veto secondary particles)
        if (isTestSample){
          h_recosliceid_pion_cuts->Fill(reco_sliceID);
          h_truesliceid_pion_cuts->Fill(true_sliceID);
        }
        else{
          uf.eff_num_Inc->Fill(true_sliceID);
          uf.pur_num_Inc->Fill(reco_sliceID);
          uf.response_SliceID_Inc.Fill(reco_sliceID, true_sliceID);
        }
      }
      else { // this beam pion event is not selected
        if (!isTestSample) uf.response_SliceID_Inc.Miss(true_sliceID);
      }
      
      if ((*evt.true_beam_endProcess) == "pi+Inelastic"){ // true pion beam interaction event (exclude elastics)
        if (isTestSample){
          h_truesliceid_pioninelastic_all->Fill(true_sliceID);
        }
        else{
          uf.eff_den_Int->Fill(true_sliceID);
        }
        if (hadana.PassPiCuts(evt) && evt.reco_beam_true_byE_matched){
          if (isTestSample){
            h_recosliceid_pioninelastic_cuts->Fill(reco_sliceID);
            h_truesliceid_pioninelastic_cuts->Fill(true_sliceID);
          }
          else{
            uf.eff_num_Int->Fill(true_sliceID);
            uf.pur_num_Int->Fill(reco_sliceID);
            uf.response_SliceID_Int.Fill(reco_sliceID, true_sliceID);
          }
        }
        else{
          if (!isTestSample) uf.response_SliceID_Int.Miss(true_sliceID);
        }
      }
    }
    if (hadana.PassPiCuts(evt)){ // the event passed full selections
      if (isTestSample){
        h_recosliceid_allevts_cuts->Fill(reco_sliceID);
      }
      else {
        uf.pur_den->Fill(reco_sliceID);
      }
    }
  }
}

void ThinSlice::FillHistograms(int cut, const anavar & evt){
  if (hadana.fAllTrackCheck) {
    if (cut>=0 && cut < pi::nCuts){
      FillHistVec1D(hreco_beam_type[cut], evt.reco_beam_type, hadana.pitype);
      FillHistVec1D(htrklen_csda_proton[cut], hadana.trklen_csda_proton, hadana.pitype);
      FillHistVec1D(hChi2_proton[cut], hadana.chi2_proton, hadana.pitype);
      FillHistVec1D(hreco_reconstructable_beam_event[cut], evt.reco_reconstructable_beam_event, hadana.pitype);
      
      FillHistVec1D(htrue_beam_endZ[cut], evt.true_beam_endZ_SCE, hadana.pitype);
      FillHistVec1D(htrue_beam_endZ_SCE[cut], evt.true_beam_endZ, hadana.pitype); // it seems SCE is reversed? and I didn't find true_beam_endZ_SCE on wiki?
      FillHistVec1D(htrue_sliceID[cut], true_sliceID, hadana.pitype);
      //    if (!evt.reco_beam_calo_wire_allTrack->empty()){
      FillHistVec1D(hreco_beam_endZ[cut], evt.reco_beam_allTrack_endZ, hadana.pitype);
      FillHistVec1D(hreco_true_beam_endZ[cut], evt.reco_beam_allTrack_endZ - evt.true_beam_endZ_SCE, hadana.pitype);
      FillHistVec2D(hreco_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_allTrack_endZ, hadana.pitype);
      FillHistVec2D(hreco_true_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_allTrack_endZ - evt.true_beam_endZ_SCE, hadana.pitype);
      
      FillHistVec1D(hreco_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ_allTrack, hadana.pitype);
      FillHistVec1D(hreco_true_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ_allTrack - evt.true_beam_endZ, hadana.pitype);
      FillHistVec2D(hreco_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ_allTrack, hadana.pitype);
      FillHistVec2D(hreco_true_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ_allTrack - evt.true_beam_endZ, hadana.pitype);
      
      FillHistVec1D(hreco_sliceID[cut], reco_sliceID, hadana.pitype);
      FillHistVec1D(hreco_true_sliceID[cut], reco_sliceID - true_sliceID, hadana.pitype);
      FillHistVec2D(hreco_vs_true_sliceID[cut], true_sliceID, reco_sliceID, hadana.pitype);
      FillHistVec2D(hreco_true_vs_true_sliceID[cut], true_sliceID, reco_sliceID - true_sliceID, hadana.pitype);
      
      // below are variables not provided by evt directly (calculated in hadana)
      FillHistVec1D(hmediandEdx[cut], hadana.median_dEdx, hadana.pitype);
      FillHistVec1D(hdaughter_michel_score[cut], hadana.daughter_michel_score, hadana.pitype);
      if (evt.reco_beam_calo_endZ_allTrack>300 && hadana.median_dEdx<2.4){ // likely to be a cosmic muon?
        if (hadana.daughter_michel_score>=0){
          FillHistVec1D(hdaughter_michel_scoreMu[cut], hadana.daughter_michel_score, hadana.pitype);
          //if (!evt.MC) cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<" "<<hadana.daughter_michel_score<<" "<<evt.reco_beam_calo_wire_allTrack->back()<<" "<<evt.reco_beam_calo_tick->back()<<" "<<evt.reco_beam_calo_wire_allTrack->front()<<" "<<evt.reco_beam_calo_tick->front()<<endl;
        }
        int nhits = 0;
        double michelscore = 0;
        for (size_t i = 0; i<evt.reco_daughter_PFP_michelScore_collection->size(); ++i){
          nhits += (*evt.reco_daughter_PFP_nHits_collection)[i];
          michelscore += (*evt.reco_daughter_PFP_michelScore_collection)[0] * (*evt.reco_daughter_PFP_nHits_collection)[i];
        }
        if (nhits && michelscore>=0){
          michelscore/=nhits;
          FillHistVec1D(hdaughter_michel_score2Mu[cut], michelscore, hadana.pitype); // what's PFP and what's difference between hdaughter_michel_scoreMu and hdaughter_michel_score2Mu?
        }
  //      if (hadana.pitype == kMuon && hadana.daughter_michel_score < 0.01){
  //        cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<endl;
  //      }
      }
      if (evt.reco_beam_calo_endZ_allTrack<100 && hadana.median_dEdx<2.4){
        if (hadana.daughter_michel_score>=0){
          //if (!evt.MC) cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<" "<<hadana.daughter_michel_score<<" "<<evt.reco_beam_calo_wire_allTrack->back()<<" "<<evt.reco_beam_calo_tick->back()<<" "<<evt.reco_beam_calo_wire_allTrack->front()<<" "<<evt.reco_beam_calo_tick->front()<<endl;
          FillHistVec1D(hdaughter_michel_scorePi[cut], hadana.daughter_michel_score, hadana.pitype);
        }
      }
      /*if (reco_sliceID>=0 && reco_sliceID<pi::nthinslices){
        FillHistVec1D(hmediandEdxSlice[reco_sliceID][cut], hadana.median_dEdx, hadana.pitype);
        FillHistVec1D(hChi2_protonSlice[reco_sliceID][cut], hadana.chi2_proton, hadana.pitype);
        FillHistVec1D(hdaughter_michel_scoreSlice[reco_sliceID][cut], hadana.daughter_michel_score, hadana.pitype);
        FillHistVec1D(hcosthetaSlice[reco_sliceID][cut], hadana.beam_costh, hadana.pitype);
      }*/

      FillHistVec1D(htrackscore[cut], evt.reco_beam_PFP_trackScore_collection, hadana.pitype);
      FillHistVec1D(hemscore[cut], evt.reco_beam_PFP_emScore_collection, hadana.pitype);
  //    if (cut == kAPA3 && evt.reco_beam_PFP_emScore_collection > 0.9){
  //      cout<<evt.run<<" "<<evt.event<<" "<<evt.reco_beam_PFP_emScore_collection<<" "<<evt.reco_beam_calo_wire_allTrack->front()<<" "<<evt.reco_beam_calo_tick->back()<<endl;
  //    }
      FillHistVec1D(hdEdx_5cm[cut], hadana.dEdx_5cm, hadana.pitype);

      FillHistVec1D(hdeltax[cut], hadana.beam_dx, hadana.pitype);
      FillHistVec1D(hdeltay[cut], hadana.beam_dy, hadana.pitype);
      FillHistVec1D(hdeltaz[cut], hadana.beam_dz, hadana.pitype);
      FillHistVec1D(hcostheta[cut], hadana.beam_costh, hadana.pitype);

      FillHistVec1D(hreco_beam_true_byE_matched[cut], evt.reco_beam_true_byE_matched, hadana.pitype);
      FillHistVec1D(hreco_trklen[cut], hadana.reco_trklen, hadana.pitype);
      FillHistVec1D(htrue_trklen[cut], hadana.true_trklen, hadana.pitype);
      FillHistVec1D(hdiff_trklen[cut], hadana.reco_trklen - hadana.true_trklen, hadana.pitype);
      FillHistVec2D(hreco_vs_true_trklen[cut], hadana.true_trklen, hadana.reco_trklen, hadana.pitype);

      FillHistVec1D(hreco_beam_startX_SCE[cut], evt.reco_beam_calo_startX_allTrack, hadana.pitype);
      FillHistVec1D(hreco_beam_startY_SCE[cut], evt.reco_beam_calo_startY_allTrack, hadana.pitype);
      FillHistVec1D(hreco_beam_startZ_SCE[cut], evt.reco_beam_calo_startZ_allTrack, hadana.pitype);

      if (!evt.reco_beam_calo_wire_allTrack->empty()){
        TVector3 pt0(evt.reco_beam_calo_startX_allTrack,
                     evt.reco_beam_calo_startY_allTrack,
                     evt.reco_beam_calo_startZ_allTrack);
        TVector3 pt1(evt.reco_beam_calo_endX_allTrack,
                     evt.reco_beam_calo_endY_allTrack,
                     evt.reco_beam_calo_endZ_allTrack);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        FillHistVec1D(hreco_beam_dcosX_SCE[cut], dir.X(), hadana.pitype);
        FillHistVec1D(hreco_beam_dcosY_SCE[cut], dir.Y(), hadana.pitype);
        FillHistVec1D(hreco_beam_dcosZ_SCE[cut], dir.Z(), hadana.pitype);
        FillHistVec1D(hreco_beam_angleX_SCE[cut], acos(dir.X())*180/TMath::Pi(), hadana.pitype);
        FillHistVec1D(hreco_beam_angleY_SCE[cut], acos(dir.Y())*180/TMath::Pi(), hadana.pitype);
        FillHistVec1D(hreco_beam_angleZ_SCE[cut], acos(dir.Z())*180/TMath::Pi(), hadana.pitype);
      }

      FillHistVec2D(hreco_beam_startXY_SCE[cut], evt.reco_beam_calo_startX_allTrack, evt.reco_beam_calo_startY_allTrack, hadana.pitype);

    }
  }
  
  else{
    if (cut>=0 && cut < pi::nCuts){
      FillHistVec1D(hreco_beam_type[cut], evt.reco_beam_type, hadana.pitype);
      FillHistVec1D(htrklen_csda_proton[cut], hadana.trklen_csda_proton, hadana.pitype);
      FillHistVec1D(hChi2_proton[cut], hadana.chi2_proton, hadana.pitype);
      FillHistVec1D(hreco_reconstructable_beam_event[cut], evt.reco_reconstructable_beam_event, hadana.pitype);
      
      FillHistVec1D(htrue_beam_endZ[cut], evt.true_beam_endZ_SCE, hadana.pitype);
      FillHistVec1D(htrue_beam_endZ_SCE[cut], evt.true_beam_endZ, hadana.pitype); // it seems SCE is reversed? and I didn't find true_beam_endZ_SCE on wiki?
      FillHistVec1D(htrue_sliceID[cut], true_sliceID, hadana.pitype);
      //    if (!evt.reco_beam_calo_wire->empty()){
      FillHistVec1D(hreco_beam_endZ[cut], evt.reco_beam_endZ, hadana.pitype);
      FillHistVec1D(hreco_true_beam_endZ[cut], evt.reco_beam_endZ - evt.true_beam_endZ_SCE, hadana.pitype);
      FillHistVec2D(hreco_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ, hadana.pitype);
      FillHistVec2D(hreco_true_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ - evt.true_beam_endZ_SCE, hadana.pitype);
      
      FillHistVec1D(hreco_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ, hadana.pitype);
      FillHistVec1D(hreco_true_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ - evt.true_beam_endZ, hadana.pitype);
      FillHistVec2D(hreco_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ, hadana.pitype);
      FillHistVec2D(hreco_true_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ - evt.true_beam_endZ, hadana.pitype);
      
      FillHistVec1D(hreco_sliceID[cut], reco_sliceID, hadana.pitype);
      FillHistVec1D(hreco_true_sliceID[cut], reco_sliceID - true_sliceID, hadana.pitype);
      FillHistVec2D(hreco_vs_true_sliceID[cut], true_sliceID, reco_sliceID, hadana.pitype);
      FillHistVec2D(hreco_true_vs_true_sliceID[cut], true_sliceID, reco_sliceID - true_sliceID, hadana.pitype);
      
      // below are variables not provided by evt directly (calculated in hadana)
      FillHistVec1D(hmediandEdx[cut], hadana.median_dEdx, hadana.pitype);
      FillHistVec1D(hdaughter_michel_score[cut], hadana.daughter_michel_score, hadana.pitype);
      FillHistVec1D(henergy_calorimetry_SCE[cut], hadana.energy_calorimetry_SCE, hadana.pitype);
      FillHistVec1D(hdEdx_SCE[cut], hadana.energy_calorimetry_SCE/hadana.reco_trklen, hadana.pitype);
      if (evt.reco_beam_calo_endZ>300 && hadana.median_dEdx<2.4){ // likely to be a cosmic muon?
        if (hadana.daughter_michel_score>=0){
          FillHistVec1D(hdaughter_michel_scoreMu[cut], hadana.daughter_michel_score, hadana.pitype);
          //if (!evt.MC) cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<" "<<hadana.daughter_michel_score<<" "<<evt.reco_beam_calo_wire->back()<<" "<<evt.reco_beam_calo_tick->back()<<" "<<evt.reco_beam_calo_wire->front()<<" "<<evt.reco_beam_calo_tick->front()<<endl;
        }
        int nhits = 0;
        double michelscore = 0;
        for (size_t i = 0; i<evt.reco_daughter_PFP_michelScore_collection->size(); ++i){
          nhits += (*evt.reco_daughter_PFP_nHits_collection)[i];
          michelscore += (*evt.reco_daughter_PFP_michelScore_collection)[0] * (*evt.reco_daughter_PFP_nHits_collection)[i];
        }
        if (nhits && michelscore>=0){
          michelscore/=nhits;
          FillHistVec1D(hdaughter_michel_score2Mu[cut], michelscore, hadana.pitype); // what's PFP and what's difference between hdaughter_michel_scoreMu and hdaughter_michel_score2Mu?
        }
  //      if (hadana.pitype == kMuon && hadana.daughter_michel_score < 0.01){
  //        cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<endl;
  //      }
      }
      if (evt.reco_beam_calo_endZ<100 && hadana.median_dEdx<2.4){
        if (hadana.daughter_michel_score>=0){
          //if (!evt.MC) cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<" "<<hadana.daughter_michel_score<<" "<<evt.reco_beam_calo_wire->back()<<" "<<evt.reco_beam_calo_tick->back()<<" "<<evt.reco_beam_calo_wire->front()<<" "<<evt.reco_beam_calo_tick->front()<<endl;
          FillHistVec1D(hdaughter_michel_scorePi[cut], hadana.daughter_michel_score, hadana.pitype);
        }
      }
      /*if (reco_sliceID>=0 && reco_sliceID<pi::nthinslices){
        FillHistVec1D(hmediandEdxSlice[reco_sliceID][cut], hadana.median_dEdx, hadana.pitype);
        FillHistVec1D(hChi2_protonSlice[reco_sliceID][cut], hadana.chi2_proton, hadana.pitype);
        FillHistVec1D(hdaughter_michel_scoreSlice[reco_sliceID][cut], hadana.daughter_michel_score, hadana.pitype);
        FillHistVec1D(hcosthetaSlice[reco_sliceID][cut], hadana.beam_costh, hadana.pitype);
      }*/

      FillHistVec1D(htrackscore[cut], evt.reco_beam_PFP_trackScore_collection, hadana.pitype);
      FillHistVec1D(hemscore[cut], evt.reco_beam_PFP_emScore_collection, hadana.pitype);
  //    if (cut == kAPA3 && evt.reco_beam_PFP_emScore_collection > 0.9){
  //      cout<<evt.run<<" "<<evt.event<<" "<<evt.reco_beam_PFP_emScore_collection<<" "<<evt.reco_beam_calo_wire->front()<<" "<<evt.reco_beam_calo_tick->back()<<endl;
  //    }
      FillHistVec1D(hdEdx_5cm[cut], hadana.dEdx_5cm, hadana.pitype);

      FillHistVec1D(hdeltax[cut], hadana.beam_dx, hadana.pitype);
      FillHistVec1D(hdeltay[cut], hadana.beam_dy, hadana.pitype);
      FillHistVec1D(hdeltaz[cut], hadana.beam_dz, hadana.pitype);
      FillHistVec1D(hcostheta[cut], hadana.beam_costh, hadana.pitype);

      FillHistVec1D(hreco_beam_true_byE_matched[cut], evt.reco_beam_true_byE_matched, hadana.pitype);
      FillHistVec1D(hreco_trklen[cut], hadana.reco_trklen, hadana.pitype);
      FillHistVec1D(htrue_trklen[cut], hadana.true_trklen, hadana.pitype);
      FillHistVec1D(hdiff_trklen[cut], hadana.reco_trklen - hadana.true_trklen, hadana.pitype);
      FillHistVec2D(hreco_vs_true_trklen[cut], hadana.true_trklen, hadana.reco_trklen, hadana.pitype);
      FillHistVec1D(hbeam_score[cut], hadana.beam_score, hadana.pitype);
      FillHistVec2D(beam_score_vs_hreco_trklen[cut], hadana.reco_trklen, hadana.beam_score, hadana.pitype);
      
      //$$$temp
     /* if ( hadana.true_trklen>20 && evt.reco_beam_alt_len>20){
        int printout = kFALSE;
        if ( hadana.true_trklen>250 && abs(evt.reco_beam_alt_len-230)<5 ){
          cout<<"$$$$$ red bar ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          printout = kTRUE;
        }
        if ( hadana.true_trklen<200 && abs(evt.reco_beam_alt_len-230)<5 ){
          cout<<"$$$$$ blue bar ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          printout = kTRUE;
        }
        if ( abs(hadana.true_trklen-evt.reco_beam_alt_len-230)<5 ){
          cout<<"$$$$$ green bar ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          printout = kTRUE;
        }
        if (printout == kTRUE){
          cout<<"Run: "<<evt.run<<";\t"
          <<"SubRun: "<<evt.subrun<<";\t"
          <<"Event: "<<evt.event<<endl;
          cout<<"True trklen: "<<hadana.true_trklen<<";\t"
          <<"Reco trklen: "<<evt.reco_beam_alt_len<<endl;
          cout<<"Start: ("<<evt.reco_beam_calo_startX<<", "<<evt.reco_beam_calo_startY<<", "<<evt.reco_beam_calo_startZ<<");\n";
          cout<<"End: ("<<evt.reco_beam_calo_endX<<", "<<evt.reco_beam_calo_endY<<", "<<evt.reco_beam_calo_endZ<<")\n";
        }
      }*/
      /*if (!evt.MC){
        if (hadana.beam_dy>4 && hadana.beam_dy<6){
          cout<<"$$$$$ delta Y [4,6] ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          cout<<"Run: "<<evt.run<<";\t"
          <<"SubRun: "<<evt.subrun<<";\t"
          <<"Event: "<<evt.event<<endl;
          cout<<"Start: ("<<evt.reco_beam_calo_startX<<", "<<evt.reco_beam_calo_startY<<", "<<evt.reco_beam_calo_startZ<<");\t trklen: "<<hadana.reco_trklen<<endl;
          cout<<"End: ("<<evt.reco_beam_calo_endX<<", "<<evt.reco_beam_calo_endY<<", "<<evt.reco_beam_calo_endZ<<")\n";
          
        }
      }*/

      FillHistVec1D(hreco_beam_startX_SCE[cut], evt.reco_beam_calo_startX, hadana.pitype);
      FillHistVec1D(hreco_beam_startY_SCE[cut], evt.reco_beam_calo_startY, hadana.pitype);
      FillHistVec1D(hreco_beam_startZ_SCE[cut], evt.reco_beam_calo_startZ, hadana.pitype);

      if (!evt.reco_beam_calo_wire->empty()){
        TVector3 pt0(evt.reco_beam_calo_startX,
                     evt.reco_beam_calo_startY,
                     evt.reco_beam_calo_startZ);
        TVector3 pt1(evt.reco_beam_calo_endX,
                     evt.reco_beam_calo_endY,
                     evt.reco_beam_calo_endZ);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        FillHistVec1D(hreco_beam_dcosX_SCE[cut], dir.X(), hadana.pitype);
        FillHistVec1D(hreco_beam_dcosY_SCE[cut], dir.Y(), hadana.pitype);
        FillHistVec1D(hreco_beam_dcosZ_SCE[cut], dir.Z(), hadana.pitype);
        FillHistVec1D(hreco_beam_angleX_SCE[cut], acos(dir.X())*180/TMath::Pi(), hadana.pitype);
        FillHistVec1D(hreco_beam_angleY_SCE[cut], acos(dir.Y())*180/TMath::Pi(), hadana.pitype);
        FillHistVec1D(hreco_beam_angleZ_SCE[cut], acos(dir.Z())*180/TMath::Pi(), hadana.pitype);
      }

      FillHistVec2D(hreco_beam_startXY_SCE[cut], evt.reco_beam_calo_startX, evt.reco_beam_calo_startY, hadana.pitype);

    }
  }
}

void ThinSlice::FillSliceHist(const anavar & evt, int constraint_type, int cut){
  if (hadana.fAllTrackCheck) {
    cout<<"Warning: AllTrackCheck hasn't been fully implemented!!!"<<endl; // Do we need AllTrackCheck in main branch?
  }
  else {
    if (constraint_type == 1) { // muon
      FillHistVec1D(hdaughter_michel_score_bkg[cut], hadana.daughter_michel_score, hadana.pitype, false, false);
    }
    else if (constraint_type == 2) { // proton
      FillHistVec1D(hmediandEdx_bkg[cut], hadana.median_dEdx, hadana.pitype, false, true);
      FillHistVec1D(hChi2_proton_bkg[cut], hadana.chi2_proton, hadana.pitype, false, true);
    }
    else if (constraint_type == 3) { // secondary pion
      FillHistVec1D(hcostheta_bkg[cut], hadana.beam_costh, hadana.pitype, true, false);
    }
    // in each slice
    if (reco_sliceID>=0 && reco_sliceID<pi::nthinslices){
      if (constraint_type == 1) { // muon
        FillHistVec1D(hdaughter_michel_scoreSlice[reco_sliceID][cut], hadana.daughter_michel_score, hadana.pitype, false, false);
      }
      else if (constraint_type == 2) { // proton
        FillHistVec1D(hmediandEdxSlice[reco_sliceID][cut], hadana.median_dEdx, hadana.pitype, false, true);
        FillHistVec1D(hChi2_protonSlice[reco_sliceID][cut], hadana.chi2_proton, hadana.pitype, false, true);
      }
      else if (constraint_type == 3) { // secondary pion
        FillHistVec1D(hcosthetaSlice[reco_sliceID][cut], hadana.beam_costh, hadana.pitype, true, false);
      }
    }
  }
}

void ThinSlice::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
  h_truesliceid_pion_uf->Write("h_truesliceid_pion_uf");
  h_truesliceid_pioninelastic_uf->Write("h_truesliceid_pioninelastic_uf");
  //response_SliceID_Pion->Write("response_SliceID_Pion");
  //response_SliceID_PionInEl->Write("response_SliceID_PionInEl");
}

void ThinSlice::CalcXS(const Unfold & uf){

  double slcid[pi::nthinslices] = {0};
  double avg_trueincE[pi::nthinslices] = {0};
  double avg_recoincE[pi::nthinslices] = {0};
  double err_trueincE[pi::nthinslices] = {0};
  double err_recoincE[pi::nthinslices] = {0};
  double reco_trueincE[pi::nthinslices] = {0};
  double err_reco_trueincE[pi::nthinslices] = {0};
  double truexs[pi::nthinslices] = {0};
  double err_truexs[pi::nthinslices] = {0};
  double true_cosangle = 1.;

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // 1.396 g/cm^3

  for (int i = 0; i<pi::nthinslices; ++i){
    
    slcid[i] = i;
    avg_trueincE[i] = true_incE[i]->GetMean();
    err_trueincE[i] = true_incE[i]->GetMeanError();
    avg_recoincE[i] = reco_incE[i]->GetMean();
    err_recoincE[i] = reco_incE[i]->GetMeanError();
    reco_trueincE[i] = avg_recoincE[i] - avg_trueincE[i];
    err_reco_trueincE[i] = sqrt(pow(err_trueincE[i],2)+pow(err_recoincE[i],2)); // is it proper to simply use root_sum_square, since the two seem not independent?
    //std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
    if (true_incidents[i] && true_interactions[i]){
      //true_cosangle = true_AngCorr->GetMean(); // no need to include angle correction
      truexs[i] = MAr/(Density*NA*pi::thinslicewidth/true_cosangle)*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs[i] = MAr/(Density*NA*pi::thinslicewidth/true_cosangle)*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
  }

  TGraphErrors *gr_trueincE = new TGraphErrors(pi::nthinslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
  TGraphErrors *gr_recoincE = new TGraphErrors(pi::nthinslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
  TGraphErrors *gr_reco_trueincE = new TGraphErrors(pi::nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

  gr_trueincE->Write("gr_trueincE");
  gr_recoincE->Write("gr_recoincE");
  gr_reco_trueincE->Write("gr_reco_trueincE");

  TGraphErrors *gr_truexs = new TGraphErrors(pi::nthinslices, &(avg_trueincE[0]), &(truexs[0]), 0, &(err_truexs[0]));
  
  gr_truexs->Write("gr_truexs");

  TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
  TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
  hinc->Multiply(uf.pur_Inc);
  hint->Multiply(uf.pur_Int);

//  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 4);
//  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 4);

  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 4);

//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, hinc, 20);   // OR
//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, hint, 20);   // OR

//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 20);   // OR
//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 20);   // OR

  h_truesliceid_pion_uf = (TH1D*) unfold_Inc.Hreco();
  h_truesliceid_pioninelastic_uf = (TH1D*) unfold_Int.Hreco();

}

void ThinSlice::Run(anavar & evt, Unfold & uf, Long64_t nentries=-1){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (selectCosmics){
      if (!hadana.isCosmics(evt)) continue;
    }
    else{
      if (!hadana.isSelectedPart(evt)) continue;
    }
    ProcessEvent(evt, uf);
    // can change order of cuts
    FillHistograms(pi::kNocut, evt);
    if (hadana.PassPandoraSliceCut(evt)){
      FillHistograms(pi::kPandoraSlice, evt);
      if (hadana.PassCaloSizeCut(evt)){
        FillHistograms(pi::kCaloSize, evt);
        if (hadana.PassBeamQualityCut()){
          FillHistograms(pi::kBeamQuality, evt);
          if (hadana.PassProtonCut()){
            FillHistograms(pi::kProtonCut, evt);
            if (hadana.PassMichelScoreCut()){
              FillHistograms(pi::kMichelScore, evt);
              if (hadana.PassAPA3Cut(evt)){
                FillHistograms(pi::kAPA3, evt);
              }
            }
          }
        }
        // for background constraints
        if (hadana.PassAPA3Cut(evt)){
          if (hadana.PassBeamQualityCut() && hadana.PassProtonCut()) { // to constrain muon
            FillSliceHist(evt, 1);
          }
          if (hadana.PassBeamQualityCut() && hadana.PassMichelScoreCut()) { // to constrain proton
            FillSliceHist(evt, 2);
          }
          if (hadana.PassBeamQualityCut(false) && hadana.PassProtonCut() && hadana.PassMichelScoreCut()) { // to constrain secondary pion
            FillSliceHist(evt, 3);
          }
        }
      }
    }
    FillSliceHist(evt, 1, 0);
    FillSliceHist(evt, 2, 0);
    FillSliceHist(evt, 3, 0);
  }
  
  uf.SaveHistograms();
  CalcXS(uf);
  SaveHistograms();
}

void ThinSlice::SetSelectCosmics(bool sc){

  selectCosmics = sc;

}
