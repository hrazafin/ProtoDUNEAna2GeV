#include "EffEval.h"
#include "HadAna.h"
#include "anavar.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "util.h"

#include <iostream>

EffEval::EffEval(){
  hadana.InitPi();
}

void EffEval::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");
  h_true_Ppi_all = new TH1D("h_true_Ppi_all","All;True p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_all->Sumw2();

  h_true_Ppi_sel = new TH1D("h_true_Ppi_sel","Selected;True p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_sel->Sumw2();

  h_reco_true_Ppi_sel = new TH2D("h_reco_true_Ppi_sel","Selected;True p_{#pi} (MeV/c); Reco/True -1", 120,0,1200,40,-1,1);
  h_res_Ppi_sel = new TH1D("h_res_Ppi_sel","Selected; Reco/True -1 (p_{#pi})", 40,-1,1);
  h_res_Ppi_sel->Sumw2();

  h_res_Ppi_michelscore = new TH2D("h_res_Ppi_michelscore",";Michel score;Reco/True -1 (p_{#pi})", 100,0,1,40,-1,1);

  h_true_Ppi_michel = new TH1D("h_true_Ppi_michel","With Michel score cut;True p_{#pi} (MeV/c);Events",120,0,1200);
  h_true_Ppi_michel->Sumw2();

  h_reco_true_Ppi_michel = new TH2D("h_reco_true_Ppi_michel","With Michel score cut;True p_{#pi} (MeV/c); Reco/True -1", 120,0,1200,40,-1,1);

  h_res_Ppi_michel = new TH1D("h_res_Ppi_michel","With Michel score cut; Reco/True -1 (p_{#pi})", 40,-1,1);
  h_res_Ppi_michel->Sumw2();

  ////////////////////////////////

  h_true_thetapi_all = new TH1D("h_true_thetapi_all","All;True #theta_{#pi};Events",180,0,180);
  h_true_thetapi_all->Sumw2();

  h_true_thetapi_sel = new TH1D("h_true_thetapi_sel","Selected;True #theta_{#pi};Events",180,0,180);
  h_true_thetapi_sel->Sumw2();

  h_reco_true_thetapi_sel = new TH2D("h_reco_true_thetapi_sel","Selected;True #theta_{#pi}; Reco-True (#theta_{#pi})", 180,0,180,60,-30,30);

  h_res_thetapi_sel = new TH1D("h_res_thetapi_sel","Selected; Reco - True (#theta_{#pi})", 60,-30,30);
  h_res_thetapi_sel->Sumw2();

  h_true_thetapi_michel = new TH1D("h_true_thetapi_michel","With Michel score cut;True #theta_{#pi};Events",180,0,180);
  h_true_thetapi_michel->Sumw2();

  h_reco_true_thetapi_michel = new TH2D("h_reco_true_thetapi_michel","With Michel score cut;True #theta_{#pi}; Reco - True (#theta_{#pi})", 180,0,180,60,-30,30);

  h_res_thetapi_michel = new TH1D("h_res_thetapi_michel","With Michel score cut; Reco-True (#theta_{#pi})", 60,-30,30);
  h_res_thetapi_michel->Sumw2();

  ////////////////////////////////

  h_true_phipi_all = new TH1D("h_true_phipi_all","All;True #phi_{#pi};Events",180,-180,180);
  h_true_phipi_all->Sumw2();

  h_true_phipi_sel = new TH1D("h_true_phipi_sel","Selected;True #phi_{#pi};Events",180,-180,180);
  h_true_phipi_sel->Sumw2();

  h_reco_true_phipi_sel = new TH2D("h_reco_true_phipi_sel","Selected;True #phi_{#pi}; Reco-True (#phi_{#pi})", 180,-180,180,60,-30,30);

  h_res_phipi_sel = new TH1D("h_res_phipi_sel","Selected; Reco - True (#phi_{#pi})", 60,-30,30);
  h_res_phipi_sel->Sumw2();

  h_true_phipi_michel = new TH1D("h_true_phipi_michel","With Michel score cut;True #phi_{#pi};Events",180,-180,180);
  h_true_phipi_michel->Sumw2();

  h_reco_true_phipi_michel = new TH2D("h_reco_true_phipi_michel","With Michel score cut;True #phi_{#pi}; Reco - True (#phi_{#pi})", 180,-180,180,60,-30,30);

  h_res_phipi_michel = new TH1D("h_res_phipi_michel","With Michel score cut; Reco-True (#phi_{#pi})", 60,-30,30);
  h_res_phipi_michel->Sumw2();

  ////////////////////////////////

  h_true_Pp_all = new TH1D("h_true_Pp_all","All;True p_{p} (MeV/c);Events",120,0,1200);
  h_true_Pp_all->Sumw2();

  h_true_Pp_sel = new TH1D("h_true_Pp_sel","Selected;True p_{p} (MeV/c);Events",120,0,1200);
  h_true_Pp_sel->Sumw2();

  h_reco_true_Pp_sel = new TH2D("h_reco_true_Pp_sel","Selected;True p_{p} (MeV/c); Reco/True -1", 120,0,1200,40,-1,1);
  h_res_Pp_sel = new TH1D("h_res_Pp_sel","Selected; Reco/True -1 (p_{p})", 40,-1,1);
  h_res_Pp_sel->Sumw2();

  ////////////////////////////////

  h_true_thetap_all = new TH1D("h_true_thetap_all","All;True #theta_{p};Events",180,0,180);
  h_true_thetap_all->Sumw2();

  h_true_thetap_sel = new TH1D("h_true_thetap_sel","Selected;True #theta_{p};Events",180,0,180);
  h_true_thetap_sel->Sumw2();

  h_reco_true_thetap_sel = new TH2D("h_reco_true_thetap_sel","Selected;True #theta_{p}; Reco-True (#theta_{p})", 180,0,180,60,-30,30);

  h_res_thetap_sel = new TH1D("h_res_thetap_sel","Selected; Reco - True (#theta_{p})", 60,-30,30);
  h_res_thetap_sel->Sumw2();

  ////////////////////////////////

  h_true_phip_all = new TH1D("h_true_phip_all","All;True #phi_{p};Events",180,-180,180);
  h_true_phip_all->Sumw2();

  h_true_phip_sel = new TH1D("h_true_phip_sel","Selected;True #phi_{p};Events",180,-180,180);
  h_true_phip_sel->Sumw2();

  h_reco_true_phip_sel = new TH2D("h_reco_true_phip_sel","Selected;True #phi_{p}; Reco-True (#phi_{p})", 180,-180,180,60,-30,30);

  h_res_phip_sel = new TH1D("h_res_phip_sel","Selected; Reco - True (#phi_{p})", 60,-30,30);
  h_res_phip_sel->Sumw2();

}

void EffEval::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void EffEval::FillHistograms(const anavar & evt){
  for (size_t i = 0; i<evt.true_beam_daughter_ID->size(); ++i){
    //Only look at pion
    if (std::abs((*evt.true_beam_daughter_PDG)[i]) != 211 &&
        (*evt.true_beam_daughter_PDG)[i] != 2212) continue;
    double true_mom = (*evt.true_beam_daughter_startP)[i]*1000;
    double true_theta = GetTheta((*evt.true_beam_daughter_startPx)[i],
                                 (*evt.true_beam_daughter_startPy)[i],
                                 (*evt.true_beam_daughter_startPz)[i]);
    double true_phi = GetPhi((*evt.true_beam_daughter_startPx)[i],
                                 (*evt.true_beam_daughter_startPy)[i],
                                 (*evt.true_beam_daughter_startPz)[i]);
    bool foundreco = false;
    double besttrack = -1;
    double tracklength = -1;
    double michelscore = -1;
    double reco_mom = -1;
    double res_mom = -1;
    double reco_theta = -1;
    double res_theta = -1;
    double reco_phi = -1;
    double res_phi = -1;
    for (size_t j = 0; j<evt.reco_daughter_PFP_ID->size(); ++j){
      if ((*evt.reco_daughter_PFP_true_byHits_ID)[j] == (*evt.true_beam_daughter_ID)[i]
          && !(*evt.reco_daughter_allTrack_dQdX_SCE)[j].empty()){
        foundreco = true;
        if ((*evt.reco_daughter_PFP_true_byHits_purity)[j]*(*evt.reco_daughter_PFP_true_byHits_completeness)[j]>besttrack){
          besttrack = (*evt.reco_daughter_PFP_true_byHits_purity)[j]*(*evt.reco_daughter_PFP_true_byHits_completeness)[j];
          tracklength = (*evt.reco_daughter_allTrack_alt_len)[j];
          if (std::abs((*evt.true_beam_daughter_PDG)[i]) == 211){
            double reco_KE = GetPionKE(tracklength);
            //reco_mom = sqrt(pow(reco_KE+139.57,2)-pow(139.57,2));
            reco_mom = (*evt.reco_daughter_allTrack_momByRange_alt_muon)[j]*1000;
            res_mom = (reco_mom - true_mom)/true_mom;
          }
          if ((*evt.true_beam_daughter_PDG)[i] == 2212){
            reco_mom = (*evt.reco_daughter_allTrack_momByRange_alt_proton)[j]*1000;
            res_mom = (reco_mom - true_mom)/true_mom;
          }
          reco_theta = (*evt.reco_daughter_allTrack_Theta)[j]*180/TMath::Pi();
          res_theta = reco_theta - true_theta;
          reco_phi = (*evt.reco_daughter_allTrack_Phi)[j]*180/TMath::Pi();
          res_phi = reco_phi - true_phi;
          if ((*evt.reco_daughter_allTrack_vertex_nHits)[j]){
            michelscore = (*evt.reco_daughter_allTrack_vertex_michel_score)[j]/(*evt.reco_daughter_allTrack_vertex_nHits)[j];
            if (michelscore<0) michelscore = 0;
            if (michelscore>=1) michelscore = 0.999999;
          }
        }
      }
    }
    
    if (std::abs((*evt.true_beam_daughter_PDG)[i]) == 211){
      h_true_Ppi_all->Fill(true_mom);
      h_true_thetapi_all->Fill(true_theta);
      h_true_phipi_all->Fill(true_phi);
      if (foundreco){
        h_true_Ppi_sel->Fill(true_mom);
        h_reco_true_Ppi_sel->Fill(true_mom, res_mom);
        h_res_Ppi_sel->Fill(res_mom);
        h_res_Ppi_michelscore->Fill(michelscore, res_mom);

        h_true_thetapi_sel->Fill(true_theta);
        h_reco_true_thetapi_sel->Fill(true_theta, res_theta);
        h_res_thetapi_sel->Fill(res_theta);

        h_true_phipi_sel->Fill(true_phi);
        h_reco_true_phipi_sel->Fill(true_phi, res_phi);
        h_res_phipi_sel->Fill(res_phi);

        if (michelscore>0.5){
          h_true_Ppi_michel->Fill(true_mom);
          h_reco_true_Ppi_michel->Fill(true_mom, res_mom);
          h_res_Ppi_michel->Fill(res_mom);

          h_true_thetapi_michel->Fill(true_theta);
          h_reco_true_thetapi_michel->Fill(true_theta, res_theta);
          h_res_thetapi_michel->Fill(res_theta);

          h_true_phipi_michel->Fill(true_phi);
          h_reco_true_phipi_michel->Fill(true_phi, res_phi);
          h_res_phipi_michel->Fill(res_phi);
        }
      }
    }

    if ((*evt.true_beam_daughter_PDG)[i] == 2212){
      h_true_Pp_all->Fill(true_mom);
      h_true_thetap_all->Fill(true_theta);
      h_true_phip_all->Fill(true_phi);
      if (foundreco){
        h_true_Pp_sel->Fill(true_mom);
        h_reco_true_Pp_sel->Fill(true_mom, res_mom);
        h_res_Pp_sel->Fill(res_mom);

        h_true_thetap_sel->Fill(true_theta);
        h_reco_true_thetap_sel->Fill(true_theta, res_theta);
        h_res_thetap_sel->Fill(res_theta);

        h_true_phip_sel->Fill(true_phi);
        h_reco_true_phip_sel->Fill(true_phi, res_phi);
        h_res_phip_sel->Fill(res_phi);
      }
    }

  }
}

void EffEval::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
//  for (int i = 0; i<h_true_Ppi_sel->GetNbinsX()+2; ++i){
//    cout<<i<<" "<<h_true_Ppi_sel->GetBinContent(i)<<" "<<h_true_Ppi_all->GetBinContent(i)<<endl;
//    if (h_true_Ppi_all->GetBinContent(i)<h_true_Ppi_sel->GetBinContent(i)) cout<<"Inconsistent"<<endl;
//  }
  if (TEfficiency::CheckConsistency(*h_true_Ppi_sel, *h_true_Ppi_all)){
    h_Eff_Ppi = new TEfficiency(*h_true_Ppi_sel, *h_true_Ppi_all);
    h_Eff_Ppi->Write("h_Eff_Ppi");
  }

  if (TEfficiency::CheckConsistency(*h_true_Ppi_michel, *h_true_Ppi_all)){
    h_Eff_Ppi_michel = new TEfficiency(*h_true_Ppi_michel, *h_true_Ppi_all);
    h_Eff_Ppi_michel->Write("h_Eff_Ppi_michel");
  }

  if (TEfficiency::CheckConsistency(*h_true_thetapi_sel, *h_true_thetapi_all)){
    h_Eff_thetapi = new TEfficiency(*h_true_thetapi_sel, *h_true_thetapi_all);
    h_Eff_thetapi->Write("h_Eff_thetapi");
  }

  if (TEfficiency::CheckConsistency(*h_true_thetapi_michel, *h_true_thetapi_all)){
    h_Eff_thetapi_michel = new TEfficiency(*h_true_thetapi_michel, *h_true_thetapi_all);
    h_Eff_thetapi_michel->Write("h_Eff_thetapi_michel");
  }

  if (TEfficiency::CheckConsistency(*h_true_phipi_sel, *h_true_phipi_all)){
    h_Eff_phipi = new TEfficiency(*h_true_phipi_sel, *h_true_phipi_all);
    h_Eff_phipi->Write("h_Eff_phipi");
  }

  if (TEfficiency::CheckConsistency(*h_true_phipi_michel, *h_true_phipi_all)){
    h_Eff_phipi_michel = new TEfficiency(*h_true_phipi_michel, *h_true_phipi_all);
    h_Eff_phipi_michel->Write("h_Eff_phipi_michel");
  }

  ///////////////////////////////////////
  if (TEfficiency::CheckConsistency(*h_true_Pp_sel, *h_true_Pp_all)){
    h_Eff_Pp = new TEfficiency(*h_true_Pp_sel, *h_true_Pp_all);
    h_Eff_Pp->Write("h_Eff_Pp");
  }

  if (TEfficiency::CheckConsistency(*h_true_thetap_sel, *h_true_thetap_all)){
    h_Eff_thetap = new TEfficiency(*h_true_thetap_sel, *h_true_thetap_all);
    h_Eff_thetap->Write("h_Eff_thetap");
  }

  if (TEfficiency::CheckConsistency(*h_true_phip_sel, *h_true_phip_all)){
    h_Eff_phip = new TEfficiency(*h_true_phip_sel, *h_true_phip_all);
    h_Eff_phip->Write("h_Eff_phip");
  }

}

void EffEval::Run(anavar & evt, Long64_t nentries=-1){

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
    if (!hadana.isSelectedPart(evt)) continue;
    ProcessEvent(evt);
    if (!hadana.PassPiCuts(evt)) continue;
    FillHistograms(evt);
  }
  SaveHistograms();
}
