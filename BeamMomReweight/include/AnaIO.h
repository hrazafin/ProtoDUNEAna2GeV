#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "PlotUtils.h"
#include "AnaFunctions.h"
// Kinematic fitting codes
#include "../../../Fitting/UserKF/UserKF.cxx"
#include "../src/BetheBloch.cxx"

namespace AnaIO
{
  //====================== Output Tree Variables ==================//
  double trueBeamMom;
  double recoFFKE;
  double instP;
  
  //====================== Reco (MC and Data)======================//
  //====================== Reco Beam ======================//
  // Declare variables
  Int_t           event = -999;
  Bool_t          reco_reconstructable_beam_event;
  Int_t           reco_beam_type = -999;
  Double_t        reco_beam_startX;
  Double_t        reco_beam_startY;
  Double_t        reco_beam_startZ;
  Double_t        reco_beam_trackDirX;
  Double_t        reco_beam_trackDirY;
  Double_t        reco_beam_trackDirZ;
  Double_t        reco_beam_endX;
  Double_t        reco_beam_endY;
  Double_t        reco_beam_endZ;
  Double_t        reco_beam_trackEndDirX;
  Double_t        reco_beam_trackEndDirY;
  Double_t        reco_beam_trackEndDirZ;
  Double_t        reco_beam_interactingEnergy;
  Double_t        reco_beam_Chi2_proton;
  Int_t           reco_beam_Chi2_ndof;
  vector<double>  *reco_beam_calibrated_dEdX_SCE = 0x0;
  vector<double>  *reco_beam_TrkPitch_SCE = 0x0;
  vector<double>  *reco_beam_TrkPitch_NoSCE = 0x0;
  vector<double>  *reco_beam_incidentEnergies = 0x0;
  vector<double>  *reco_beam_new_incidentEnergies = 0x0;


  // Beam new variables with possible SCE correction
  Double_t        reco_beam_calo_startX;
  Double_t        reco_beam_calo_startY;
  Double_t        reco_beam_calo_startZ;
  Double_t        reco_beam_calo_endX;
  Double_t        reco_beam_calo_endY;
  Double_t        reco_beam_calo_endZ;
  vector<double>  *reco_beam_calo_X;
  vector<double>  *reco_beam_calo_Y;
  vector<double>  *reco_beam_calo_Z;


  Double_t        reco_beam_momByRange_alt_muon;
  Double_t        reco_beam_momByRange_muon;
  vector<double>  *reco_beam_calo_startDirX;
  vector<double>  *reco_beam_calo_startDirY;
  vector<double>  *reco_beam_calo_startDirZ;
  vector<double>  *reco_beam_calo_endDirX;
  vector<double>  *reco_beam_calo_endDirY;
  vector<double>  *reco_beam_calo_endDirZ;
  vector<double>  *reco_beam_calo_wire;

  Double_t        reco_beam_vertex_michel_score;
  Int_t           reco_beam_vertex_nHits;

  Int_t           reco_beam_PFP_nHits;
  Double_t        reco_beam_PFP_trackScore_collection;
  Double_t        reco_beam_PFP_trackScore_collection_weight_by_charge;
  Double_t        reco_beam_PFP_emScore_collection;
  Double_t        reco_beam_PFP_emScore_collection_weight_by_charge;
  Double_t        reco_beam_vertex_michel_score_weight_by_charge;

  Int_t           reco_beam_true_byHits_PDG;
  Int_t           reco_beam_true_byHits_origin;
  Bool_t          reco_beam_true_byHits_matched;
  Double_t        reco_beam_true_byHits_startPx;
  Double_t        reco_beam_true_byHits_startPy;
  Double_t        reco_beam_true_byHits_startPz;
  Double_t        reco_beam_true_byHits_startE;
  Double_t        reco_beam_true_byHits_startP;
  Double_t        reco_beam_true_byHits_endPx;
  Double_t        reco_beam_true_byHits_endPy;
  Double_t        reco_beam_true_byHits_endPz;
  Double_t        reco_beam_true_byHits_endE;
  Double_t        reco_beam_true_byHits_endP;

  string         *reco_beam_true_byHits_process;
  string         *reco_beam_true_byHits_endProcess;

  Bool_t          reco_beam_true_byE_matched;
  Int_t           reco_beam_true_byE_PDG;
  
  //====================== Reco Daughter ======================//

  vector<int>             *reco_daughter_PFP_ID = 0x0;
  vector<int>             *reco_daughter_PFP_nHits = 0x0; 
  vector<int>             *reco_daughter_PFP_nHits_collection = 0x0;
  vector<double>          *reco_daughter_PFP_trackScore = 0x0;
  vector<double>          *reco_daughter_PFP_trackScore_collection = 0x0;
  vector<double>          *reco_daughter_PFP_emScore = 0x0;
  vector<double>          *reco_daughter_PFP_emScore_collection = 0x0;
  vector<double>          *reco_daughter_PFP_michelScore = 0x0;
  vector<double>          *reco_daughter_PFP_michelScore_collection =0x0;
  vector<int>             *reco_daughter_allTrack_ID = 0x0;
  vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX_SCE = 0x0;
  vector<double>          *reco_daughter_allTrack_Chi2_proton = 0x0;
  vector<int>             *reco_daughter_allTrack_Chi2_ndof = 0x0;
  vector<double>          *reco_daughter_allTrack_momByRange_proton = 0x0;
  vector<double>          *reco_daughter_allTrack_momByRange_muon = 0x0;
  vector<double>          *reco_daughter_allTrack_momByRange_alt_proton = 0x0;
  vector<double>          *reco_daughter_allTrack_momByRange_alt_muon = 0x0;
  vector<double>          *reco_daughter_allTrack_Theta = 0x0;
  vector<double>          *reco_daughter_allTrack_Phi = 0x0;
  vector<double>          *reco_daughter_allTrack_startX = 0x0;
  vector<double>          *reco_daughter_allTrack_startY = 0x0;
  vector<double>          *reco_daughter_allTrack_startZ = 0x0;

 
  vector<int>             *reco_daughter_allShower_ID = 0x0; 
  vector<double>          *reco_daughter_allShower_startX = 0x0;
  vector<double>          *reco_daughter_allShower_startY = 0x0;
  vector<double>          *reco_daughter_allShower_startZ = 0x0;
  vector<double>          *reco_daughter_allShower_len = 0x0;
  vector<double>          *reco_daughter_allShower_dirX=0x0;
  vector<double>          *reco_daughter_allShower_dirY=0x0;
  vector<double>          *reco_daughter_allShower_dirZ=0x0;
  vector<double>          *reco_daughter_allShower_energy=0x0;

  vector<int>             *reco_daughter_pandora_type = 0x0;

  // Declare histograms
  
  
  //====================== Reco (Data and MC)======================//
  // Declare variables
  Double_t        beam_inst_X;
  Double_t        beam_inst_Y;
  Double_t        beam_inst_dirX;
  Double_t        beam_inst_dirY;
  Double_t        beam_inst_dirZ;
  Int_t           beam_inst_trigger;
  Int_t           beam_inst_nMomenta;
  Int_t           beam_inst_nTracks;
  vector<int> *   beam_inst_PDG_candidates = 0x0;
  Double_t        beam_inst_P; 
  
  //====================== Truth (MC only)======================//
  // Declare variables
  bool Signal = false;
  Int_t           true_beam_PDG;
  string          *true_beam_endProcess;
  Int_t           true_daughter_nPi0;
  Int_t           true_daughter_nPiPlus;
  Int_t           true_daughter_nProton;
  Int_t           true_daughter_nNeutron;
  Int_t           true_daughter_nPiMinus;
  Int_t           true_daughter_nNucleus;
  Double_t        true_beam_startX;
  Double_t        true_beam_startY;
  Double_t        true_beam_startZ;  
  Double_t        true_beam_startDirX;
  Double_t        true_beam_startDirY;
  Double_t        true_beam_startDirZ;
  Double_t        true_beam_endX;
  Double_t        true_beam_endY;
  Double_t        true_beam_endZ;
  Double_t        true_beam_endPx;
  Double_t        true_beam_endPy;
  Double_t        true_beam_endPz;
  Double_t        true_beam_startP;
  Double_t        true_beam_endP;

  vector<double>  *true_beam_traj_X;
  vector<double>  *true_beam_traj_Y;
  vector<double>  *true_beam_traj_Z;
  vector<double>  *true_beam_traj_X_SCE;
  vector<double>  *true_beam_traj_Y_SCE;
  vector<double>  *true_beam_traj_Z_SCE;
  
  vector<double>  *true_beam_traj_KE;
  vector<double>  *true_beam_slices_deltaE;
  vector<double>  *true_beam_incidentEnergies;

  vector<int>    *true_beam_daughter_ID = 0x0;  
  vector<int>    *true_beam_daughter_PDG = 0x0;
  vector<string> *true_beam_daughter_Process = 0x0;  
  vector<string> *true_beam_daughter_endProcess = 0x0;  
  vector<double> *true_beam_daughter_startP = 0x0;
  vector<double> *true_beam_daughter_startPx = 0x0;
  vector<double> *true_beam_daughter_startPy = 0x0;
  vector<double> *true_beam_daughter_startPz = 0x0;
  vector<double> *true_beam_daughter_startX = 0x0;
  vector<double> *true_beam_daughter_startY = 0x0;
  vector<double> *true_beam_daughter_startZ = 0x0;
  vector<double> *true_beam_daughter_endX = 0x0;
  vector<double> *true_beam_daughter_endY = 0x0;
  vector<double> *true_beam_daughter_endZ = 0x0;

  vector<int>     *true_beam_Pi0_decay_ID = 0x0;
  vector<int>     *true_beam_Pi0_decay_parID = 0x0;
  vector<int>     *true_beam_Pi0_decay_PDG = 0x0;
  vector<double>  *true_beam_Pi0_decay_startP = 0x0;
  vector<double>  *true_beam_Pi0_decay_startPx = 0x0;
  vector<double>  *true_beam_Pi0_decay_startPy = 0x0;
  vector<double>  *true_beam_Pi0_decay_startPz = 0x0;
  vector<double>  *true_beam_Pi0_decay_startX = 0x0;
  vector<double>  *true_beam_Pi0_decay_startY = 0x0;
  vector<double>  *true_beam_Pi0_decay_startZ = 0x0;
  vector<double>  *true_beam_Pi0_decay_len = 0x0;
  vector<int>     *true_beam_Pi0_decay_nHits = 0x0;
  
  vector<int> *true_beam_grand_daughter_ID= 0x0;
  vector<int> *true_beam_grand_daughter_PDG = 0x0;
  vector<string> *true_beam_grand_daughter_Process = 0x0;  
  vector<string> *true_beam_grand_daughter_endProcess = 0x0; 
 
  vector<int> *reco_daughter_PFP_true_byHits_PDG = 0x0;
  vector<int> *reco_daughter_PFP_true_byHits_origin = 0x0;
  vector<int> *reco_daughter_PFP_true_byHits_ID = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_startPx = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_startPy = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_startPz = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_purity = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_completeness = 0x0;
  vector<double>  *reco_daughter_PFP_true_byE_completeness = 0x0;

  // Declare histograms
  
  TH2D * hBeamTrackLength = 0x0;
  TH1D * hTrackLenRatio = 0x0;
  TH1D * hTrackLenRatio_Reco = 0x0;
  TH2D * hTrackLenRatio_Reco_STK = 0x0;
  TH2D * hStopMuonInstP_STK = 0x0;
  TH2D * hStopMuonKEff_STK = 0x0;
  TH1D * hTrueBeamMom = 0x0;
  TH1D * hTrueBeamMom_matched = 0x0;

  TH2D * hBeamTrackLength_STK = 0x0;
  TH2D * hBeamLongTrackLength_STK = 0x0;

  TH1D * hBeamLongTrackLength_Muon = 0x0;
  TH1D * hBeamLongTrackLength_Other = 0x0;

  // Get input tree
  TTree * GetInputTree(TFile * fin, const TString tname, const TString tag)
  {
    TTree * tree=(TTree*) fin->Get(tname);
    if(!tree){
      cout<<"no tree!"<<endl;
      gDirectory->ls();
      exit(1);
    }
    else cout << "Successuflly get " << tag << " tree with name: " << tname << endl;
    
    // Set Branch values
    //====================== Reco (MC and Data)======================// 
    tree->SetBranchAddress("event", &event); 
    tree->SetBranchAddress("reco_reconstructable_beam_event", &reco_reconstructable_beam_event);
    tree->SetBranchAddress("reco_beam_type", &reco_beam_type);  
    tree->SetBranchAddress("reco_beam_startX", &reco_beam_startX);
    tree->SetBranchAddress("reco_beam_startY", &reco_beam_startY);
    tree->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ);
    tree->SetBranchAddress("reco_beam_trackDirX", &reco_beam_trackDirX);
    tree->SetBranchAddress("reco_beam_trackDirY", &reco_beam_trackDirY);
    tree->SetBranchAddress("reco_beam_trackDirZ", &reco_beam_trackDirZ);
    tree->SetBranchAddress("reco_beam_endX", &reco_beam_endX);
    tree->SetBranchAddress("reco_beam_endY", &reco_beam_endY);
    tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
    tree->SetBranchAddress("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
    tree->SetBranchAddress("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
    tree->SetBranchAddress("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);
    tree->SetBranchAddress("reco_beam_interactingEnergy", &reco_beam_interactingEnergy);
    tree->SetBranchAddress("reco_beam_calibrated_dEdX_SCE", &reco_beam_calibrated_dEdX_SCE);
    tree->SetBranchAddress("reco_beam_TrkPitch_SCE", &reco_beam_TrkPitch_SCE);
    tree->SetBranchAddress("reco_beam_TrkPitch_NoSCE", &reco_beam_TrkPitch_NoSCE);
    tree->SetBranchAddress("reco_beam_incidentEnergies", &reco_beam_incidentEnergies);
    tree->SetBranchAddress("reco_beam_incidentEnergies", &reco_beam_new_incidentEnergies);


    

    tree->SetBranchAddress("reco_beam_Chi2_proton",&reco_beam_Chi2_proton);
    tree->SetBranchAddress("reco_beam_Chi2_ndof",&reco_beam_Chi2_ndof);

    tree->SetBranchAddress("reco_beam_calo_startX", &reco_beam_calo_startX);
    tree->SetBranchAddress("reco_beam_calo_startY", &reco_beam_calo_startY);
    tree->SetBranchAddress("reco_beam_calo_startZ", &reco_beam_calo_startZ);
    tree->SetBranchAddress("reco_beam_calo_endX", &reco_beam_calo_endX);
    tree->SetBranchAddress("reco_beam_calo_endY", &reco_beam_calo_endY);
    tree->SetBranchAddress("reco_beam_calo_endZ", &reco_beam_calo_endZ);
    tree->SetBranchAddress("reco_beam_calo_X", &reco_beam_calo_X);
    tree->SetBranchAddress("reco_beam_calo_Y", &reco_beam_calo_Y);
    tree->SetBranchAddress("reco_beam_calo_Z", &reco_beam_calo_Z);

    
    tree->SetBranchAddress("reco_beam_momByRange_muon", &reco_beam_momByRange_muon);
    tree->SetBranchAddress("reco_beam_momByRange_alt_muon", &reco_beam_momByRange_alt_muon);

    tree->SetBranchAddress("reco_beam_calo_startDirX", &reco_beam_calo_startDirX);
    tree->SetBranchAddress("reco_beam_calo_startDirY", &reco_beam_calo_startDirY);
    tree->SetBranchAddress("reco_beam_calo_startDirZ", &reco_beam_calo_startDirZ);
    tree->SetBranchAddress("reco_beam_calo_endDirX", &reco_beam_calo_endDirX);
    tree->SetBranchAddress("reco_beam_calo_endDirY", &reco_beam_calo_endDirY);
    tree->SetBranchAddress("reco_beam_calo_endDirZ", &reco_beam_calo_endDirZ);
    tree->SetBranchAddress("reco_beam_calo_wire", &reco_beam_calo_wire);
    tree->SetBranchAddress("reco_beam_vertex_michel_score", &reco_beam_vertex_michel_score);
    tree->SetBranchAddress("reco_beam_vertex_nHits", &reco_beam_vertex_nHits);
    tree->SetBranchAddress("reco_beam_PFP_nHits", &reco_beam_PFP_nHits);  
    tree->SetBranchAddress("reco_beam_PFP_trackScore_collection",&reco_beam_PFP_trackScore_collection);
    tree->SetBranchAddress("reco_beam_PFP_trackScore_collection_weight_by_charge",&reco_beam_PFP_trackScore_collection_weight_by_charge);
    tree->SetBranchAddress("reco_beam_PFP_emScore_collection",&reco_beam_PFP_emScore_collection);
    tree->SetBranchAddress("reco_beam_PFP_emScore_collection_weight_by_charge",&reco_beam_PFP_emScore_collection_weight_by_charge);
    tree->SetBranchAddress("reco_beam_vertex_michel_score_weight_by_charge",&reco_beam_vertex_michel_score_weight_by_charge);

    tree->SetBranchAddress("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG);
    tree->SetBranchAddress("reco_beam_true_byHits_origin", &reco_beam_true_byHits_origin);
    tree->SetBranchAddress("reco_beam_true_byHits_matched", &reco_beam_true_byHits_matched);
    tree->SetBranchAddress("reco_beam_true_byHits_endPx", &reco_beam_true_byHits_endPx);
    tree->SetBranchAddress("reco_beam_true_byHits_endPy", &reco_beam_true_byHits_endPy);
    tree->SetBranchAddress("reco_beam_true_byHits_endPz", &reco_beam_true_byHits_endPz);
    tree->SetBranchAddress("reco_beam_true_byHits_endE", &reco_beam_true_byHits_endE);
    tree->SetBranchAddress("reco_beam_true_byHits_endP", &reco_beam_true_byHits_endP);
    tree->SetBranchAddress("reco_beam_true_byHits_startPx", &reco_beam_true_byHits_startPx);
    tree->SetBranchAddress("reco_beam_true_byHits_startPy", &reco_beam_true_byHits_startPy);
    tree->SetBranchAddress("reco_beam_true_byHits_startPz", &reco_beam_true_byHits_startPz);
    tree->SetBranchAddress("reco_beam_true_byHits_startE", &reco_beam_true_byHits_startE);
    tree->SetBranchAddress("reco_beam_true_byHits_startP", &reco_beam_true_byHits_startP);

    tree->SetBranchAddress("reco_beam_true_byHits_process",&reco_beam_true_byHits_process);
    tree->SetBranchAddress("reco_beam_true_byHits_endProcess",&reco_beam_true_byHits_endProcess);
    tree->SetBranchAddress("reco_beam_true_byE_matched",&reco_beam_true_byE_matched);
    tree->SetBranchAddress("reco_beam_true_byE_PDG",&reco_beam_true_byE_PDG);

    tree->SetBranchAddress("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
    tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
    tree->SetBranchAddress("reco_daughter_PFP_nHits_collection", &reco_daughter_PFP_nHits_collection);
    tree->SetBranchAddress("reco_daughter_PFP_trackScore", &reco_daughter_PFP_trackScore);
    tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
    tree->SetBranchAddress("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore);
    tree->SetBranchAddress("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
    tree->SetBranchAddress("reco_daughter_PFP_michelScore", &reco_daughter_PFP_michelScore);
    tree->SetBranchAddress("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);
    tree->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
    tree->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);
    tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
    tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);
    tree->SetBranchAddress("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
    tree->SetBranchAddress("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);
    tree->SetBranchAddress("reco_daughter_allTrack_momByRange_alt_proton", &reco_daughter_allTrack_momByRange_alt_proton);
    tree->SetBranchAddress("reco_daughter_allTrack_momByRange_alt_muon", &reco_daughter_allTrack_momByRange_alt_muon);
    tree->SetBranchAddress("reco_daughter_allTrack_Theta", &reco_daughter_allTrack_Theta);
    tree->SetBranchAddress("reco_daughter_allTrack_Phi", &reco_daughter_allTrack_Phi);
    tree->SetBranchAddress("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX);
    tree->SetBranchAddress("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY);
    tree->SetBranchAddress("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ);

    tree->SetBranchAddress("reco_daughter_allShower_ID", &reco_daughter_allShower_ID);
    tree->SetBranchAddress("reco_daughter_allShower_startX", &reco_daughter_allShower_startX);
    tree->SetBranchAddress("reco_daughter_allShower_startY", &reco_daughter_allShower_startY);
    tree->SetBranchAddress("reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ);
    tree->SetBranchAddress("reco_daughter_allShower_len", &reco_daughter_allShower_len);
    tree->SetBranchAddress("reco_daughter_allShower_dirX", &reco_daughter_allShower_dirX);
    tree->SetBranchAddress("reco_daughter_allShower_dirY", &reco_daughter_allShower_dirY);
    tree->SetBranchAddress("reco_daughter_allShower_dirZ", &reco_daughter_allShower_dirZ);
    tree->SetBranchAddress("reco_daughter_allShower_energy", &reco_daughter_allShower_energy);

    tree->SetBranchAddress("reco_daughter_pandora_type", &reco_daughter_pandora_type);

    //====================== Reco (Data and MC)======================//
    tree->SetBranchAddress("beam_inst_PDG_candidates", &beam_inst_PDG_candidates);
    tree->SetBranchAddress("beam_inst_X", &beam_inst_X);
    tree->SetBranchAddress("beam_inst_Y", &beam_inst_Y);
    tree->SetBranchAddress("beam_inst_dirX", &beam_inst_dirX);
    tree->SetBranchAddress("beam_inst_dirY", &beam_inst_dirY);
    tree->SetBranchAddress("beam_inst_dirZ", &beam_inst_dirZ);
    tree->SetBranchAddress("beam_inst_trigger", &beam_inst_trigger); 
    tree->SetBranchAddress("beam_inst_nMomenta", &beam_inst_nMomenta);
    tree->SetBranchAddress("beam_inst_nTracks", &beam_inst_nTracks);
    tree->SetBranchAddress("beam_inst_P", &beam_inst_P);

    

    //====================== Truth (MC only)======================//
    tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
    tree->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess);
    tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
    tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);
    tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);
    tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
    tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
    tree->SetBranchAddress("true_daughter_nNucleus", &true_daughter_nNucleus);
    tree->SetBranchAddress("true_beam_startX", &true_beam_startX);
    tree->SetBranchAddress("true_beam_startY", &true_beam_startY);
    tree->SetBranchAddress("true_beam_startZ", &true_beam_startZ);
    tree->SetBranchAddress("true_beam_startDirX", &true_beam_startDirX);
    tree->SetBranchAddress("true_beam_startDirY", &true_beam_startDirY);
    tree->SetBranchAddress("true_beam_startDirZ", &true_beam_startDirZ);  
    tree->SetBranchAddress("true_beam_endX", &true_beam_endX);
    tree->SetBranchAddress("true_beam_endY", &true_beam_endY);
    tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ);
    tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
    tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
    tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);
    tree->SetBranchAddress("true_beam_startP", &true_beam_startP);
    tree->SetBranchAddress("true_beam_endP", &true_beam_endP);

    tree->SetBranchAddress("true_beam_traj_X", &true_beam_traj_X);
    tree->SetBranchAddress("true_beam_traj_Y", &true_beam_traj_Y);
    tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
    tree->SetBranchAddress("true_beam_traj_X_SCE", &true_beam_traj_X_SCE);
    tree->SetBranchAddress("true_beam_traj_Y_SCE", &true_beam_traj_Y_SCE);
    tree->SetBranchAddress("true_beam_traj_Z_SCE", &true_beam_traj_Z_SCE);
    tree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE);
    tree->SetBranchAddress("true_beam_slices_deltaE", &true_beam_slices_deltaE);
    tree->SetBranchAddress("true_beam_incidentEnergies", &true_beam_incidentEnergies);

    tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID); 
    tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
    tree->SetBranchAddress("true_beam_daughter_Process", &true_beam_daughter_Process); 
    tree->SetBranchAddress("true_beam_daughter_endProcess", &true_beam_daughter_endProcess);
    tree->SetBranchAddress("true_beam_daughter_startP", &true_beam_daughter_startP);
    tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
    tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
    tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

    tree->SetBranchAddress("true_beam_daughter_startX", &true_beam_daughter_startX);
    tree->SetBranchAddress("true_beam_daughter_startY", &true_beam_daughter_startY);
    tree->SetBranchAddress("true_beam_daughter_startZ", &true_beam_daughter_startZ);
    tree->SetBranchAddress("true_beam_daughter_endX", &true_beam_daughter_endX);
    tree->SetBranchAddress("true_beam_daughter_endY", &true_beam_daughter_endY);
    tree->SetBranchAddress("true_beam_daughter_endZ", &true_beam_daughter_endZ);
   
    tree->SetBranchAddress("true_beam_Pi0_decay_ID", &true_beam_Pi0_decay_ID);
    tree->SetBranchAddress("true_beam_Pi0_decay_parID", &true_beam_Pi0_decay_parID);
    tree->SetBranchAddress("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG);
    tree->SetBranchAddress("true_beam_Pi0_decay_startP", &true_beam_Pi0_decay_startP);
    tree->SetBranchAddress("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
    tree->SetBranchAddress("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
    tree->SetBranchAddress("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
    tree->SetBranchAddress("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
    tree->SetBranchAddress("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
    tree->SetBranchAddress("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);
    tree->SetBranchAddress("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);
    tree->SetBranchAddress("true_beam_Pi0_decay_nHits", &true_beam_Pi0_decay_nHits);
    tree->SetBranchAddress("true_beam_grand_daughter_ID", &true_beam_grand_daughter_ID);
    tree->SetBranchAddress("true_beam_grand_daughter_PDG", &true_beam_grand_daughter_PDG);
    tree->SetBranchAddress("true_beam_grand_daughter_Process", &true_beam_grand_daughter_Process); 
    tree->SetBranchAddress("true_beam_grand_daughter_endProcess", &true_beam_grand_daughter_endProcess);
 
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_origin", &reco_daughter_PFP_true_byHits_origin);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_startPx", &reco_daughter_PFP_true_byHits_startPx);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_startPy", &reco_daughter_PFP_true_byHits_startPy);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_startPz", &reco_daughter_PFP_true_byHits_startPz);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_purity",&reco_daughter_PFP_true_byHits_purity);
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_completeness",&reco_daughter_PFP_true_byHits_completeness);
    tree->SetBranchAddress("reco_daughter_PFP_true_byE_completeness",&reco_daughter_PFP_true_byE_completeness);
    return tree;
  } // End of GetInputTree
  
  
  TTree * GetOutputTree(TList * lout, const TString tag)
  {
    TTree * tout = new TTree("tree", tag); lout->Add(tout);

    // Definitely need this to avoid memory-resident Tree
    tout->SetDirectory(0);
    tout->Branch("trueBeamMom",&trueBeamMom);
    tout->Branch("recoFFKE",&recoFFKE);
    tout->Branch("instP",&instP);
    
    return tout;
  }

  // Initialise histograms
  void IniHist(TList * lout, const bool kMC)
  {  
    const int nbeamType = 8;
    const double beamTypemin = -0.5;
    const double beamTypemax = 7.5;

    const int nx = 51;
    const int xmin = -10;
    const int xmax = 500;

    //const int nx1 = 50;
    //const int xmin1 = 0;
    //const int xmax1 = 500;

    // Define the bin edges
    /*Double_t xbins[nx+1];
    for (Int_t i = 0; i < nx; i++) {
      xbins[i] = xmin + i*10;
    }
    xbins[nx] = xmax;*/

    hBeamTrackLength = new TH2D("a001hBeamTrackLength_STK",";Beam Track Length (cm);Candidates", nx, xmin, xmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamTrackLength);


    hTrackLenRatio_Reco = new TH1D("a002hTrackLenRatio_Reco",";Track Length / RangebyKE; Candidates", 120, 0, 1.2); 
    lout->Add(hTrackLenRatio_Reco);

    hTrackLenRatio_Reco_STK = new TH2D("a003hBeamTrackLenRatio_Reco_STK",";Track Length / RangebyKE; Candidates", 100, -0.3, 1.2, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hTrackLenRatio_Reco_STK);

    hStopMuonInstP_STK = new TH2D("a003hBeamStopMuonInstP_STK",";Inst. P (GeV/c); Candidates", 20, 0.6, 1.3, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hStopMuonInstP_STK);

    hStopMuonKEff_STK = new TH2D("a004hBeamStopMuonKEff_STK",";#mu^{+} Front-Face Energy (MeV/c); Candidates", 20, 600, 1300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hStopMuonKEff_STK);

    hBeamTrackLength_STK = new TH2D("a005hBeamTrackLength_STK",";Track Length (cm); Candidates", nx, xmin, xmax, nbeamType, beamTypemin, beamTypemax);
    lout->Add(hBeamTrackLength_STK);

    const double LongTrkBin[] = {-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 500};

    //hBeamLongTrackLength_STK = new TH2D("a006hBeamLongTrackLength_STK",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin, nbeamType, beamTypemin, beamTypemax);
    hBeamLongTrackLength_STK = new TH2D("a006hBeamLongTrackLength_STK",";Long Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin, nbeamType, beamTypemin, beamTypemax);
    lout->Add(hBeamLongTrackLength_STK);


    hBeamLongTrackLength_Muon = new TH1D("a007hBeamLongTrackLength_Muon",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    lout->Add(hBeamLongTrackLength_Muon);

    hBeamLongTrackLength_Other = new TH1D("a008hBeamLongTrackLength_Other",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    lout->Add(hBeamLongTrackLength_Other);

    if(kMC){
      hTrackLenRatio = new TH1D("b001hTrackLenByRange",";Track Length / RangebyKE; Candidates", 100, -0.3, 1.2); 
      lout->Add(hTrackLenRatio);
      hTrueBeamMom = new TH1D("b002hTrueBeamMom",";True Beam Start Mom. (GeV/c); Candidates", 30, 0.7, 1.3);
      lout->Add(hTrueBeamMom);
      hTrueBeamMom_matched = new TH1D("b003hTrueBeamMom_matched",";True Beam Start Mom. (GeV/c); Candidates", 30, 0.7, 1.3);
      lout->Add(hTrueBeamMom_matched);
    }
  }// End of IniHist

} // End of namespace
#endif
