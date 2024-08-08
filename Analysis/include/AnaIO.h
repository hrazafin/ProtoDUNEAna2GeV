#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "PlotUtils.h"
#include "AnaFunctions.h"
// Kinematic fitting codes
#include "../src/UserKF.cxx"
#include "../src/BetheBloch.cxx"

namespace AnaIO
{
  //====================== Output Tree Variables ==================//
  
  double LeadingShowerEnergy;
  double SubLeadingShowerEnergy;
  double LeadingShowerEnergyRaw;
  double SubLeadingShowerEnergyRaw;
  double OpeningAngle;
  double LeadingShowerEnergyTruth;
  double SubLeadingShowerEnergyTruth;
  double OpeningAngleTruth;

  TVector3 LeadingShowerEnergyUnitDir;
  TVector3 SubLeadingShowerEnergyUnitDir;
  TVector3 LeadingShowerEnergyUnitDirTruth;
  TVector3 SubLeadingShowerEnergyUnitDirTruth;

  double dalphat;
  double dphit;
  double dpt;
  double pn;
  double iniPimomentum;
  double iniPitheta;
  double finPimomentum;
  double finPitheta;
  double finProtonmomentum;
  double finProtontheta;
  double fin2Pmom;

  //====================== Reco (MC and Data)======================//
  //====================== Reco Beam ======================//
  // Declare variables
  Int_t           run = -999;
  Int_t           subrun= -999;
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
  vector<double>  *reco_beam_dEdX_SCE = 0x0;
  vector<double>  *reco_beam_dEdX_NoSCE = 0x0;
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

  // Declare histograms Herilala
  TH2D * hxstartystart = 0x0;
  TH2D * hxstartystartBeam = 0x0;
  TH1D * hBeamDeltaXYSigma = 0x0;
  TH1D * hrstartX = 0x0;
  TH1D * hrstartY = 0x0;
  TH1D * hrstartZ = 0x0;

  TH1D * hreco_beam_trackDirX = 0x0;
  TH1D * hreco_beam_trackDirY = 0x0;
  TH1D * hreco_beam_trackDirZ = 0x0;

  TH1D * hrendX = 0x0;
  TH1D * hrendY = 0x0;
  TH1D * hrendZ = 0x0;

  // Class A - beam cut related
  TH1I * hCutBeamPDGPass = 0x0;
  TH1I * hCutPandoraSlicePass = 0x0;
  TH1I * hCutCaloSizePass = 0x0;
  TH1I * hCutBeamQualityPass = 0x0;
  TH1I * hCutAPA3EndZPass = 0x0;
  TH1I * hCutMichelScorePass = 0x0;
  TH1I * hCutProtonChi2Pass = 0x0;
  TH1I * hCutBeamScraperPass = 0x0;

  TH1D * hRecBeamStartX = 0x0;
  TH1D * hRecBeamStartY = 0x0;
  TH1D * hRecBeamStartZ = 0x0;
  TH1D * hRecBeamThetaX = 0x0;
  TH1D * hRecBeamThetaY = 0x0;
  TH1D * hRecBeamThetaZ = 0x0;
  TH1D * hCutBeamDeltaXYSigma = 0x0;
  TH1D * hCutBeamDeltaZSigma = 0x0;
  TH1D * hCutBeamCosTheta = 0x0;
  TH2D * hCutBeamMichelScoreNhits = 0x0;
  TH2D * hCutBeamProtonChi2DOF = 0x0;
  TH1D * hCutBeamProtonChi2DOFType = 0x0;
  TH2D * hCutBeamQualityXY = 0x0;
  TH2D * hCutBeamQualityX = 0x0;
  TH2D * hCutBeamQualityY = 0x0;
  TH2D * hCutBeamQualityZ = 0x0;
  TH2D * hCutBeamQualityTheta = 0x0;
  TH2D * hCutAPA3EndZ = 0x0;

  TH2D * hBeamtrackScore = 0x0;
  TH2D * hBeamtrackScore_wbc = 0x0;
  TH2D * hBeamemScore = 0x0;
  TH2D * hBeamemScore_wbc = 0x0;
  TH2D * hBeamDaughterNumbers = 0x0;
  TH2D * hBeamStartX = 0x0;
  TH2D * hBeamStartCaloX = 0x0;
  TH2D * hBeamStartY = 0x0;
  TH2D * hBeamStartCaloY = 0x0;
  TH2D * hBeamStartZ = 0x0;
  TH2D * hBeamStartCaloZ = 0x0;
  TH2D * hBeamEndX = 0x0;
  TH2D * hBeamEndCaloX = 0x0;
  TH2D * hBeamEndY = 0x0;
  TH2D * hBeamEndCaloY = 0x0;
  TH2D * hBeamEndZ = 0x0;
  TH2D * hBeamEndCaloZ = 0x0;

  TH2D * hCutBeamQualityInstXY = 0x0;
  TH2D * hCutBeamQualityInstX = 0x0;
  TH2D * hCutBeamQualityInstY = 0x0;
  TH1D * hRecBeamInstX = 0x0;
  TH1D * hRecBeamInstY = 0x0;
  
  TH2D * hBeamDirDiff = 0x0;
  TH1D * hRecBeamTrackPitch_SCE = 0x0;
  TH1D * hRecBeamTrackPitch_NoSCE = 0x0;
  TH1D * hRecBeamdEdx_SCE = 0x0;
  TH1D * hRecBeamdEdx_NoSCE = 0x0;  
  TH1D * hRecBeamDeltaE_SCE = 0x0;
  TH1D * hRecBeamDeltaE_NoSCE = 0x0;
  TH1D * hRecIntEdiff = 0x0;
  // Class B - reconstructed beam/FS particles
  TH2D * hRecBeamPhi = 0x0;
  TH2D * hRecBeamTheta = 0x0;
  TH2D * hRecBeamMomentum = 0x0;
  TH2D * hRecBeamPhi_XsEvt = 0x0;
  TH2D * hRecBeamTheta_XsEvt = 0x0;
  TH2D * hRecBeamMomentum_XsEvt = 0x0;
  TH2D * hRecProtonTheta = 0x0;
  TH2D * hRecProtonMomentum = 0x0;
  TH2D * hRecPiPlusTheta = 0x0;
  TH2D * hRecPiPlusMomentum = 0x0;
  TH2D * hRecShowerEnergy = 0x0;
  TH2D * hRecShowerEnergyRaw = 0x0;
  TH2D * hRecShowerTheta = 0x0; 
  TH2D * hRecShowerLength = 0x0;
  TH2D * hRecPi0Nshower = 0x0;
  TH2D * hRecLeadingShowerEnergy = 0x0;
  TH2D * hRecSubLeadingShowerEnergy = 0x0;
  TH2D * hRecLeadingShowerEnergyRaw = 0x0;
  TH2D * hRecSubLeadingShowerEnergyRaw = 0x0;
  TH2D * hRecShowerOpenAngle = 0x0;
  TH2D * hRecPi0ShowerSep_OVERLAY = 0x0;
  TH2D * hRecPi0Mass = 0x0;
  TH2D * hRecPi0MassRaw = 0x0;
  TH2D * hRecPi0Momentum = 0x0;
  TH2D * hRecPi0MomentumRaw = 0x0;

  TH2D * hRecPi0Mass_OVERLAY = 0x0;
  TH2D * hRecPi0Mass_OVERLAY_EVT = 0x0;

  TH2D * testhRecPi0Mass_OVERLAY = 0x0;
  TH2D * test1hRecPi0Mass_OVERLAY = 0x0;

  TH2D * testhRecPi0Mass_OVERLAY_EVT = 0x0;
  TH2D * test1hRecPi0Mass_OVERLAY_EVT = 0x0;
  
  TH2D * hRecPi0MassRaw_OVERLAY = 0x0;
  TH2D * hRecPi0Momentum_OVERLAY = 0x0;
  TH2D * hRecPi0MomentumRaw_OVERLAY = 0x0;
  TH2D * hRecPi0OA_OVERLAY = 0x0;
  TH2D * testhRecPi0OA_OVERLAY = 0x0;
  TH2D * test1hRecPi0OA_OVERLAY = 0x0;

  TH2D * hRecPi0OARaw_OVERLAY = 0x0;
  TH2D * hRecPi0Energy_OVERLAY = 0x0;

  TH2D * hRecPi0Theta_OVERLAY = 0x0;
  TH2D * hRecPi0Phi_OVERLAY = 0x0;

  TH2D * hRecPi0ThetaRaw_OVERLAY = 0x0;
  TH2D * hRecPi0PhiRaw_OVERLAY = 0x0;

  TH2D * hRecPi0Energy_OVERLAY_After = 0x0;
  TH2D * hRecPi0Theta_OVERLAY_After = 0x0;
  TH2D * hRecPi0Phi_OVERLAY_After = 0x0;


  TH2D * hRecPi0Energy_OVERLAY_After_Kinetic = 0x0;
  TH2D * hRecPi0Energy_OVERLAY_After_EVT = 0x0;
  TH2D * hRecPi0Energy_OVERLAY_After_EVT_High = 0x0;
  TH2D * hRecPi0Energy_OVERLAY_After_EVT_Medium = 0x0;
  TH2D * hRecPi0Energy_OVERLAY_After_EVT_Low = 0x0;

  TH2D * hRecPi0TotalE_Default = 0x0;
  TH2D * hRecPi0TotalE_Fitted = 0x0;
  TH2D * hRecPi0TotalEafterSel_Default = 0x0;
  TH2D * hRecPi0TotalEafterSel_Fitted = 0x0;

  TH2D * hRecPi0KineticEnergyVSPi0Mass = 0x0;
  


  TH2D * hRecPi0Energy_OVERLAY_AfterTOP_EVT = 0x0;

  TH1D * hPi0MassLowE1 = 0x0;
  TH1D * hPi0MassHighE1 = 0x0;

  // Common Reco TKI
  TH2D * hRecdalphat = 0x0;
  TH2D * hRecdphit = 0x0;
  TH2D * hRecdpt = 0x0;
  TH2D * hRecpn = 0x0;

  TH2D * hRecdalphat_truth = 0x0;
  TH2D * hRecdphit_truth = 0x0;
  TH2D * hRecdpt_truth = 0x0;
  TH2D * hRecpn_truth = 0x0;

  TH2D * hdalphat_REG = 0x0;
  TH2D * hdphit_REG = 0x0;
  TH2D * hdpt_REG = 0x0;
  TH2D * hpn_REG = 0x0;

  TH2D * hdalphat_RES = 0x0;
  TH2D * hdphit_RES = 0x0;
  TH2D * hdpt_RES = 0x0;
  TH2D * hpn_RES = 0x0;

  TH1D * hLDShower_PreFit = 0x0;
  TH1D * hLDShower_PostFit = 0x0;

  TH1D * hSLShower_PreFit = 0x0;
  TH1D * hSLShower_PostFit = 0x0;

  TH1D * hOAShower_PreFit = 0x0;
  TH1D * hOAShower_PostFit = 0x0;

  TH1D * hPi0EnergyE1E2_PreFit = 0x0;
  TH1D * hPi0EnergyE1E2_PostFit = 0x0;

  TH1D * hPi0EnergyE1OA_PreFit = 0x0;
  TH1D * hPi0EnergyE1OA_PostFit = 0x0;

  TH1D * hPi0EnergyAsym_PreFit = 0x0;
  TH1D * hPi0EnergyAsym_PostFit = 0x0;

  TH2D * hShowerE1Compare_REG = 0x0;
  TH2D * hShowerE2Compare_REG = 0x0;
  TH2D * hShowerOACompare_REG = 0x0;
  TH2D * hPi0EnergyComparem1_REG = 0x0;
  TH2D * hPi0EnergyComparem2_REG = 0x0;


  // FS particle cut
  TH1I * hCutDaughterPandoraShowerPass = 0x0;
  TH1I * hCutDaughterShowerScorePass = 0x0;
  TH1I * hCutDaughterShowernHitsPass = 0x0;
  TH1I * hCutDaughterShowerNonEmptyEPass = 0x0;
  TH1I * hCutDaughterShowerStartZPass = 0x0;
  TH1I * hCutDaughterShowerDistPass = 0x0;
  TH1I * hCutDaughterShowerIPPass = 0x0;

  TH1I * hCutDaughterTrackScorePass = 0x0;
  TH1I * hCutDaughterTracknHitsPass = 0x0;
  TH1I * hCutDaughterProtonSubPIDPass = 0x0;

  TH1I * hCutDaughterPionTrackScorePass = 0x0;
  TH1I * hCutDaughterPionTracknHitsPass = 0x0;
  TH1I * hCutDaughterPionSubPIDPass = 0x0;

  TH1I * hCutDaughterPi0NOCutsPass = 0x0;
  TH1I * hCutDaughterPi0MassPass = 0x0;
  TH1I * hCutDaughterPi0OAPass = 0x0;

  TH2D * hCutTracknHits = 0x0;
  TH2D * hCutTrackScore = 0x0;
  TH2D * hCutlastTME = 0x0;
  TH2D * hCutChi2NDF = 0x0;
  TH2D * hCutemScore = 0x0;
  // Check SCE affect shower score
  TH2D * hCutemScore_R1 = 0x0;
  TH2D * hCutemScore_R2 = 0x0;
  TH2D * hCutemScore_R3 = 0x0;
  TH2D * hCutemScore_R4 = 0x0;

  TH2D * hCutemScoreVSStartX = 0x0;
  TH2D * hCutemScoreVSStartY = 0x0;
  TH2D * hCutemScoreVSStartZ = 0x0;
  TH2D * hCutemScoreVSEnergy = 0x0;

  TH2D * hCutStartXVSStartY = 0x0;
  TH2D * hCutStartZVSStartX = 0x0;

  TH2D * hAllCutemScoreVSStartX = 0x0;
  TH2D * hAllCutemScoreVSStartY = 0x0;
  TH2D * hAllCutemScoreVSStartZ = 0x0;

  TH2D * hAllCutStartXVSStartY = 0x0;
  TH2D * hAllCutStartZVSStartX = 0x0;

  TH2D * hAllTrackCutStartXVSStartY = 0x0;
  TH2D * hAllTrackCutStartZVSStartX = 0x0;

  TH2D * hBeamEndZ_CutPDG = 0x0;
  TH2D * hBeamEndZ_CutPandora = 0x0;
  TH2D * hBeamEndZ_CutCaloSize = 0x0;
  TH2D * hBeamEndZ_CutBeamQuality = 0x0;
  TH2D * hBeamEndZ_CutAPA3 = 0x0;
  TH2D * hBeamEndZ_CutMichelScore = 0x0;
  TH2D * hBeamEndZ_CutChi2DOF = 0x0;
  TH2D * hBeamEndZ_CutBeamScraper = 0x0;

  TH2D * hBeamEndZ_Channels = 0x0;
  TH2D * hBeamEndZ_Channels_afterEvtTop = 0x0;
  TH1D * hBeamEndZ_ChannelsTrueSignal = 0x0;
  TH2D * hBeamEndZ_Channels_BeamPDG = 0x0;
  TH2D * hBeamEndZ_Channels_BeamPandora = 0x0;
  TH2D * hBeamEndZ_Channels_BeamCalo = 0x0;
  TH2D * hBeamEndZ_Channels_BeamQuality = 0x0;
  TH2D * hBeamEndZ_Channels_BeamAPA3 = 0x0;
  TH2D * hBeamEndZ_Channels_BeamMichel = 0x0;
  TH2D * hBeamEndZ_Channels_BeamChi2 = 0x0;
  TH2D * hBeamEndZ_Channels_BeamScraper = 0x0;
  TH2D * hBeamEndZ_Channels_TwoShowers = 0x0;
  TH2D * hBeamEndZ_Channels_NoPiPlus = 0x0;
  TH2D * hBeamEndZ_Channels_NoMichel = 0x0;
  TH2D * hBeamEndZ_Channels_OnePi0 = 0x0;


  TH2D * hBeamEndZ_XsEvt_BeamPDG = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamPandora = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamCalo = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamQuality = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamAPA3 = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamMichel = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamChi2 = 0x0;
  TH2D * hBeamEndZ_XsEvt_BeamScraper = 0x0;
  TH2D * hBeamEndZ_XsEvt_TwoShowers = 0x0;
  TH2D * hBeamEndZ_XsEvt_NoPiPlus = 0x0;
  TH2D * hBeamEndZ_XsEvt_NoMichel = 0x0;
  TH2D * hBeamEndZ_XsEvt_OnePi0 = 0x0;
  

  TH1D * hPiMom_InelasticChannels = 0x0;
  TH2D * hBeamInstXY = 0x0;
  TH2D * hBeamInstXY_misIDproton = 0x0;
  TH2D * hBeamInstXY_misIDpion = 0x0;

  TH1D * hBeamEndZ_TrueAvailable = 0x0;
  TH1D * hBeamEndZ_ChargeExcChannel = 0x0;

  TH2D * hShowerEnergy_NoCut = 0x0;
  TH2D * hShowerEnergy_CutEMscore = 0x0;
  TH2D * hShowerEnergy_CutnHits = 0x0;
  TH2D * hShowerEnergy_CutStartZ = 0x0;
  TH2D * hShowerEnergy_CutDistance = 0x0;
  TH2D * hShowerEnergy_CutIP = 0x0;

  TH2D * hPi0Energy_NoCut = 0x0;
  TH2D * hPi0Energy_CutMass = 0x0;
  TH2D * hPi0Energy_CutOA = 0x0;

  TH2D * hProtonMom_NoCut = 0x0;
  TH2D * hProtonMom_CutTrackScore = 0x0;
  TH2D * hProtonMom_CutnHits = 0x0;
  TH2D * hProtonMom_CutSubPID = 0x0;

  TH2D * hPiPlusMom_NoCut = 0x0;
  TH2D * hPiPlusMom_CutTrackScore = 0x0;
  TH2D * hPiPlusMom_CutnHits = 0x0;
  TH2D * hPiPlusMom_CutSubPID = 0x0;

  TH1D * hRecoIncidentHistData = 0x0;
  TH1D * hRecoInteractingHistData = 0x0;
  TH1D * hRecoBeamInitialHistData = 0x0;
  TH1D * hRecoBeamInteractingHistData = 0x0;

  TH1D * hUnFoldedInteractingHistData = 0x0;
  TH1D * hUnFoldedIncidentHistData = 0x0;
  TH1D * hUnFoldedBeamInitialHistData = 0x0;
  TH1D * hUnFoldedBeamInteractingHistData = 0x0;
  TH1D * hUnfoldedBeamIncidentHistData = 0x0;



  TH1D * hRecoPi0KEHistData = 0x0;
  TH1D * hUnFoldedPi0KEHistData = 0x0;

  TH1D * hRecoPi0CosThetaHistData = 0x0;
  TH1D * hUnFoldedPi0CosThetaHistData = 0x0;

  TH1D * hRecoPi0ThetaHistData = 0x0;
  TH1D * hUnFoldedPi0ThetaHistData = 0x0;

 



  TH2D * hCutemScore_AfterCut = 0x0;
  
  TH2D * hCutmichelScore = 0x0;
  TH2D * hCutShowerDist = 0x0;
  TH2D * hCutShowerIP = 0x0;
  TH2D * hCutShowerStartZ = 0x0;

  TH2D * hCutnproton = 0x0;
  TH2D * hCutnpiplus = 0x0;
  TH2D * hCutnshower = 0x0;
  TH2D * hCutnmichel = 0x0;
  TH2D * hCutnpi0 = 0x0;
  TH2D * hCutKFPass = 0x0;

  TH2D * hShowerStartX = 0x0;
  TH2D * hShowerStartY = 0x0;
  TH2D * hShowerStartZ = 0x0;
  TH2D * hShowernHitsColl = 0x0;

  TH2D * hNhitsCollection = 0x0;

  TH2D * htrackScoreCollection = 0x0;
  TH2D * hemScoreCollection = 0x0;
  TH2D * hmichelScoreCollection = 0x0;

  TH2D * hShowerAngleOffset = 0x0;

  TH2D * hShowerMichelScore = 0x0;
  
  //====================== Reco (Data and MC)======================//
  // Declare variables
  Double_t        beam_inst_X;
  Double_t        beam_inst_Y;
  Double_t        beam_inst_Z;
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
 
  //vector<vector<double>>  *g4rw_full_grid_piplus_weights = 0x0;
  vector<vector<double>>  *g4rw_full_grid_piplus_coeffs = 0x0;

  // Declare histograms
  TH1D * hTruthBeamMomentum = 0x0;
  TH1I * hTruthBeamType = 0x0;
  TH1I * hTruthSignal = 0x0;
  TH1I * hTruthFSParticleNumber = 0x0;  
  TH1I * hTruthFSParticleType = 0x0;
  TH1I * hTruthFSPi0Number = 0x0;
  TH1I * hTruthFSMultiPi0 = 0x0;
  TH1D * hTruthProtonP = 0x0;
  TH1D * hTruthLeadingProtonP = 0x0;
  TH1D * hTruthSubLeadingProtonP = 0x0;
  TH1D * hTruthLeadingPiZeroP = 0x0;
  TH1D * hTruthLeadingPiZeroE = 0x0;
  TH1D * hTruthSubLeadingPiZeroP = 0x0;
  TH1D * hTruthGammaMaxE = 0x0;
  TH1D * hTruthPi0DecayParticleNumber = 0x0;
  TH1D * hTruthLeadingPi0GammaP = 0x0;
  TH1D * hTruthLeadingPi0GammaPBin1 = 0x0;
  TH1D * hTruthSubLeadingPi0GammaP = 0x0;
  TH1D * hTruthPi0OA = 0x0;
  TH1D * hTruthLeadingPi0GammaOA = 0x0;
  TH1D * hTruthSubLeadingPi0GammaOA = 0x0;

  TH1D * hTruthLeadingPiZeroGammaDist = 0x0;
  TH1D * hTruthSubLeadingPiZeroGammaDist = 0x0;

  TH2D * hTruthPiZeroGammaE1E2 = 0x0;
  TH2D * hTruthPiZeroGammaE1OA = 0x0;

  TH1D * hTruthPiPlusP = 0x0;
  TH1D * hTruthLeadingPiPlusP = 0x0;
  TH1D * hTruthSubLeadingPiPlusP = 0x0;

  TH1D * hTruthPiMinusP = 0x0;

  TH1D * hTruthPDGCode_PionMatched = 0x0;

  TH1D * hTruthMatchedBeamSample = 0x0;
  TH1D * hTruthMatchedProtonSample = 0x0;
  TH1D * hTruthMatchedPiPlusSample = 0x0;
  TH1D * hTruthMatchedShowerSample = 0x0;
  TH1D * hTruthMatchedPi0Sample = 0x0;
  TH1D * hTruthMatchedChannelCEXSample = 0x0;
  TH1D * hTruthMatchedXsEvtCEXSample = 0x0;

  TH2D * hTruthPi0GammaEnergy = 0x0;
  TH1D * hTruthRarePi0GammaP = 0x0;
  TH1D * hTruthRarePi0ElectronP = 0x0;
  TH1D * hTruthRarePi0PositronP = 0x0;

  TH1I * hTruthSignalFSParticleNumber = 0x0;  
  TH1I * hTruthSignalFSParticleType = 0x0;

  TH1D * hTruthIncidentHist = 0x0;
  TH1D * hTruthInteractingHist = 0x0;
  TH1D * hNewTruthInteractingHist = 0x0;
  TH1D * hNewTruthInteractingHistTest = 0x0;

  TH1D * hTruthSingleIncidentHist = 0x0;
  TH1D * hTruthCalcIncidentHist = 0x0;
  TH1D * hTruthSingleInteractingHist = 0x0;
  TH1D * hTruthBeamInteractingHist = 0x0;
  TH1D * hNewTruthBeamInteractingHist = 0x0;

  TH1D * hTruthBeamInteractingHist_50MeVbin = 0x0;

  TH1D * hTruthBeamCalcIncidentHist = 0x0;
  TH1D * hNewTruthBeamCalcIncidentHist = 0x0;

  TH1D * hTruthBeamCalcIncidentHist_50MeVbin = 0x0;

  TH1D * hTruthInitialHist = 0x0;
  TH1D * hTruthBeamInitialHist = 0x0;
  TH1D * hNewTruthBeamInitialHist = 0x0;
  TH1D * hTruthBeamIncidentHist = 0x0;
  TH1D * hNewTruthBeamIncidentHist = 0x0;

  TH1D * hTruthBeamIncidentHistOldM = 0x0;


  TH1D * hTruthBeamInitialHist_50MeVbin = 0x0;
  TH1D * hTruthBeamIncidentHist_50MeVbin = 0x0;

  TH1D * hTruthTestSameBin = 0x0;
  TH1D * hTruthTestFFEnergyM1 = 0x0;
  TH1D * hTruthTestFFEnergyM2 = 0x0;
  TH1D * hTruthTestFFEnergyDiff = 0x0;

  TH1D * hTruthTestIntEnergyM1 = 0x0;
  TH1D * hTruthTestIntEnergyM2 = 0x0;
  TH1D * hTruthTestIntEnergyDiff = 0x0;


  TH1D * hTruthTotalXSecHist = 0x0;
  TH1D * hTruthCEXInteractingHist = 0x0;
  TH1D * hNewTruthCEXInteractingHist = 0x0;
  TH1D * hNewTruthCEXInteractingHistTest = 0x0;

  TH1D * hTruthDiffCEXInteractingHist = 0x0;

  TH1D * hTruthDiffCEXInteractingHist_700MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHist_800MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHist_650to800MeV = 0x0;
  TH1D * hTruthSingleDiffCEXInteractingHist_800MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHist_900MeV = 0x0;

  TH1D * hTruthDiffCEXInteractingHistTheta_700MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHistTheta_800MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHistTheta_650to800MeV = 0x0;

  TH1D * hTruthDiffCEXInteractingHistTheta_900MeV = 0x0;

  TH1D * pi0theta = 0x0;

  TH1D * hTruthDiffCEXInteractingHistCosTheta_700MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHistCosTheta_800MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHistCosTheta_650to800MeV = 0x0;

  TH1D * hTruthDiffCEXInteractingHistCosTheta_900MeV = 0x0;

  TH1D * hTruthDiffCEXInteractingHist_500MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHist_300MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHist_200MeV = 0x0;
  TH1D * hTruthDiffCEXInteractingHist_100MeV = 0x0;

  TH2D * hTruthOutgoingKEcosTheta_Low = 0x0;
  TH2D * hTruthOutgoingKEcosTheta_Mid = 0x0;
  TH2D * hTruthOutgoingKEcosTheta_High = 0x0;

  TH2D * hTruthOutgoingKEcosTheta = 0x0;

  

  TH1D * hTruthIncidentZ = 0x0;
  TH1D * hTruthIncidentZ_SCE = 0x0;

  TH2D * hTruthCEXDaughters_Low = 0x0;
  TH2D * hTruthCEXDaughters_Middle = 0x0;
  TH2D * hTruthCEXDaughters_High = 0x0;

  TH1D * hTruthPi0DaughtersCosTheta_Low = 0x0;
  TH1D * hTruthPi0DaughtersCosTheta_Middle = 0x0;
  TH1D * hTruthPi0DaughtersCosTheta_High = 0x0;



  TH1D * hRecoInitialHist = 0x0;
  TH1D * hRecoBeamInteractingHist = 0x0;
  TH1D * hRecoIncidentHist = 0x0;
  TH1D * hRecoInteractingHist = 0x0;

  TH1D * hNewRecoInitialHist = 0x0;
  TH1D * hNewRecoBeamInteractingHist = 0x0;
  TH1D * hNewRecoIncidentHist = 0x0;
  TH1D * hNewRecoInteractingHist = 0x0;

  TH1D * hTruthMatchedIncidentHist = 0x0;
  TH1D * hTruthMatchedInteractingHist = 0x0;

  TH1D * hRecoBckSubIncidentHist = 0x0;
  TH1D * hRecoBckSubInteractingHist = 0x0;

  TH1D * hNonPionBeamIncidentHist = 0x0;
  TH1D * hNonPionBeamInteractingHist = 0x0;

  TH1D * hUnFoldedInitialHist = 0x0;
  TH1D * hUnFoldedBeamInteractingHist = 0x0;
  TH1D * hUnFoldedInteractingHist = 0x0;
  TH1D * hUnFoldedIncidentHist = 0x0;
  TH1D * hUnfoldedBeamIncidentHist = 0x0;

  TH1D * hRecoPi0KEHist = 0x0;
  TH1D * hUnFoldedPi0KEHist = 0x0;

  TH1D * hRecoPi0CosThetaHist = 0x0;
  TH1D * hUnFoldedPi0CosThetaHist = 0x0;

  TH1D * hRecoPi0ThetaHist = 0x0;
  TH1D * hUnFoldedPi0ThetaHist = 0x0;

  TH2D * hRecPiPlusEnergy_OVERLAY_After_EVTXS = 0x0;
  TH2D * hRecPi0Energy_OVERLAY_After_EVTXS = 0x0;

  TH2D * hRecPiPlusInteractingEnergy = 0x0;
  TH2D * hRecPiPlusInstMomentum = 0x0;
  TH2D * hRecPiPlusInstMomentumNoSmearing = 0x0;
  TH2D * hRecPiPlusInstEnergyNoSmearing = 0x0;
  TH2D * hRecPiPlusFrontFaceEnergy = 0x0;
  TH2D * hRecPiPlusIncidentEnergy = 0x0;
  TH2D * hRecPiPlusInitialEnergy = 0x0;
  TH2D * hRecPiPlusInitialEnergyPosZCut = 0x0;
  TH2D * hRecPiPlusIncidentEnergyNew = 0x0;
  TH2D * hRecPiPlusInteractingEnergyPar = 0x0;
  TH2D * hRecPiPlusInteractingEnergyPar_NoTune = 0x0;
  TH2D * hRecPiPlusTrackLength = 0x0;
  TH2D * hRecPiPlusCaloZ = 0x0;

  TH2D * hRawInstMomentumRaw = 0x0;
  TH2D * hRawInstEnergyRaw = 0x0;
  TH2D * hRawFrontFaceERaw = 0x0;
  TH2D * hRecPiPlusFrontFaceENoSmearing = 0x0;
  TH2D * hRecPiPlusFrontFaceE = 0x0;


  TH2D * hRecPiPlusInteractingEnergyBckSub = 0x0;
  TH2D * hRecPiPlusInteractingEnergyBckSubCheck = 0x0;
  
  TH2D * hRecPiPlusInteractingEnergyEvt = 0x0;
  TH2D * hRecPiPlusInteractingEnergyEvt_NoTune = 0x0;
  TH2D * hRecPiZeroKineticEnergyEvtNoWeight = 0x0;
  TH2D * hRecPiZeroKineticEnergyEvt = 0x0;
  TH2D * hRecPiZeroSliceKineticEnergyEvt = 0x0;
  TH2D * hRecPiZeroSlice2KineticEnergyEvt = 0x0;
  TH2D * hRecPiZeroSlice3KineticEnergyEvt = 0x0;

  TH2D * hRecPiZeroRangeKineticEnergyEvt = 0x0;
  TH2D * hRecPiZeroRangeKineticEnergyEvtNoWeight = 0x0;
  TH2D * hRecPiZeroRangeKineticEnergyEvtOneWeight = 0x0;
  TH2D * hRecPiZeroRangeKineticEnergyEvtAnotherWeight = 0x0;
  TH2D * hRecPiZeroRangeKineticEnergyEvtNoWeight_LargeBin = 0x0;
  TH2D * hRecPiZeroRangeKineticEnergyEvtRaw_LargeBin = 0x0;
  
  TH2D * hRecPiZeroRangeCosThetaEvt = 0x0;
  TH2D * hRecPiZeroRangeCosThetaEvtNoWeight = 0x0;
  TH2D * hRecPiZeroRangeCosThetaEvtOneWeight = 0x0;
  TH2D * hRecPiZeroRangeCosThetaEvtAnotherWeight = 0x0;

  TH2D * hRecPiZeroRangeThetaEvt = 0x0;
  TH2D * hRecPiZeroRangeThetaEvtNoWeight = 0x0;
  TH2D * hRecPiZeroRangeThetaEvtOneWeight = 0x0;
  TH2D * hRecPiZeroRangeThetaEvtAnotherWeight = 0x0;

  

  TH1D * hRecPiZeroSliceKineticEnergyEvtMC = 0x0;
  TH1D * hRecPiZeroSliceKineticEnergyEvtData = 0x0;

  // E-slice method histogram
  TH1D * hOldIncSizeDiff = 0x0;
  TH1D * hNewIncIntervalDiff = 0x0;

  // Beam Scraper cuts
  TH2D * hBeamInstVSTruthKEffNoCuts = 0x0;
  TH2D * hBeamInstVSTruthKEffAfterCuts = 0x0;
  TH1D * hUpStreamELoss700MeV = 0x0;
  TH1D * hUpStreamELoss800MeV = 0x0;
  TH1D * hUpStreamELoss900MeV = 0x0;
  TH1D * hUpStreamELoss1000MeV = 0x0;

  TH2D * hBeamInstXVSBeamInstYBeam700MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYScraper700MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYBeam800MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYScraper800MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYBeam900MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYScraper900MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYBeam1000MeV = 0x0;
  TH2D * hBeamInstXVSBeamInstYScraper1000MeV = 0x0;

  TH2D * hBeamInstVSTruthKEffAfterScraperCuts = 0x0;
  TH1D * hUpStreamELoss700MeVAfterScraperCuts = 0x0;
  TH1D * hUpStreamELoss800MeVAfterScraperCuts = 0x0;
  TH1D * hUpStreamELoss900MeVAfterScraperCuts = 0x0;
  TH1D * hUpStreamELoss1000MeVAfterScraperCuts = 0x0;

  TH2D * hUpStreamELossAfterSmearingAndWeight = 0x0;
  TH2D * hUpStreamELossRes = 0x0;

  // Truth TKI related variables
  TH1I * hTruthNproton = 0x0;
  TH1I * hTruthNneutron = 0x0;
  TH1I * hTruthNPiZero = 0x0;

  TH1D * hTruthMomIniPi = 0x0;
  TH1D * hTruthThetaIniPi = 0x0;
  TH1D * hTruthMomFinPi = 0x0;
  TH1D * hTruthThetaFinPi = 0x0;

  TH1D * hTruthMomFinProton = 0x0;
  TH1D * hTruthThetaFinProton = 0x0;

  TH1D * hTruthMomFin2Proton = 0x0;

  TH1D * hTruthDalphat = 0x0;
  TH1D * hTruthDphit = 0x0;
  TH1D * hTruthDpt = 0x0;
  TH1D * hTruthPn = 0x0;

  THStack * stkTruthDalphat = 0x0;
  TH1D * hTruthDalphat1p0n = 0x0;
  TH1D * hTruthDalphatNp0n = 0x0;
  TH1D * hTruthDalphat1pMn = 0x0; 
  TH1D * hTruthDalphatNpMn = 0x0;

  THStack * stkTruthDphit = 0x0;
  TH1D * hTruthDphit1p0n = 0x0;
  TH1D * hTruthDphitNp0n = 0x0;
  TH1D * hTruthDphit1pMn = 0x0; 
  TH1D * hTruthDphitNpMn = 0x0;

  THStack * stkTruthDpt = 0x0;
  TH1D * hTruthDpt1p0n = 0x0;
  TH1D * hTruthDptNp0n = 0x0;
  TH1D * hTruthDpt1pMn = 0x0; 
  TH1D * hTruthDptNpMn = 0x0;

  THStack * stkTruthPn = 0x0;
  TH1D * hTruthPn1p0n = 0x0;
  TH1D * hTruthPnNp0n = 0x0;
  TH1D * hTruthPn1pMn = 0x0; 
  TH1D * hTruthPnNpMn = 0x0;

  THStack * stkTruthMomIniPi = 0x0;
  TH1D * hTruthMomIniPi1p0n = 0x0;
  TH1D * hTruthMomIniPiNp0n = 0x0;
  TH1D * hTruthMomIniPi1pMn = 0x0; 
  TH1D * hTruthMomIniPiNpMn = 0x0;

  THStack * stkTruthThetaIniPi = 0x0;
  TH1D * hTruthThetaIniPi1p0n = 0x0;
  TH1D * hTruthThetaIniPiNp0n = 0x0;
  TH1D * hTruthThetaIniPi1pMn = 0x0; 
  TH1D * hTruthThetaIniPiNpMn = 0x0;

  THStack * stkTruthMomFinPi = 0x0;
  TH1D * hTruthMomFinPi1p0n = 0x0;
  TH1D * hTruthMomFinPiNp0n = 0x0;
  TH1D * hTruthMomFinPi1pMn = 0x0; 
  TH1D * hTruthMomFinPiNpMn = 0x0;

  THStack * stkTruthThetaFinPi = 0x0;
  TH1D * hTruthThetaFinPi1p0n = 0x0;
  TH1D * hTruthThetaFinPiNp0n = 0x0;
  TH1D * hTruthThetaFinPi1pMn = 0x0; 
  TH1D * hTruthThetaFinPiNpMn = 0x0;

  THStack * stkTruthMomFinProton = 0x0;
  TH1D * hTruthMomFinProton1p0n = 0x0;
  TH1D * hTruthMomFinProtonNp0n = 0x0;
  TH1D * hTruthMomFinProton1pMn = 0x0; 
  TH1D * hTruthMomFinProtonNpMn = 0x0;

  THStack * stkTruthThetaFinProton = 0x0;
  TH1D * hTruthThetaFinProton1p0n = 0x0;
  TH1D * hTruthThetaFinProtonNp0n = 0x0;
  TH1D * hTruthThetaFinProton1pMn = 0x0; 
  TH1D * hTruthThetaFinProtonNpMn = 0x0;

  THStack * stkTruthMomFin2Proton = 0x0;
  TH1D * hTruthMomFin2Proton1p0n = 0x0;
  TH1D * hTruthMomFin2ProtonNp0n = 0x0;
  TH1D * hTruthMomFin2Proton1pMn = 0x0; 
  TH1D * hTruthMomFin2ProtonNpMn = 0x0;

  // Resolution histograms need to use truth info
  TH2D * hBeamTrackLengthRes = 0x0;
  TH2D * hBeamPhiRes = 0x0;
  TH2D * hBeamThetaRes = 0x0;
  TH2D * hBeamMomentumRes = 0x0;
  TH2D * hProtonThetaRes = 0x0;
  TH2D * hProtonThetaLabRawRes = 0x0;
  TH2D * hProtonPhiLabRawRes = 0x0;
  TH2D * hProtonThetaLabRes = 0x0;
  TH2D * hProtonPhiLabRes = 0x0;
  TH2D * hProtonMomentumRes = 0x0;
  TH2D * hProtonMomentumRawRes = 0x0;
  TH2D * hPiPlusThetaRes = 0x0;
  TH2D * hPiPlusMomentumRes = 0x0;
  TH2D * hShowerEnergyRes = 0x0;
  TH2D * hShowerEnergyResRaw = 0x0;
  TH2D * hShowerThetaRes = 0x0; 
  TH2D * hShowerThetaResRaw = 0x0;
  TH2D * hShowerPhiRes = 0x0; 
  TH2D * hShowerPhiResRaw = 0x0;

  TH2D * hLeadingPhotonAngleRes = 0x0;
  TH2D * hSubLeadingPhotonAngleRes = 0x0;

  TH2D * hOpeningAngleRes = 0x0;

  TH2D * hProtonMomentumRawRecVSTruth_REG = 0x0;
  TH2D * hProtonMomentumRecVSTruth_REG = 0x0;
  TH2D * hProtonTransverseMomentumRawRecVSTruth_REG = 0x0;
  TH2D * hProtonTransverseMomentumRecVSTruth_REG = 0x0;

  TH2D * hProtonMomentumRecVSTruth_REG_Correction = 0x0;
  TH2D * hProtonMomentumRecVSTruth_REG_AfterCor = 0x0;
  TH2D * hProtonTransverseMomentumRecVSTruth_REG_Correction = 0x0;
  TH2D * hShowerEnergyRecVSTruth_REG_Correction = 0x0;
  TH2D * hShowerEnergyRecVSTruth_REG_AfterCor = 0x0;
  TH2D * hShowerThetaRecVSTruth_REG_AfterCor = 0x0;
  TH2D * hShowerPhiRecVSTruth_REG_AfterCor = 0x0;
  TH2D * hProtonThetaRecVSTruth_REG_AfterCor = 0x0;
  TH2D * hProtonPhiRecVSTruth_REG_AfterCor = 0x0;

  TH2D * hLeadingShowerEnergyRecVSTruth_REG_Correction = 0x0;
  TH2D * hSubLeadingShowerEnergyRecVSTruth_REG_Correction = 0x0;
  TH2D * hOpeningAngleRecVSTruth_REG_Correction = 0x0;
  TH2D * hShowerThetaRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanShowerTheta = 0x0;
  TH2D * hShowerPhiRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanShowerPhi = 0x0;
  TH2D * hProtonThetaRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanProtonTheta = 0x0;
  TH2D * hProtonPhiRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanProtonPhi = 0x0;

  TH2D * hBeamMomentumRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanBeamPMom = 0x0;
  TH2D * hBeamThetaRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanBeamTheta = 0x0;
  TH2D * hBeamPhiRecVSTruth_REG_Correction = 0x0;
  TH1D * hMeanBeamPhi = 0x0;

  TH1D * hMeanPMom = 0x0;
  TH1D * hMeanPMomT = 0x0;
  TH1D * hMeanShowerE = 0x0;
  TH1D * hMeanLDShowerE = 0x0;
  TH1D * hMeanSLShowerE = 0x0;

  TH2D * hShowerEnergyRecVSTruth_REG = 0x0;
  TH2D * hShowerEnergyRawRecVSTruth_REG = 0x0;
  TH2D * hShowerEnergyResVSnHits_REG = 0x0;
  TH2D * hShowerEnergyResVSIP_REG = 0x0;
 
  TH2D * hLeadingShowerEnergyRes = 0x0;
  TH2D * hSubLeadingShowerEnergyRes = 0x0;
  TH2D * hLeadingShowerEnergyResRaw = 0x0;
  TH2D * hSubLeadingShowerEnergyResRaw = 0x0;

  TH2D * hShowerOpenAngleRes = 0x0;
  TH2D * hShowerOpenAngleResFit = 0x0;
  TH2D * hPi0MomentumRes = 0x0;
  TH2D * hPi0MassRes = 0x0;
  TH2D * hPi0MomentumResRaw = 0x0;
  TH2D * hPi0MassResRaw = 0x0;
  TH2D * hPi0MomentumResFitKF = 0x0;
  TH2D * hPi0MomentumResFitKF_M2 = 0x0;
  TH2D * hPi0MassResFit = 0x0;

  TH2D * hPi0ThetaRes = 0x0;
  TH2D * hPi0ThetaResRaw = 0x0;
  TH2D * hPi0PhiRes = 0x0;
  TH2D * hPi0PhiResRaw = 0x0;

  TH1D * hMatchedTruthBeamMomentum = 0x0;
  TH1D * hMatchedTruthPi0Energy = 0x0;
  TH1D * hMatchedTruthldShowerTheta = 0x0;
  TH1D * hMatchedTruthslShowerTheta = 0x0;
  TH1D * hMatchedTruthPi0Mass = 0x0;
  TH1D * hMatchedTruthldShowerEnergy = 0x0;
  TH1D * hMatchedTruthslShowerEnergy = 0x0;

  TH1D * hMatchedTruthProtonMomentum = 0x0;
  TH1D * hMatchedTruthPiPlusMomentum = 0x0;
  TH1D * hMatchedTruthLeadingPiPlusMomentum = 0x0;
  TH1D * hMatchedTruthLeadingShowerEnergy = 0x0;
  TH1D * hMatchedTruthLeadingProtonMomentum = 0x0;

  TH1D * hMatchedTruthPi0OA = 0x0;
  TH1D * hMatchedTruthldShowerOA = 0x0;
  TH1D * hMatchedTruthslShowerOA = 0x0;

  TH2D * hPi0ThetaResFit = 0x0;
  TH2D * hPi0PhiResFit = 0x0;
  TH2D * hPi0MomRes = 0x0;
  TH2D * hPi0MomPreRes = 0x0;

  // CVM bin histograms
  TH2D * hV11_1x1 = 0x0;
  TH2D * hV12_1x1 = 0x0;
  TH2D * hV22_1x1 = 0x0;

  TH2D * hV11_2x2 = 0x0;
  TH2D * hV12_2x2 = 0x0;
  TH2D * hV22_2x2 = 0x0;

  TH2D * hV11_3x3 = 0x0;
  TH2D * hV12_3x3 = 0x0;
  TH2D * hV22_3x3 = 0x0;

  TH2D * hV11_4x4 = 0x0;
  TH2D * hV12_4x4 = 0x0;
  TH2D * hV22_4x4 = 0x0;

  TH2D * hBinSize_1x1 = 0x0;
  TH2D * hBinSize_2x2 = 0x0;
  TH2D * hBinSize_3x3 = 0x0;
  TH2D * hBinSize_4x4 = 0x0;

  TH2D * hCVM = 0x0;

  TH2D * hsigma11 = 0x0;
  TH1D * hsigma11Mean = 0x0;

  TH2D * hsigma22 = 0x0;
  TH1D * hsigma22Mean = 0x0;

  TH2D * hShowerE1PreFitRes = 0x0;
  TH2D * hShowerE1PostFitRes = 0x0;

  TH2D * hShowerE2PreFitRes = 0x0;
  TH2D * hShowerE2PostFitRes = 0x0;

  TH2D * hShowerOAPreFitRes = 0x0;
  TH2D * hShowerOAPostFitRes = 0x0;

  TH2D * hPi0MomPreFitRes = 0x0;
  TH2D * hPi0MomPostFitRes = 0x0;


  TH2D * hPi0MomPreFitResPaper = 0x0;
  TH2D * hPi0MomPostFitResPaper = 0x0;

  TH2D * hPi0ThetaFitRes = 0x0;
  TH2D * hLDShowerThetaFitRes = 0x0;

  TH1D * hShowerE1Compare = 0x0;
  TH1D * hShowerE2Compare = 0x0;
  TH1D * hShowerOACompare = 0x0;

  TH1D * hShowerE1ComparePost = 0x0;
  TH1D * hShowerE2ComparePost = 0x0;
  TH1D * hShowerOAComparePost = 0x0;

  TH1D * hPi0MassCompare = 0x0;
  TH1D * hPi0MassComparePost = 0x0;

  TH1D * hPi0MomCompare = 0x0;
  TH1D * hPi0MomComparePost = 0x0;
  TH1D * hPi0MomComparePre = 0x0;
  TH1D * hPi0MomNorm = 0x0;
  TH1D * hPi0MomEOA = 0x0;
  TH1D * hPi0MomAsym = 0x0;

  TH1D * hPrePi0MomNorm = 0x0;
  TH1D * hPrePi0MomEOA = 0x0;
  TH1D * hPrePi0MomAsym = 0x0;

  TH1D * hPi0MomComparePaper = 0x0;
  TH1D * hPi0MomComparePostPaper = 0x0;
  TH1D * hPi0MomComparePrePaper = 0x0;

  TH1D * hPi0MomCompareStandard = 0x0;
  TH1D * hPi0MomComparePostStandard = 0x0;
  TH1D * hPi0MomComparePreStandard = 0x0;

  TH2D * hPi0ShowerE1Compare_REG = 0x0;
  TH2D * hPi0ShowerE2Compare_REG = 0x0;
  TH2D * hPi0ShowerOACompare_REG = 0x0;
  TH2D * hPi0MomCompare_REG = 0x0;

  TH2D * hPi0MomCompare_REGPaper = 0x0;
  TH2D * hPi0MomCompare_REGStand = 0x0;

  TH1D * hLDShowerE = 0x0;
  TH1D * hSLShowerE = 0x0;
  TH1D * hOAShower = 0x0;

  TH1D * hKFPassRate = 0x0;
  TH1D * hPi0Total = 0x0;
  TH1D * hPi0Selected = 0x0;

  // Completeness and Purity
  TH2D * hTrackPurityVSnHits = 0x0;
  TH2D * hTrackCompletenessVSnHits = 0x0;
  TH2D * hTrackPurityVSnHitsColl = 0x0;
  TH2D * hTrackCompletenessVSnHitsColl = 0x0;
  TH2D * hShowerPurityVSnHits = 0x0;
  TH2D * hShowerCompletenessVSnHits = 0x0;
  TH2D * hShowerPurityVSnHitsColl = 0x0;
  TH2D * hShowerCompletenessVSnHitsColl = 0x0;
  TH2D * hShowerPurityVSIP = 0x0;
  TH2D * hShowerCompletenessVSIP = 0x0;
  TH2D * hShowerPurityVSEnergy = 0x0;
  TH2D * hShowerCompletenessVSEnergy = 0x0;
  TH2D * hShowernHitsVSEnergy = 0x0;

  TH2D * hTrackCompleteness = 0x0;
  TH2D * hShowerCompleteness = 0x0;

  TH1D * hTruthRecRatioProtonMom = 0x0;
  

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
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
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
    tree->SetBranchAddress("reco_beam_dEdX_SCE", &reco_beam_dEdX_SCE);
    tree->SetBranchAddress("reco_beam_dEdX_NoSCE", &reco_beam_dEdX_NoSCE);

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
    tree->SetBranchAddress("beam_inst_Z", &beam_inst_Z);
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
    tree->SetBranchAddress("g4rw_full_grid_piplus_coeffs",&g4rw_full_grid_piplus_coeffs);
    //tree->SetBranchAddress("g4rw_full_grid_piplus_weights",&g4rw_full_grid_piplus_weights);
    
    return tree;
  } // End of GetInputTree
  
  TTree * GetOutputTree(TList * lout, const TString tag)
  {
    TTree * tout = new TTree("tree", tag); lout->Add(tout);

    // Definitely need this to avoid memory-resident Tree
    tout->SetDirectory(0);
    tout->Branch("LeadingShowerEnergy",&LeadingShowerEnergy);
    tout->Branch("SubLeadingShowerEnergy",&SubLeadingShowerEnergy);
    tout->Branch("LeadingShowerEnergyRaw",&LeadingShowerEnergyRaw);
    tout->Branch("SubLeadingShowerEnergyRaw",&SubLeadingShowerEnergyRaw);
    tout->Branch("OpeningAngle",&OpeningAngle);
    tout->Branch("LeadingShowerEnergyTruth",&LeadingShowerEnergyTruth);
    tout->Branch("SubLeadingShowerEnergyTruth",&SubLeadingShowerEnergyTruth);
    tout->Branch("OpeningAngleTruth",&OpeningAngleTruth);

    tout->Branch("LeadingShowerEnergyUnitDir",&LeadingShowerEnergyUnitDir);
    tout->Branch("SubLeadingShowerEnergyUnitDir",&SubLeadingShowerEnergyUnitDir);
    tout->Branch("LeadingShowerEnergyUnitDirTruth",&LeadingShowerEnergyUnitDirTruth);
    tout->Branch("SubLeadingShowerEnergyUnitDirTruth",&SubLeadingShowerEnergyUnitDirTruth);
    return tout;
  }

  // Initialise histograms
  void IniHist(TList * lout, const bool kMC)
  {
    // Binning 
    const int nbeamType = 8;
    const double beamTypemin = -0.5;
    const double beamTypemax = 7.5;

    const int nchannelType = 7;
    const double channelTypemin = -0.5;
    const double channelTypemax = 6.5;

    const int nparType = 11;
    const double parTypemin = 0.5;
    const double parTypemax = 11.5;

    const int nevtXSType = 6;
    const double evtXSTypemin = -0.5;
    const double evtXSTypemax = 5.5;

    const int nevtType = 3;
    const double evtTypemin = -0.5;
    const double evtTypemax = 2.5;

    const int nevtTKIType = 5;
    const double evtTKITypemin = -0.5;
    const double evtTKITypemax = 4.5;

    const int ncounter = 10;
    const double countermin = -0.5;
    const double countermax = 9.5;

    const int nPass = 2;
    const double Passmin = -0.5;
    const double Passmax = 1.5;

    // Beam cuts EndZ binning
    const int nEndZ = 100;
    const double EndZmin = -100;
    const double EndZmax = 600;

    
    //=== Herilala mod ====//
    hxstartystart = new TH2D("hxstartystart",";StartX (cm);StartY (cm)", 400, -200, 200, 300, 300, 600); 
    lout->Add(hxstartystart);
    hxstartystartBeam = new TH2D("hxstartystartBeam",";StartX (cm);StartY (cm)", 400, -200, 200, 300, 300, 600);
    lout->Add(hxstartystartBeam); 
    hBeamDeltaXYSigma = new TH1D("hBeamDeltaXYSigma","",100, 0, 20);
    lout->Add(hBeamDeltaXYSigma);

    hrstartX = new TH1D("hrstartX",";StartX (cm)",100, -80, 20);
    lout->Add(hrstartX);
    hrstartY = new TH1D("hrstartY",";StartY (cm)",100, 360, 500);
    lout->Add(hrstartY);
    hrstartZ = new TH1D("hrstartZ",";StartZ (cm)",100, 20, 40);
    lout->Add(hrstartZ);

    hrendX = new TH1D("hrendX",";EndX (cm)",180, -160, 20);
    lout->Add(hrendX);
    hrendY = new TH1D("hrendY",";EndY (cm)",200, 240, 500);
    lout->Add(hrendY);
    hrendZ = new TH1D("hrendZ",";EndZ (cm)",120, -20, 100);
    lout->Add(hrendZ);

    hreco_beam_trackDirX = new TH1D("hreco_beam_trackDirX"," ",200, -4, 4);
    lout->Add(hreco_beam_trackDirX);
    hreco_beam_trackDirY = new TH1D("hreco_beam_trackDirY"," ",200, -4, 4);
    lout->Add(hreco_beam_trackDirY);
    hreco_beam_trackDirZ = new TH1D("hreco_beam_trackDirZ"," ",200, -4, 4);
    lout->Add(hreco_beam_trackDirZ);
    //====================== Reco (MC and Data)======================//

    // Class A - beam cut related
    hCutBeamPDGPass = new TH1I("a001hCutBeamPDGPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutBeamPDGPass);
    hCutPandoraSlicePass = new TH1I("a002hCutPandoraSlicePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutPandoraSlicePass); 
    hCutCaloSizePass = new TH1I("a003hCutCaloSizePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutCaloSizePass);
    hCutBeamQualityPass = new TH1I("a004hCutBeamQualityPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutBeamQualityPass);
    hCutAPA3EndZPass = new TH1I("a005hCutAPA3EndZPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutAPA3EndZPass);
    hCutMichelScorePass = new TH1I("a006hCutMichelScorePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutMichelScorePass);
    hCutProtonChi2Pass = new TH1I("a007hCutProtonChi2Pass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutProtonChi2Pass);
    hCutBeamScraperPass = new TH1I("a008hCutBeamScraperPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax);
    lout->Add(hCutBeamScraperPass);
    hRecBeamStartX = new TH1D("a009hRecBeamStartX_FCN",";StartX (cm);Area Normalised",100, -80, 20);
    lout->Add(hRecBeamStartX);
    hRecBeamStartY = new TH1D("a010hRecBeamStartY_FCN",";StartY (cm);Area Normalised",100, 350, 500);
    lout->Add(hRecBeamStartY);
    hRecBeamStartZ = new TH1D("a011hRecBeamStartZ_FCN",";StartZ (cm);Area Normalised",100, 20,40);
    //hRecBeamStartZ = new TH1D("a011hRecBeamStartZ_FCN",";StartZ (cm);Area Normalised",100, -20,80);
    lout->Add(hRecBeamStartZ);
    hRecBeamThetaX = new TH1D("a012hRecBeamThetaX_FCN",";Angle #theta_{x} (deg);Area Normalised",120, 60, 140);
    //hRecBeamThetaX = new TH1D("a012hRecBeamThetaX_FCN",";Angle #theta_{x} (deg);Area Normalised",217.5, -5, 140);
    lout->Add(hRecBeamThetaX);
    hRecBeamThetaY = new TH1D("a013hRecBeamThetaY_FCN",";Angle #theta_{y} (deg);Area Normalised",120, 60, 140);
    //hRecBeamThetaY = new TH1D("a013hRecBeamThetaY_FCN",";Angle #theta_{y} (deg);Area Normalised",217.5, -5, 140);
    lout->Add(hRecBeamThetaY);
    hRecBeamThetaZ = new TH1D("a014hRecBeamThetaZ_FCN",";Angle #theta_{z} (deg);Area Normalised",100, 5, 45);
    //hRecBeamThetaZ = new TH1D("a014hRecBeamThetaZ_FCN",";Angle #theta_{z} (deg);Area Normalised",285.29, -10, 180);
    lout->Add(hRecBeamThetaZ);
    hCutBeamDeltaXYSigma = new TH1D("a015hCutBeamDeltaXYSigma","",100, 0, 20);
    lout->Add(hCutBeamDeltaXYSigma);
    hCutBeamDeltaZSigma = new TH1D("a016hCutBeamDeltaZSigma","",100, -10, 10);
    lout->Add(hCutBeamDeltaZSigma);
    hCutBeamCosTheta = new TH1D("a017hCutBeamCosTheta","",100, 0.9, 1);
    lout->Add(hCutBeamCosTheta);
    //hCutBeamMichelScoreNhits = new TH2D("a018hCutBeamMichelScoreNhits_STK","",100, 0, 1, nbeamType, beamTypemin, beamTypemax);
    hCutBeamMichelScoreNhits = new TH2D("a018hCutBeamMichelScoreNhits_STK","",50, 0, 1, nbeamType, beamTypemin, beamTypemax);
    lout->Add(hCutBeamMichelScoreNhits);
    //hCutBeamProtonChi2DOF = new TH2D("a019hCutBeamProtonChi2DOF_STK","", 100, 0, 350, nbeamType, beamTypemin, beamTypemax); 
    hCutBeamProtonChi2DOF = new TH2D("a019hCutBeamProtonChi2DOF_STK","", 50, 0, 300, nbeamType, beamTypemin, beamTypemax);
    lout->Add(hCutBeamProtonChi2DOF);
    hCutBeamProtonChi2DOFType = new TH1D("a020hCutBeamProtonChi2DOFType","",23, -0.5, 22.5);
    lout->Add(hCutBeamProtonChi2DOFType);
    hCutBeamQualityXY = new TH2D("a021hCutBeamQualityXY_STK",";#sqrt{(#Deltax/#sigma_{x})^{2}+(#Deltay/#sigma_{y})^{2}};Candidates", 100, 0, 5, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityXY);
    hCutBeamQualityX = new TH2D("a022hCutBeamQualityX_STK",";#Deltax/#sigma_{x};Candidates", 100, -10, 10, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityX);
    hCutBeamQualityY = new TH2D("a023hCutBeamQualityY_STK",";#Deltay/#sigma_{y};Candidates", 100, -10, 10, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityY);
    //hCutBeamQualityZ = new TH2D("a024hCutBeamQualityZ_STK",";#Deltaz/#sigma_{z};Candidates", 100, -10, 10, nbeamType, beamTypemin, beamTypemax);
    hCutBeamQualityZ = new TH2D("a024hCutBeamQualityZ_STK",";#Deltaz/#sigma_{z};Candidates", 80, -8, 8, nbeamType, beamTypemin, beamTypemax);
    lout->Add(hCutBeamQualityZ);
    hCutBeamQualityTheta = new TH2D("a025hCutBeamQualityTheta_STK",";cos(#alpha);Candidates", 100, 0.9, 1, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityTheta);

    hCutAPA3EndZ = new TH2D("a026hCutAPA3EndZ_STK",";End Z (cm);Candidates", 100, 0, 500, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutAPA3EndZ);

    hBeamemScore = new TH2D("a027hBeamemScore_STK",";Beam EM Shower Score;Candidates", 50, 0, 1, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamemScore);
    hBeamemScore_wbc = new TH2D("a028hBeamemScore_wbc_STK",";Beam EM Shower Score;Candidates", 50, 0, 1, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamemScore_wbc);
    hBeamDaughterNumbers = new TH2D("a029hBeamDaughterNumbers_STK",";Number of Daughter Particles;Candidates", 15, 0, 15, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamDaughterNumbers);

    hBeamtrackScore = new TH2D("a030hBeamtrackScore_STK",";Beam Track Shower Score;Candidates", 50, 0, 1, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamtrackScore);
    hBeamtrackScore_wbc = new TH2D("a031hBeamtrackScore_wbc_STK",";Beam Track Shower Score;Candidates", 50, 0, 1, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamtrackScore_wbc);

    hBeamStartX = new TH2D("a032hBeamStartX_STK",";StartX (cm);Candidates",100, -80, 20, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamStartX);
    hBeamStartY = new TH2D("a033hBeamStartY_STK",";StartY (cm);Candidates",100, 350, 500, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamStartY);
    hBeamStartZ = new TH2D("a034hBeamStartZ_STK",";StartZ (cm);Candidates",100, -5, 10, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamStartZ);
    hBeamStartCaloX = new TH2D("a035hBeamStartCaloX_STK",";StartCaloX (cm);Candidates",100, -80, 20, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamStartCaloX);
    hBeamStartCaloY = new TH2D("a036hBeamStartCaloY_STK",";StartCaloY (cm);Candidates",100, 350, 500, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamStartCaloY);
    hBeamStartCaloZ = new TH2D("a037hBeamStartCaloZ_STK",";StartCaloZ (cm);Candidates",100, -5, 10, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamStartCaloZ);

    hCutBeamQualityInstXY = new TH2D("a038hCutBeamQualityInstXY_STK","", 80, 0, 5, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityInstXY);
    hCutBeamQualityInstX = new TH2D("a039hCutBeamQualityInstX_STK","", 80, -5, 5, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityInstX);
    hCutBeamQualityInstY = new TH2D("a040hCutBeamQualityInstY_STK","", 80, -5, 5, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hCutBeamQualityInstY);

    hRecBeamInstX = new TH1D("a041hRecBeamInstX_FCN",";InstX (cm);Area Normalised",100, -50, 0);
    lout->Add(hRecBeamInstX);
    hRecBeamInstY = new TH1D("a042hRecBeamInstY_FCN",";InstY (cm);Area Normalised",100, 400, 440);
    lout->Add(hRecBeamInstY);


    hBeamEndX = new TH2D("a043hBeamEndX_STK",";EndX (cm);Candidates",100, -140, 40, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndX);
    hBeamEndY = new TH2D("a044hBeamEndY_STK",";EndY (cm);Candidates",100, 300, 500, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndY);
    hBeamEndZ = new TH2D("a045hBeamEndZ_STK",";EndZ (cm);Candidates",100, 0, 300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ);
    hBeamEndCaloX = new TH2D("a046hBeamEndCaloX_STK",";EndCaloX (cm);Candidates",100, -140, 40, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndCaloX);
    hBeamEndCaloY = new TH2D("a047hBeamEndCaloY_STK",";EndCaloY (cm);Candidates",100, 300, 500, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndCaloY);
    hBeamEndCaloZ = new TH2D("a048hBeamEndCaloZ_STK",";EndCaloZ (cm);Candidates",100, 0, 300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndCaloZ);

    hBeamDirDiff = new TH2D("a049hBeamDirDiff_STK",";Angle Diff (deg.);Candidates",30, -10, 50, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamDirDiff);

    hRecBeamTrackPitch_SCE = new TH1D("a050hRecBeamTrackPitch_SCE_FCN",";Track Pitch (cm);Area Normalised",80, 0.2, 1);
    lout->Add(hRecBeamTrackPitch_SCE);
    hRecBeamTrackPitch_NoSCE = new TH1D("a051hRecBeamTrackPitch_NoSCE_FCN",";Track Pitch (cm);Area Normalised",80, 0.2, 1);
    lout->Add(hRecBeamTrackPitch_NoSCE);

    hRecBeamdEdx_SCE = new TH1D("a052hRecBeamdEdx_SCE_FCN",";dEdx (MeV/cm);Area Normalised",80, 0.0, 7);
    lout->Add(hRecBeamdEdx_SCE);
    hRecBeamdEdx_NoSCE = new TH1D("a053hRecBeamdEdx_NoSCE_FCN",";dEdx (MeV/cm);Area Normalised",80, 0.0, 7);
    lout->Add(hRecBeamdEdx_NoSCE);

    hRecBeamDeltaE_SCE = new TH1D("a054hRecBeamDeltaE_SCE_FCN",";DeltaE (MeV);Area Normalised",80, 0.0, 3.5);
    lout->Add(hRecBeamDeltaE_SCE);
    hRecBeamDeltaE_NoSCE = new TH1D("a055hRecBeamDeltaE_NoSCE_FCN",";DeltaE (MeV);Area Normalised",80, 0.0, 3.5);
    lout->Add(hRecBeamDeltaE_NoSCE);

    hRecIntEdiff = new TH1D("a056hRecIntEdiff",";DeltaE (MeV);Area Normalised",80, -200, 200);
    lout->Add(hRecIntEdiff);

    // Class B - reconstructed beam/FS particles 
    hRecBeamPhi = new TH2D("b000hRecBeamPhi_STK",";Beam #Phi (deg);Candidates", 80 , -180, 0, 3, -0.5, 2.5); 
    lout->Add(hRecBeamPhi);
    hRecBeamPhi_XsEvt = new TH2D("b000hRecBeamPhi_COMPOSE_XsEvt",";Beam #Phi (deg);Candidates", 80 , -180, 0, nevtXSType, evtXSTypemin, evtXSTypemax);
    lout->Add(hRecBeamPhi_XsEvt);
    hRecBeamTheta = new TH2D("b001hRecBeamTheta_STK",";Beam #theta (deg);Candidates", 80 , 0, 60, 3, -0.5, 2.5); 
    lout->Add(hRecBeamTheta);
    hRecBeamTheta_XsEvt = new TH2D("b001hRecBeamTheta_COMPOSE_XsEvt",";Beam #theta (deg);Candidates", 80 , 0, 60, nevtXSType, evtXSTypemin, evtXSTypemax);
    lout->Add(hRecBeamTheta_XsEvt);
    //hRecBeamMomentum = new TH2D("b002hRecBeamMomentum_STK",";Beam Momentum (GeV/c);Candidates", 50, 0, 2, 3, -0.5, 2.5); 
    hRecBeamMomentum = new TH2D("b002hRecBeamMomentum_STK",";Beam Momentum (GeV/c);Candidates", 50, 0, 2, 3, -0.5, 2.5);
    lout->Add(hRecBeamMomentum);
    hRecBeamMomentum_XsEvt = new TH2D("b002hRecBeamMomentum_COMPOSE_XsEvt",";Beam Momentum (GeV/c);Candidates", 50, 1, 2.5, nevtXSType, evtXSTypemin, evtXSTypemax);
    lout->Add(hRecBeamMomentum_XsEvt);

    hRecProtonTheta = new TH2D("b003hRecProtonTheta_STK",";Proton #theta (deg);Candidates",  20, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecProtonTheta);
    hRecProtonMomentum = new TH2D("b004hRecProtonMomentum_STK",";Proton Momentum (GeV/c);Candidates",  20, 0.1, 1.3, nparType, parTypemin, parTypemax); 
    lout->Add(hRecProtonMomentum);
    hRecPiPlusTheta = new TH2D("b005hRecPiPlusTheta_STK",";#pi^{+} #theta (deg);Candidates", 30, 0, 180, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusTheta);
    hRecPiPlusMomentum = new TH2D("b006hRecPiPlusMomentum_STK",";#pi^{+} Momentum (GeV/c);Candidates",  24, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusMomentum);
    hRecShowerEnergy = new TH2D("b007hRecShowerEnergy_STK",";Shower Energy (GeV);Candidates", 20, 0, 0.5, nparType, parTypemin, parTypemax); 
    lout->Add(hRecShowerEnergy);
    hRecShowerEnergyRaw = new TH2D("b008hRecShowerEnergy_STK_RAW",";Shower Energy Raw (GeV);Candidates", 20, 0, 0.5, nparType, parTypemin, parTypemax);
    lout->Add(hRecShowerEnergyRaw);
    hRecShowerTheta = new TH2D("b009hRecShowerTheta_STK",";Shower #theta (deg);Candidates", 20, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecShowerTheta);
    hRecShowerLength = new TH2D("b010hRecShowerLength_STK",";Shower Length (cm);Candidates", 20, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hRecShowerLength);
    hRecPi0Nshower = new TH2D("b011hRecPi0Nshower_STK",";Number of Showers;Candidates", 7, -0.5, 6.5, 3, -0.5, 2.5); 
    lout->Add(hRecPi0Nshower);
    hRecLeadingShowerEnergy = new TH2D("b012hRecLeadingShowerEnergy_STK",";Leading Shower Energy (GeV);Candidates", 20, 0, 0.8, nparType, parTypemin, parTypemax);
    lout->Add(hRecLeadingShowerEnergy);
    hRecSubLeadingShowerEnergy = new TH2D("b013hRecSubLeadingShowerEnergy_STK",";SubLeading Shower Energy (GeV);Candidates", 20, 0, 0.4, nparType, parTypemin, parTypemax);
    lout->Add(hRecSubLeadingShowerEnergy); 
    hRecLeadingShowerEnergyRaw = new TH2D("b014hRecLeadingShowerEnergy_STK_RAW",";Leading Shower Energy Raw (GeV);Candidates", 20, 0, 0.8, nparType, parTypemin, parTypemax);
    lout->Add(hRecLeadingShowerEnergyRaw);
    hRecSubLeadingShowerEnergyRaw = new TH2D("b015hRecSubLeadingShowerEnergy_STK_RAW",";SubLeading Shower Energy Raw(GeV);Candidates", 20, 0, 0.4, nparType, parTypemin, parTypemax);
    lout->Add(hRecSubLeadingShowerEnergyRaw);
    hRecShowerOpenAngle = new TH2D("b016hRecShowerOpenAngle_STK",";Shower Opening Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 2.5); 
    lout->Add(hRecShowerOpenAngle);
    hRecPi0Mass = new TH2D("b018hRecPi0Mass_STK",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Mass);
    hRecPi0MassRaw = new TH2D("b019hRecPi0Mass_STK_RAW",";#pi^{0} Mass Raw(GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5);
    lout->Add(hRecPi0MassRaw);
    hRecPi0Momentum = new TH2D("b020hRecPi0Momentum_STK", ";#pi^{0} Momentum (GeV/c);Candidates", 40, 0, 2, 3, -0.5, 3.5);
    lout->Add(hRecPi0Momentum);
    hRecPi0MomentumRaw = new TH2D("b021hRecPi0Momentum_STK_RAW", ";#pi^{0} Momentum Raw (GeV/c);Candidates", 40, 0, 2, 3, -0.5, 3.5);
    lout->Add(hRecPi0MomentumRaw);

    hPi0MassLowE1 = new TH1D("b020hPi0MassLowE1","", 20, 0, 0.5);
    lout->Add(hPi0MassLowE1);
    hPi0MassHighE1 = new TH1D("b021hPi0MassHighE1","", 20, 0, 0.5);
    lout->Add(hPi0MassHighE1);

    hRecPi0Mass_OVERLAY = new TH2D("b022hRecPi0Mass_COMPOSE",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Mass_OVERLAY);
    hRecPi0Mass_OVERLAY_EVT = new TH2D("b022hRecPi0Mass_COMPOSE_EVT",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, nevtXSType, evtXSTypemin, evtXSTypemax);
    lout->Add(hRecPi0Mass_OVERLAY_EVT);
 
    testhRecPi0Mass_OVERLAY = new TH2D("testb022hRecPi0Mass_COMPOSE",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(testhRecPi0Mass_OVERLAY);
    test1hRecPi0Mass_OVERLAY = new TH2D("test1b022hRecPi0Mass_COMPOSE",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(test1hRecPi0Mass_OVERLAY);

    testhRecPi0Mass_OVERLAY_EVT = new TH2D("testb022hRecPi0Mass_COMPOSE_EVT",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, nevtXSType, evtXSTypemin, evtXSTypemax);
    lout->Add(testhRecPi0Mass_OVERLAY_EVT);
    test1hRecPi0Mass_OVERLAY_EVT = new TH2D("test1b022hRecPi0Mass_COMPOSE_EVT",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, nevtXSType, evtXSTypemin, evtXSTypemax);
    lout->Add(test1hRecPi0Mass_OVERLAY_EVT);

    hRecPi0Momentum_OVERLAY = new TH2D("b023hRecPi0Momentum_COMPOSE",";#pi^{0} Momentum (GeV/c);Candidates", 40, 0, 2, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Momentum_OVERLAY);

    hRecPi0OA_OVERLAY = new TH2D("b024hRecPi0OA_COMPOSE",";Shower Opening Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0OA_OVERLAY);

    testhRecPi0OA_OVERLAY = new TH2D("testb024hRecPi0OA_COMPOSE",";Shower Opening Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(testhRecPi0OA_OVERLAY);
    test1hRecPi0OA_OVERLAY = new TH2D("test1b024hRecPi0OA_COMPOSE",";Shower Opening Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(test1hRecPi0OA_OVERLAY);

    hRecPi0Energy_OVERLAY = new TH2D("b023hRecPi0Energy_COMPOSE",";#pi^{0} Energy (GeV);Candidates", 40, 0, 2, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Energy_OVERLAY);

    hRecPi0ShowerSep_OVERLAY = new TH2D("b024hRecPi0ShowerSep_COMPOSE",";Shower Separation (cm);Candidates", 20, 0, 150, 3, -0.5, 3.5);
    lout->Add(hRecPi0ShowerSep_OVERLAY);

    hRecPi0MassRaw_OVERLAY = new TH2D("b025hRecPi0MassRaw_COMPOSE",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(hRecPi0MassRaw_OVERLAY);

    hRecPi0MomentumRaw_OVERLAY = new TH2D("b026hRecPi0MomentumRaw_COMPOSE",";#pi^{0} Momentum (GeV/c);Candidates", 40, 0, 2, 3, -0.5, 3.5); 
    lout->Add(hRecPi0MomentumRaw_OVERLAY);

    hRecPi0OARaw_OVERLAY = new TH2D("b024hRecPi0OARaw_COMPOSE",";Shower Opening Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0OARaw_OVERLAY);

    hRecPi0Theta_OVERLAY = new TH2D("b027hRecPi0Theta_COMPOSE",";#pi^{0} #theta Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Theta_OVERLAY);

    hRecPi0Phi_OVERLAY = new TH2D("b028hRecPi0Phi_COMPOSE",";#pi^{0} #phi Angle (deg);Candidates", 20, -180, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Phi_OVERLAY);

    hRecPi0ThetaRaw_OVERLAY = new TH2D("b029hRecPi0Theta_COMPOSE",";#pi^{0} #theta Angle Raw (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0ThetaRaw_OVERLAY);

    hRecPi0PhiRaw_OVERLAY = new TH2D("b030hRecPi0Phi_COMPOSE",";#pi^{0} #phi Angle Raw (deg);Candidates", 20, -180, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0PhiRaw_OVERLAY);


    hRecPi0Energy_OVERLAY_After = new TH2D("b0233hRecPi0Energy_COMPOSE_After",";#pi^{0} Energy (GeV);Candidates", 40, 0, 2, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After);
    hRecPi0Theta_OVERLAY_After = new TH2D("b0277hRecPi0Theta_COMPOSE_After",";#pi^{0} #theta Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Theta_OVERLAY_After);
    hRecPi0Phi_OVERLAY_After = new TH2D("b0277hRecPi0Phi_COMPOSE_After",";#pi^{0} #phi Angle (deg);Candidates", 20, -180, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Phi_OVERLAY_After);

    hRecPi0Energy_OVERLAY_After_Kinetic = new TH2D("b0244hRecPi0Energy_COMPOSE_After_Kinetic",";#pi^{0} Kinetic Energy (GeV);Candidates", 40, 0, 2, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After_Kinetic);

    hRecPi0Energy_OVERLAY_After_EVT = new TH2D("b0255hRecPi0KineticEnergy_COMPOSE_After_EVT",";#pi^{0} Kinetic Energy (GeV);Candidates", 40, 0, 2, 7, -0.5, 6.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After_EVT);

    hRecPi0Energy_OVERLAY_AfterTOP_EVT = new TH2D("b0266hRecPi0Energy_COMPOSE_AfterTOP_EVT",";#pi^{0} Kinetic Energy (GeV);Candidates", 40, 0, 4, 7, -0.5, 6.5); 
    lout->Add(hRecPi0Energy_OVERLAY_AfterTOP_EVT);

    hRecPi0Energy_OVERLAY_After_EVT_High = new TH2D("b0277hRecPi0Energy_COMPOSE_After_EVT_High",";#pi^{0} Kinetic Energy (GeV);Candidates", 40, 0, 2000, 7, -0.5, 6.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After_EVT_High);
    hRecPi0Energy_OVERLAY_After_EVT_Medium = new TH2D("b0288hRecPi0Energy_COMPOSE_After_EVT_Medium",";#pi^{0} Kinetic Energy (GeV);Candidates", 40, 0, 2000, 7, -0.5, 6.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After_EVT_Medium);
    hRecPi0Energy_OVERLAY_After_EVT_Low = new TH2D("b0299hRecPi0Energy_COMPOSE_After_EVT_Low",";#pi^{0} Kinetic Energy (GeV);Candidates", 40, 0, 2000, 7, -0.5, 6.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After_EVT_Low);

    hRecPi0TotalE_Default = new TH2D("test001hpi0TotalE_Default_COMPOSE_EVT",";Total E;Candidates", 40, 0, 2, 7, -0.5, 6.5); 
    lout->Add(hRecPi0TotalE_Default);
    hRecPi0TotalE_Fitted = new TH2D("test002hpi0TotalE_Fitted_COMPOSE_EVT",";Total E;Candidates", 40, 0, 2, 7, -0.5, 6.5); 
    lout->Add(hRecPi0TotalE_Fitted);

    hRecPi0TotalEafterSel_Default = new TH2D("test003hpi0TotalEafterSel_Default_COMPOSE_EVT",";Total E;Candidates", 40, 0, 2, 7, -0.5, 6.5); 
    lout->Add(hRecPi0TotalEafterSel_Default);
    hRecPi0TotalEafterSel_Fitted = new TH2D("test004hpi0TotalEafterSel_Fitted_COMPOSE_EVT",";Total E;Candidates", 40, 0, 2, 7, -0.5, 6.5); 
    lout->Add(hRecPi0TotalEafterSel_Fitted);

    hRecPi0KineticEnergyVSPi0Mass = new TH2D("b099hRecPi0KineticEnergyVSPi0Mass",";Kinetic Energy (GeV);Mass (GeV/c^{2})", 100, 0, 2, 50, 0, 0.5); 
    lout->Add(hRecPi0KineticEnergyVSPi0Mass);

    // Class C - event topoplogy cut related
    hCutDaughterPandoraShowerPass = new TH1I("c000hCutDaughterPandoraShowerPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPandoraShowerPass);
    hCutDaughterShowerScorePass = new TH1I("c001hCutDaughterShowerScorePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerScorePass);
    hCutDaughterShowernHitsPass = new TH1I("c002hCutDaughterShowernHitsPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowernHitsPass);
    hCutDaughterShowerNonEmptyEPass = new TH1I("c003hCutDaughterShowerNonEmptyEPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerNonEmptyEPass);
    hCutDaughterShowerStartZPass = new TH1I("c004hCutDaughterShowerStartZPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerStartZPass);
    hCutDaughterShowerDistPass = new TH1I("c005hCutDaughterShowerDistPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerDistPass);
    hCutDaughterShowerIPPass = new TH1I("c006hCutDaughterShowerIPPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerIPPass);

    hCutDaughterTrackScorePass = new TH1I("c007hCutDaughterTrackScorePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterTrackScorePass);
    hCutDaughterTracknHitsPass = new TH1I("c008hCutDaughterTracknHitsPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterTracknHitsPass);
    hCutDaughterProtonSubPIDPass = new TH1I("c009hCutDaughterProtonSubPIDPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterProtonSubPIDPass);


    hCutDaughterPionTrackScorePass = new TH1I("c007hCutDaughterPionTrackScorePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPionTrackScorePass);
    hCutDaughterPionTracknHitsPass = new TH1I("c008hCutDaughterPionTracknHitsPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPionTracknHitsPass);
    hCutDaughterPionSubPIDPass = new TH1I("c009hCutDaughterPionSubPIDPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPionSubPIDPass);

    hCutDaughterPi0NOCutsPass = new TH1I("c009hCutDaughterPi0NOCutsPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPi0NOCutsPass);
    hCutDaughterPi0MassPass = new TH1I("c009hCutDaughterPi0MassPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPi0MassPass);
    hCutDaughterPi0OAPass = new TH1I("c010hCutDaughterPi0OAPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPi0OAPass);


    hCutTracknHits = new TH2D("c001hCutTracknHits_STK",";Number of Hits - total;Candidates", 50, 0, 500, nparType, parTypemin, parTypemax); 
    lout->Add(hCutTracknHits);
    hCutTrackScore = new TH2D("c002hCutTrackScore_STK",";Track Score;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax);
    lout->Add(hCutTrackScore);
    hCutlastTME = new TH2D("c003hCutlastTME_STK",";Truncated Mean dEdx (MeV/cm);Candidates", 60, 0, 10, nparType, parTypemin, parTypemax); 
    lout->Add(hCutlastTME);
    hCutChi2NDF = new TH2D("c004hCutChi2NDF_STK",";#chi^{2}/NDF;Candidates", 30, 0, 500, nparType, parTypemin, parTypemax);   
    lout->Add(hCutChi2NDF);
    hCutemScore = new TH2D("c005hCutemScore_STK",";EM Shower Score;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore);
    // Check SCE affect shower score
    hCutemScore_R1 = new TH2D("c005hCutemScore_R1_STK",";EM Shower Score;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore_R1);
    hCutemScore_R2 = new TH2D("c005hCutemScore_R2_STK",";EM Shower Score;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore_R2);
    hCutemScore_R3 = new TH2D("c005hCutemScore_R3_STK",";EM Shower Score;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore_R3);
    hCutemScore_R4 = new TH2D("c005hCutemScore_R4_STK",";EM Shower Score;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore_R4);
    hCutmichelScore = new TH2D("c006hCutmichelScore_STK",";Michel Score;Candidates",50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutmichelScore);
    hCutShowerDist = new TH2D("c007hCutShowerDist_STK",";Shower Distance (cm);Candidates", 31, 0, 93, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerDist);
    hCutShowerIP = new TH2D("c008hCutShowerIP_STK",";Shower Impact Parameter (cm);Candidates", 30, 0, 30, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerIP);
    hCutShowerStartZ = new TH2D("c009hCutShowerStartZ_STK",";Shower Energy (GeV);Candidates", 30, 0, 240, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerStartZ);
    hCutnproton = new TH2D("c010hCutnproton_STK",";Number of Protons;Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnproton);
    hCutnshower = new TH2D("c011hCutnshower_STK",";Number of Showers;Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnshower);
    hCutnpiplus = new TH2D("c012hCutnpiplus_STK",";Number of Pions (#pi^{+});Candidates",  ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnpiplus);
    hCutnmichel = new TH2D("c013hCutnmichel_STK",";Number of Michel Electron;Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnmichel);
    hNhitsCollection = new TH2D("c014NhitsCollection_COMPOSE_LOG",";Number of Hits -coll.;Candidates", 50, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hNhitsCollection);
    hCutnpi0 = new TH2D("c015hCutnpi0_STK",";Number of #pi^{0};Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnpi0);
    hCutKFPass = new TH2D("c015hCutKFPass_STK",";KF Pass;Candidates",  nPass, Passmin, Passmax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutKFPass);

    hShowerStartX = new TH2D("c016hShowerStartX_STK",";StartX (cm);Candidates", 30, -100, 20, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerStartX);
    hShowerStartY = new TH2D("c017hShowerStartY_STK",";StartY (cm);Candidates", 30, 340, 460, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerStartY);
    hShowerStartZ = new TH2D("c018hShowerStartZ_STK",";StartZ (cm);Candidates", 30, 0, 240, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerStartZ);
    hShowernHitsColl = new TH2D("c019hShowernHitsColl_STK",";Number of Hits -coll.;Candidates", 50, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hShowernHitsColl);

    htrackScoreCollection = new TH2D("c020htrackScoreCollection_STK",";Track Score -coll.;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(htrackScoreCollection);
    hemScoreCollection = new TH2D("c021hemScoreCollection_STK",";EM Score -coll.;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hemScoreCollection);
    hmichelScoreCollection = new TH2D("c022hmichelScoreCollection_STK",";Michel Score -coll.;Candidates", 50, 0, 0.3, nparType, parTypemin, parTypemax); 
    lout->Add(hmichelScoreCollection);

    hShowerAngleOffset = new TH2D("c023hShowerAngleOffset_STK",";Shower Angle Offset;Candidates", 30, 0, 180, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerAngleOffset);

    hShowerMichelScore = new TH2D("c024hShowerMichelScore_STK",";Michel Score;Candidates",50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerMichelScore);

    hKFPassRate = new TH1D("k001hKFPassRate","",  nPass, Passmin, Passmax); 
    lout->Add(hKFPassRate); 
    hPi0Total = new TH1D("k002hPi0Total","",  30, 0, 3); 
    lout->Add(hPi0Total);
    hPi0Selected = new TH1D("k003hPi0Selected","",  30, 0, 3);  
    lout->Add(hPi0Selected);

    // Class D - reconstructed TKI variables
    hRecdalphat = new TH2D("d001Recdalphat_STK",";#delta#alpha_{T} (deg);Candidates", 8, 0, 180,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecdalphat);
    hRecdphit = new TH2D("d002Recphit_STK","", 8, 0, 180,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecdphit);
    hRecdpt = new TH2D("d003Recdpt_STK","", 8, 0, 1,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecdpt);
    hRecpn = new TH2D("d004Recpn_STK",";p_{n} (GeV/c);Candidates", 8, 0, 1,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecpn);

    hRecdalphat_truth = new TH2D("d005Recdalphat_STK_TM","", 8, 0, 180,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecdalphat_truth);
    hRecdphit_truth = new TH2D("d006Recphit_STK_TM","", 8, 0, 180,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecdphit_truth);
    hRecdpt_truth = new TH2D("d007Recdpt_STK_TM","", 8, 0, 1,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecdpt_truth);
    hRecpn_truth = new TH2D("d008Recpn_STK_TM","", 8, 0, 1,nevtTKIType, evtTKITypemin, evtTKITypemax); 
    lout->Add(hRecpn_truth);


    hCutemScoreVSStartX = new TH2D("s001hCutemScoreVSStartX_MAP",";StartX (cm);EM Shower Score", 20, -100, 20, 10, 0.5, 1); 
    lout->Add(hCutemScoreVSStartX);
    hCutemScoreVSStartY = new TH2D("s002hCutemScoreVSStartY_MAP",";StartY (cm);EM Shower Score", 20, 340, 460, 10, 0.5, 1); 
    lout->Add(hCutemScoreVSStartY);
    hCutemScoreVSStartZ = new TH2D("s003hCutemScoreVSStartZ_MAP",";StartZ (cm);EM Shower Score", 20, 0, 240, 10, 0.5, 1); 
    lout->Add(hCutemScoreVSStartZ);
    hCutStartXVSStartY = new TH2D("s004hCutStartXVSStartY_MAP",";StartX (cm);StartY (cm)", 15, -100, 20, 15, 340, 460); 
    lout->Add(hCutStartXVSStartY);
    hCutStartZVSStartX = new TH2D("s005hCutStartZVSStartX_MAP",";StartZ (cm);StartX (cm)", 15, 0, 240, 15, -100, 20); 
    lout->Add(hCutStartZVSStartX);
    hCutemScoreVSEnergy = new TH2D("s006hCutemScoreVSEnergy_MAP",";Energy (GeV);EM Shower Score", 20, 0, 1, 10, 0, 1); 
    lout->Add(hCutemScoreVSEnergy);
    hCutemScore_AfterCut = new TH2D("s007hCutemScore_AfterCut_STK",";EM Shower Score;Candidates", 25, 0.5, 1, nparType, parTypemin, parTypemax);
    lout->Add(hCutemScore_AfterCut);

    hAllCutStartXVSStartY = new TH2D("s008hAllCutStartXVSStartY_MAP",";StartX (cm);StartY (cm)", 20, -100, 20, 20, 340, 460); 
    lout->Add(hAllCutStartXVSStartY);
    hAllCutStartZVSStartX = new TH2D("s009hAllCutStartZVSStartX_MAP",";StartZ (cm);StartX (cm)", 20, 0, 240, 20, -100, 20); 
    lout->Add(hAllCutStartZVSStartX);
    hAllCutemScoreVSStartX = new TH2D("s010hAllCutemScoreVSStartX_MAP",";StartX (cm);EM Shower Score", 20, -100, 20, 20, 0, 1); 
    lout->Add(hAllCutemScoreVSStartX);
    hAllCutemScoreVSStartY = new TH2D("s011hAllCutemScoreVSStartY_MAP",";StartY (cm);EM Shower Score", 20, 340, 460, 20, 0, 1); 
    lout->Add(hAllCutemScoreVSStartY);
    hAllCutemScoreVSStartZ = new TH2D("s012hAllCutemScoreVSStartZ_MAP",";StartZ (cm);EM Shower Score", 20, 0, 240, 20, 0, 1); 
    lout->Add(hAllCutemScoreVSStartZ);

    hAllTrackCutStartXVSStartY = new TH2D("s013hAllTrackCutStartXVSStartY_MAP",";StartX (cm);StartY (cm)", 20, -100, 20, 20, 340, 460); 
    lout->Add(hAllTrackCutStartXVSStartY);
    hAllTrackCutStartZVSStartX = new TH2D("s014hAllTrackCutStartZVSStartX_MAP",";StartZ (cm);StartX (cm)", 20, 0, 240, 20, -100, 20); 
    lout->Add(hAllTrackCutStartZVSStartX);

    hBeamEndZ_CutPDG = new TH2D("t001hBeamEndZ_CutPDG_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutPDG);
    hBeamEndZ_CutPandora = new TH2D("t002hBeamEndZ_CutPandora_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutPandora);
    hBeamEndZ_CutCaloSize = new TH2D("t003hBeamEndZ_CutCaloSize_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutCaloSize);
    hBeamEndZ_CutBeamQuality = new TH2D("t004hBeamEndZ_CutBeamQuality_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutBeamQuality);
    hBeamEndZ_CutAPA3 = new TH2D("t005hBeamEndZ_CutAPA3_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutAPA3);
    hBeamEndZ_CutMichelScore = new TH2D("t006hBeamEndZ_CutMichelScore_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutMichelScore);
    hBeamEndZ_CutChi2DOF = new TH2D("t007hBeamEndZ_CutChi2DOF_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutChi2DOF);
    hBeamEndZ_CutBeamScraper = new TH2D("t008hBeamEndZ_CutBeamScraper_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hBeamEndZ_CutBeamScraper);

    hBeamEndZ_Channels = new TH2D("t009hBeamEndZ_Channels_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels);
    hBeamEndZ_Channels_afterEvtTop = new TH2D("t010hBeamEndZ_Channels_afterEvtTop_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_afterEvtTop);
    hBeamEndZ_ChannelsTrueSignal = new TH1D("t011hBeamEndZ_ChannelsTrueSignal",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax); 
    lout->Add(hBeamEndZ_ChannelsTrueSignal);
    
    hBeamEndZ_Channels_BeamPDG = new TH2D("t012hBeamEndZ_Channels_BeamPDG_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamPDG);
    hBeamEndZ_Channels_BeamPandora = new TH2D("t013hBeamEndZ_Channels_BeamPandora_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamPandora);
    hBeamEndZ_Channels_BeamCalo = new TH2D("t014hBeamEndZ_Channels_BeamCalo_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamCalo);
    hBeamEndZ_Channels_BeamQuality = new TH2D("t015hBeamEndZ_Channels_BeamQuality_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamQuality);
    hBeamEndZ_Channels_BeamAPA3 = new TH2D("t016hBeamEndZ_Channels_BeamAPA3_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamAPA3);
    hBeamEndZ_Channels_BeamMichel = new TH2D("t017hBeamEndZ_Channels_BeamMichel_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamMichel);
    hBeamEndZ_Channels_BeamChi2 = new TH2D("t018hBeamEndZ_Channels_BeamChi2_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamChi2);
    hBeamEndZ_Channels_BeamScraper = new TH2D("t019hBeamEndZ_Channels_BeamScraper_STK",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_BeamScraper);
    hBeamEndZ_Channels_TwoShowers = new TH2D("t020hBeamEndZ_Channels_TwoShowers_STK",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_TwoShowers);
    hBeamEndZ_Channels_NoPiPlus = new TH2D("t021hBeamEndZ_Channels_NoPiPlus_STK",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_NoPiPlus);
    hBeamEndZ_Channels_NoMichel = new TH2D("t021hBeamEndZ_Channels_NoMichel_STK",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_NoMichel);
    hBeamEndZ_Channels_OnePi0 = new TH2D("t022hBeamEndZ_Channels_OnePi0_STK",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nchannelType, channelTypemin, channelTypemax); 
    lout->Add(hBeamEndZ_Channels_OnePi0);
    // Truth vars - Eff calculation
    hPiMom_InelasticChannels = new TH1D("t023hPiMom_InelasticChannels",";Momentum (GeV);Candidates", 20, 0, 1); 
    lout->Add(hPiMom_InelasticChannels);
    hBeamEndZ_TrueAvailable = new TH1D("t024hBeamEndZ_TrueAvailable",";True Beam End Z (cm);Candidates", 20, 0, 1); 
    lout->Add(hBeamEndZ_TrueAvailable);
    hBeamEndZ_ChargeExcChannel = new TH1D("t025hBeamEndZ_ChargeExcChannel",";Momentum (GeV);Candidates", 20, 0, 1); 
    lout->Add(hBeamEndZ_ChargeExcChannel);

    hBeamEndZ_XsEvt_BeamPDG = new TH2D("t026hBeamEndZ_XsEvt_BeamPDG_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamPDG);
    hBeamEndZ_XsEvt_BeamPandora = new TH2D("t027hBeamEndZ_XsEvt_BeamPandora_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamPandora);
    hBeamEndZ_XsEvt_BeamCalo = new TH2D("t028hBeamEndZ_XsEvt_BeamCalo_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamCalo);
    hBeamEndZ_XsEvt_BeamQuality = new TH2D("t029hBeamEndZ_XsEvt_BeamQuality_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamQuality);
    hBeamEndZ_XsEvt_BeamAPA3 = new TH2D("t030hBeamEndZ_XsEvt_BeamAPA3_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamAPA3);
    hBeamEndZ_XsEvt_BeamMichel = new TH2D("t031hBeamEndZ_XsEvt_BeamMichel_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamMichel);
    hBeamEndZ_XsEvt_BeamChi2 = new TH2D("t032hBeamEndZ_XsEvt_BeamChi2_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamChi2);
    hBeamEndZ_XsEvt_BeamScraper = new TH2D("t033hBeamEndZ_XsEvt_BeamScraper_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_BeamScraper);
    hBeamEndZ_XsEvt_TwoShowers = new TH2D("t034hBeamEndZ_XsEvt_TwoShowers_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_TwoShowers);
    hBeamEndZ_XsEvt_NoPiPlus = new TH2D("t035hBeamEndZ_XsEvt_NoPiPlus_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_NoPiPlus);
    hBeamEndZ_XsEvt_NoMichel = new TH2D("t036hBeamEndZ_XsEvt_NoMichel_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_NoMichel);
    hBeamEndZ_XsEvt_OnePi0 = new TH2D("t037hBeamEndZ_XsEvt_OnePi0_COMPOSE",";Beam EndZ (cm);Candidates", nEndZ/2, EndZmin, EndZmax, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hBeamEndZ_XsEvt_OnePi0);


    //hBeamInstXY = new TH2D("u001hBeamInstXY_REG",";Beam InstX (cm);Beam InstY (cm)", 40, -50, 0, 40, 405, 450); 
    hBeamInstXY = new TH2D("u001hBeamInstXY_REG",";Beam InstX (cm);Beam InstY (cm)", 40, -55, 5, 40, 390, 435);
    lout->Add(hBeamInstXY);
    //hBeamInstXY_misIDproton = new TH2D("u002hBeamInstXY_misIDproton_REG",";Beam InstX (cm);Beam InstY (cm)", 40, -50, 0, 40, 405, 450); 
    hBeamInstXY_misIDproton = new TH2D("u002hBeamInstXY_misIDproton_REG",";Beam InstX (cm);Beam InstY (cm)", 40, -55, 5, 40, 390, 435);
    lout->Add(hBeamInstXY_misIDproton);
    //hBeamInstXY_misIDpion = new TH2D("u003hBeamInstXY_misIDpion_REG",";Beam InstX (cm);Beam InstY (cm)", 40, -50, 0, 40, 405, 450); 
    hBeamInstXY_misIDpion = new TH2D("u003hBeamInstXY_misIDpion_REG",";Beam InstX (cm);Beam InstY (cm)", 40, -55, 5, 40, 390, 435);
    lout->Add(hBeamInstXY_misIDpion);

    hShowerEnergy_NoCut = new TH2D("r001hShowerEnergy_NoCut_STK",";Shower Energy (GeV);Candidates", 20, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerEnergy_NoCut);
    hShowerEnergy_CutEMscore = new TH2D("r002hShowerEnergy_CutEMscore_STK",";Shower Energy (GeV);Candidates", 20, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerEnergy_CutEMscore);
    hShowerEnergy_CutnHits = new TH2D("r003hShowerEnergy_CutnHits_STK",";Shower Energy (GeV);Candidates", 20, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerEnergy_CutnHits);
    hShowerEnergy_CutStartZ = new TH2D("r004hShowerEnergy_CutStartZ_STK",";Shower Energy (GeV);Candidates", 20, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerEnergy_CutStartZ);
    hShowerEnergy_CutDistance = new TH2D("r005hShowerEnergy_CutDistance_STK",";Shower Energy (GeV);Candidates", 20, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerEnergy_CutDistance);
    hShowerEnergy_CutIP = new TH2D("r006hShowerEnergy_CutIP_STK",";Shower Energy (GeV);Candidates", 20, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerEnergy_CutIP);

    hPi0Energy_NoCut = new TH2D("r007hPi0Energy_NoCut_STK",";Pi0 Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 3.5);
    lout->Add(hPi0Energy_NoCut);
    hPi0Energy_CutMass = new TH2D("r008hPi0Energy_CutMass_STK",";Pi0 Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 3.5);
    lout->Add(hPi0Energy_CutMass);
    hPi0Energy_CutOA = new TH2D("r009hPi0Energy_CutOA_STK",";Pi0 Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 3.5);
    lout->Add(hPi0Energy_CutOA);

    hProtonMom_NoCut = new TH2D("r010hProtonMom_NoCut_STK",";Proton Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hProtonMom_NoCut);
    hProtonMom_CutTrackScore = new TH2D("r011hProtonMom_CutTrackScore_STK",";Proton Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hProtonMom_CutTrackScore);
    hProtonMom_CutnHits = new TH2D("r012hProtonMom_CutnHits_STK",";Proton Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hProtonMom_CutnHits);
    hProtonMom_CutSubPID = new TH2D("r013hProtonMom_CutSubPID_STK",";Proton Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hProtonMom_CutSubPID);

    hPiPlusMom_NoCut = new TH2D("r014hPiPlusMom_NoCut_STK",";PiPlus Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hPiPlusMom_NoCut);
    hPiPlusMom_CutTrackScore = new TH2D("r015hPiPlusMom_CutTrackScore_STK",";PiPlus Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hPiPlusMom_CutTrackScore);
    hPiPlusMom_CutnHits = new TH2D("r016hPiPlusMom_CutnHits_STK",";PiPlus Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hPiPlusMom_CutnHits);
    hPiPlusMom_CutSubPID = new TH2D("r017hPiPlusMom_CutSubPID_STK",";PiPlus Mom. (GeV/c);Candidates", 20, 0, 2, nparType, parTypemin, parTypemax); 
    lout->Add(hPiPlusMom_CutSubPID);

    hRecPiPlusEnergy_OVERLAY_After_EVTXS = new TH2D("i040hRecPiPlusEnergy_COMPOSE_After_EVTXS",";#pi^{+} Kinetic Energy (GeV);Candidates", 20, 0, 1000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusEnergy_OVERLAY_After_EVTXS);

    hRecPi0Energy_OVERLAY_After_EVTXS = new TH2D("i040hRecPi0Energy_COMPOSE_After_EVTXS",";#pi^{0} Kinetic Energy (GeV);Candidates", 20, 0, 1000, 7, -0.5, 6.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After_EVTXS);


    hRecPiPlusInteractingEnergy = new TH2D("i050hRecPiPlusInteractingEnergy_STK",";#pi^{+} Beam Interacting Energy (MeV);Candidates", 40, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInteractingEnergy);
    hRecPiPlusInstMomentum = new TH2D("i051hRecPiPlusInstMomentum_STK",";#pi^{+} Inst. Momentum (MeV);Candidates", 60, 600, 1200, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInstMomentum);
    hRecPiPlusFrontFaceEnergy = new TH2D("i052hRecPiPlusFrontFaceEnergy_STK",";#pi^{+} Front-Face Energy (MeV);Candidates", 70, 500, 1200, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusFrontFaceEnergy);
    hRecPiPlusIncidentEnergy = new TH2D("i053hRecPiPlusIncidentEnergy_STK",";#pi^{+} Incident Energy (MeV);Candidates", 100, 0, 1000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusIncidentEnergy);
    //hRecPiPlusInitialEnergy = new TH2D("i054hRecPiPlusInitialEnergy_STK",";#pi^{+} Initial Energy (MeV);Candidates", 70, 500, 1200, nbeamType, beamTypemin, beamTypemax); 
    //lout->Add(hRecPiPlusInitialEnergy);
    hRecPiPlusInitialEnergy = new TH2D("i054hRecPiPlusInitialEnergy_STK",";#pi^{+} Initial Energy (MeV);Candidates", 40, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInitialEnergy);
    hRecPiPlusInitialEnergyPosZCut = new TH2D("i054hRecPiPlusInitialEnergyPosZCut_STK",";#pi^{+} Initial Energy (MeV);Candidates", 40, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInitialEnergyPosZCut);
    hRecPiPlusIncidentEnergyNew = new TH2D("i055hRecPiPlusIncidentEnergyNew_STK",";#pi^{+} Incident Energy (MeV);Candidates", 200, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusIncidentEnergyNew);
    hRecPiPlusInteractingEnergyPar  = new TH2D("i056hRecPiPlusInteractingEnergyPar_STK",";#pi^{+} Interacting Energy (MeV);Candidates", 40, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInteractingEnergyPar);
    hRecPiPlusInstMomentumNoSmearing = new TH2D("i057hRecPiPlusInstMomentumNoSmearing_STK",";#pi^{+} Inst. Momentum (MeV);Candidates", 60, 1600, 2200, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInstMomentumNoSmearing);
    hRecPiPlusInstEnergyNoSmearing = new TH2D("i058hRecPiPlusInstEnergyNoSmearing_STK",";#pi^{+} Inst. Energy (MeV);Candidates", 60, 1600, 2200, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInstEnergyNoSmearing);
    
    // Raw 
    hRawInstMomentumRaw = new TH2D("i058hRawInstMomentumRaw_STK",";#pi^{+} Inst. Momentum (MeV/c);Candidates", 60, 600, 1200, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRawInstMomentumRaw);
    hRawInstEnergyRaw = new TH2D("i069hRawInstEnergyRaw_STK",";#pi^{+} Inst. Energy (MeV);Candidates", 60, 700, 1300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRawInstEnergyRaw);
    hRawFrontFaceERaw = new TH2D("i060hRawFrontFaceERaw_STK",";#pi^{+} Front Face Energy (MeV);Candidates", 60, 700, 1300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRawFrontFaceERaw);
    hRecPiPlusFrontFaceENoSmearing = new TH2D("i061hRecPiPlusFrontFaceENoSmearing_STK",";#pi^{+} Front Face Energy (MeV);Candidates", 60, 700, 1300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusFrontFaceENoSmearing);
    hRecPiPlusFrontFaceE = new TH2D("i062hRecPiPlusFrontFaceE_STK",";#pi^{+} Front Face Energy (MeV);Candidates", 60, 700, 1300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusFrontFaceE);
    // Tuned case see i056
    hRecPiPlusInteractingEnergyPar_NoTune  = new TH2D("i053hRecPiPlusInteractingEnergyPar_NoTune_STK",";#pi^{+} Interacting Energy (MeV);Candidates", 20, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusInteractingEnergyPar_NoTune);

    hRecPiPlusInteractingEnergyEvt = new TH2D("i076hRecPiPlusInteractingEnergyEvt_COMPOSE_XsEvt",";#pi^{+} Interacting Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiPlusInteractingEnergyEvt);
    hRecPiPlusInteractingEnergyEvt_NoTune = new TH2D("i077hRecPiPlusInteractingEnergyEvt_NoTune_COMPOSE_XsEvt",";#pi^{+} Interacting Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiPlusInteractingEnergyEvt_NoTune);
    
    hRecPiZeroKineticEnergyEvtNoWeight = new TH2D("i086hRecPiZeroKineticEnergyEvtNoWeight_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroKineticEnergyEvtNoWeight);
    hRecPiZeroKineticEnergyEvt = new TH2D("i087hRecPiZeroKineticEnergyEvt_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroKineticEnergyEvt);
    hRecPiZeroSliceKineticEnergyEvt = new TH2D("i088hRecPiZeroSliceKineticEnergyEvt_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroSliceKineticEnergyEvt);
    hRecPiZeroSlice2KineticEnergyEvt = new TH2D("i089hRecPiZeroSlice2KineticEnergyEvt_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroSlice2KineticEnergyEvt);
    hRecPiZeroSlice3KineticEnergyEvt = new TH2D("i086hRecPiZeroSlice3KineticEnergyEvt_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroSlice3KineticEnergyEvt);

    hRecPiZeroSliceKineticEnergyEvtMC = new TH1D("i088hRecPiZeroSliceKineticEnergyEvtMC_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000); 
    lout->Add(hRecPiZeroSliceKineticEnergyEvtMC);
    hRecPiZeroSliceKineticEnergyEvtData = new TH1D("i088hRecPiZeroSliceKineticEnergyEvtData_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000); 
    lout->Add(hRecPiZeroSliceKineticEnergyEvtData);

    hRecPiPlusCaloZ = new TH2D("i060hRecPiPlusCaloZ_STK",";#pi^{+} Calo Z (cm);Candidates", 100, -100, 300, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusCaloZ);
    hRecPiPlusTrackLength = new TH2D("i061hRecPiPlusTrackLength_STK",";#pi^{+} Track Length (cm);Candidates", 100, 0, 2000, nbeamType, beamTypemin, beamTypemax); 
    lout->Add(hRecPiPlusTrackLength);

    hRecPiPlusInteractingEnergyBckSub = new TH2D("i099hRecPiPlusInteractingEnergyBckSub_COMPOSE_XsEvt",";#pi^{+} Interacting Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiPlusInteractingEnergyBckSub);

    hRecPiPlusInteractingEnergyBckSubCheck = new TH2D("i100hRecPiPlusInteractingEnergyBckSubCheck_COMPOSE_XsEvt",";#pi^{+} Interacting Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiPlusInteractingEnergyBckSubCheck);

    // =========== Pi0 KE =========== //
    hRecPiZeroRangeKineticEnergyEvtNoWeight = new TH2D("i300hRecPiZeroRangeKineticEnergyEvtNoWeight_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeKineticEnergyEvtNoWeight);
    hRecPiZeroRangeKineticEnergyEvtOneWeight = new TH2D("i301hRecPiZeroRangeKineticEnergyEvtOneWeight_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeKineticEnergyEvtOneWeight);
    hRecPiZeroRangeKineticEnergyEvtAnotherWeight = new TH2D("i302hRecPiZeroRangeKineticEnergyEvtAnotherWeight_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeKineticEnergyEvtAnotherWeight);
    hRecPiZeroRangeKineticEnergyEvt = new TH2D("i303hRecPiZeroRangeKineticEnergyEvt_COMPOSE_XsEvt",";#pi^{0} Kinetic Energy (MeV);Candidates", 20, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeKineticEnergyEvt);
    hRecPiZeroRangeKineticEnergyEvtNoWeight_LargeBin = new TH2D("i304hRecPiZeroRangeKineticEnergyEvtNoWeight_COMPOSE_XsEvt_LargeBin",";#pi^{0} Kinetic Energy (MeV);Candidates", 10, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeKineticEnergyEvtNoWeight_LargeBin);
    hRecPiZeroRangeKineticEnergyEvtRaw_LargeBin = new TH2D("i305hRecPiZeroRangeKineticEnergyEvtRaw_COMPOSE_XsEvt_LargeBin",";#pi^{0} Kinetic Energy (MeV);Candidates", 10, 0, 2000, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeKineticEnergyEvtRaw_LargeBin);
    // =========== Pi0 Costheta =========== //
    hRecPiZeroRangeCosThetaEvtNoWeight = new TH2D("i400hRecPiZeroRangeCosThetaEvtNoWeight_COMPOSE_XsEvt",";#pi^{0} CosTheta;Candidates", 10, -1, 1, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeCosThetaEvtNoWeight);
    hRecPiZeroRangeCosThetaEvtOneWeight = new TH2D("i401hRecPiZeroRangeCosThetaEvtOneWeight_COMPOSE_XsEvt",";#pi^{0} CosTheta;Candidates", 10, -1, 1, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeCosThetaEvtOneWeight);
    hRecPiZeroRangeCosThetaEvtAnotherWeight = new TH2D("i402hRecPiZeroRangeCosThetaEvtAnotherWeight_COMPOSE_XsEvt",";#pi^{0} CosTheta;Candidates", 10, -1, 1, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeCosThetaEvtAnotherWeight);
    hRecPiZeroRangeCosThetaEvt = new TH2D("i403hRecPiZeroRangeCosThetaEvt_COMPOSE_XsEvt",";#pi^{0} CosTheta;Candidates", 10, -1, 1, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeCosThetaEvt);
    // =========== Pi0 Theta =========== //
    hRecPiZeroRangeThetaEvtNoWeight = new TH2D("i500hRecPiZeroRangeThetaEvtNoWeight_COMPOSE_XsEvt",";#pi^{0} Theta (deg.);Candidates", 10, 0, 180, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeThetaEvtNoWeight);
    hRecPiZeroRangeThetaEvtOneWeight = new TH2D("i501hRecPiZeroRangeThetaEvtOneWeight_COMPOSE_XsEvt",";#pi^{0} Theta (deg.);Candidates", 10, 0, 180, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeThetaEvtOneWeight);
    hRecPiZeroRangeThetaEvtAnotherWeight = new TH2D("i502hRecPiZeroRangeThetaEvtAnotherWeight_COMPOSE_XsEvt",";#pi^{0} Theta (deg.);Candidates", 10, 0, 180, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeThetaEvtAnotherWeight);
    hRecPiZeroRangeThetaEvt = new TH2D("i503hRecPiZeroRangeThetaEvt_COMPOSE_XsEvt",";#pi^{0} Theta (deg.);Candidates", 10, 0, 180, nevtXSType, evtXSTypemin, evtXSTypemax); 
    lout->Add(hRecPiZeroRangeThetaEvt);

    int xsecmin = 0;
    int xsecmax = 2000;
    int xsecbin = 20; //40 not sure

    int xsecthetamin = 0;
    int xsecthetamax = 180;
    int xsecthetabin = 10;

    int xseccosthetamin = -1;
    int xseccosthetamax = 1;
    int xseccosthetabin = 10;

    if(!kMC){
      //hRecoBeamInitialHistData = new TH1D("i019hRecoBeamInitialHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 1000, 0, 1000);
      hRecoBeamInitialHistData = new TH1D("i019hRecoBeamInitialHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hRecoBeamInitialHistData);
      //hRecoBeamInteractingHistData = new TH1D("i020hRecoBeamInteractingHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 1000, 0, 1000);
      hRecoBeamInteractingHistData = new TH1D("i020hRecoBeamInteractingHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hRecoBeamInteractingHistData);
      hRecoInteractingHistData = new TH1D("i018hRecoInteractingHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoInteractingHistData);
      hRecoPi0KEHistData = new TH1D("i030hRecoPi0KEHistData", ";Reco #pi^{0} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoPi0KEHistData);
      hRecoPi0CosThetaHistData = new TH1D("i032hRecoPi0CosThetaHistData", ";#pi^{0} Cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hRecoPi0CosThetaHistData);
      hRecoPi0ThetaHistData = new TH1D("i034hRecoPi0ThetaHistData", ";#pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hRecoPi0ThetaHistData);
    }

    //====================== Truth (MC only)======================//
    if(kMC){
      // Class E - truth beam/FS particles
      hTruthBeamMomentum = new TH1D("e000hTruthBeamMomentum", ";Truth Beam Mom. (GeV/c);Candidates", 20, 0.4, 1.0); 
      lout->Add(hTruthBeamMomentum);
      hTruthBeamType = new TH1I("e001hTruthBeamType",  ";Truth Beam Type;Candidates", 20, -0.5, 19.5); 
      lout->Add(hTruthBeamType);
      hTruthFSParticleNumber = new TH1I("e002hTruthFSParticleNumber", ";Truth FS Particle Number;Candidates", 36, -0.5, 35.5); 
      lout->Add(hTruthFSParticleNumber);
      hTruthFSParticleType = new TH1I("e003hTruthFSParticleType", ";Truth FS Particle Type;Candidates", 23, -0.5, 22.5);
      lout->Add(hTruthFSParticleType);
      hTruthFSPi0Number = new TH1I("e004hTruthFSPi0Number", ";Truth FS #pi^{0} Number;Candidates", 5, 0, 5);
      lout->Add(hTruthFSPi0Number);
      hTruthFSMultiPi0 = new TH1I("e005hTruthFSMultiPi0", ";Truth FS Multi #pi^{0};Candidates", 5, 0, 5);
      lout->Add(hTruthFSMultiPi0);
      hTruthProtonP = new TH1D("e006hTruthProtonP", ";Truth Proton Mom.(GeV/c);Candidates", 20, 0.2, 1.0);
      lout->Add(hTruthProtonP);
      hTruthLeadingProtonP = new TH1D("e006hTruthLeadingProtonP", ";Truth Leading Proton Mom.(GeV/c);Candidates", 20, 0.2, 1.0);
      lout->Add(hTruthLeadingProtonP);
      hTruthSubLeadingProtonP = new TH1D("e007hTruthSubLeadingProtonP", ";Truth SubLeading Proton Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingProtonP);
      hTruthLeadingPiZeroP = new TH1D("e008hTruthLeadingPiZeroP", ";Truth Leading #pi^{0} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPiZeroP);
      hTruthLeadingPiZeroE = new TH1D("e008hTruthLeadingPiZeroE", ";Truth Leading #pi^{0} Energy(GeV);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPiZeroE);
      hTruthSubLeadingPiZeroP = new TH1D("e009hTruthSubLeadingPiZeroP", ";Truth SubLeading #pi^{0} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingPiZeroP);
      hTruthGammaMaxE = new TH1D("e010hTruthGammaMaxE", ";Truth Gamma Max. Energy (MeV);Candidates", 50, 0, 15);
      lout->Add(hTruthGammaMaxE);
      hTruthPi0DecayParticleNumber = new TH1D("e011hTruthPi0DecayParticleNumber", ";Truth #pi^{0} Decay Particle Number;Candidates", 10, -0.5, 9.5);
      lout->Add(hTruthPi0DecayParticleNumber);
      hTruthLeadingPi0GammaP = new TH1D("e012hTruthLeadingPi0GammaP", ";Truth Leading #gamma Energy (GeV);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPi0GammaP);
      hTruthLeadingPi0GammaPBin1 = new TH1D("e012hTruthLeadingPi0GammaPBin1", ";Truth Leading #gamma Energy (GeV);Candidates", 18, 0.1, 1);
      lout->Add(hTruthLeadingPi0GammaPBin1);
      hTruthSubLeadingPi0GammaP = new TH1D("e013hTruthSubLeadingPi0GammaP", ";Truth SubLeading #gamma Energy (GeV);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingPi0GammaP);
      hTruthPi0GammaEnergy = new TH2D("e014hTruthPi0GammaEnergy", ";Truth #gamma Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 2.5);
      lout->Add(hTruthPi0GammaEnergy);
      hTruthRarePi0GammaP = new TH1D("e015hTruthRarePi0GammaP", ";Truth Rare #gamma Energy (GeV);Candidates", 20, 0, 0.8);
      lout->Add(hTruthRarePi0GammaP);
      hTruthRarePi0ElectronP = new TH1D("e016hTruthRarePi0ElectronP", ";Truth Rare #e^{-} Energy (GeV);Candidates", 20, 0, 0.8);
      lout->Add(hTruthRarePi0ElectronP);
      hTruthRarePi0PositronP = new TH1D("e017hTruthRarePi0PositronP", ";Truth Rare #e^{+} Energy (GeV);Candidates", 20, 0, 0.8);
      lout->Add(hTruthRarePi0PositronP);
      hTruthSignal = new TH1I("e018hTruthSignal", ";Truth Signal/Background;Candidates",  nPass, Passmin, Passmax); 
      lout->Add(hTruthSignal); 
      hTruthSignalFSParticleNumber = new TH1I("e019hTruthSignalFSParticleNumber", ";Truth Signal FS Particle Number;Candidates", 36, -0.5, 35.5); 
      lout->Add(hTruthSignalFSParticleNumber);
      hTruthSignalFSParticleType = new TH1I("e020hTruthSignalFSParticleType", ";Truth Signal FSParticle Type;Candidates", 23, -0.5, 22.5);
      lout->Add(hTruthSignalFSParticleType);
      hTruthPi0OA = new TH1D("e021hTruthPi0OA", ";Truth Opening Angle;Candidates", 20, 0, 180);
      lout->Add(hTruthPi0OA);
      hTruthLeadingPi0GammaOA = new TH1D("e022hTruthLeadingPi0GammaOA", ";Truth LD Opening Angle;Candidates", 20, 0, 180);
      lout->Add(hTruthLeadingPi0GammaOA);
      hTruthSubLeadingPi0GammaOA = new TH1D("e023hTruthSubLeadingPi0GammaOA", ";Truth SL Opening Angle;Candidates", 20, 0, 180);
      lout->Add(hTruthSubLeadingPi0GammaOA);

      hTruthLeadingPiZeroGammaDist = new TH1D("e024hTruthLeadingPiZeroGammaDist", ";Truth LD #pi^{0} #gamma Distance (cm);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPiZeroGammaDist);
      hTruthSubLeadingPiZeroGammaDist = new TH1D("e025hTruthSubLeadingPiZeroGammaDist", ";Truth LD #pi^{0} #gamma Distance (cm);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingPiZeroGammaDist);

      // 2D Pi0 daughter photons E and OA
      hTruthPiZeroGammaE1E2 = new TH2D("e026hTruthPiZeroGammaE1E2_REG", ";Ld #gamma Energy (GeV);Sl #gamma Energy (GeV)", 20, 0, 1, 12, 0, 0.6);
      lout->Add(hTruthPiZeroGammaE1E2); 
      hTruthPiZeroGammaE1OA = new TH2D("e027hTruthPiZeroGammaE1OA_REG", ";Ld #gamma Energy (GeV);Opening Angle (deg.)", 20, 0, 1, 20, 0, 180);
      lout->Add(hTruthPiZeroGammaE1OA);

      hTruthPiPlusP = new TH1D("e028hTruthPiPlusP", ";Truth #pi^{+} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthPiPlusP);
      hTruthPiMinusP = new TH1D("e029hTruthPiMinusP", ";Truth #pi^{+} Mom.(GeV/c);Candidates", 48, 0, 1.2);
      lout->Add(hTruthPiMinusP);

      hTruthPDGCode_PionMatched = new TH1D("e030hTruthPDGCode_PionMatched", ";PDG Code;Candidates", 3000, 0, 3000);
      lout->Add(hTruthPDGCode_PionMatched);

      hTruthMatchedBeamSample = new TH1D("e031hTruthMatchedBeamSample", "", nEndZ, EndZmin, EndZmax);
      lout->Add(hTruthMatchedBeamSample);
      hTruthMatchedProtonSample = new TH1D("e032hTruthMatchedProtonSample", "", 20, 0, 2);
      lout->Add(hTruthMatchedProtonSample);
      hTruthMatchedShowerSample = new TH1D("e033hTruthMatchedShowerSample", "", 20, 0, 1);
      lout->Add(hTruthMatchedShowerSample);
      hTruthMatchedPi0Sample = new TH1D("e034hTruthMatchedPi0Sample", "", 20, 0, 1);
      lout->Add(hTruthMatchedPi0Sample); 

      hTruthMatchedChannelCEXSample = new TH1D("e035hTruthMatchedChannelCEXSample", "", nEndZ, EndZmin, EndZmax);
      lout->Add(hTruthMatchedChannelCEXSample);  
      hTruthMatchedXsEvtCEXSample = new TH1D("e036hTruthMatchedXsEvtCEXSample", "", nEndZ, EndZmin, EndZmax);
      lout->Add(hTruthMatchedXsEvtCEXSample);

      hTruthLeadingPiPlusP = new TH1D("e037hTruthLeadingPiPlusP", ";Truth #pi^{+} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPiPlusP);
      hTruthSubLeadingPiPlusP = new TH1D("e038hTruthSubLeadingPiPlusP", ";Truth #pi^{+} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingPiPlusP);

      hTruthMatchedPiPlusSample = new TH1D("e039hTruthMatchedPiPlusSample", "", 20, 0, 2);
      lout->Add(hTruthMatchedPiPlusSample);  

      // Class F - truth TKI related 
      hTruthNproton = new TH1I("f001hTruthNproton","",11, -0.5, 10.5); 
      lout->Add(hTruthNproton);
      hTruthNneutron = new TH1I("f002hTruthNneutron","",11, -0.5, 10.5); 
      lout->Add(hTruthNneutron);
      hTruthNPiZero = new TH1I("f003hTruthNPiZero","",11, -0.5, 10.5); 
      lout->Add(hTruthNPiZero);
      hTruthMomIniPi = new TH1D("f004hTruthMomIniPi","", 30, 0.2, 1.5); 
      lout->Add(hTruthMomIniPi);
      hTruthThetaIniPi = new TH1D("f005hTruthThetaIniPi","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPi);
      hTruthMomFinPi = new TH1D("f006hTruthMomFinPi","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPi);
      hTruthThetaFinPi = new TH1D("f007hTruthThetaFinPi","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPi);
      hTruthMomFinProton = new TH1D("f008hTruthMomFinProton","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProton);
      hTruthThetaFinProton = new TH1D("f009hTruthThetaFinProton","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProton);
      hTruthMomFin2Proton = new TH1D("f010hTruthMomFin2Proton","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2Proton);
      hTruthDalphat = new TH1D("f011hTruthDalphat","", 18, 0, 180); 
      lout->Add(hTruthDalphat);
      hTruthDphit = new TH1D("f012hTruthDphit","", 18, 0, 180); 
      lout->Add(hTruthDphit);
      hTruthDpt = new TH1D("f013hTruthDpt","", 30, 0, 1.2); 
      lout->Add(hTruthDpt);
      hTruthPn = new TH1D("f014hTruthPn","", 30, 0, 1.2); 
      lout->Add(hTruthPn);

      stkTruthDalphat = new THStack("f100stkTruthDalphat",";#delta#alpha_{T} (deg);Candidates"); 
      lout->Add(stkTruthDalphat);
      hTruthDalphat1p0n = new TH1D("f100hTruthDalphat1p0n",";#delta#alpha_{T} (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphat1p0n);
      hTruthDalphatNp0n = new TH1D("f100hTruthDalphatNp0n",";#delta#alpha_{T} (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphatNp0n);
      hTruthDalphat1pMn = new TH1D("f100hTruthDalphat1pMn",";#delta#alpha_{T} (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphat1pMn);
      hTruthDalphatNpMn = new TH1D("f100hTruthDalphatNpMn",";#delta#alpha_{T} (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphatNpMn);
      stkTruthDphit = new THStack("f101stkTruthDphit",""); 
      lout->Add(stkTruthDphit);
      hTruthDphit1p0n = new TH1D("f101hTruthDphit1p0n","", 18, 0, 180); 
      lout->Add(hTruthDphit1p0n);
      hTruthDphitNp0n = new TH1D("f101hTruthDphitNp0n","", 18, 0, 180); 
      lout->Add(hTruthDphitNp0n);
      hTruthDphit1pMn = new TH1D("f101hTruthDphit1pMn","", 18, 0, 180); 
      lout->Add(hTruthDphit1pMn);
      hTruthDphitNpMn = new TH1D("f101hTruthDphitNpMn","", 18, 0, 180); 
      lout->Add(hTruthDphitNpMn);
      stkTruthDpt = new THStack("f102stkTruthDpt",""); 
      lout->Add(stkTruthDpt);
      hTruthDpt1p0n = new TH1D("f102hTruthDpt1p0n","", 30, 0, 1.2); 
      lout->Add(hTruthDpt1p0n);
      hTruthDptNp0n = new TH1D("f102hTruthDptNp0n","", 30, 0, 1.2); 
      lout->Add(hTruthDptNp0n);
      hTruthDpt1pMn = new TH1D("f102hTruthDpt1pMn","", 30, 0, 1.2); 
      lout->Add(hTruthDpt1pMn);
      hTruthDptNpMn = new TH1D("f102hTruthDptNpMn","", 30, 0, 1.2); 
      lout->Add(hTruthDptNpMn);
      stkTruthPn = new THStack("f103stkTruthPn",";p_n (GeV/c);Candidates"); 
      lout->Add(stkTruthPn);
      hTruthPn1p0n = new TH1D("f103hTruthPn1p0n",";p_n (GeV/c);Candidates", 24, 0, 1.2); 
      lout->Add(hTruthPn1p0n);
      hTruthPnNp0n = new TH1D("f103hTruthPnNp0n",";p_n (GeV/c);Candidates", 24, 0, 1.2); 
      lout->Add(hTruthPnNp0n);
      hTruthPn1pMn = new TH1D("f103hTruthPn1pMn",";p_n (GeV/c);Candidates", 24, 0, 1.2); 
      lout->Add(hTruthPn1pMn);
      hTruthPnNpMn = new TH1D("f103hTruthPnNpMn",";p_n (GeV/c);Candidates", 24, 0, 1.2); 
      lout->Add(hTruthPnNpMn);
      stkTruthMomIniPi = new THStack("f200stkTruthMomIniPi",""); 
      lout->Add(stkTruthMomIniPi);
      hTruthMomIniPi1p0n = new TH1D("f200hTruthMomIniPi1p0n","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPi1p0n);
      hTruthMomIniPiNp0n = new TH1D("f200hTruthMomIniPiNp0n","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPiNp0n);
      hTruthMomIniPi1pMn = new TH1D("f200hTruthMomIniPi1pMn","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPi1pMn);
      hTruthMomIniPiNpMn = new TH1D("f200hTruthMomIniPiNpMn","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPiNpMn);
      stkTruthThetaIniPi = new THStack("f201stkTruthThetaIniPi",""); 
      lout->Add(stkTruthThetaIniPi);
      hTruthThetaIniPi1p0n = new TH1D("f201hTruthThetaIniPi1p0n","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPi1p0n);
      hTruthThetaIniPiNp0n = new TH1D("f201hTruthThetaIniPiNp0n","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPiNp0n);
      hTruthThetaIniPi1pMn = new TH1D("f201hTruthThetaIniPi1pMn","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPi1pMn);
      hTruthThetaIniPiNpMn = new TH1D("f201hTruthThetaIniPiNpMn","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPiNpMn);
      stkTruthMomFinPi = new THStack("f202stkTruthMomFinPi",""); 
      lout->Add(stkTruthMomFinPi);
      hTruthMomFinPi1p0n = new TH1D("f202hTruthMomFinPi1p0n","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPi1p0n);
      hTruthMomFinPiNp0n = new TH1D("f202hTruthMomFinPiNp0n","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPiNp0n);
      hTruthMomFinPi1pMn = new TH1D("f202hTruthMomFinPi1pMn","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPi1pMn);
      hTruthMomFinPiNpMn = new TH1D("f202hTruthMomFinPiNpMn","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPiNpMn);
      stkTruthThetaFinPi = new THStack("f203stkTruthThetaFinPi",""); 
      lout->Add(stkTruthThetaFinPi);
      hTruthThetaFinPi1p0n = new TH1D("f203hTruthThetaFinPi1p0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPi1p0n);
      hTruthThetaFinPiNp0n = new TH1D("f203hTruthThetaFinPiNp0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPiNp0n);
      hTruthThetaFinPi1pMn = new TH1D("f203hTruthThetaFinPi1pMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPi1pMn);
      hTruthThetaFinPiNpMn = new TH1D("f203hTruthThetaFinPiNpMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPiNpMn);
      stkTruthMomFinProton = new THStack("f204stkTruthMomFinProton",""); 
      lout->Add(stkTruthMomFinProton);
      hTruthMomFinProton1p0n = new TH1D("f204hTruthMomFinProton1p0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProton1p0n);
      hTruthMomFinProtonNp0n = new TH1D("f204hTruthMomFinProtonNp0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProtonNp0n);
      hTruthMomFinProton1pMn = new TH1D("f204hTruthMomFinProton1pMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProton1pMn);
      hTruthMomFinProtonNpMn = new TH1D("f204hTruthMomFinProtonNpMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProtonNpMn);
      stkTruthThetaFinProton = new THStack("f205stkTruthThetaFinProton",""); 
      lout->Add(stkTruthThetaFinProton);
      hTruthThetaFinProton1p0n = new TH1D("f205hTruthThetaFinProton1p0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProton1p0n);
      hTruthThetaFinProtonNp0n = new TH1D("f205hTruthThetaFinProtonNp0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProtonNp0n);
      hTruthThetaFinProton1pMn = new TH1D("f205hTruthThetaFinProton1pMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProton1pMn);
      hTruthThetaFinProtonNpMn = new TH1D("f205hTruthThetaFinProtonNpMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProtonNpMn);
      stkTruthMomFin2Proton = new THStack("f206stkTruthMomFin2Proton",""); 
      lout->Add(stkTruthMomFin2Proton);
      hTruthMomFin2Proton1p0n = new TH1D("f206hTruthMomFin2Proton1p0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2Proton1p0n);
      hTruthMomFin2ProtonNp0n = new TH1D("f206hTruthMomFin2ProtonNp0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2ProtonNp0n);
      hTruthMomFin2Proton1pMn = new TH1D("f206hTruthMomFin2Proton1pMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2Proton1pMn);
      hTruthMomFin2ProtonNpMn = new TH1D("f206hTruthMomFin2ProtonNpMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2ProtonNpMn); 

      // Class G - Resolution related (truth/reco - 1)
      hBeamTrackLengthRes = new TH2D("g000hBeamTrackLengthFitRes_RES",";Rec. Beam Length (cm);Rec. - Truth / Rec.Beam Length", 80 , 0, 300, 30, -0.5, 0.5);
      lout->Add(hBeamTrackLengthRes);
      hBeamPhiRes = new TH2D("g000hBeamPhiFitRes_RES",";Truth Beam #phi (deg);Rec. - Truth (deg)", 80 , -180, 0, 30, -20, 20);
      lout->Add(hBeamPhiRes);
      hBeamThetaRes = new TH2D("g001hBeamThetaFitRes_RES",";Truth Beam #theta (deg);Rec. - Truth (deg)", 80 , 0, 60, 30, -20, 20);
      lout->Add(hBeamThetaRes);
      hBeamMomentumRes  = new TH2D("g002hBeamMomentumFitRes_RES",";Truth Beam Momentum (GeV/c);Rec./Truth - 1", 50, 0.2, 1.6, 30, -0.5, 0.5);
      lout->Add(hBeamMomentumRes);
      hProtonThetaRes = new TH2D("g003hProtonThetaFitRes_RES",";Truth Proton #theta (deg);Rec. - Truth (deg)", 40, 0, 180, 30, -20, 30); 
      lout->Add(hProtonThetaRes);
      hProtonThetaLabRawRes = new TH2D("g004hProtonThetaLabFitRawRes_RES",";Truth Proton #theta Lab (deg);Rec. - Truth (deg)", 40, 0, 180, 30, -20, 30); 
      lout->Add(hProtonThetaLabRawRes);
      hProtonPhiLabRawRes = new TH2D("g005hProtonPhiLabFitRawRes_RES",";Truth Proton #phi Lab (deg);Rec. - Truth (deg)", 80, -180, 180, 30, -20, 30); 
      lout->Add(hProtonPhiLabRawRes);
      hProtonThetaLabRes = new TH2D("g044hProtonThetaLabFitRes_RES",";Truth Proton #theta Lab (deg);Rec. - Truth (deg)", 40, 0, 180, 30, -20, 30); 
      lout->Add(hProtonThetaLabRes);
      hProtonPhiLabRes = new TH2D("g055hProtonPhiLabFitRes_RES",";Truth Proton #phi Lab (deg);Rec. - Truth (deg)", 80, -180, 180, 30, -20, 30); 
      lout->Add(hProtonPhiLabRes);

      hProtonMomentumRes = new TH2D("g006hProtonMomentumFitRes_RES",";Truth Proton Momentum (GeV/c);Rec./Truth - 1",40, 0.3, 1.2, 20, -0.2, 0.2); 
      lout->Add(hProtonMomentumRes);
      hProtonMomentumRawRes = new TH2D("g007hProtonMomentumFitRes_RES_RAW",";Truth Proton Momentum (GeV/c);Rec./Truth - 1",40, 0.3, 1.2, 20, -0.2, 0.2); 
      lout->Add(hProtonMomentumRawRes);
      hPiPlusThetaRes = new TH2D("g008hPiPlusThetaFitRes_RES",";Truth #pi^{+} #theta (deg);Rec. - Truth (deg)", 15, 0, 180, 25, -20, 30); 
      lout->Add(hPiPlusThetaRes);
      hPiPlusMomentumRes = new TH2D("g009hPiPlusMomentumFitRes_RES",";Truth #pi^{+} Momentum (GeV/c);Rec./Truth - 1", 24, 0, 1.2, 20, -1.0, 1.0); 
      lout->Add(hPiPlusMomentumRes);
      hShowerThetaRes = new TH2D("g010hShowerThetaFitRes_RES",";Truth Photon #theta (deg);Rec. - Truth (deg)", 15, 0, 180, 50, -20, 30);
      lout->Add(hShowerThetaRes);
      hShowerThetaResRaw = new TH2D("g011hShowerThetaFitRes_RES_RAW",";Truth Photon #theta Raw (deg);Rec. - Truth (deg)", 15, 0, 180, 50, -20, 30);
      lout->Add(hShowerThetaResRaw);
      hShowerEnergyRes = new TH2D("g012hShowerEnergyFitRes_RES",";Truth Photon Energy (GeV);Rec./Truth - 1", 20, 0, 0.8, 50, -1.1, 1.1);
      lout->Add(hShowerEnergyRes);
      hShowerEnergyResRaw = new TH2D("g013hShowerEnergyFitRawRes_RES_RAW",";Truth Photon Energy(GeV);Rec./Truth - 1", 20, 0, 0.8, 50, -1.1, 1.1);
      lout->Add(hShowerEnergyResRaw);
      hShowerPhiRes = new TH2D("g014hShowerPhiFitRes_RES",";Truth Photon #phi (deg);Rec. - Truth (deg)", 15, -180, 180, 50, -20, 30);
      lout->Add(hShowerPhiRes);
      hShowerPhiResRaw = new TH2D("g015hShowerPhiFitRawRes_RES_RAW",";Truth Photon #phi Raw (deg);Rec. - Truth (deg)", 15, -180, 180, 50, -20, 30);
      lout->Add(hShowerPhiResRaw);

      hLeadingPhotonAngleRes = new TH2D("g016hLeadingPhotonAngleRes_RES",";Truth Leading Photon Energy (GeV);Photon/Shower Angle (deg)", 20, 0, 1, 50, 0, 50);
      lout->Add(hLeadingPhotonAngleRes);

      hSubLeadingPhotonAngleRes = new TH2D("g017hSubLeadingPhotonAngleRes_RES",";Truth SubLeading Photon Energy (GeV);Photon/Shower Angle (deg)", 20, 0, 0.5, 50, 0, 50);
      lout->Add(hSubLeadingPhotonAngleRes);

      hOpeningAngleRes = new TH2D("g018hOpeningAngleRes_RES",";Truth #pi^{0} OA (deg.);Photon/Shower Angle (deg)", 20, 0, 1, 50, -50, 50);
      lout->Add(hOpeningAngleRes);


      
      // Class H - Rec VS truth 
      hProtonMomentumRawRecVSTruth_REG  = new TH2D("h001hProtonMomentumRawRecVSTruth_REG_DIAG",";Rec. Proton Raw Mom. (GeV/c);Truth Proton Mom. (GeV/c)",20, 0.2, 1.2, 20, 0.2, 1.2);
      lout->Add(hProtonMomentumRawRecVSTruth_REG);
      hProtonMomentumRecVSTruth_REG  = new TH2D("h002hProtonMomentumRecVSTruth_REG_DIAG",";Rec. Proton Mom. (GeV/c);Truth Proton Mom. (GeV/c)",20, 0.2, 1.2, 20, 0.2, 1.2);
      lout->Add(hProtonMomentumRecVSTruth_REG);
      hProtonTransverseMomentumRawRecVSTruth_REG  = new TH2D("h003hProtonTransverseMomentumRawRecVSTruth_REG_DIAG","",20, 0, 1, 20, 0, 1);
      lout->Add(hProtonTransverseMomentumRawRecVSTruth_REG);
      hProtonTransverseMomentumRecVSTruth_REG  = new TH2D("h004hProtonTransverseMomentumRecVSTruth_REG_DIAG","",20, 0, 1, 20, 0, 1);
      lout->Add(hProtonTransverseMomentumRecVSTruth_REG);
      hShowerEnergyRecVSTruth_REG = new TH2D("h005hShowerEnergyRecVSTruth_REG_DIAG","", 20, 0, 1, 20, 0, 1);
      lout->Add(hShowerEnergyRecVSTruth_REG);
      hShowerEnergyRawRecVSTruth_REG = new TH2D("h006hShowerEnergyRawRecVSTruth_REG_DIAG","", 20, 0, 1, 20, 0, 1);
      lout->Add(hShowerEnergyRawRecVSTruth_REG);
      hShowerEnergyResVSnHits_REG = new TH2D("h007ShowerEnergyResVSnHits_REG","", 50, 0, 1000, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResVSnHits_REG); 
      hShowerEnergyResVSIP_REG = new TH2D("h008ShowerEnergyResVSIP_REG","", 20, 0, 30, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResVSIP_REG);
      const double ProtonBin[] = {0.45, 0.48, 0.51, 0.54, 0.57, 0.6, 0.63, 0.68, 0.73, 0.78, 0.85, 0.95, 1.1};
      const double ShowerBin[] = {0, 0.05, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.25, 0.3, 0.36, 0.48, 0.7};
      const double ThetaBin[] = {0, 15, 23, 30, 37, 46, 55, 65, 75, 90, 105, 120, 180};
      hProtonMomentumRecVSTruth_REG_AfterCor = new TH2D("h009hProtonMomentumRecVSTruth_REG_AfterCor","",sizeof(ProtonBin)/sizeof(double)-1, ProtonBin, sizeof(ProtonBin)/sizeof(double)-1, ProtonBin);
      lout->Add(hProtonMomentumRecVSTruth_REG_AfterCor);
      hShowerEnergyRecVSTruth_REG_AfterCor  = new TH2D("h010hShowerEnergyRecVSTruth_REG_AfterCor","",sizeof(ShowerBin)/sizeof(double)-1, ShowerBin, sizeof(ShowerBin)/sizeof(double)-1, ShowerBin);
      lout->Add(hShowerEnergyRecVSTruth_REG_AfterCor);
      hShowerThetaRecVSTruth_REG_AfterCor = new TH2D("h011hShowerThetaRecVSTruth_REG_AfterCor","",sizeof(ThetaBin)/sizeof(double)-1, ThetaBin, sizeof(ThetaBin)/sizeof(double)-1, ThetaBin);
      lout->Add(hShowerThetaRecVSTruth_REG_AfterCor);
      hShowerPhiRecVSTruth_REG_AfterCor = new TH2D("h012hShowerPhiRecVSTruth_REG_AfterCor","",12,-180,180,12,-180,180);
      lout->Add(hShowerPhiRecVSTruth_REG_AfterCor);
      hProtonThetaRecVSTruth_REG_AfterCor = new TH2D("h013hProtonThetaRecVSTruth_REG_AfterCor","",sizeof(ThetaBin)/sizeof(double)-1, ThetaBin, sizeof(ThetaBin)/sizeof(double)-1, ThetaBin);
      lout->Add(hProtonThetaRecVSTruth_REG_AfterCor);
      hProtonPhiRecVSTruth_REG_AfterCor = new TH2D("h014hProtonPhiRecVSTruth_REG_AfterCor","",12,-180,180,12,-180,180);
      lout->Add(hProtonPhiRecVSTruth_REG_AfterCor);

      // Class L Energy/momentum correction
      hProtonMomentumRecVSTruth_REG_Correction  = new TH2D("l001hProtonMomentumRecVSTruth_REG_Correction_protonPAll","",sizeof(ProtonBin)/sizeof(double)-1, ProtonBin, 100, -0.3, 0.3);
      lout->Add(hProtonMomentumRecVSTruth_REG_Correction);
      hProtonTransverseMomentumRecVSTruth_REG_Correction  = new TH2D("l002hProtonTransverseMomentumRecVSTruth_REG_Correction_protonPT","",12, 0.2, 0.8, 100, -0.6, 0.6);
      lout->Add(hProtonTransverseMomentumRecVSTruth_REG_Correction);
      hShowerEnergyRecVSTruth_REG_Correction  = new TH2D("l003hShowerEnergyRecVSTruth_REG_Correction_showerE","",sizeof(ShowerBin)/sizeof(double)-1, ShowerBin, 50, -1, 1);
      lout->Add(hShowerEnergyRecVSTruth_REG_Correction);  
      hMeanPMom = new TH1D("l004hProtonMomentum_CorMean","", sizeof(ProtonBin)/sizeof(double)-1, ProtonBin);
      lout->Add(hMeanPMom);
      hMeanPMomT = new TH1D("l005hProtonMomentumT_CorMean","",12, 0.2, 0.8);
      lout->Add(hMeanPMomT);
      hMeanShowerE = new TH1D("l006hShowerEnergy_CorMean",";Mean Rec./Truth -1;Rec. Shower Energy (GeV)", sizeof(ShowerBin)/sizeof(double)-1, ShowerBin);
      lout->Add(hMeanShowerE);

      const double LdShowerBin[] = {0, 0.09, 0.12, 0.14, 0.17, 0.2, 0.23, 0.27, 0.32, 0.38, 0.45, 0.58, 0.9};
      const double SlShowerBin[] = {0, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.19, 0.25, 0.9};
      hLeadingShowerEnergyRecVSTruth_REG_Correction  = new TH2D("l007hLeadingShowerEnergyRecVSTruth_REG_Correction_showerLDE","",sizeof(LdShowerBin)/sizeof(double)-1, LdShowerBin, 50, -1, 1);
      lout->Add(hLeadingShowerEnergyRecVSTruth_REG_Correction);
      hSubLeadingShowerEnergyRecVSTruth_REG_Correction  = new TH2D("l008hSubLeadingShowerEnergyRecVSTruth_REG_Correction_showerSLE","",sizeof(SlShowerBin)/sizeof(double)-1, SlShowerBin, 50, -1, 1);
      lout->Add(hSubLeadingShowerEnergyRecVSTruth_REG_Correction);
      hOpeningAngleRecVSTruth_REG_Correction  = new TH2D("l009hOpeningAngleRecVSTruth_REG_Correction","",12, 0, 180, 50, -50, 50);
      lout->Add(hOpeningAngleRecVSTruth_REG_Correction);
      hShowerThetaRecVSTruth_REG_Correction  = new TH2D("l010hShowerThetaRecVSTruth_REG_Correction_showerTheta","",sizeof(ThetaBin)/sizeof(double)-1, ThetaBin, 50, -50, 50);
      lout->Add(hShowerThetaRecVSTruth_REG_Correction);  
      hMeanShowerTheta = new TH1D("l011hMeanShowerTheta_CorMean","", sizeof(ThetaBin)/sizeof(double)-1, ThetaBin);
      lout->Add(hMeanShowerTheta);

      hShowerPhiRecVSTruth_REG_Correction  = new TH2D("l012hShowerPhiRecVSTruth_REG_Correction_showerPhi","",12, -180, 180, 50, -50, 50);
      lout->Add(hShowerPhiRecVSTruth_REG_Correction);  
      hMeanShowerPhi = new TH1D("l013hMeanShowerPhi_CorMean","", 12, -180, 180);
      lout->Add(hMeanShowerPhi);

      hProtonThetaRecVSTruth_REG_Correction  = new TH2D("l014hProtonThetaRecVSTruth_REG_Correction_ProtonTheta","",sizeof(ThetaBin)/sizeof(double)-1, ThetaBin, 50, -50, 50);
      lout->Add(hProtonThetaRecVSTruth_REG_Correction);  
      hMeanProtonTheta = new TH1D("l015hMeanProtonTheta_CorMean","", sizeof(ThetaBin)/sizeof(double)-1, ThetaBin);
      lout->Add(hMeanProtonTheta);

      hProtonPhiRecVSTruth_REG_Correction  = new TH2D("l016hProtonPhiRecVSTruth_REG_Correction_ProtonPhi","",12, -180, 180, 50, -50, 50);
      lout->Add(hProtonPhiRecVSTruth_REG_Correction);  
      hMeanProtonPhi = new TH1D("l017hMeanProtonPhi_CorMean","", 12, -180, 180);
      lout->Add(hMeanProtonPhi);

      hMeanLDShowerE = new TH1D("l018hShowerLDEnergy_CorMean","", sizeof(LdShowerBin)/sizeof(double)-1, LdShowerBin);
      lout->Add(hMeanLDShowerE);
      hMeanSLShowerE = new TH1D("l019hShowerSLEnergy_CorMean","", sizeof(SlShowerBin)/sizeof(double)-1, SlShowerBin);
      lout->Add(hMeanSLShowerE);


      // Beam correction
      const double BeamMomBin[] = {0.5, 0.62, 0.72, 0.78, 0.81, 0.84, 0.87, 0.9, 0.94, 0.98, 1.03, 1.2};
      const double BeamThetaBin[] = {0, 10, 13, 15, 16, 17, 18, 19, 20, 22, 25, 35, 80};
      const double BeamPhiBin[] = {-180, -160, -150, -145, -140, -135, -130, -124, -119, -100, -50};

      hBeamThetaRecVSTruth_REG_Correction  = new TH2D("l020hBeamThetaRecVSTruth_REG_Correction_BeamTheta","",sizeof(BeamThetaBin)/sizeof(double)-1, BeamThetaBin, 50, -50, 50);
      lout->Add(hBeamThetaRecVSTruth_REG_Correction);  
      hMeanBeamTheta = new TH1D("l021hMeanBeamTheta_CorMean","", sizeof(BeamThetaBin)/sizeof(double)-1, BeamThetaBin);
      lout->Add(hMeanBeamTheta);

      hBeamPhiRecVSTruth_REG_Correction  = new TH2D("l022hBeamPhiRecVSTruth_REG_Correction_BeamPhi","",sizeof(BeamPhiBin)/sizeof(double)-1, BeamPhiBin, 50, -50, 50);
      lout->Add(hBeamPhiRecVSTruth_REG_Correction);  
      hMeanBeamPhi = new TH1D("l023hMeanBeamPhi_CorMean","", sizeof(BeamPhiBin)/sizeof(double)-1, BeamPhiBin);
      lout->Add(hMeanBeamPhi);

      hBeamMomentumRecVSTruth_REG_Correction  = new TH2D("l024hBeamMomentumRecVSTruth_REG_Correction_BeamPAll","",sizeof(BeamMomBin)/sizeof(double)-1, BeamMomBin, 100, -0.3, 0.3);
      lout->Add(hBeamMomentumRecVSTruth_REG_Correction);  
      hMeanBeamPMom = new TH1D("l025hBeamMomentum_CorMean","", sizeof(BeamMomBin)/sizeof(double)-1, BeamMomBin);
      lout->Add(hMeanBeamPMom);


      // Treat Leading and SubLeading showers differently
      
      // Class M - Pi0 shower related (rec/trurh matched)
      hLeadingShowerEnergyRes = new TH2D("m001hLeadingShowerEnergy_RES","", 20, 0, 0.8, 20, -1.1, 1.1);
      lout->Add(hLeadingShowerEnergyRes);
      hSubLeadingShowerEnergyRes = new TH2D("m002hSubLeadingShowerEnergy_RES","", 20, 0, 0.5, 20, -1.1, 1.1);
      lout->Add(hSubLeadingShowerEnergyRes);  
      hLeadingShowerEnergyResRaw = new TH2D("m003hLeadingShowerEnergy_RES_RAW","", 20, 0, 0.8, 20, -1.1, 1.1);
      lout->Add(hLeadingShowerEnergyResRaw);
      hSubLeadingShowerEnergyResRaw = new TH2D("m004hSubLeadingShowerEnergy_RES_RAW","", 20, 0, 0.5, 20, -1.1, 1.1);
      lout->Add(hSubLeadingShowerEnergyResRaw);
      hShowerOpenAngleRes = new TH2D("m005hShowerOpenAngle_RES","", 20, 0, 180, 20, -50, 50);
      lout->Add(hShowerOpenAngleRes); 
      hPi0MomentumResRaw = new TH2D("m006Pi0MomentumFitRes_RES_RAW",";Truth #pi^{0} Energy Raw (GeV);Rec./Truth - 1", 25, 0, 1, 20, -1, 1);
      lout->Add(hPi0MomentumResRaw);
      hPi0MomentumRes = new TH2D("m007hPi0MomentumFitRes_RES",";Truth #pi^{0} Energy (GeV);Rec./Truth - 1", 25, 0, 1, 20, -1, 1);
      lout->Add(hPi0MomentumRes);
      hPi0MassResRaw = new TH2D("m008Pi0Mass_RES_RAW","", 20, 0, 0.5, 20, -0.5, 0.5);
      lout->Add(hPi0MassResRaw);
      hPi0MassRes = new TH2D("m009Pi0Mass_RES","", 20, 0, 0.5, 20, -0.5, 0.5); 
      lout->Add(hPi0MassRes);
      hPi0MomentumResFitKF = new TH2D("m099Pi0MomentumFitRes_RES_FITKF",";Truth #pi^{0} Energy (GeV);Rec./Truth - 1", 25, 0, 1, 20, -1, 1);
      lout->Add(hPi0MomentumResFitKF);
      hPi0MomentumResFitKF_M2 = new TH2D("m100Pi0MomentumFitRes_RES_FITKF_M2",";Truth #pi^{0} Energy (GeV);Rec./Truth - 1", 25, 0, 1, 20, -1, 1);
      lout->Add(hPi0MomentumResFitKF_M2);

      hMatchedTruthPi0Energy = new TH1D("m010hMatchedTruthPi0Energy",";Truth Energy (GeV);Efficiency", 20, 0, 1);
      lout->Add(hMatchedTruthPi0Energy);
      hMatchedTruthldShowerTheta = new TH1D("m011hMatchedTruthldShowerTheta","", 30, 0, 180);
      lout->Add(hMatchedTruthldShowerTheta);
      hMatchedTruthslShowerTheta = new TH1D("m012hMatchedTruthslShowerTheta","", 30, 0, 180);
      lout->Add(hMatchedTruthslShowerTheta);
      hMatchedTruthPi0Mass = new TH1D("m013hMatchedTruthPi0Mass","", 20, 0, 0.4);
      lout->Add(hMatchedTruthPi0Mass);

      hMatchedTruthProtonMomentum = new TH1D("m010hMatchedTruthProtonMomentum",";Truth Proton Mom. (GeV/c);Efficiency", 20, 0.2, 1.0);
      lout->Add(hMatchedTruthProtonMomentum);
      hMatchedTruthPiPlusMomentum = new TH1D("m010hMatchedTruthPiPlusMomentum",";Truth #pi^{+} Mom. (GeV/c);Efficiency", 20, 0, 1);
      lout->Add(hMatchedTruthPiPlusMomentum);
      
      hMatchedTruthLeadingShowerEnergy = new TH1D("m010hMatchedTruthLeadingShowerEnergy",";Truth Leading #gamma Energy (GeV);Efficiency", 18, 0.1, 1);
      lout->Add(hMatchedTruthLeadingShowerEnergy);
      hMatchedTruthLeadingProtonMomentum = new TH1D("m010hMatchedTruthLeadingProtonMomentum",";Truth Leading Proton Mom. (GeV/c);Efficiency", 20, 0.2, 1.0);
      lout->Add(hMatchedTruthLeadingProtonMomentum);
      hMatchedTruthLeadingPiPlusMomentum = new TH1D("m010hMatchedTruthLeadingPiPlusMomentum",";Truth Leading #pi^{+} Mom. (GeV/c);Efficiency", 20, 0, 1);
      lout->Add(hMatchedTruthLeadingPiPlusMomentum);

      hPi0ThetaResFit = new TH2D("m014hPi0ThetaFitRes_RES",";Truth #pi^{0} #theta Angle (deg);Rec.- Truth (deg)", 30, 0, 180, 40, -50, 50);
      lout->Add(hPi0ThetaResFit);
      hPi0PhiResFit = new TH2D("m015hPi0PhiFitRes_RES",";Truth #pi^{0} #phi Angle (deg);Rec.- Truth (deg)", 30, -180, 180, 40, -50, 50);
      lout->Add(hPi0PhiResFit);
      hPi0MomRes = new TH2D("m016hPi0MomFitRes_RES",";Truth #pi^{0} Energy (GeV);Rec.(post fit)/Truth - 1", 25, 0, 1, 20, -1, 1);
      lout->Add(hPi0MomRes);
      hPi0MomPreRes = new TH2D("m017hPi0MomPreFitRes_RES",";Truth #pi^{0} Energy (GeV);Rec.(pre fit)/Truth - 1", 25, 0, 1, 20, -1, 1);
      lout->Add(hPi0MomPreRes);

      hMatchedTruthldShowerEnergy = new TH1D("m011hMatchedTruthldShowerEnergy","", 20, 0, 1);
      lout->Add(hMatchedTruthldShowerEnergy);
      hMatchedTruthslShowerEnergy = new TH1D("m012hMatchedTruthslShowerEnergy","", 20, 0, 1);
      lout->Add(hMatchedTruthslShowerEnergy);

      hMatchedTruthPi0OA = new TH1D("m013hMatchedTruthPi0OA","", 20, 0, 180);
      lout->Add(hMatchedTruthPi0OA);
      hMatchedTruthldShowerOA = new TH1D("m014hMatchedTruthldShowerOA","", 20, 0, 180);
      lout->Add(hMatchedTruthldShowerOA);
      hMatchedTruthslShowerOA = new TH1D("m015hMatchedTruthslShowerOA","", 20, 0, 180);
      lout->Add(hMatchedTruthslShowerOA);

      hPi0ThetaRes = new TH2D("m018hPi0ThetaFitRes_RES",";Truth #pi^{0} #theta Angle (deg);Rec.- Truth (deg)", 30, 0, 180, 40, -50, 50);
      lout->Add(hPi0ThetaRes);
      hPi0PhiRes = new TH2D("m019hPi0PhiFitRes_RES",";Truth #pi^{0} #phi Angle (deg);Rec.- Truth (deg)", 30, -180, 180, 40, -50, 50);
      lout->Add(hPi0PhiRes);

      hPi0ThetaResRaw = new TH2D("m018hPi0ThetaFitRes_RES_RAW",";Truth #pi^{0} #theta Angle (deg);Rec.(raw) - Truth (deg)", 30, 0, 180, 40, -50, 50);
      lout->Add(hPi0ThetaResRaw);
      hPi0PhiResRaw = new TH2D("m019hPi0PhiFitRes_RES_RAW",";Truth #pi^{0} #phi Angle (deg);Rec.(raw) - Truth (deg)", 30, -180, 180, 40, -50, 50);
      lout->Add(hPi0PhiResRaw);

      
      hMatchedTruthBeamMomentum = new TH1D("m020hMatchedTruthBeamMomentum",";Truth Beam #pi^{+} Mom. (GeV/c);Efficiency", 20, 0.4, 1.0);
      lout->Add(hMatchedTruthBeamMomentum);

      
      //hPi0MassResFit = new TH2D("i009Pi0Mass_RES_FIT","", 20, 0, 0.5, 20, -0.5, 0.5);
      //lout->Add(hPi0MassResFit);
      //hShowerOpenAngleResFit = new TH2D("i010ShowerOpenAngle_RES_FIT","", 20, 0, 180, 20, -1.1, 1.1);
      //lout->Add(hShowerOpenAngleResFit);

      

      hTruthTestFFEnergyM1 = new TH1D("i000hTruthTestFFEnergyM1", ";Truth #pi^{+} FFE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin+400, xsecmax+400);
      lout->Add(hTruthTestFFEnergyM1);
      hTruthTestFFEnergyM2 = new TH1D("i000hTruthTestFFEnergyM2", ";Truth #pi^{+} FFE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin+400, xsecmax+400);
      lout->Add(hTruthTestFFEnergyM2);
      hTruthTestFFEnergyDiff = new TH1D("i000hTruthTestFFEnergyDiff", ";Truth #pi^{+} FFE (MeV);#sigma_{CEX} (mb)", 50, -10, 5);
      lout->Add(hTruthTestFFEnergyDiff);

      hTruthTestIntEnergyM1 = new TH1D("i100hTruthTestIntEnergyM1", ";Truth #pi^{+} IntE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthTestIntEnergyM1);
      hTruthTestIntEnergyM2 = new TH1D("i100hTruthTestIntEnergyM2", ";Truth #pi^{+} IntE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthTestIntEnergyM2);
      hTruthTestIntEnergyDiff = new TH1D("i100hTruthTestIntEnergyDiff", ";Truth #pi^{+} IntE (MeV);#sigma_{CEX} (mb)", 100, -50, 5);
      lout->Add(hTruthTestIntEnergyDiff);

      hTruthTestSameBin = new TH1D("i000hTruthTestSameBin", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, 0, 20);
      lout->Add(hTruthTestSameBin);
      

      hTruthInitialHist = new TH1D("i000hTruthInitialHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthInitialHist);
      hTruthBeamInitialHist = new TH1D("i000hTruthBeamInitialHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2050, 0, 2050);
      lout->Add(hTruthBeamInitialHist);
      hNewTruthBeamInitialHist = new TH1D("i000hNewTruthBeamInitialHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 20, 0, 2000);
      lout->Add(hNewTruthBeamInitialHist);
      hTruthBeamIncidentHist = new TH1D("i000hTruthBeamIncidentHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hTruthBeamIncidentHist);
      hNewTruthBeamIncidentHist = new TH1D("i000hNewTruthBeamIncidentHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 20, 0, 2000);
      lout->Add(hNewTruthBeamIncidentHist);

      hTruthBeamInitialHist_50MeVbin = new TH1D("i000hTruthBeamInitialHist_50MeVbin", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 20, 0, 2000);
      lout->Add(hTruthBeamInitialHist_50MeVbin);
      hTruthBeamIncidentHist_50MeVbin = new TH1D("i000hTruthBeamIncidentHist_50MeVbin", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 20, 0, 2000);
      lout->Add(hTruthBeamIncidentHist_50MeVbin);

      hTruthBeamIncidentHistOldM = new TH1D("i000hTruthBeamIncidentHistOldM", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hTruthBeamIncidentHistOldM);

      hTruthIncidentHist = new TH1D("i001hTruthIncidentHist", ";Indident #pi^{+} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthIncidentHist);
      hTruthSingleIncidentHist = new TH1D("i001hTruthSingleIncidentHist", ";Indident #pi^{+} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthSingleIncidentHist);
      hTruthCalcIncidentHist = new TH1D("i001hTruthCalcIncidentHist", ";Indident #pi^{+} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthCalcIncidentHist);

      hTruthBeamCalcIncidentHist = new TH1D("i001hTruthBeamCalcIncidentHist", ";Indident #pi^{+} KE (MeV);Candidates", 2000, 0, 2000);
      lout->Add(hTruthBeamCalcIncidentHist);
      hNewTruthBeamCalcIncidentHist = new TH1D("i001hNewTruthBeamCalcIncidentHist", ";Indident #pi^{+} KE (MeV);Candidates", 20, 0, 2000);  
      lout->Add(hNewTruthBeamCalcIncidentHist);

      hTruthBeamCalcIncidentHist_50MeVbin = new TH1D("i001hTruthBeamCalcIncidentHist_50MeVbin", ";Indident #pi^{+} KE (MeV);Candidates", 20, 0, 2000);
      lout->Add(hTruthBeamCalcIncidentHist_50MeVbin);

      hTruthInteractingHist = new TH1D("i002hTruthInteractingHist", ";Interacting #pi^{+} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthInteractingHist);
      hNewTruthInteractingHist = new TH1D("i002hNewTruthInteractingHist", ";Interacting #pi^{+} KE (MeV);Candidates", 20, 0, 2000);
      lout->Add(hNewTruthInteractingHist);
      hNewTruthInteractingHistTest = new TH1D("testi002hNewTruthInteractingHist", ";Interacting #pi^{+} KE (MeV);Candidates", 20, 0, 1);
      lout->Add(hNewTruthInteractingHistTest);
      hTruthSingleInteractingHist = new TH1D("i002hTruthSingleInteractingHist", ";Interacting #pi^{+} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthSingleInteractingHist);
      hTruthBeamInteractingHist = new TH1D("i002hTruthBeamInteractingHist", ";Interacting #pi^{+} KE (MeV);Candidates", 2050, 0, 2050);
      lout->Add(hTruthBeamInteractingHist);
      hNewTruthBeamInteractingHist = new TH1D("i002hNewTruthBeamInteractingHist", ";Interacting #pi^{+} KE (MeV);Candidates", 20, 0, 2000);
      lout->Add(hNewTruthBeamInteractingHist);

      hTruthBeamInteractingHist_50MeVbin = new TH1D("i002hTruthBeamInteractingHist_50MeVbin", ";Interacting #pi^{+} KE (MeV);Candidates", 20, 0, 2000);
      lout->Add(hTruthBeamInteractingHist_50MeVbin);

      hTruthTotalXSecHist = new TH1D("i003hTruthTotalXSecHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthTotalXSecHist);

      hTruthCEXInteractingHist = new TH1D("i004hTruthCEXInteractingHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthCEXInteractingHist);
      hNewTruthCEXInteractingHist = new TH1D("i004hNewTruthCEXInteractingHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 20, 0, 2000);
      lout->Add(hNewTruthCEXInteractingHist);
      hNewTruthCEXInteractingHistTest = new TH1D("testi004hNewTruthCEXInteractingHist", ";Truth #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 20, 0, 1);
      lout->Add(hNewTruthCEXInteractingHistTest);

      // KE
      hTruthDiffCEXInteractingHist_700MeV = new TH1D("i005hTruthDiffCEXInteractingHist_700MeV", ";Outgoing #pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_700MeV);

      hTruthDiffCEXInteractingHist_800MeV = new TH1D("i006hTruthDiffCEXInteractingHist_800MeV", ";Outgoing #pi^{0} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_800MeV);
      hTruthDiffCEXInteractingHist_650to800MeV = new TH1D("i006hTruthDiffCEXInteractingHist_650to800MeV", ";Outgoing #pi^{0} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_650to800MeV);
      hTruthSingleDiffCEXInteractingHist_800MeV = new TH1D("i006hTruthSingleDiffCEXInteractingHist_800MeV", ";Outgoing #pi^{0} KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthSingleDiffCEXInteractingHist_800MeV);

      hTruthDiffCEXInteractingHist_900MeV = new TH1D("i007hTruthDiffCEXInteractingHist_900MeV", ";Outgoing #pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_900MeV);

      // Theta
      hTruthDiffCEXInteractingHistTheta_700MeV = new TH1D("i055hTruthDiffCEXInteractingHistTheta_700MeV", ";Outgoing #pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hTruthDiffCEXInteractingHistTheta_700MeV);

      hTruthDiffCEXInteractingHistTheta_800MeV = new TH1D("i066hTruthDiffCEXInteractingHistTheta_800MeV", ";Outgoing #pi^{0} #theta (degree);Candidates", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hTruthDiffCEXInteractingHistTheta_800MeV);
      hTruthDiffCEXInteractingHistTheta_650to800MeV = new TH1D("i066hTruthDiffCEXInteractingHistTheta_650to800MeV", ";Outgoing #pi^{0} #theta (degree);Candidates", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hTruthDiffCEXInteractingHistTheta_650to800MeV);
      
      hTruthDiffCEXInteractingHistTheta_900MeV = new TH1D("i077hTruthDiffCEXInteractingHistTheta_900MeV", ";Outgoing #pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hTruthDiffCEXInteractingHistTheta_900MeV);

      pi0theta = new TH1D("pi0theta", ";Outgoing #pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(pi0theta);

      // Cos Theta
      hTruthDiffCEXInteractingHistCosTheta_700MeV = new TH1D("i555hTruthDiffCEXInteractingHistCosTheta_700MeV", ";Outgoing #pi^{0} cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hTruthDiffCEXInteractingHistCosTheta_700MeV);

      hTruthDiffCEXInteractingHistCosTheta_800MeV = new TH1D("i666hTruthDiffCEXInteractingHistCosTheta_800MeV", ";Outgoing #pi^{0} cos#theta;Candidates", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hTruthDiffCEXInteractingHistCosTheta_800MeV);
      hTruthDiffCEXInteractingHistCosTheta_650to800MeV = new TH1D("i666hTruthDiffCEXInteractingHistCosTheta_650to800MeV", ";Outgoing #pi^{0} cos#theta;Candidates", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hTruthDiffCEXInteractingHistCosTheta_650to800MeV);
     
      hTruthDiffCEXInteractingHistCosTheta_900MeV = new TH1D("i777hTruthDiffCEXInteractingHistCosTheta_900MeV", ";Outgoing #pi^{0} cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hTruthDiffCEXInteractingHistCosTheta_900MeV);


      hTruthDiffCEXInteractingHist_500MeV = new TH1D("i008hTruthDiffCEXInteractingHist_500MeV", ";Outgoing #pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_500MeV);

      hTruthDiffCEXInteractingHist_300MeV = new TH1D("i009hTruthDiffCEXInteractingHist_300MeV", ";Outgoing #pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_300MeV);

      hTruthDiffCEXInteractingHist_200MeV = new TH1D("i010hTruthDiffCEXInteractingHist_200MeV", ";Outgoing #pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_200MeV);

      hTruthDiffCEXInteractingHist_100MeV = new TH1D("i011hTruthDiffCEXInteractingHist_100MeV", ";Outgoing #pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist_100MeV);

      hTruthDiffCEXInteractingHist = new TH1D("i012hTruthDiffCEXInteractingHist", ";Outgoing Pion KE (MeV);Candidates", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthDiffCEXInteractingHist);


      hTruthOutgoingKEcosTheta_Low = new TH2D("i013hTruthOutgoingKEcosTheta_Low_REG", ";Outgoing Pion KE (MeV);Outgoing Pion cos #theta", xsecbin, xsecmin, xsecmax, 20, -1, 1);
      lout->Add(hTruthOutgoingKEcosTheta_Low);
      hTruthOutgoingKEcosTheta_Mid = new TH2D("i014hTruthOutgoingKEcosTheta_Mid_REG", ";Outgoing Pion KE (MeV);Outgoing Pion cos #theta", xsecbin, xsecmin, xsecmax, 20, -1, 1);
      lout->Add(hTruthOutgoingKEcosTheta_Mid);
      hTruthOutgoingKEcosTheta_High = new TH2D("i015hTruthOutgoingKEcosTheta_High_REG", ";Outgoing Pion KE (MeV);Outgoing Pion cos #theta", xsecbin, xsecmin, xsecmax, 20, -1, 1);
      lout->Add(hTruthOutgoingKEcosTheta_High);

      hTruthOutgoingKEcosTheta = new TH2D("i016hTruthOutgoingKEcosTheta_REG", ";Outgoing Pion KE (MeV);Outgoing Pion cos #theta", xsecbin, xsecmin, xsecmax, 20, -1, 1);
      lout->Add(hTruthOutgoingKEcosTheta);

      //hRecoInitialHist = new TH1D("i016hRecoInitialHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 1000, 0, 1000);
      hRecoInitialHist = new TH1D("i016hRecoInitialHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hRecoInitialHist);
      hRecoBeamInteractingHist = new TH1D("i016hRecoBeamInteractingHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hRecoBeamInteractingHist);

      hRecoIncidentHist = new TH1D("i017hRecoIncidentHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoIncidentHist);
      hRecoInteractingHist = new TH1D("i018hRecoInteractingHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoInteractingHist);

      hTruthMatchedIncidentHist = new TH1D("i019hTruthMatchedIncidentHist", ";TruthMatched #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthMatchedIncidentHist);
      hTruthMatchedInteractingHist = new TH1D("i020hTruthMatchedInteractingHist", ";TruthMatched #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hTruthMatchedInteractingHist);

      hRecoBckSubIncidentHist = new TH1D("i021hRecoBckSubIncidentHist", ";RecoBckSub #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoBckSubIncidentHist);
      hRecoBckSubInteractingHist = new TH1D("i022hRecoBckSubInteractingHist", ";RecoBckSub #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoBckSubInteractingHist);

      hNonPionBeamIncidentHist = new TH1D("i023hNonPionBeamIncidentHist", ";NonPionBeam #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hNonPionBeamIncidentHist);
      hNonPionBeamInteractingHist = new TH1D("i024hNonPionBeamInteractingHist", ";NonPionBeam #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hNonPionBeamInteractingHist);

      hUnFoldedInitialHist = new TH1D("i024hUnFoldedInitialHist", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hUnFoldedInitialHist);
      hUnFoldedBeamInteractingHist = new TH1D("i024hUnFoldedBeamInteractingHist", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hUnFoldedBeamInteractingHist);

      hUnFoldedInteractingHist = new TH1D("i025hUnFoldedInteractingHist", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hUnFoldedInteractingHist);
      hUnFoldedIncidentHist = new TH1D("i026hUnFoldedIncidentHist", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hUnFoldedIncidentHist);
      hUnfoldedBeamIncidentHist = new TH1D("i027hUnFoldedBeamIncidentHist", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hUnfoldedBeamIncidentHist);

      hRecoPi0KEHist = new TH1D("i030hRecoPi0KEHist", ";#pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoPi0KEHist);
      hUnFoldedPi0KEHist = new TH1D("i031hUnFoldedPi0KEHist", ";#pi^{0} KE (MeV);d#sigma/dT_{#pi^{0}} (mb/MeV)", xsecbin, xsecmin, xsecmax);
      lout->Add(hUnFoldedPi0KEHist);

      hRecoPi0CosThetaHist = new TH1D("i032hRecoPi0CosThetaHist", ";#pi^{0} Cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hRecoPi0CosThetaHist);
      hUnFoldedPi0CosThetaHist = new TH1D("i033hUnFoldedPi0CosThetaHist", ";#pi^{0} Cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hUnFoldedPi0CosThetaHist);

      hRecoPi0ThetaHist = new TH1D("i034hRecoPi0ThetaHist", ";#pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hRecoPi0ThetaHist);
      hUnFoldedPi0ThetaHist = new TH1D("i035hUnFoldedPi0ThetaHist", ";#pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hUnFoldedPi0ThetaHist);


      hTruthIncidentZ = new TH1D("i003hTruthIncidentZ", ";Truth Position Z (cm);Candidates", 50, -1, 500);
      lout->Add(hTruthIncidentZ);
      hTruthIncidentZ_SCE = new TH1D("i004hTruthIncidentZ_SCE", ";Truth Position Z SCE (cm);Candidates", 50, -1, 500);
      lout->Add(hTruthIncidentZ_SCE);
      
      //removing overlay
      hTruthCEXDaughters_Low = new TH2D("i101hTruthCEXDaughters_Low", ";Number of Particles;Candidates", 15, 0, 15, 3, -0.5, 2.5);
      lout->Add(hTruthCEXDaughters_Low);
      hTruthCEXDaughters_Middle = new TH2D("i102hTruthCEXDaughters_Middle", ";Number of Particles;Candidates", 15, 0, 15, 3, -0.5, 2.5);
      lout->Add(hTruthCEXDaughters_Middle);
      hTruthCEXDaughters_High = new TH2D("i103hTruthCEXDaughters_High", ";Number of Particles;Candidates", 15, 0, 15, 3, -0.5, 2.5);
      lout->Add(hTruthCEXDaughters_High);

      hTruthPi0DaughtersCosTheta_Low = new TH1D("i104hTruthPi0DaughtersCosTheta_Low", ";Cos#theta;Candidates", 20, -1, 1);
      lout->Add(hTruthPi0DaughtersCosTheta_Low);
      hTruthPi0DaughtersCosTheta_Middle = new TH1D("i105hTruthPi0DaughtersCosTheta_Middle", ";Cos#theta;Candidates", 20, -1, 1);
      lout->Add(hTruthPi0DaughtersCosTheta_Middle);
      hTruthPi0DaughtersCosTheta_High = new TH1D("i104hTruthPi0DaughtersCosTheta_High", ";Cos#theta;Candidates", 20, -1, 1);
      lout->Add(hTruthPi0DaughtersCosTheta_High);


      hRecoIncidentHistData = new TH1D("i017hRecoIncidentHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hRecoIncidentHistData);
      //hRecoInteractingHistData = new TH1D("i018hRecoInteractingHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      //lout->Add(hRecoInteractingHistData);
      //hRecoBeamInitialHistData = new TH1D("i019hRecoBeamInitialHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 1000, 0, 1000);
      //lout->Add(hRecoBeamInitialHistData);
      //hRecoBeamInteractingHistData = new TH1D("i020hRecoBeamInteractingHistData", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", 1000, 0, 1000);
      //lout->Add(hRecoBeamInteractingHistData);

      hUnFoldedInteractingHistData = new TH1D("i025hUnFoldedInteractingHistData", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hUnFoldedInteractingHistData);
      hUnFoldedIncidentHistData = new TH1D("i026hUnFoldedIncidentHistData", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hUnFoldedIncidentHistData);

      hUnFoldedBeamInitialHistData = new TH1D("i027hUnFoldedBeamInitialHistData", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hUnFoldedBeamInitialHistData);
      hUnFoldedBeamInteractingHistData = new TH1D("i028hUnFoldedBeamInteractingHistData", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hUnFoldedBeamInteractingHistData);
      hUnfoldedBeamIncidentHistData = new TH1D("i029hUnFoldedBeamIncidentHistData", ";#pi^{+} KE (MeV);#sigma_{CEX} (mb)", 2000, 0, 2000);
      lout->Add(hUnfoldedBeamIncidentHistData);

      //hRecoPi0KEHistData = new TH1D("i030hRecoPi0KEHistData", ";Reco #pi^{0} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      //lout->Add(hRecoPi0KEHistData);
      hUnFoldedPi0KEHistData = new TH1D("i031hUnFoldedPi0KEHistData", ";Reco #pi^{0} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hUnFoldedPi0KEHistData);

      //hRecoPi0CosThetaHistData = new TH1D("i032hRecoPi0CosThetaHistData", ";#pi^{0} Cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      //lout->Add(hRecoPi0CosThetaHistData);
      hUnFoldedPi0CosThetaHistData = new TH1D("i033hUnFoldedPi0CosThetaHistData", ";#pi^{0} Cos#theta;d#sigma/dcos#theta_{#pi^{0}} (mb)", xseccosthetabin, xseccosthetamin, xseccosthetamax);
      lout->Add(hUnFoldedPi0CosThetaHistData);

      //hRecoPi0ThetaHistData = new TH1D("i034hRecoPi0ThetaHistData", ";#pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      //lout->Add(hRecoPi0ThetaHistData);
      hUnFoldedPi0ThetaHistData = new TH1D("i035hUnFoldedPi0ThetaHistData", ";#pi^{0} #theta (degree);d#sigma/d#theta_{#pi^{0}} (mb/degree)", xsecthetabin, xsecthetamin, xsecthetamax);
      lout->Add(hUnFoldedPi0ThetaHistData);

      
      hOldIncSizeDiff = new TH1D("i032hOldIncSizeDiff", ";Old Size diff;Candidates", 40, -20, 20);
      lout->Add(hOldIncSizeDiff);
      hNewIncIntervalDiff = new TH1D("i033hNewIncIntervalDiff", ";New Interval diff;Candidates", 40, -100, 100);
      lout->Add(hNewIncIntervalDiff);

      // Beam scraper cuts study
      hBeamInstVSTruthKEffNoCuts = new TH2D("j001hBeamInstVSTruthKEffNoCuts_RES",";KE_{Beam  Inst.}^{reco.} (MeV);KE_{ff}^{true} (MeV)", 300, 0, 1500, 300, 0, 1500);
      lout->Add(hBeamInstVSTruthKEffNoCuts);
      hBeamInstVSTruthKEffAfterCuts = new TH2D("j002hBeamInstVSTruthKEffAfterCuts_RES",";KE_{Beam  Inst.}^{reco.} (MeV);KE_{ff}^{true} (MeV)", 300, 0, 1500, 300, 0, 1500);
      lout->Add(hBeamInstVSTruthKEffAfterCuts);
      hUpStreamELoss700MeV = new TH1D("j003hUpStreamELoss700MeV", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss700MeV);
      hUpStreamELoss800MeV = new TH1D("j003hUpStreamELoss800MeV", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss800MeV);
      hUpStreamELoss900MeV = new TH1D("j003hUpStreamELoss900MeV", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss900MeV);
      hUpStreamELoss1000MeV = new TH1D("j003hUpStreamELoss1000MeV", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss1000MeV);

      hBeamInstXVSBeamInstYBeam700MeV = new TH2D("j004hBeamInstXVSBeamInstYBeam700MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYBeam700MeV);
      hBeamInstXVSBeamInstYScraper700MeV = new TH2D("j005hBeamInstXVSBeamInstYScraper700MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYScraper700MeV);
      
      hBeamInstXVSBeamInstYBeam800MeV = new TH2D("j004hBeamInstXVSBeamInstYBeam800MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYBeam800MeV);
      hBeamInstXVSBeamInstYScraper800MeV = new TH2D("j005hBeamInstXVSBeamInstYScraper800MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYScraper800MeV);

      hBeamInstXVSBeamInstYBeam900MeV = new TH2D("j004hBeamInstXVSBeamInstYBeam900MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYBeam900MeV);
      hBeamInstXVSBeamInstYScraper900MeV = new TH2D("j005hBeamInstXVSBeamInstYScraper900MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYScraper900MeV);

      hBeamInstXVSBeamInstYBeam1000MeV = new TH2D("j004hBeamInstXVSBeamInstYBeam1000MeV_RES",";X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYBeam1000MeV);
      hBeamInstXVSBeamInstYScraper1000MeV = new TH2D("j005hBeamInstXVSBeamInstYScraper1000MeV_RES",";;X_{ff}^{reco.} (cm);Y_{ff}^{reco.} (cm)", 100, -50, 10, 100, 400, 460);
      lout->Add(hBeamInstXVSBeamInstYScraper1000MeV);

      hBeamInstVSTruthKEffAfterScraperCuts = new TH2D("j006hBeamInstVSTruthKEffAfterScraperCuts_RES",";KE_{Beam  Inst.}^{reco.} (MeV);KE_{ff}^{true} (MeV)", 300, 0, 1500, 300, 0, 1500);
      lout->Add(hBeamInstVSTruthKEffAfterScraperCuts);
      hUpStreamELoss700MeVAfterScraperCuts = new TH1D("j007hUpStreamELoss700MeVAfterScraperCuts", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss700MeVAfterScraperCuts);
      hUpStreamELoss800MeVAfterScraperCuts = new TH1D("j007hUpStreamELoss800MeVAfterScraperCuts", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss800MeVAfterScraperCuts);
      hUpStreamELoss900MeVAfterScraperCuts = new TH1D("j007hUpStreamELoss900MeVAfterScraperCuts", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss900MeVAfterScraperCuts);
      hUpStreamELoss1000MeVAfterScraperCuts = new TH1D("j007hUpStreamELoss1000MeVAfterScraperCuts", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300);
      lout->Add(hUpStreamELoss1000MeVAfterScraperCuts);

      hUpStreamELossAfterSmearingAndWeight = new TH2D("j008hUpStreamELossAfterSmearingAndWeight", ";KE_{Beam  Inst.}^{reco.} - KE_{ff}^{true} (MeV);Candidates", 120, -300, 300, 8, -0.5, 7.5);
      lout->Add(hUpStreamELossAfterSmearingAndWeight);

      hUpStreamELossRes = new TH2D("j009hUpStreamELossFitGaus_RES", ";KE_{Beam  Inst.}^{reco.} (MeV); Upstream E Loss/KE_{Beam  Inst.}^{reco.}", 60, 600, 1200, 30, -0.3, 0.3);
      lout->Add(hUpStreamELossRes);

      // Sungbin's new method
      hNewRecoInitialHist = new TH1D("n001hNewRecoInitialHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hNewRecoInitialHist);
      hNewRecoBeamInteractingHist = new TH1D("n002hNewRecoBeamInteractingHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hNewRecoBeamInteractingHist);
      hNewRecoIncidentHist = new TH1D("n003hNewRecoIncidentHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hNewRecoIncidentHist);
      hNewRecoInteractingHist = new TH1D("n004hNewRecoInteractingHist", ";Reco #pi^{+} KE (MeV);#sigma_{CEX} (mb)", xsecbin, xsecmin, xsecmax);
      lout->Add(hNewRecoInteractingHist);

      // Class N - CVM bin histograms
      //const double BinE1_1x1[] = {0,1.2};
      const double BinE2_1x1[] = {0,0.8};
      const double BinE1_1x1[] = {0,180*TMath::DegToRad()};

      //const double BinE1_2x2[] = {0,0.3,1.2};
      const double BinE2_2x2[] = {0,0.15,0.8};
      const double BinE1_2x2[] = {0,40*TMath::DegToRad(),180*TMath::DegToRad()};

      //const double BinE1_3x3[] = {0,0.25,0.35,1.2};
      const double BinE2_3x3[] = {0,0.12,0.2,0.8};
      const double BinE1_3x3[] = {0,25*TMath::DegToRad(),60*TMath::DegToRad(),180*TMath::DegToRad()};

      //const double BinE1_4x4[] = {0,0.2,0.3,0.4,1.2};
      const double BinE2_4x4[] = {0,0.1,0.16,0.26,0.8};
      const double BinE1_4x4[] = {0,20*TMath::DegToRad(),40*TMath::DegToRad(),80*TMath::DegToRad(),180*TMath::DegToRad()};


      hV11_1x1 = new TH2D("hV11_1x1Bin", "", sizeof(BinE1_1x1)/sizeof(double)-1, BinE1_1x1, sizeof(BinE2_1x1)/sizeof(double)-1, BinE2_1x1);
      lout->Add(hV11_1x1);
      hV12_1x1 = new TH2D("hV12_1x1Bin", "", sizeof(BinE1_1x1)/sizeof(double)-1, BinE1_1x1, sizeof(BinE2_1x1)/sizeof(double)-1, BinE2_1x1);
      lout->Add(hV12_1x1);
      hV22_1x1 = new TH2D("hV22_1x1Bin", "", sizeof(BinE1_1x1)/sizeof(double)-1, BinE1_1x1, sizeof(BinE2_1x1)/sizeof(double)-1, BinE2_1x1);
      lout->Add(hV22_1x1);

      hV11_2x2 = new TH2D("hV11_2x2Bin", "", sizeof(BinE1_2x2)/sizeof(double)-1, BinE1_2x2, sizeof(BinE2_2x2)/sizeof(double)-1, BinE2_2x2);
      lout->Add(hV11_2x2);
      hV12_2x2 = new TH2D("hV12_2x2Bin", "", sizeof(BinE1_2x2)/sizeof(double)-1, BinE1_2x2, sizeof(BinE2_2x2)/sizeof(double)-1, BinE2_2x2);
      lout->Add(hV12_2x2);
      hV22_2x2 = new TH2D("hV22_2x2Bin", "", sizeof(BinE1_2x2)/sizeof(double)-1, BinE1_2x2, sizeof(BinE2_2x2)/sizeof(double)-1, BinE2_2x2);
      lout->Add(hV22_2x2);

      hV11_3x3 = new TH2D("hV11_3x3Bin", "", sizeof(BinE1_3x3)/sizeof(double)-1, BinE1_3x3, sizeof(BinE2_3x3)/sizeof(double)-1, BinE2_3x3);
      lout->Add(hV11_3x3);
      hV12_3x3 = new TH2D("hV12_3x3Bin", "", sizeof(BinE1_3x3)/sizeof(double)-1, BinE1_3x3, sizeof(BinE2_3x3)/sizeof(double)-1, BinE2_3x3);
      lout->Add(hV12_3x3);
      hV22_3x3 = new TH2D("hV22_3x3Bin", "", sizeof(BinE1_3x3)/sizeof(double)-1, BinE1_3x3, sizeof(BinE2_3x3)/sizeof(double)-1, BinE2_3x3);
      lout->Add(hV22_3x3);

      hV11_4x4 = new TH2D("hV11_4x4Bin", "", sizeof(BinE1_4x4)/sizeof(double)-1, BinE1_4x4, sizeof(BinE2_4x4)/sizeof(double)-1, BinE2_4x4);
      lout->Add(hV11_4x4);
      hV12_4x4 = new TH2D("hV12_4x4Bin", "", sizeof(BinE1_4x4)/sizeof(double)-1, BinE1_4x4, sizeof(BinE2_4x4)/sizeof(double)-1, BinE2_4x4);
      lout->Add(hV12_4x4);
      hV22_4x4 = new TH2D("hV22_4x4Bin", "", sizeof(BinE1_4x4)/sizeof(double)-1, BinE1_4x4, sizeof(BinE2_4x4)/sizeof(double)-1, BinE2_4x4);
      lout->Add(hV22_4x4);

      hBinSize_1x1 = new TH2D("hBinSize_1x1Bin", "", sizeof(BinE1_1x1)/sizeof(double)-1, BinE1_1x1, sizeof(BinE2_1x1)/sizeof(double)-1, BinE2_1x1);
      lout->Add(hBinSize_1x1);
      hBinSize_2x2 = new TH2D("hBinSize_2x2Bin", "", sizeof(BinE1_2x2)/sizeof(double)-1, BinE1_2x2, sizeof(BinE2_2x2)/sizeof(double)-1, BinE2_2x2);
      lout->Add(hBinSize_2x2);
      hBinSize_3x3 = new TH2D("hBinSize_3x3Bin", "", sizeof(BinE1_3x3)/sizeof(double)-1, BinE1_3x3, sizeof(BinE2_3x3)/sizeof(double)-1, BinE2_3x3);
      lout->Add(hBinSize_3x3);
      hBinSize_4x4 = new TH2D("hBinSize_4x4Bin", "", sizeof(BinE1_4x4)/sizeof(double)-1, BinE1_4x4, sizeof(BinE2_4x4)/sizeof(double)-1, BinE2_4x4);
      lout->Add(hBinSize_4x4);

      hCVM = new TH2D("hCVM_Bin", "", 3, 0, 3, 3, 0, 3);
      lout->Add(hCVM);

      hsigma11 = new TH2D("hsigma11", "", 6, 0.1, 0.9, 10, 0, 1);
      lout->Add(hsigma11);
      hsigma11Mean = new TH1D("hsigma11Mean", "", 6, 0.1, 0.9);
      lout->Add(hsigma11Mean);

      hsigma22 = new TH2D("hsigma22", "", 6, 0.05, 0.4, 10, 0, 1);
      lout->Add(hsigma22);
      hsigma22Mean = new TH1D("hsigma22Mean", "", 6, 0.05, 0.4);
      lout->Add(hsigma22Mean);


      hShowerE1PreFitRes = new TH2D("x001hShowerE1Pre_RES","",20, 0, 0.8, 20, -1.1, 1.1);
      lout->Add(hShowerE1PreFitRes);
      hShowerE1PostFitRes = new TH2D("x002hShowerE1Post_RES"," ",20, 0, 0.8, 20, -1.1, 1.1);
      lout->Add(hShowerE1PostFitRes);

      hShowerE2PreFitRes = new TH2D("x003hShowerE2Pre_RES"," ",20, 0, 0.5, 20, -1.1, 1.1);
      lout->Add(hShowerE2PreFitRes);
      hShowerE2PostFitRes = new TH2D("x004hShowerE2Post_RES"," ",20, 0, 0.5, 20, -1.1, 1.1);
      lout->Add(hShowerE2PostFitRes);

      hShowerOAPreFitRes = new TH2D("x005hShowerOAPre_RES"," ",20, 0, 180, 20, -30, 30);
      lout->Add(hShowerOAPreFitRes);
      hShowerOAPostFitRes = new TH2D("x006hShowerOAPost_RES"," ",20, 0, 180, 20, -30, 30);
      lout->Add(hShowerOAPostFitRes);

      hPi0MomPreFitRes = new TH2D("x007hPi0MomPreFit_REG"," ",20, 0, 1, 20, 0, 1);
      lout->Add(hPi0MomPreFitRes);
      hPi0MomPostFitRes = new TH2D("x008hPi0MomPostFit_REG"," ",20, 0, 1, 20, 0, 1);
      lout->Add(hPi0MomPostFitRes);

      hPi0MomPreFitResPaper = new TH2D("x007hPi0MomPreFit_REG_Paper"," ",20, 0, 1, 20, 0, 1);
      lout->Add(hPi0MomPreFitResPaper);
      hPi0MomPostFitResPaper = new TH2D("x008hPi0MomPostFit_REG_Paper"," ",20, 0, 1, 20, 0, 1);
      lout->Add(hPi0MomPostFitResPaper);

      hPi0ThetaFitRes = new TH2D("x009hPi0ThetaFitRes_RES"," ",30, 0, 180, 25, -50, 50);
      lout->Add(hPi0ThetaFitRes);
      hLDShowerThetaFitRes = new TH2D("x010hLDShowerThetaFitRes_RES"," ",20, 0, 180, 20, -50, 50);
      lout->Add(hLDShowerThetaFitRes);

      hShowerE1Compare = new TH1D("y001hShowerE1Compare",";Rec./Truth - 1;Candidates", 40, -1.1, 1.1);
      lout->Add(hShowerE1Compare);
      hShowerE2Compare = new TH1D("y002hShowerE2Compare",";Rec./Truth - 1;Candidates", 40, -1.1, 1.1);
      lout->Add(hShowerE2Compare);
      hShowerOACompare = new TH1D("y003hShowerOACompare",";Rec. - Truth (deg);Candidates", 40, -50, 50);
      lout->Add(hShowerOACompare);
      hPi0MassCompare = new TH1D("y004hPi0MassCompare",";Rec./Truth - 1;Candidates", 100, 0, 0.3);
      lout->Add(hPi0MassCompare);
      hPi0MomCompare = new TH1D("y005hPi0MomCompare",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomCompare);

      hPi0MomComparePre = new TH1D("y005hPi0MomComparePre",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePre);

      hPi0MomComparePaper = new TH1D("y005hPi0MomComparePaper",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePaper);

      hPi0MomComparePrePaper = new TH1D("y005hPi0MomComparePrePaper",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePrePaper);

      hPi0MomComparePostPaper = new TH1D("y005hPi0MomComparePostPaper",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePostPaper);

      hPi0MomCompareStandard = new TH1D("y005hPi0MomCompareStandard",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomCompareStandard);

      hPi0MomComparePreStandard = new TH1D("y005hPi0MomComparePreStandard",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePreStandard);

      hPi0MomComparePostStandard = new TH1D("y005hPi0MomComparePostStandard",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePostStandard);


      hLDShowerE = new TH1D("y006hLDShowerE","", 40, -1, 1);
      lout->Add(hLDShowerE);
      hSLShowerE = new TH1D("y006hSLShowerE","", 40, -1, 1);
      lout->Add(hSLShowerE);
      hOAShower = new TH1D("y006hOAShower","", 40, -1, 1);
      lout->Add(hOAShower);

      hShowerE1ComparePost = new TH1D("y001hShowerE1ComparePost",";Rec./Truth - 1;Candidates", 40, -1.1, 1.1);
      lout->Add(hShowerE1ComparePost);
      hShowerE2ComparePost = new TH1D("y002hShowerE2ComparePost",";Rec./Truth - 1;Candidates", 40, -1.1, 1.1);
      lout->Add(hShowerE2ComparePost);
      hShowerOAComparePost = new TH1D("y003hShowerOAComparePost",";Rec. - Truth (deg);Candidates", 40, -50, 50);
      lout->Add(hShowerOAComparePost);
      hPi0MassComparePost = new TH1D("y004hPi0MassComparePost",";Rec./Truth - 1;Candidates", 100, 0, 0.3);
      lout->Add(hPi0MassComparePost);
      hPi0MomComparePost = new TH1D("y005hPi0MomComparePost",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomComparePost);

      

      hPi0MomNorm = new TH1D("y006hPi0MomNorm",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomNorm);
      hPi0MomEOA = new TH1D("y006hPi0MomEOA",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomEOA);
      hPi0MomAsym = new TH1D("y006hPi0MomAsym",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPi0MomAsym);

      hPrePi0MomNorm = new TH1D("y007hPrePi0MomNorm",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPrePi0MomNorm);
      hPrePi0MomEOA = new TH1D("y007hPrePi0MomEOA",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPrePi0MomEOA);
      hPrePi0MomAsym = new TH1D("y007hPrePi0MomAsym",";Rec./Truth - 1;Candidates", 40, -1, 1);
      lout->Add(hPrePi0MomAsym);
      

      hPi0ShowerE1Compare_REG = new TH2D("y010hPi0ShowerE1Compare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1.1, 1.1, 40, -1.1, 1.1);
      lout->Add(hPi0ShowerE1Compare_REG);
      hPi0ShowerE2Compare_REG = new TH2D("y011hPi0ShowerE2Compare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1.1, 1.1, 40, -1.1, 1.1);
      lout->Add(hPi0ShowerE2Compare_REG);
      hPi0ShowerOACompare_REG = new TH2D("y012hPi0ShowerOACompare_REG",";Rec. - Truth (Pre Fit);Rec. - Truth (Post Fit)", 40, -50, 50, 40, -50, 50);
      lout->Add(hPi0ShowerOACompare_REG);
      hPi0MomCompare_REG = new TH2D("y013hPi0MomCompare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1, 1, 40, -1, 1);
      lout->Add(hPi0MomCompare_REG);

      hPi0MomCompare_REGPaper = new TH2D("y013hPi0MomCompare_REG_Paper",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1, 1, 40, -1, 1);
      lout->Add(hPi0MomCompare_REGPaper);

      hPi0MomCompare_REGStand = new TH2D("y013hPi0MomCompare_REG_Stand",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1, 1, 40, -1, 1);
      lout->Add(hPi0MomCompare_REGStand);

      hTrackPurityVSnHits = new TH2D("z001hTrackPurityVSnHits_REG",";Number of Hits -total;Track Purity", 25, 0, 1000, 20, 0, 1.2);
      lout->Add(hTrackPurityVSnHits);
      hTrackCompletenessVSnHits = new TH2D("z001hTrackCompletenessVSnHits_REG",";Number of Hits -total;Track Completeness", 25, 0, 1000, 20, 0, 1.2);
      lout->Add(hTrackCompletenessVSnHits);

      hTrackPurityVSnHitsColl = new TH2D("z002hTrackPurityVSnHitsColl_REG",";Number of Hits -coll.;Track Purity", 20, 0, 500, 20, 0, 1.2);
      lout->Add(hTrackPurityVSnHitsColl);
      hTrackCompletenessVSnHitsColl = new TH2D("z002hTrackCompletenessVSnHitsColl_REG",";Number of Hits -coll.;Track Completeness", 20, 0, 500, 20, 0, 1.2);
      lout->Add(hTrackCompletenessVSnHitsColl);

      hShowerPurityVSnHits = new TH2D("z003hShowerPurityVSnHits_REG",";Number of Hits -total;Shower Purity", 25, 0, 1000, 20, 0, 1.2);
      lout->Add(hShowerPurityVSnHits);
      hShowerCompletenessVSnHits = new TH2D("z003hShowerCompletenessVSnHits_REG",";Number of Hits -total;Shower Completeness", 25, 0, 1000, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSnHits);

      hShowerPurityVSnHitsColl = new TH2D("z004hShowerPurityVSnHitsColl_REG",";Number of Hits -coll.;Shower Purity", 20, 0, 400, 20, 0, 1.2);
      lout->Add(hShowerPurityVSnHitsColl);
      hShowerCompletenessVSnHitsColl = new TH2D("z004hShowerCompletenessVSnHitsColl_REG",";Number of Hits -coll.;Shower Completeness", 20, 0, 400, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSnHitsColl);

      hTrackCompleteness = new TH2D("z005hTrackCompleteness_STK","", 20, 0, 1, nparType, parTypemin, parTypemax );
      lout->Add(hTrackCompleteness);
      hShowerCompleteness = new TH2D("z006hShowerCompleteness_STK","", 20, 0, 1, nparType, parTypemin, parTypemax );
      lout->Add(hShowerCompleteness);

      hShowerPurityVSIP = new TH2D("z007hShowerPurityVSIP_REG",";Impact Parameter (cm);Shower Purity", 25, 0, 50, 20, 0, 1.2);
      lout->Add(hShowerPurityVSIP);
      hShowerCompletenessVSIP = new TH2D("z007hShowerCompletenessVSIP_REG",";Impact Parameter (cm);Shower Completeness", 25, 0, 50, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSIP);

      hShowerPurityVSEnergy = new TH2D("z009hShowerPurityVSEnergy_REG",";Shower Energy (GeV);Shower Purity", 25, 0, 0.8, 20, 0, 1.2);
      lout->Add(hShowerPurityVSEnergy);
      hShowerCompletenessVSEnergy = new TH2D("z009hShowerCompletenessVSEnergy_REG",";Shower Energy (GeV);Shower Completeness", 25, 0, 0.8, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSEnergy);

      hShowernHitsVSEnergy = new TH2D("z011hShowernHitsVSEnergy_REG",";Number of Hits -total;Shower Energy (GeV)", 50, 0, 1000, 50, 0, 1);
      lout->Add(hShowernHitsVSEnergy);


      //hTruthRecRatioProtonMom = new TH1D("zz001hTruthRecRatioProtonMom","", 60, 0.5, 1.5);
      hTruthRecRatioProtonMom = new TH1D("zz001hTruthRecRatioProtonMom","", 200, -0.6, 0.6);
      lout->Add(hTruthRecRatioProtonMom);

      hdalphat_REG = new TH2D("dd001hdalphat_REG",";Rec. #delta#alpha_{T} (deg);Truth #delta#alpha_{T} (deg)", 8, 0, 180, 8, 0, 180); 
      lout->Add(hdalphat_REG);
      hdphit_REG = new TH2D("dd002hdphit_REG","", 8, 0, 180, 8, 0, 180); 
      lout->Add(hdphit_REG);
      hdpt_REG = new TH2D("dd003hdpt_REG","", 8, 0, 1, 8, 0, 1); 
      lout->Add(hdpt_REG);
      hpn_REG = new TH2D("dd004hpn_REG",";Rec. p_{n} (GeV/c);Truth p_{n} (GeV/c)", 8, 0, 1, 8, 0, 1); 
      lout->Add(hpn_REG);

      hdalphat_RES = new TH2D("dd005hdalphat_RES","", 8, 0, 180, 8, -30, 30); 
      lout->Add(hdalphat_RES);
      hdphit_RES = new TH2D("dd006hdphit_RES","", 8, 0, 180, 8, -30, 30); 
      lout->Add(hdphit_RES);
      hdpt_RES = new TH2D("dd007hdpt_RES","", 8, 0, 1, 8, -1, 1); 
      lout->Add(hdpt_RES);
      hpn_RES = new TH2D("dd008hpn_RES","", 8, 0, 1, 8, -1, 1); 
      lout->Add(hpn_RES);

      hLDShower_PreFit = new TH1D("d001hLDShower_PreFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hLDShower_PreFit);
      hLDShower_PostFit = new TH1D("d001hLDShower_PostFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hLDShower_PostFit);

      hSLShower_PreFit = new TH1D("d002hSLShower_PreFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hSLShower_PreFit);
      hSLShower_PostFit = new TH1D("d002hSLShower_PostFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hSLShower_PostFit);

      hOAShower_PreFit = new TH1D("d003hOAShower_PreFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hOAShower_PreFit);
      hOAShower_PostFit = new TH1D("d003hOAShower_PostFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hOAShower_PostFit);

      hPi0EnergyE1E2_PreFit = new TH1D("d004hPi0EnergyE1E2_PreFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hPi0EnergyE1E2_PreFit);
      hPi0EnergyE1E2_PostFit = new TH1D("d004hPi0EnergyE1E2_PostFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hPi0EnergyE1E2_PostFit);

      hPi0EnergyE1OA_PreFit = new TH1D("d005hPi0EnergyE1OA_PreFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hPi0EnergyE1OA_PreFit);
      hPi0EnergyE1OA_PostFit = new TH1D("d005hPi0EnergyE1OA_PostFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hPi0EnergyE1OA_PostFit);

      hPi0EnergyAsym_PreFit = new TH1D("d006hPi0EnergyAsym_PreFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hPi0EnergyAsym_PreFit);
      hPi0EnergyAsym_PostFit = new TH1D("d006hPi0EnergyAsym_PostFit", ";Rec./Truth - 1;Candidates", 35, -1, 1);
      lout->Add(hPi0EnergyAsym_PostFit);

      hShowerE1Compare_REG = new TH2D("d006hShowerE1Compare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 35, -1.1, 1.1, 35, -1.1, 1.1);
      lout->Add(hShowerE1Compare_REG);
      hShowerE2Compare_REG = new TH2D("d007hShowerE2Compare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 35, -1.1, 1.1, 35, -1.1, 1.1);
      lout->Add(hShowerE2Compare_REG);
      hShowerOACompare_REG = new TH2D("d008hShowerOACompare_REG",";Rec. - Truth (Pre Fit);Rec. - Truth (Post Fit)", 35, -50, 50, 35, -50, 50);
      lout->Add(hShowerOACompare_REG);
      hPi0EnergyComparem1_REG = new TH2D("d009hPi0EnergyComparem1_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 35, -1, 1, 35, -1, 1);
      lout->Add(hPi0EnergyComparem1_REG);
      hPi0EnergyComparem2_REG = new TH2D("d010hPi0EnergyComparem2_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 20, -1, 1, 20, -1, 1);
      lout->Add(hPi0EnergyComparem2_REG);

    }
  }// End of IniHist

} // End of namespace
#endif
