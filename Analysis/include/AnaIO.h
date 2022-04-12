#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "PlotUtils.h"
#include "AnaFunctions.h"
// Kinematic fitting codes
#include "../../Fitting/UserKF/KF/src/UserKF.cxx"

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

  // Beam new variables with possible SCE correction
  Double_t        reco_beam_calo_startX;
  Double_t        reco_beam_calo_startY;
  Double_t        reco_beam_calo_startZ;
  Double_t        reco_beam_calo_endX;
  Double_t        reco_beam_calo_endY;
  Double_t        reco_beam_calo_endZ;
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

  Int_t           reco_beam_PFP_nHits;;
  Double_t        reco_beam_PFP_emScore_collection;
  Double_t        reco_beam_PFP_emScore_collection_weight_by_charge;
  Double_t        reco_beam_vertex_michel_score_weight_by_charge;

  Int_t           reco_beam_true_byHits_PDG;
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
  
  //====================== Reco Daughter ======================//

  vector<int>             *reco_daughter_PFP_ID = 0x0;
  vector<int>             *reco_daughter_PFP_nHits = 0x0; 
  vector<int>             *reco_daughter_PFP_nHits_collection = 0x0; 
  vector<double>          *reco_daughter_PFP_trackScore_collection = 0x0;
  vector<double>          *reco_daughter_PFP_emScore = 0x0;
  vector<double>          *reco_daughter_PFP_emScore_collection = 0x0;
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
  // Class A - beam cut related
  TH1I * hCutBeamPDGPass = 0x0;
  TH1I * hCutPandoraSlicePass = 0x0;
  TH1I * hCutCaloSizePass = 0x0;
  TH1I * hCutBeamQualityPass = 0x0;
  TH1I * hCutAPA3EndZPass = 0x0;
  TH1I * hCutMichelScorePass = 0x0;
  TH1I * hCutProtonChi2Pass = 0x0;

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

  TH2D * hBeamemScore = 0x0;
  TH2D * hBeamemScore_wbc = 0x0;
  TH2D * hBeamDaughterNumbers = 0x0;
  TH2D * hBeamGrandDaughterNumbers = 0x0;


  // Class B - reconstructed beam/FS particles
  TH2D * hRecBeamTheta = 0x0;
  TH2D * hRecBeamMomentum = 0x0;
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
  TH2D * hRecPi0MassRaw_OVERLAY = 0x0;
  TH2D * hRecPi0Momentum_OVERLAY = 0x0;
  TH2D * hRecPi0MomentumRaw_OVERLAY = 0x0;
  TH2D * hRecPi0OA_OVERLAY = 0x0;
  TH2D * hRecPi0OARaw_OVERLAY = 0x0;
  TH2D * hRecPi0Energy_OVERLAY = 0x0;

  TH2D * hRecPi0Theta_OVERLAY = 0x0;
  TH2D * hRecPi0Phi_OVERLAY = 0x0;

  TH2D * hRecPi0ThetaRaw_OVERLAY = 0x0;
  TH2D * hRecPi0PhiRaw_OVERLAY = 0x0;

  TH2D * hRecPi0Energy_OVERLAY_After = 0x0;

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

  // FS particle cut
  TH1I * hCutDaughterPandoraShowerPass = 0x0;
  TH1I * hCutDaughterShowerScorePass = 0x0;
  TH1I * hCutDaughterShowernHitsPass = 0x0;
  TH1I * hCutDaughterShowerNonEmptyEPass = 0x0;
  TH1I * hCutDaughterShowerDistPass = 0x0;
  TH1I * hCutDaughterShowerIPPass = 0x0;

  TH2D * hCutTracknHits = 0x0;
  TH2D * hCutTrackScore = 0x0;
  TH2D * hCutlastTME = 0x0;
  TH2D * hCutChi2NDF = 0x0;
  TH2D * hCutemScore = 0x0;
  TH2D * hCutmichelScore = 0x0;
  TH2D * hCutShowerDist = 0x0;
  TH2D * hCutShowerIP = 0x0;
  TH2D * hCutShowerEnergy = 0x0;

  TH2D * hCutnproton = 0x0;
  TH2D * hCutnpiplus = 0x0;
  TH2D * hCutnshower = 0x0;
  TH2D * hCutnmichel = 0x0;
  TH2D * hCutnpi0 = 0x0;

  TH2D * hCutShowerStartX = 0x0;
  TH2D * hCutShowerStartY = 0x0;
  TH2D * hCutShowerStartZ = 0x0;
  TH2D * hShowernHitsColl = 0x0;

  TH2D * hNhitsCollection = 0x0;

  TH2D * htrackScoreCollection = 0x0;
  TH2D * hemScoreCollection = 0x0;
  TH2D * hmichelScoreCollection = 0x0;

  TH2D * hShowerAngleOffset = 0x0;
  
  //====================== Reco (Data only)======================//
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
  
  //====================== Truth (MC only)======================//
  // Declare variables
  bool Signal = false;
  Int_t           true_beam_PDG;
  string          *true_beam_endProcess;
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

  vector<double> *true_beam_daughter_startPx = 0x0;
  vector<double> *true_beam_daughter_startPy = 0x0;
  vector<double> *true_beam_daughter_startPz = 0x0;
  vector<double> *true_beam_daughter_startX = 0x0;
  vector<double> *true_beam_daughter_startY = 0x0;
  vector<double> *true_beam_daughter_startZ = 0x0;
  vector<double> *true_beam_daughter_endX = 0x0;
  vector<double> *true_beam_daughter_endY = 0x0;
  vector<double> *true_beam_daughter_endZ = 0x0;

  vector<int> *true_beam_daughter_ID = 0x0;  
  vector<int> *true_beam_daughter_PDG = 0x0;
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
 
  vector<int> *reco_daughter_PFP_true_byHits_PDG = 0x0;
  vector<int> *reco_daughter_PFP_true_byHits_ID = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_startPx = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_startPy = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_startPz = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_purity = 0x0;
  vector<double>  *reco_daughter_PFP_true_byHits_completeness = 0x0;
  vector<double>  *reco_daughter_PFP_true_byE_completeness = 0x0;

  // Declare histograms
  TH1I * hTruthBeamType = 0x0;
  TH1I * hTruthSignal = 0x0;
  TH1I * hTruthFSParticleNumber = 0x0;  
  TH1I * hTruthFSParticleType = 0x0;
  TH1I * hTruthFSPi0Number = 0x0;
  TH1I * hTruthFSMultiPi0 = 0x0;
  TH1D * hTruthLeadingProtonP = 0x0;
  TH1D * hTruthSubLeadingProtonP = 0x0;
  TH1D * hTruthLeadingPiZeroP = 0x0;
  TH1D * hTruthLeadingPiZeroE = 0x0;
  TH1D * hTruthSubLeadingPiZeroP = 0x0;
  TH1D * hTruthGammaMaxE = 0x0;
  TH1D * hTruthPi0DecayParticleNumber = 0x0;
  TH1D * hTruthLeadingPi0GammaP = 0x0;
  TH1D * hTruthSubLeadingPi0GammaP = 0x0;
  TH1D * hTruthPi0OA = 0x0;
  TH1D * hTruthLeadingPi0GammaOA = 0x0;
  TH1D * hTruthSubLeadingPi0GammaOA = 0x0;

  TH1D * hTruthLeadingPiZeroGammaDist = 0x0;
  TH1D * hTruthSubLeadingPiZeroGammaDist = 0x0;

  TH2D * hTruthPi0GammaEnergy = 0x0;
  TH1D * hTruthRarePi0GammaP = 0x0;
  TH1D * hTruthRarePi0ElectronP = 0x0;
  TH1D * hTruthRarePi0PositronP = 0x0;

  TH1I * hTruthSignalFSParticleNumber = 0x0;  
  TH1I * hTruthSignalFSParticleType = 0x0;

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
  TH2D * hBeamThetaRes = 0x0;
  TH2D * hBeamMomentumRes = 0x0;
  TH2D * hProtonThetaRes = 0x0;
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
  TH2D * hPi0MomentumResFit = 0x0;
  TH2D * hPi0MassResFit = 0x0;

  TH2D * hPi0ThetaRes = 0x0;
  TH2D * hPi0ThetaResRaw = 0x0;
  TH2D * hPi0PhiRes = 0x0;
  TH2D * hPi0PhiResRaw = 0x0;

  TH1D * hMatchedTruthPi0Momentum = 0x0;
  TH1D * hMatchedTruthldShowerTheta = 0x0;
  TH1D * hMatchedTruthslShowerTheta = 0x0;
  TH1D * hMatchedTruthPi0Mass = 0x0;
  TH1D * hMatchedTruthldShowerEnergy = 0x0;
  TH1D * hMatchedTruthslShowerEnergy = 0x0;

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

  TH2D * hShowerE1Compare_REG = 0x0;
  TH2D * hShowerE2Compare_REG = 0x0;
  TH2D * hShowerOACompare_REG = 0x0;
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
    tree->SetBranchAddress("reco_beam_Chi2_proton",&reco_beam_Chi2_proton);
    tree->SetBranchAddress("reco_beam_Chi2_ndof",&reco_beam_Chi2_ndof);

    tree->SetBranchAddress("reco_beam_calo_startX", &reco_beam_calo_startX);
    tree->SetBranchAddress("reco_beam_calo_startY", &reco_beam_calo_startY);
    tree->SetBranchAddress("reco_beam_calo_startZ", &reco_beam_calo_startZ);
    tree->SetBranchAddress("reco_beam_calo_endX", &reco_beam_calo_endX);
    tree->SetBranchAddress("reco_beam_calo_endY", &reco_beam_calo_endY);
    tree->SetBranchAddress("reco_beam_calo_endZ", &reco_beam_calo_endZ);
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
    tree->SetBranchAddress("reco_beam_PFP_emScore_collection",&reco_beam_PFP_emScore_collection);
    tree->SetBranchAddress("reco_beam_PFP_emScore_collection_weight_by_charge",&reco_beam_PFP_emScore_collection_weight_by_charge);
    tree->SetBranchAddress("reco_beam_vertex_michel_score_weight_by_charge",&reco_beam_vertex_michel_score_weight_by_charge);

    tree->SetBranchAddress("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG);
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


    tree->SetBranchAddress("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
    tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
    tree->SetBranchAddress("reco_daughter_PFP_nHits_collection", &reco_daughter_PFP_nHits_collection);
    tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
    tree->SetBranchAddress("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore);
    tree->SetBranchAddress("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
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

    //====================== Reco (Data only)======================//
    tree->SetBranchAddress("beam_inst_PDG_candidates", &beam_inst_PDG_candidates);
    tree->SetBranchAddress("beam_inst_X", &beam_inst_X);
    tree->SetBranchAddress("beam_inst_Y", &beam_inst_Y);
    tree->SetBranchAddress("beam_inst_dirX", &beam_inst_dirX);
    tree->SetBranchAddress("beam_inst_dirY", &beam_inst_dirY);
    tree->SetBranchAddress("beam_inst_dirZ", &beam_inst_dirZ);
    tree->SetBranchAddress("beam_inst_trigger", &beam_inst_trigger); 
    tree->SetBranchAddress("beam_inst_nMomenta", &beam_inst_nMomenta);
    tree->SetBranchAddress("beam_inst_nTracks", &beam_inst_nTracks);

    //====================== Truth (MC only)======================//
    tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
    tree->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess);
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

    tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
    tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
    tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

    tree->SetBranchAddress("true_beam_daughter_startX", &true_beam_daughter_startX);
    tree->SetBranchAddress("true_beam_daughter_startY", &true_beam_daughter_startY);
    tree->SetBranchAddress("true_beam_daughter_startZ", &true_beam_daughter_startZ);
    tree->SetBranchAddress("true_beam_daughter_endX", &true_beam_daughter_endX);
    tree->SetBranchAddress("true_beam_daughter_endY", &true_beam_daughter_endY);
    tree->SetBranchAddress("true_beam_daughter_endZ", &true_beam_daughter_endZ);
   
    tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID); 
    tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
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
 
    tree->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG);
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
    const int nparType = 11;
    const double parTypemin = 0.5;
    const double parTypemax = 11.5;

    const int nevtType = 3;
    const double evtTypemin = -0.5;
    const double evtTypemax = 2.5;

    const int ncounter = 10;
    const double countermin = -0.5;
    const double countermax = 9.5;

    const int nPass = 2;
    const double Passmin = -0.5;
    const double Passmax = 1.5;


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
    hRecBeamStartX = new TH1D("a008hRecBeamStartX_FCN",";StartX (cm);Area Normalised",100, -80, 20);
    lout->Add(hRecBeamStartX);
    hRecBeamStartY = new TH1D("a009hRecBeamStartY_FCN",";StartY (cm);Area Normalised",100, 350, 500);
    lout->Add(hRecBeamStartY);
    hRecBeamStartZ = new TH1D("a010hRecBeamStartZ_FCN",";StartZ (cm);Area Normalised",100, -5, 10);
    lout->Add(hRecBeamStartZ);
    hRecBeamThetaX = new TH1D("a011hRecBeamThetaX_FCN",";Angle #theta_{x} (deg);Area Normalised",120, 60, 140);
    lout->Add(hRecBeamThetaX);
    hRecBeamThetaY = new TH1D("a012hRecBeamThetaY_FCN",";Angle #theta_{y} (deg);Area Normalised",120, 60, 140);
    lout->Add(hRecBeamThetaY);
    hRecBeamThetaZ = new TH1D("a013hRecBeamThetaZ_FCN",";Angle #theta_{z} (deg);Area Normalised",100, 5, 45);
    lout->Add(hRecBeamThetaZ);
    hCutBeamDeltaXYSigma = new TH1D("a014hCutBeamDeltaXYSigma","",100, 0, 20);
    lout->Add(hCutBeamDeltaXYSigma);
    hCutBeamDeltaZSigma = new TH1D("a015hCutBeamDeltaZSigma","",100, -10, 10);
    lout->Add(hCutBeamDeltaZSigma);
    hCutBeamCosTheta = new TH1D("a016hCutBeamCosTheta","",100, 0.9, 1);
    lout->Add(hCutBeamCosTheta);
    hCutBeamMichelScoreNhits = new TH2D("a017hCutBeamMichelScoreNhits_STK","",100, 0, 1, 8, 0.5, 8.5);
    lout->Add(hCutBeamMichelScoreNhits);
    hCutBeamProtonChi2DOF = new TH2D("a018hCutBeamProtonChi2DOF_STK","", 100, 0, 350, 8, 0.5, 8.5); 
    lout->Add(hCutBeamProtonChi2DOF);
    hCutBeamProtonChi2DOFType = new TH1D("a019hCutBeamProtonChi2DOFType","",23, -0.5, 22.5);
    lout->Add(hCutBeamProtonChi2DOFType);
    hCutBeamQualityXY = new TH2D("a020hCutBeamQualityXY_STK","", 100, 0, 20, 8, 0.5, 8.5); 
    lout->Add(hCutBeamQualityXY);
    hCutBeamQualityX = new TH2D("a021hCutBeamQualityX_STK","", 100, -10, 10, 8, 0.5, 8.5); 
    lout->Add(hCutBeamQualityX);
    hCutBeamQualityY = new TH2D("a022hCutBeamQualityY_STK","", 100, -10, 10, 8, 0.5, 8.5); 
    lout->Add(hCutBeamQualityY);
    hCutBeamQualityZ = new TH2D("a023hCutBeamQualityZ_STK","", 100, -10, 10, 8, 0.5, 8.5); 
    lout->Add(hCutBeamQualityZ);
    hCutBeamQualityTheta = new TH2D("a024hCutBeamQualityTheta_STK","", 100, 0.9, 1, 8, 0.5, 8.5); 
    lout->Add(hCutBeamQualityTheta);

    hCutAPA3EndZ = new TH2D("a025hCutAPA3EndZ_STK",";End Z (cm);Candidates", 100, 0, 500, 8, 0.5, 8.5); 
    lout->Add(hCutAPA3EndZ);

    hBeamemScore = new TH2D("a025hBeamemScore_STK",";Beam EM Shower Score;Candidates", 50, 0, 1, 8, 0.5, 8.5); 
    lout->Add(hBeamemScore);

    hBeamemScore_wbc = new TH2D("a026hBeamemScore_wbc_STK",";Beam EM Shower Score;Candidates", 50, 0, 1, 8, 0.5, 8.5); 
    lout->Add(hBeamemScore_wbc);

    hBeamDaughterNumbers = new TH2D("a027hBeamDaughterNumbers_STK",";Number of Daughter Particles;Candidates", 15, 0, 15, 8, 0.5, 8.5); 
    lout->Add(hBeamDaughterNumbers);
    hBeamGrandDaughterNumbers = new TH2D("a028hBeamGrandDaughterNumbers_STK",";Number of Grand Daughter Particles;Candidates", 15, 0, 15, 8, 0.5, 8.5); 
    lout->Add(hBeamGrandDaughterNumbers);

    // Class B - reconstructed beam/FS particles
    hRecBeamTheta = new TH2D("b001hRecBeamTheta_STK",";Beam #theta (deg);Candidates", 80 , 0, 60, 3, -0.5, 2.5); 
    lout->Add(hRecBeamTheta);
    hRecBeamMomentum = new TH2D("b002hRecBeamMomentum_STK",";Beam Momentum (GeV/c);Candidates", 50, 0, 2, 3, -0.5, 2.5); 
    lout->Add(hRecBeamMomentum);
    hRecProtonTheta = new TH2D("b003hRecProtonTheta_STK",";Proton #theta (deg);Candidates",  20, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecProtonTheta);
    hRecProtonMomentum = new TH2D("b004hRecProtonMomentum_STK",";Proton Momentum (GeV/c);Candidates",  20, 0.1, 1.3, nparType, parTypemin, parTypemax); 
    lout->Add(hRecProtonMomentum);
    hRecPiPlusTheta = new TH2D("b005hRecPiPlusTheta_STK",";#pi^{+} #theta (deg);Candidates", 30, 0, 180, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusTheta);
    hRecPiPlusMomentum = new TH2D("b006hRecPiPlusMomentum_STK",";#pi^{+} Momentum (GeV/c);Candidates",  20, 0, 1.2, nparType, parTypemin, parTypemax); 
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
    hRecPi0Momentum = new TH2D("b020hRecPi0Momentum_STK", ";#pi^{0} Momentum (GeV/c);Candidates", 20, 0, 1, 3, -0.5, 3.5);
    lout->Add(hRecPi0Momentum);
    hRecPi0MomentumRaw = new TH2D("b021hRecPi0Momentum_STK_RAW", ";#pi^{0} Momentum Raw (GeV/c);Candidates", 20, 0, 1, 3, -0.5, 3.5);
    lout->Add(hRecPi0MomentumRaw);

    hPi0MassLowE1 = new TH1D("b020hPi0MassLowE1","", 20, 0, 0.5);
    lout->Add(hPi0MassLowE1);
    hPi0MassHighE1 = new TH1D("b021hPi0MassHighE1","", 20, 0, 0.5);
    lout->Add(hPi0MassHighE1);

    hRecPi0Mass_OVERLAY = new TH2D("b022hRecPi0Mass_COMPOSE",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Mass_OVERLAY);

    hRecPi0Momentum_OVERLAY = new TH2D("b023hRecPi0Momentum_COMPOSE",";#pi^{0} Momentum (GeV/c);Candidates", 20, 0, 1, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Momentum_OVERLAY);

    hRecPi0OA_OVERLAY = new TH2D("b024hRecPi0OA_COMPOSE",";Shower Opening Angle (deg);Candidates", 20, 0, 180, 3, -0.5, 3.5); 
    lout->Add(hRecPi0OA_OVERLAY);

    hRecPi0Energy_OVERLAY = new TH2D("b023hRecPi0Energy_COMPOSE",";#pi^{0} Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Energy_OVERLAY);

    hRecPi0ShowerSep_OVERLAY = new TH2D("b024hRecPi0ShowerSep_COMPOSE",";Shower Separation (cm);Candidates", 20, 0, 150, 3, -0.5, 3.5);
    lout->Add(hRecPi0ShowerSep_OVERLAY);

    hRecPi0MassRaw_OVERLAY = new TH2D("b025hRecPi0MassRaw_COMPOSE",";#pi^{0} Mass (GeV/c^{2});Candidates", 20, 0, 0.5, 3, -0.5, 3.5); 
    lout->Add(hRecPi0MassRaw_OVERLAY);

    hRecPi0MomentumRaw_OVERLAY = new TH2D("b026hRecPi0MomentumRaw_COMPOSE",";#pi^{0} Momentum (GeV/c);Candidates", 20, 0, 1, 3, -0.5, 3.5); 
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


    hRecPi0Energy_OVERLAY_After = new TH2D("b0233hRecPi0Energy_COMPOSE_After",";#pi^{0} Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 3.5); 
    lout->Add(hRecPi0Energy_OVERLAY_After);

    // Class C - event topoplogy cut related
    hCutDaughterPandoraShowerPass = new TH1I("c000hCutDaughterPandoraShowerPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterPandoraShowerPass);
    hCutDaughterShowerScorePass = new TH1I("c001hCutDaughterShowerScorePass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerScorePass);
    hCutDaughterShowernHitsPass = new TH1I("c002hCutDaughterShowernHitsPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowernHitsPass);
    hCutDaughterShowerNonEmptyEPass = new TH1I("c003hCutDaughterShowerNonEmptyEPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerNonEmptyEPass);
    hCutDaughterShowerDistPass = new TH1I("c004hCutDaughterShowerDistPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerDistPass);
    hCutDaughterShowerIPPass = new TH1I("c005hCutDaughterShowerIPPass",";Fail/Pass;Candidates", nPass, Passmin, Passmax); 
    lout->Add(hCutDaughterShowerIPPass);

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
    hCutmichelScore = new TH2D("c006hCutmichelScore_STK",";Michel Score;Candidates",50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutmichelScore);
    hCutShowerDist = new TH2D("c007hCutShowerDist_STK",";Shower Distance (cm);Candidates", 31, 0, 93, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerDist);
    hCutShowerIP = new TH2D("c008hCutShowerIP_STK",";Shower Impact Parameter (cm);Candidates", 30, 0, 30, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerIP);
    hCutShowerEnergy = new TH2D("c009hCutShowerEnergy_STK",";Shower Energy (GeV);Candidates", 50, 0, 0.5, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerEnergy);
    hCutnproton = new TH2D("c009hCutnproton_STK",";Number of Protons;Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnproton);
    hCutnshower = new TH2D("c010hCutnshower_STK",";Number of Showers;Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnshower);
    hCutnpiplus = new TH2D("c011hCutnpiplus_STK",";Number of Pions (#pi^{+});Candidates",  ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnpiplus);
    hCutnmichel = new TH2D("c012hCutnmichel_STK",";Number of Michel Electron;Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnmichel);
    hNhitsCollection = new TH2D("c013NhitsCollection_COMPOSE_LOG",";Number of Hits -coll.;Candidates", 50, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hNhitsCollection);
    hCutnpi0 = new TH2D("c012hCutnpi0_STK",";Number of #pi^{0};Candidates", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnpi0);

    hCutShowerStartX = new TH2D("c015hCutShowerStartX_STK",";StartX (cm);Candidates", 30, -100, 20, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerStartX);
    hCutShowerStartY = new TH2D("c016hCutShowerStartY_STK",";StartY (cm);Candidates", 30, 340, 460, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerStartY);
    hCutShowerStartZ = new TH2D("c017hCutShowerStartZ_STK",";StartZ (cm);Candidates", 30, 0, 240, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerStartZ);
    hShowernHitsColl = new TH2D("c018hShowernHitsColl_STK",";Number of Hits -coll.;Candidates", 50, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hShowernHitsColl);

    htrackScoreCollection = new TH2D("c019htrackScoreCollection_STK",";Track Score -coll.;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(htrackScoreCollection);
    hemScoreCollection = new TH2D("c019hemScoreCollection_STK",";EM Score -coll.;Candidates", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hemScoreCollection);
    hmichelScoreCollection = new TH2D("c019hmichelScoreCollection_STK",";Michel Score -coll.;Candidates", 50, 0, 0.3, nparType, parTypemin, parTypemax); 
    lout->Add(hmichelScoreCollection);

    hShowerAngleOffset = new TH2D("c020hShowerAngleOffset_STK",";Shower Angle Offset;Candidates", 30, 0, 180, nparType, parTypemin, parTypemax); 
    lout->Add(hShowerAngleOffset);

    hKFPassRate = new TH1D("k001hKFPassRate","",  nPass, Passmin, Passmax); 
    lout->Add(hKFPassRate); 
    hPi0Total = new TH1D("k002hPi0Total","",  30, 0, 3); 
    lout->Add(hPi0Total);
    hPi0Selected = new TH1D("k003hPi0Selected","",  30, 0, 3);  
    lout->Add(hPi0Selected);

    // Class D - reconstructed TKI variables
    hRecdalphat = new TH2D("d001Recdalphat_STK","", 9, 0, 180,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdalphat);
    hRecdphit = new TH2D("d002Recphit_STK","", 9, 0, 180,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdphit);
    hRecdpt = new TH2D("d003Recdpt_STK","", 10, 0, 1,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdpt);
    hRecpn = new TH2D("d004Recpn_STK","", 10, 0, 1,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecpn);

    hRecdalphat_truth = new TH2D("d005Recdalphat_STK_TM","", 9, 0, 180,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdalphat_truth);
    hRecdphit_truth = new TH2D("d006Recphit_STK_TM","", 9, 0, 180,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdphit_truth);
    hRecdpt_truth = new TH2D("d007Recdpt_STK_TM","", 10, 0, 1,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdpt_truth);
    hRecpn_truth = new TH2D("d008Recpn_STK_TM","", 10, 0, 1,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecpn_truth);



    //====================== Truth (MC only)======================//
    if(kMC){
      // Class E - truth beam/FS particles
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
      hTruthLeadingProtonP = new TH1D("e006hTruthLeadingProtonP", ";Truth Leading Proton Mom.(GeV/c);Candidates", 20, 0, 2);
      lout->Add(hTruthLeadingProtonP);
      hTruthSubLeadingProtonP = new TH1D("e007hTruthSubLeadingProtonP", ";Truth SubLeading Proton Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingProtonP);
      hTruthLeadingPiZeroP = new TH1D("e008hTruthLeadingPiZeroP", ";Truth Leading #pi^{0} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPiZeroP);
      hTruthLeadingPiZeroE = new TH1D("e008hTruthLeadingPiZeroE", ";Truth Leading #pi^{0} Energy(GeV);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPiZeroE);
      hTruthSubLeadingPiZeroP = new TH1D("e009hTruthSubLeadingPiZeroP", ";Truth SubLeading #pi^{0} Mom.(GeV/c);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingPiZeroP);
      hTruthGammaMaxE = new TH1D("e010hTruthGammaMaxE", ";Truth Gamma Max. Energy (GeV);Candidates", 20, 0, 2);
      lout->Add(hTruthGammaMaxE);
      hTruthPi0DecayParticleNumber = new TH1D("e011hTruthPi0DecayParticleNumber", ";Truth #pi^{0} Decay Particle Number;Candidates", 10, -0.5, 9.5);
      lout->Add(hTruthPi0DecayParticleNumber);
      hTruthLeadingPi0GammaP = new TH1D("e012hTruthLeadingPi0GammaP", ";Truth Leading #gamma Energy (GeV);Candidates", 20, 0, 1);
      lout->Add(hTruthLeadingPi0GammaP);
      hTruthSubLeadingPi0GammaP = new TH1D("e013hTruthSubLeadingPi0GammaP", ";Truth SubLeading #gamma Energy (GeV);Candidates", 20, 0, 1);
      lout->Add(hTruthSubLeadingPi0GammaP);
      hTruthPi0GammaEnergy = new TH2D("e014hTruthPi0GammaEnergy_OVERLAY", ";Truth #gamma Energy (GeV);Candidates", 20, 0, 1, 3, -0.5, 2.5);
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

      stkTruthDalphat = new THStack("f100stkTruthDalphat",";#delta#alpha_T (deg);Candidates"); 
      lout->Add(stkTruthDalphat);
      hTruthDalphat1p0n = new TH1D("f100hTruthDalphat1p0n",";#delta#alpha_T (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphat1p0n);
      hTruthDalphatNp0n = new TH1D("f100hTruthDalphatNp0n",";#delta#alpha_T (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphatNp0n);
      hTruthDalphat1pMn = new TH1D("f100hTruthDalphat1pMn",";#delta#alpha_T (deg);Candidates", 18, 0, 180); 
      lout->Add(hTruthDalphat1pMn);
      hTruthDalphatNpMn = new TH1D("f100hTruthDalphatNpMn",";#delta#alpha_T (deg);Candidates", 18, 0, 180); 
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
      stkTruthPn = new THStack("f103stkTruthPn",";#p_n (GeV/c);Candidates"); 
      lout->Add(stkTruthPn);
      hTruthPn1p0n = new TH1D("f103hTruthPn1p0n",";#p_n (GeV/c);Candidates", 30, 0, 1.2); 
      lout->Add(hTruthPn1p0n);
      hTruthPnNp0n = new TH1D("f103hTruthPnNp0n",";#p_n (GeV/c);Candidates", 30, 0, 1.2); 
      lout->Add(hTruthPnNp0n);
      hTruthPn1pMn = new TH1D("f103hTruthPn1pMn",";#p_n (GeV/c);Candidates", 30, 0, 1.2); 
      lout->Add(hTruthPn1pMn);
      hTruthPnNpMn = new TH1D("f103hTruthPnNpMn",";#p_n (GeV/c);Candidates", 30, 0, 1.2); 
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
      hBeamThetaRes = new TH2D("g001hBeamTheta_RES",";Beam Theta (deg);Rec. - Truth (deg)", 80 , 0, 60, 30, -20, 20);
      lout->Add(hBeamThetaRes);
      hBeamMomentumRes  = new TH2D("g002hBeamMomentum_RES",";Beam Momentum (GeV/c);Rec./Truth - 1", 50, 0.2, 1.6, 30, -0.5, 0.5);
      lout->Add(hBeamMomentumRes);
      hProtonThetaRes = new TH2D("g003hProtonThetaFitRes_RES",";Proton Theta (deg);Rec. - Truth (deg)", 40, 0, 180, 30, -20, 30); 
      lout->Add(hProtonThetaRes);
      hProtonThetaLabRes = new TH2D("g004hProtonThetaLabFitRes_RES",";Proton Theta Lab (deg);Rec. - Truth (deg)", 40, 0, 180, 30, -20, 30); 
      lout->Add(hProtonThetaLabRes);
      hProtonPhiLabRes = new TH2D("g005hProtonPhiLabFitRes_RES",";Proton Phi Lab (deg);Rec. - Truth (deg)", 80, -180, 180, 30, -20, 30); 
      lout->Add(hProtonPhiLabRes);
      hProtonMomentumRes = new TH2D("g006hProtonMomentumFitRes_RES",";Proton Momentum (GeV/c);Rec./Truth - 1",40, 0.3, 1.2, 20, -0.2, 0.2); 
      lout->Add(hProtonMomentumRes);
      hProtonMomentumRawRes = new TH2D("g007hProtonMomentumFitRes_RES_RAW",";Proton Momentum Raw (GeV/c);Rec./Truth - 1",40, 0.3, 1.2, 20, -0.2, 0.2); 
      lout->Add(hProtonMomentumRawRes);
      hPiPlusThetaRes = new TH2D("g008hPiPlusTheta_RES","", 15, 0, 180, 25, -20, 30); 
      lout->Add(hPiPlusThetaRes);
      hPiPlusMomentumRes = new TH2D("g009hPiPlusMomentum_RES","", 20, 0, 1.2, 20, -0.2, 0.2); 
      lout->Add(hPiPlusMomentumRes);
      hShowerThetaRes = new TH2D("g010hShowerThetaRes_RES",";Shower Theta (deg);Rec. - Truth (deg)", 15, 0, 180, 50, -20, 30);
      lout->Add(hShowerThetaRes);
      hShowerThetaResRaw = new TH2D("g011hShowerThetaRes_RES_RAW",";Shower Theta Raw (deg);Rec. - Truth (deg)", 15, 0, 180, 50, -20, 30);
      lout->Add(hShowerThetaResRaw);
      hShowerEnergyRes = new TH2D("g012hShowerEnergy_RES",";Truth Photon Energy (GeV);Rec./Truth - 1", 20, 0, 0.8, 50, -1.1, 1.1);
      lout->Add(hShowerEnergyRes);
      hShowerEnergyResRaw = new TH2D("g013hShowerEnergy_RES_RAW",";Truth Photon Energy(GeV);Rec./Truth - 1", 20, 0, 0.8, 50, -1.1, 1.1);
      lout->Add(hShowerEnergyResRaw);
      hShowerPhiRes = new TH2D("g014hShowerPhiRes_RES",";Shower Phi (deg);Rec. - Truth (deg)", 15, -180, 180, 50, -20, 30);
      lout->Add(hShowerPhiRes);
      hShowerPhiResRaw = new TH2D("g015hShowerPhiRes_RES_RAW",";Shower Phi Raw (deg);Rec. - Truth (deg)", 15, -180, 180, 50, -20, 30);
      lout->Add(hShowerPhiResRaw);

      hLeadingPhotonAngleRes = new TH2D("g016hLeadingPhotonAngleRes_RES",";Leading Photon Energy (GeV);Photon/Shower Angle (deg)", 20, 0, 1, 50, 0, 50);
      lout->Add(hLeadingPhotonAngleRes);

      hSubLeadingPhotonAngleRes = new TH2D("g017hSubLeadingPhotonAngleRes_RES",";SubLeading Photon Energy (GeV);Photon/Shower Angle (deg)", 20, 0, 0.5, 50, 0, 50);
      lout->Add(hSubLeadingPhotonAngleRes);

      hOpeningAngleRes = new TH2D("g018hOpeningAngleRes_RES",";#pi^{0} Energy (GeV);Photon/Shower Angle (deg)", 20, 0, 1, 50, -50, 50);
      lout->Add(hOpeningAngleRes);


      
      // Class H - Rec VS truth 
      hProtonMomentumRawRecVSTruth_REG  = new TH2D("h001hProtonMomentumRawRecVSTruth_REG_DIAG","",20, 0.2, 1.2, 20, 0.2, 1.2);
      lout->Add(hProtonMomentumRawRecVSTruth_REG);
      hProtonMomentumRecVSTruth_REG  = new TH2D("h002hProtonMomentumRecVSTruth_REG_DIAG","",20, 0.2, 1.2, 20, 0.2, 1.2);
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
      hMeanShowerE = new TH1D("l006hShowerEnergy_CorMean","", sizeof(ShowerBin)/sizeof(double)-1, ShowerBin);
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
      hPi0MomentumResRaw = new TH2D("m006Pi0MomentumFitRes_RES_RAW",";Truth #pi^{0} Energy Raw (GeV);Rec./Truth - 1", 25, 0, 1, 40, -1, 1);
      lout->Add(hPi0MomentumResRaw);
      hPi0MomentumRes = new TH2D("m007hPi0MomentumFitRes_RES",";Truth #pi^{0} Energy (GeV);Rec./Truth - 1", 25, 0, 1, 40, -1, 1);
      lout->Add(hPi0MomentumRes);
      hPi0MassResRaw = new TH2D("m008Pi0Mass_RES_RAW","", 20, 0, 0.5, 20, -0.5, 0.5);
      lout->Add(hPi0MassResRaw);
      hPi0MassRes = new TH2D("m009Pi0Mass_RES","", 20, 0, 0.5, 20, -0.5, 0.5); 
      lout->Add(hPi0MassRes);

      hMatchedTruthPi0Momentum = new TH1D("m010hMatchedTruthPi0Momentum","", 20, 0, 1);
      lout->Add(hMatchedTruthPi0Momentum);
      hMatchedTruthldShowerTheta = new TH1D("m011hMatchedTruthldShowerTheta","", 30, 0, 180);
      lout->Add(hMatchedTruthldShowerTheta);
      hMatchedTruthslShowerTheta = new TH1D("m012hMatchedTruthslShowerTheta","", 30, 0, 180);
      lout->Add(hMatchedTruthslShowerTheta);
      hMatchedTruthPi0Mass = new TH1D("m013hMatchedTruthPi0Mass","", 20, 0, 0.4);
      lout->Add(hMatchedTruthPi0Mass);

      hPi0ThetaResFit = new TH2D("m014hPi0ThetaFitRes_RES",";Truth #pi^{0} #theta Angle (deg);Rec.- Truth (deg)", 30, 0, 180, 40, -50, 50);
      lout->Add(hPi0ThetaResFit);
      hPi0PhiResFit = new TH2D("m015hPi0PhiFitRes_RES",";Truth #pi^{0} #phi Angle (deg);Rec.- Truth (deg)", 30, -180, 180, 40, -50, 50);
      lout->Add(hPi0PhiResFit);
      hPi0MomRes = new TH2D("m016hPi0MomFitRes_RES",";Truth #pi^{0} Energy (GeV);Rec.(post fit)/Truth - 1", 25, 0, 1, 40, -1, 1);
      lout->Add(hPi0MomRes);
      hPi0MomPreRes = new TH2D("m017hPi0MomPreFitRes_RES",";Truth #pi^{0} Energy (GeV);Rec.(pre fit)/Truth - 1", 25, 0, 1, 40, -1, 1);
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

      //hPi0MomentumResFit = new TH2D("i008Pi0Momentum_RES_FIT","", 20, 0, 1, 10, -0.5, 0.5);
      //lout->Add(hPi0MomentumResFit);
      //hPi0MassResFit = new TH2D("i009Pi0Mass_RES_FIT","", 20, 0, 0.5, 20, -0.5, 0.5);
      //lout->Add(hPi0MassResFit);
      //hShowerOpenAngleResFit = new TH2D("i010ShowerOpenAngle_RES_FIT","", 20, 0, 180, 20, -1.1, 1.1);
      //lout->Add(hShowerOpenAngleResFit);

      

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
      

      hShowerE1Compare_REG = new TH2D("y010hShowerE1Compare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1.1, 1.1, 40, -1.1, 1.1);
      lout->Add(hShowerE1Compare_REG);
      hShowerE2Compare_REG = new TH2D("y011hShowerE2Compare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1.1, 1.1, 40, -1.1, 1.1);
      lout->Add(hShowerE2Compare_REG);
      hShowerOACompare_REG = new TH2D("y012hShowerOACompare_REG",";Rec. - Truth (Pre Fit);Rec. - Truth (Post Fit)", 40, -50, 50, 40, -50, 50);
      lout->Add(hShowerOACompare_REG);
      hPi0MomCompare_REG = new TH2D("y013hPi0MomCompare_REG",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1, 1, 40, -1, 1);
      lout->Add(hPi0MomCompare_REG);

      hPi0MomCompare_REGPaper = new TH2D("y013hPi0MomCompare_REG_Paper",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1, 1, 40, -1, 1);
      lout->Add(hPi0MomCompare_REGPaper);

      hPi0MomCompare_REGStand = new TH2D("y013hPi0MomCompare_REG_Stand",";Rec./Truth - 1 (Pre Fit);Rec./Truth - 1 (Post Fit)", 40, -1, 1, 40, -1, 1);
      lout->Add(hPi0MomCompare_REGStand);

      hTrackPurityVSnHits = new TH2D("z001hTrackPurityVSnHits_REG",";Number of Hits -total;Track Purity", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hTrackPurityVSnHits);
      hTrackCompletenessVSnHits = new TH2D("z001hTrackCompletenessVSnHits_REG",";Number of Hits -total;Track Completeness", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hTrackCompletenessVSnHits);

      hTrackPurityVSnHitsColl = new TH2D("z002hTrackPurityVSnHitsColl_REG",";Number of Hits -coll.;Track Purity", 25, 0, 500, 20, 0, 1.2);
      lout->Add(hTrackPurityVSnHitsColl);
      hTrackCompletenessVSnHitsColl = new TH2D("z002hTrackCompletenessVSnHitsColl_REG",";Number of Hits -coll.;Track Completeness", 25, 0, 500, 20, 0, 1.2);
      lout->Add(hTrackCompletenessVSnHitsColl);

      hShowerPurityVSnHits = new TH2D("z003hShowerPurityVSnHits_REG",";Number of Hits -total;Shower Purity", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hShowerPurityVSnHits);
      hShowerCompletenessVSnHits = new TH2D("z003hShowerCompletenessVSnHits_REG",";Number of Hits -total;Shower Completeness", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSnHits);

      hShowerPurityVSnHitsColl = new TH2D("z004hShowerPurityVSnHitsColl_REG",";Number of Hits -coll.;Shower Purity", 20, 0, 400, 20, 0, 1.2);
      lout->Add(hShowerPurityVSnHitsColl);
      hShowerCompletenessVSnHitsColl = new TH2D("z004hShowerCompletenessVSnHitsColl_REG",";Number of Hits -coll.;Shower Completeness", 20, 0, 400, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSnHitsColl);

      hTrackCompleteness = new TH2D("z005hTrackCompleteness_STK","", 20, 0, 1, nparType, parTypemin, parTypemax );
      lout->Add(hTrackCompleteness);
      hShowerCompleteness = new TH2D("z006hShowerCompleteness_STK","", 20, 0, 1, nparType, parTypemin, parTypemax );
      lout->Add(hShowerCompleteness);

      hShowerPurityVSIP = new TH2D("z007hShowerPurityVSIP_REG",";Impact Parameter (cm);Shower Purity", 50, 0, 50, 20, 0, 1.2);
      lout->Add(hShowerPurityVSIP);
      hShowerCompletenessVSIP = new TH2D("z007hShowerCompletenessVSIP_REG",";Impact Parameter (cm);Shower Completeness", 50, 0, 50, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSIP);

      hShowerPurityVSEnergy = new TH2D("z009hShowerPurityVSEnergy_REG",";Shower Energy (GeV);Shower Purity", 50, 0, 1, 20, 0, 1.2);
      lout->Add(hShowerPurityVSEnergy);
      hShowerCompletenessVSEnergy = new TH2D("z009hShowerCompletenessVSEnergy_REG",";Shower Energy (GeV);Shower Completeness", 50, 0, 1, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSEnergy);

      hShowernHitsVSEnergy = new TH2D("z011hShowernHitsVSEnergy_REG",";Number of Hits -total;Shower Energy (GeV)", 50, 0, 1000, 50, 0, 1);
      lout->Add(hShowernHitsVSEnergy);


      //hTruthRecRatioProtonMom = new TH1D("zz001hTruthRecRatioProtonMom","", 60, 0.5, 1.5);
      hTruthRecRatioProtonMom = new TH1D("zz001hTruthRecRatioProtonMom","", 200, -0.6, 0.6);
      lout->Add(hTruthRecRatioProtonMom);

      hdalphat_REG = new TH2D("dd001hdalphat_REG","", 18, 0, 180, 9, 0, 180); 
      lout->Add(hdalphat_REG);
      hdphit_REG = new TH2D("dd002hdphit_REG","", 18, 0, 180, 9, 0, 180); 
      lout->Add(hdphit_REG);
      hdpt_REG = new TH2D("dd003hdpt_REG","", 20, 0, 0.8, 10, 0, 0.8); 
      lout->Add(hdpt_REG);
      hpn_REG = new TH2D("dd004hpn_REG","", 20, 0, 1, 10, 0, 1); 
      lout->Add(hpn_REG);

      hdalphat_RES = new TH2D("dd005hdalphat_RES","", 18, 0, 180, 10, -30, 30); 
      lout->Add(hdalphat_RES);
      hdphit_RES = new TH2D("dd006hdphit_RES","", 18, 0, 180, 10, -30, 30); 
      lout->Add(hdphit_RES);
      hdpt_RES = new TH2D("dd007hdpt_RES","", 20, 0, 0.8, 25, -1, 1); 
      lout->Add(hdpt_RES);
      hpn_RES = new TH2D("dd008hpn_RES","", 20, 0, 1, 25, -1, 1); 
      lout->Add(hpn_RES);


    }
  }// End of IniHist

} // End of namespace
#endif
