#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "PlotUtils.h"
#include "AnaFunctions.h"

namespace AnaIO
{
  //====================== Output Tree Variables ==================//
  double LeadingShowerEnergy;
  double SubLeadingShowerEnergy;
  double LeadingShowerEnergyRaw;
  double SubLeadingShowerEnergyRaw;
  double OpeningAngle;

  vector<double> * LeadingShowerEnergyUnitDir = 0x0;
  vector<double> * SubLeadingShowerEnergyUnitDir = 0x0;

  double LeadingShowerEnergyTruth;
  double SubLeadingShowerEnergyTruth;
  double OpeningAngleTruth;

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

  vector<double> * LeadingShowerEnergyUnitDirTruth = 0x0;
  vector<double> * SubLeadingShowerEnergyUnitDirTruth = 0x0;

  //====================== Reco (MC and Data)======================//
  // Declare variables
  int event = -999;
  int reco_beam_type = -999;
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
  vector<double>  *reco_beam_calibrated_dEdX_SCE = 0x0;

  vector<int>             *reco_daughter_PFP_ID = 0x0;
  vector<int>             *reco_daughter_PFP_nHits = 0x0; 
  vector<double>          *reco_daughter_PFP_trackScore_collection = 0x0;
  vector<double>          *reco_daughter_PFP_emScore_collection = 0x0;
  vector<double>          *reco_daughter_PFP_michelScore_collection =0x0;
  vector<int>             *reco_daughter_allTrack_ID = 0x0;
  vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX_SCE = 0x0;
  vector<double>          *reco_daughter_allTrack_Chi2_proton = 0x0;
  vector<int>             *reco_daughter_allTrack_Chi2_ndof = 0x0;
  vector<double>          *reco_daughter_allTrack_momByRange_proton = 0x0;
  vector<double>          *reco_daughter_allTrack_momByRange_muon = 0x0;
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

  // Declare histograms
  TH1I * hCutBeamIDPass = 0x0;
  TH1I * hCutBeamType = 0x0;
  TH1I * hCutBeamPosPass = 0x0;
  TH1I * hCutBeamEndZPass = 0x0;
  TH1D * hCutBeamEndZ = 0x0; 
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
  TH2D * hRecShowerOpenAngle = 0x0;
  TH2D * hRecPi0Mass = 0x0;
  TH2D * hRecPi0MassRaw = 0x0;
  TH2D * hRecPi0MassFit = 0x0; 
  TH2D * hRecPi0Momentum = 0x0;
  TH2D * hRecPi0MomentumRaw = 0x0;
  TH2D * hRecPi0ShowerSep = 0x0;
  TH2D * hRecLeadingShowerEnergy = 0x0;
  TH2D * hRecSubLeadingShowerEnergy = 0x0;
  TH2D * hRecLeadingShowerEnergyRaw = 0x0;
  TH2D * hRecSubLeadingShowerEnergyRaw = 0x0;

  // Common Reco TKI
  TH2D * hRecdalphat = 0x0;
  TH2D * hRecdphit = 0x0;
  TH2D * hRecdpt = 0x0;
  TH2D * hRecpn = 0x0;
 
  TH2D * hCutTracknHits = 0x0;
  TH2D * hCutTrackScore = 0x0;
  TH2D * hCutlastTME = 0x0;
  TH2D * hCutChi2NDF = 0x0;
  TH2D * hCutemScore = 0x0;
  TH2D * hCutmichelScore = 0x0;
  TH2D * hCutShowerDist = 0x0;
  TH2D * hCutShowerIP = 0x0;

  TH2D * hCutnproton = 0x0;
  TH2D * hCutnpiplus = 0x0;
  TH2D * hCutnshower = 0x0;
  TH2D * hCutnmichel = 0x0;
  
  //====================== Reco (Data only)======================//
  // Declare variables
  Double_t        beam_inst_X;
  Double_t        beam_inst_Y;
  Double_t        beam_inst_dirX;
  Double_t        beam_inst_dirY;
  Double_t        beam_inst_dirZ;
  Int_t           beam_inst_nMomenta;
  Int_t           beam_inst_nTracks;
  vector<int> *beam_inst_PDG_candidates = 0x0;
  
  //====================== Truth (MC only)======================//
  // Declare variables
  bool Signal = false;
  int true_beam_PDG = -999;
  Double_t        true_beam_startX;
  Double_t        true_beam_startY;
  Double_t        true_beam_startZ;  
  Double_t        true_beam_startDirX;
  Double_t        true_beam_startDirY;
  Double_t        true_beam_startDirZ;
  Double_t        true_beam_endX;
  Double_t        true_beam_endY;
  Double_t        true_beam_endZ;
  double          true_beam_endPx = -999;
  double          true_beam_endPy = -999;
  double          true_beam_endPz = -999;

  vector<double> *true_beam_daughter_startPx = 0x0;
  vector<double> *true_beam_daughter_startPy = 0x0;
  vector<double> *true_beam_daughter_startPz = 0x0;
  vector<double> *true_beam_daughter_startX = 0x0;
  vector<double> *true_beam_daughter_startY = 0x0;
  vector<double> *true_beam_daughter_startZ = 0x0;

  vector<int> *true_beam_daughter_ID = 0x0;  
  vector<int> *true_beam_daughter_PDG = 0x0;
  vector<int> *true_beam_Pi0_decay_ID = 0x0;
  vector<int> *true_beam_Pi0_decay_PDG = 0x0;
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
  TH1D * hTruthGammaMaxE = 0x0;
  TH2D * hTruthPi0ShowerEnergy = 0x0;
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
  TH2D * hProtonMomentumRes = 0x0;
  TH2D * hPiPlusThetaRes = 0x0;
  TH2D * hPiPlusMomentumRes = 0x0;
  TH2D * hShowerEnergyRes = 0x0;
  TH2D * hShowerEnergyResRaw = 0x0;
  TH2D * hShowerThetaRes = 0x0; 

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

  TH1D * hTruthPi0Momentum = 0x0;
  TH1D * hTruthldShowerTheta = 0x0;
  TH1D * hTruthslShowerTheta = 0x0;
  TH1D * hTruthPi0Mass = 0x0;

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

  TH2D * hShowerE1PreFitRes = 0x0;
  TH2D * hShowerE1PostFitRes = 0x0;

  TH2D * hShowerE2PreFitRes = 0x0;
  TH2D * hShowerE2PostFitRes = 0x0;

  TH2D * hShowerOAPreFitRes = 0x0;
  TH2D * hShowerOAPostFitRes = 0x0;

  TH1D * hShowerE1Compare = 0x0;
  TH1D * hShowerE2Compare = 0x0;
  TH1D * hShowerOACompare = 0x0;

  TH1D * hShowerE1ComparePost = 0x0;
  TH1D * hShowerE2ComparePost = 0x0;
  TH1D * hShowerOAComparePost = 0x0;

  TH1D * hPi0MassCompare = 0x0;
  TH1D * hPi0MassComparePost = 0x0;

  // Completeness and Purity
  TH2D * hTrackPurityVSnHits = 0x0;
  TH2D * hTrackCompletenessVSnHits = 0x0;
  TH2D * hShowerPurityVSnHits = 0x0;
  TH2D * hShowerCompletenessVSnHits = 0x0;

  TH2D * hTrackCompleteness = 0x0;
  TH2D * hShowerCompleteness = 0x0;
  

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

    tree->SetBranchAddress("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
    tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
    tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
    tree->SetBranchAddress("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
    tree->SetBranchAddress("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);
    tree->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
    tree->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);
    tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
    tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);
    tree->SetBranchAddress("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
    tree->SetBranchAddress("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);
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

    //====================== Reco (Data only)======================//
    tree->SetBranchAddress("beam_inst_PDG_candidates", &beam_inst_PDG_candidates);
    tree->SetBranchAddress("beam_inst_X", &beam_inst_X);
    tree->SetBranchAddress("beam_inst_Y", &beam_inst_Y);
    tree->SetBranchAddress("beam_inst_dirX", &beam_inst_dirX);
    tree->SetBranchAddress("beam_inst_dirY", &beam_inst_dirY);
    tree->SetBranchAddress("beam_inst_dirZ", &beam_inst_dirZ);
    tree->SetBranchAddress("beam_inst_nMomenta", &beam_inst_nMomenta);
    tree->SetBranchAddress("beam_inst_nTracks", &beam_inst_nTracks);

    //====================== Truth (MC only)======================//
    tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
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
   
    tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID); 
    tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
    tree->SetBranchAddress("true_beam_Pi0_decay_ID", &true_beam_Pi0_decay_ID);
    tree->SetBranchAddress("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG);
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

    //tout->Branch("LeadingShowerEnergyUnitDir",&LeadingShowerEnergyUnitDir);
    //tout->Branch("SubLeadingShowerEnergyUnitDir",&SubLeadingShowerEnergyUnitDir);
    //tout->Branch("LeadingShowerEnergyUnitDirTruth",&LeadingShowerEnergyUnitDirTruth);
    //tout->Branch("SubLeadingShowerEnergyUnitDirTruth",&SubLeadingShowerEnergyUnitDirTruth);
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

    hCutBeamIDPass = new TH1I("a001CutBeamIDPass","", 2, -0.5, 1.5); 
    lout->Add(hCutBeamIDPass);
    hCutBeamType = new TH1I("a002CutBeamType","",30, -4.5, 25.5);
    lout->Add(hCutBeamType); 
    hCutBeamPosPass = new TH1I("a003CutBeamPosPass", "", 4, -0.5, 3.5); 
    lout->Add(hCutBeamPosPass);
    hCutBeamEndZPass = new TH1I("a004CutBeamEndZPass", "", 4, -0.5, 3.5);
    lout->Add(hCutBeamEndZPass);
    hCutBeamEndZ = new TH1D("a005CutBeamEndZ","",50, 0, 500);
    lout->Add(hCutBeamEndZ);
    hRecBeamTheta = new TH2D("b001RecBeamTheta_STK","", 80 , 0, 60, 3, -0.5, 2.5); 
    lout->Add(hRecBeamTheta);
    hRecBeamMomentum = new TH2D("b002RecBeamMomentum_STK","", 50, 0, 2, 3, -0.5, 2.5); 
    lout->Add(hRecBeamMomentum);
    hRecProtonTheta = new TH2D("b003RecProtonTheta_STK","",  20, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecProtonTheta);
    hRecProtonMomentum = new TH2D("b004RecProtonMomentum_STK","",  20, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecProtonMomentum);
    hRecPiPlusTheta = new TH2D("b005RecPiPlusTheta_STK","", 15, 0, 180, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusTheta);
    hRecPiPlusMomentum = new TH2D("b006RecPiPlusMomentum_STK","",  20, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusMomentum);
    hRecShowerEnergy = new TH2D("b007RecShowerEnergy_STK","", 15, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecShowerEnergy);
    hRecShowerEnergyRaw = new TH2D("b008RecShowerEnergy_STK_RAW","", 15, 0, 1.2, nparType, parTypemin, parTypemax);
    lout->Add(hRecShowerEnergyRaw);

    hRecShowerTheta = new TH2D("b009RecShowerTheta_STK","", 15, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecShowerTheta);
    hRecShowerLength = new TH2D("b010RecShowerLength_STK","", 40, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hRecShowerLength);

    hRecPi0Nshower = new TH2D("b011RecPi0Nshower_STK","", 10, -0.5, 9.5, 3, -0.5, 2.5); 
    lout->Add(hRecPi0Nshower);
    hRecShowerOpenAngle = new TH2D("b012RecShowerOpenAngle_STK","", 20, 0, 180, 3, -0.5, 2.5); 
    lout->Add(hRecShowerOpenAngle);
    hRecPi0Mass = new TH2D("b013RecPi0Mass_STK","", 15, 0, 0.5, 3, -0.5, 2.5); 
    lout->Add(hRecPi0Mass);
    hRecPi0MassRaw = new TH2D("b014RecPi0Mass_STK_RAW","", 15, 0, 0.5, 3, -0.5, 2.5);
    lout->Add(hRecPi0MassRaw);
    //hRecPi0MassFit = new TH2D("b015RecPi0Mass_STK_FIT","", 15, 0, 0.5, 3, -0.5, 2.5);
    //lout->Add(hRecPi0MassFit);
    hRecPi0Momentum = new TH2D("b016RecPi0Momentum_STK", "", 15, 0, 1, 3, -0.5, 2.5);
    lout->Add(hRecPi0Momentum);
    hRecPi0MomentumRaw = new TH2D("b017RecPi0Momentum_STK_RAW", "", 15, 0, 1, 3, -0.5, 2.5);
    lout->Add(hRecPi0MomentumRaw);

    hRecPi0ShowerSep = new TH2D("b018RecPi0ShowerSep_STK","", 50, 0, 200, 3, -0.5, 2.5);
    lout->Add(hRecPi0ShowerSep);
    hRecLeadingShowerEnergy = new TH2D("c001RecLeadingShowerEnergy_STK","", 15, 0, 0.8, nparType, parTypemin, parTypemax);
    lout->Add(hRecLeadingShowerEnergy);
    hRecSubLeadingShowerEnergy = new TH2D("c002RecSubLeadingShowerEnergy_STK","", 15, 0, 0.5, nparType, parTypemin, parTypemax);
    lout->Add(hRecSubLeadingShowerEnergy); 
    hRecLeadingShowerEnergyRaw = new TH2D("c003RecLeadingShowerEnergy_STK_RAW","", 15, 0, 0.8, nparType, parTypemin, parTypemax);
    lout->Add(hRecLeadingShowerEnergyRaw);
    hRecSubLeadingShowerEnergyRaw = new TH2D("c004RecSubLeadingShowerEnergy_STK_RAW","", 15, 0, 0.5, nparType, parTypemin, parTypemax);
    lout->Add(hRecSubLeadingShowerEnergyRaw);

    // Commom reco TKI
    hRecdalphat = new TH2D("d001Recdalphat_STK","", 9, 0, 360,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdalphat);
    hRecdphit = new TH2D("d002Recphit_STK","", 9, 0, 360,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdphit);
    hRecdpt = new TH2D("d003Recdpt_STK","", 8, 0, 1.8,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdpt);
    hRecpn = new TH2D("d004Recpn_STK","", 8, 0, 1.8,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecpn);

    hCutTracknHits = new TH2D("e001CutTracknHits_STK","", 50, 0, 500, nparType, parTypemin, parTypemax); 
    lout->Add(hCutTracknHits);
    hCutTrackScore = new TH2D("e002CutTrackScore_STK","", 50, 0, 1, nparType, parTypemin, parTypemax);
    lout->Add(hCutTrackScore);
    hCutlastTME = new TH2D("e003CutlastTME_STK","", 60, 0, 10, nparType, parTypemin, parTypemax); 
    lout->Add(hCutlastTME);
    hCutChi2NDF = new TH2D("e004CutChi2NDF_STK","", 30, 0, 500, nparType, parTypemin, parTypemax);   
    lout->Add(hCutChi2NDF);
    hCutemScore = new TH2D("e005CutemScore_STK","", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore);
    hCutmichelScore = new TH2D("e006CutmichelScore_STK","",50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutmichelScore);
    hCutShowerDist = new TH2D("e007CutShowerDist_STK","", 31, 0, 93, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerDist);
    hCutShowerIP = new TH2D("e008CutShowerIP_STK","", 30, 0, 30, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerIP);

    hCutnproton = new TH2D("e009Cutnproton_STK","", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnproton);
    hCutnshower = new TH2D("e010Cutnshower_STK","", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnshower);
    hCutnpiplus = new TH2D("e011Cutnpiplus_STK","",  ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnpiplus);
    hCutnmichel = new TH2D("e012Cutnmichel_STK","", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); 
    lout->Add(hCutnmichel);


    //====================== Truth (MC only)======================//
    if(kMC){
      hTruthBeamType = new TH1I("f001TruthBeamType",  "", 20, -0.5, 19.5); 
      lout->Add(hTruthBeamType);
      hTruthFSParticleNumber = new TH1I("f002TruthFSParticleNumber", "", 36, -0.5, 35.5); 
      lout->Add(hTruthFSParticleNumber);
      hTruthFSParticleType = new TH1I("f003TruthFSParticleType", "", 23, -0.5, 22.5);
      lout->Add(hTruthFSParticleType);
      hTruthFSPi0Number = new TH1I("f004TruthFSPi0Number", "", 5, 0, 5);
      lout->Add(hTruthFSPi0Number);
      hTruthFSMultiPi0 = new TH1I("f005TruthFSMultiPi0", "", 5, 0, 5);
      lout->Add(hTruthFSMultiPi0);
      hTruthLeadingProtonP = new TH1D("f006TruthLeadingProtonP", "", 20, 0, 2);
      lout->Add(hTruthLeadingProtonP);
      hTruthSubLeadingProtonP = new TH1D("f007TruthSubLeadingProtonP", "", 20, 0, 2);
      lout->Add(hTruthSubLeadingProtonP);
      hTruthGammaMaxE = new TH1D("f008TruthGammaMaxE", "", 20, 0, 2);
      lout->Add(hTruthGammaMaxE);
      hTruthPi0ShowerEnergy = new TH2D("f009TruthPi0ShowerEnergy_OVERLAY", "", 20, 0, 1.2, 3, -0.5, 2.5);
      lout->Add(hTruthPi0ShowerEnergy);
      hTruthSignal = new TH1I("f010TruthSignal",  "",  nPass, Passmin, Passmax); 
      lout->Add(hTruthSignal); 
      hTruthSignalFSParticleNumber = new TH1I("f011TruthSignalFSParticleNumber", "", 36, -0.5, 35.5); 
      lout->Add(hTruthSignalFSParticleNumber);
      hTruthSignalFSParticleType = new TH1I("f012TruthSignalFSParticleType", "", 23, -0.5, 22.5);
      lout->Add(hTruthSignalFSParticleType);

      hTruthNproton = new TH1I("f100TruthNproton","",11, -0.5, 10.5); 
      lout->Add(hTruthNproton);
      hTruthNneutron = new TH1I("f101TruthNneutron","",11, -0.5, 10.5); 
      lout->Add(hTruthNneutron);
      hTruthNPiZero = new TH1I("f102TruthNPiZero","",11, -0.5, 10.5); 
      lout->Add(hTruthNPiZero);

      hTruthMomIniPi = new TH1D("f103TruthMomIniPi","", 30, 0, 1.3); 
      lout->Add(hTruthMomIniPi);
      hTruthThetaIniPi = new TH1D("f104TruthThetaIniPi","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPi);
      hTruthMomFinPi = new TH1D("f105TruthMomFinPi","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPi);
      hTruthThetaFinPi = new TH1D("f106TruthThetaFinPi","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPi);
      hTruthMomFinProton = new TH1D("f107TruthMomFinProton","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProton);
      hTruthThetaFinProton = new TH1D("f108TruthThetaFinProton","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProton);
      hTruthMomFin2Proton = new TH1D("f109TruthMomFin2Proton","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2Proton);

      hTruthDalphat = new TH1D("f200TruthDalphat","", 18, 0, 180); 
      lout->Add(hTruthDalphat);
      hTruthDphit = new TH1D("f201TruthDphit","", 18, 0, 180); 
      lout->Add(hTruthDphit);
      hTruthDpt = new TH1D("f202TruthDpt","", 30, 0, 1.2); 
      lout->Add(hTruthDpt);
      hTruthPn = new TH1D("f203TruthPn","", 30, 0, 1.2); 
      lout->Add(hTruthPn);

      stkTruthDalphat = new THStack("f300stkTruthDalphat",""); 
      lout->Add(stkTruthDalphat);
      hTruthDalphat1p0n = new TH1D("f300TruthDalphat1p0n","", 18, 0, 180); 
      lout->Add(hTruthDalphat1p0n);
      hTruthDalphatNp0n = new TH1D("f300TruthDalphatNp0n","", 18, 0, 180); 
      lout->Add(hTruthDalphatNp0n);
      hTruthDalphat1pMn = new TH1D("f300TruthDalphat1pMn","", 18, 0, 180); 
      lout->Add(hTruthDalphat1pMn);
      hTruthDalphatNpMn = new TH1D("f300TruthDalphatNpMn","", 18, 0, 180); 
      lout->Add(hTruthDalphatNpMn);

      stkTruthDphit = new THStack("f301stkTruthDphit",""); 
      lout->Add(stkTruthDphit);
      hTruthDphit1p0n = new TH1D("f301TruthDphit1p0n","", 18, 0, 180); 
      lout->Add(hTruthDphit1p0n);
      hTruthDphitNp0n = new TH1D("f301TruthDphitNp0n","", 18, 0, 180); 
      lout->Add(hTruthDphitNp0n);
      hTruthDphit1pMn = new TH1D("f301TruthDphit1pMn","", 18, 0, 180); 
      lout->Add(hTruthDphit1pMn);
      hTruthDphitNpMn = new TH1D("f301TruthDphitNpMn","", 18, 0, 180); 
      lout->Add(hTruthDphitNpMn);

      stkTruthDpt = new THStack("f302stkTruthDpt",""); 
      lout->Add(stkTruthDpt);
      hTruthDpt1p0n = new TH1D("f302TruthDpt1p0n","", 30, 0, 1.2); 
      lout->Add(hTruthDpt1p0n);
      hTruthDptNp0n = new TH1D("f302TruthDptNp0n","", 30, 0, 1.2); 
      lout->Add(hTruthDptNp0n);
      hTruthDpt1pMn = new TH1D("f302TruthDpt1pMn","", 30, 0, 1.2); 
      lout->Add(hTruthDpt1pMn);
      hTruthDptNpMn = new TH1D("f302TruthDptNpMn","", 30, 0, 1.2); 
      lout->Add(hTruthDptNpMn);

      stkTruthPn = new THStack("f303stkTruthPn",""); 
      lout->Add(stkTruthPn);
      hTruthPn1p0n = new TH1D("f303TruthPn1p0n","", 30, 0, 1.2); 
      lout->Add(hTruthPn1p0n);
      hTruthPnNp0n = new TH1D("f303TruthPnNp0n","", 30, 0, 1.2); 
      lout->Add(hTruthPnNp0n);
      hTruthPn1pMn = new TH1D("f303TruthPn1pMn","", 30, 0, 1.2); 
      lout->Add(hTruthPn1pMn);
      hTruthPnNpMn = new TH1D("f303TruthPnNpMn","", 30, 0, 1.2); 
      lout->Add(hTruthPnNpMn);

      stkTruthMomIniPi = new THStack("f400stkTruthMomIniPi",""); 
      lout->Add(stkTruthMomIniPi);
      hTruthMomIniPi1p0n = new TH1D("f400TruthMomIniPi1p0n","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPi1p0n);
      hTruthMomIniPiNp0n = new TH1D("f400TruthMomIniPiNp0n","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPiNp0n);
      hTruthMomIniPi1pMn = new TH1D("f400TruthMomIniPi1pMn","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPi1pMn);
      hTruthMomIniPiNpMn = new TH1D("f400TruthMomIniPiNpMn","", 30, 0, 1.2); 
      lout->Add(hTruthMomIniPiNpMn);

      stkTruthThetaIniPi = new THStack("f401stkTruthThetaIniPi",""); 
      lout->Add(stkTruthThetaIniPi);
      hTruthThetaIniPi1p0n = new TH1D("f401TruthThetaIniPi1p0n","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPi1p0n);
      hTruthThetaIniPiNp0n = new TH1D("f401TruthThetaIniPiNp0n","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPiNp0n);
      hTruthThetaIniPi1pMn = new TH1D("f401TruthThetaIniPi1pMn","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPi1pMn);
      hTruthThetaIniPiNpMn = new TH1D("f401TruthThetaIniPiNpMn","", 30, 0, 30); 
      lout->Add(hTruthThetaIniPiNpMn);

      stkTruthMomFinPi = new THStack("f402stkTruthMomFinPi",""); 
      lout->Add(stkTruthMomFinPi);
      hTruthMomFinPi1p0n = new TH1D("f402TruthMomFinPi1p0n","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPi1p0n);
      hTruthMomFinPiNp0n = new TH1D("f402TruthMomFinPiNp0n","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPiNp0n);
      hTruthMomFinPi1pMn = new TH1D("f402TruthMomFinPi1pMn","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPi1pMn);
      hTruthMomFinPiNpMn = new TH1D("f402TruthMomFinPiNpMn","", 20, 0, 1.1); 
      lout->Add(hTruthMomFinPiNpMn);

      stkTruthThetaFinPi = new THStack("f403stkTruthThetaFinPi",""); 
      lout->Add(stkTruthThetaFinPi);
      hTruthThetaFinPi1p0n = new TH1D("f403TruthThetaFinPi1p0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPi1p0n);
      hTruthThetaFinPiNp0n = new TH1D("f403TruthThetaFinPiNp0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPiNp0n);
      hTruthThetaFinPi1pMn = new TH1D("f403TruthThetaFinPi1pMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPi1pMn);
      hTruthThetaFinPiNpMn = new TH1D("f403TruthThetaFinPiNpMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinPiNpMn);

      stkTruthMomFinProton = new THStack("f404stkTruthMomFinProton",""); 
      lout->Add(stkTruthMomFinProton);
      hTruthMomFinProton1p0n = new TH1D("f404TruthMomFinProton1p0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProton1p0n);
      hTruthMomFinProtonNp0n = new TH1D("f404TruthMomFinProtonNp0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProtonNp0n);
      hTruthMomFinProton1pMn = new TH1D("f404TruthMomFinProton1pMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProton1pMn);
      hTruthMomFinProtonNpMn = new TH1D("f404TruthMomFinProtonNpMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFinProtonNpMn);

      stkTruthThetaFinProton = new THStack("f405stkTruthThetaFinProton",""); 
      lout->Add(stkTruthThetaFinProton);
      hTruthThetaFinProton1p0n = new TH1D("f405TruthThetaFinProton1p0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProton1p0n);
      hTruthThetaFinProtonNp0n = new TH1D("f405TruthThetaFinProtonNp0n","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProtonNp0n);
      hTruthThetaFinProton1pMn = new TH1D("f405TruthThetaFinProton1pMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProton1pMn);
      hTruthThetaFinProtonNpMn = new TH1D("f405TruthThetaFinProtonNpMn","", 20, 0, 180); 
      lout->Add(hTruthThetaFinProtonNpMn);

      stkTruthMomFin2Proton = new THStack("f406stkTruthMomFin2Proton",""); 
      lout->Add(stkTruthMomFin2Proton);
      hTruthMomFin2Proton1p0n = new TH1D("f406TruthMomFin2Proton1p0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2Proton1p0n);
      hTruthMomFin2ProtonNp0n = new TH1D("f406TruthMomFin2ProtonNp0n","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2ProtonNp0n);
      hTruthMomFin2Proton1pMn = new TH1D("f406TruthMomFin2Proton1pMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2Proton1pMn);
      hTruthMomFin2ProtonNpMn = new TH1D("f406TruthMomFin2ProtonNpMn","", 30, 0, 1.5); 
      lout->Add(hTruthMomFin2ProtonNpMn); 

      hBeamThetaRes = new TH2D("g001BeamTheta_RES","", 80 , 0, 60, 25, -20, 30);
      lout->Add(hBeamThetaRes);
      hBeamMomentumRes  = new TH2D("g002BeamMomentum_RES","", 50, 0, 2, 20, -0.5, 0.5);
      hBeamMomentumRes  = new TH2D("g003BeamMomentum_RES","", 50, 0, 2, 20, -0.5, 0.5);
      lout->Add(hBeamMomentumRes);
      hProtonThetaRes = new TH2D("g004ProtonTheta_RES","", 20, 0, 180, 25, -20, 30); 
      lout->Add(hProtonThetaRes);
      hProtonMomentumRes = new TH2D("g005ProtonMomentum_RES","",20, 0, 1.2, 20, -0.2, 0.2); 
      lout->Add(hProtonMomentumRes);
      hPiPlusThetaRes = new TH2D("g006PiPlusTheta_RES","", 15, 0, 180, 25, -20, 30); 
      lout->Add(hPiPlusThetaRes);
      hPiPlusMomentumRes = new TH2D("g007PiPlusMomentum_RES","", 20, 0, 1.2, 20, -0.2, 0.2); 
      lout->Add(hPiPlusMomentumRes);

      hShowerEnergyRes = new TH2D("h001ShowerEnergy_RES","", 20, 0, 1.2, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyRes);
      hShowerEnergyResRaw = new TH2D("h002ShowerEnergy_RES_RAW","", 20, 0, 1.2, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResRaw);
      hShowerThetaRes = new TH2D("h003ShowerThetaRes_RES","", 15, 0, 180, 25, -20, 30);
      lout->Add(hShowerThetaRes);
      hShowerEnergyRecVSTruth_REG = new TH2D("h004ShowerEnergyRecVSTruth_REG","", 20, 0, 1, 20, 0, 1);
      lout->Add(hShowerEnergyRecVSTruth_REG);
      hShowerEnergyRawRecVSTruth_REG = new TH2D("h005ShowerEnergyRawRecVSTruth_REG","", 20, 0, 1, 20, 0, 1);
      lout->Add(hShowerEnergyRawRecVSTruth_REG);
      hShowerEnergyResVSnHits_REG = new TH2D("h006ShowerEnergyResVSnHits_REG","", 50, 0, 1000, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResVSnHits_REG); 
      hShowerEnergyResVSIP_REG = new TH2D("h007ShowerEnergyResVSIP_REG","", 20, 0, 30, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResVSIP_REG);
      hLeadingShowerEnergyRes = new TH2D("h008LeadingShowerEnergy_RES","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hLeadingShowerEnergyRes);
      hSubLeadingShowerEnergyRes = new TH2D("h009SubLeadingShowerEnergy_RES","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hSubLeadingShowerEnergyRes);  
      hLeadingShowerEnergyResRaw = new TH2D("h010LeadingShowerEnergy_RES_RAW","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hLeadingShowerEnergyResRaw);
      hSubLeadingShowerEnergyResRaw = new TH2D("h011SubLeadingShowerEnergy_RES_RAW","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hSubLeadingShowerEnergyResRaw);
      hShowerOpenAngleRes = new TH2D("h012ShowerOpenAngle_RES","", 20, 0, 180, 20, -50, 50);
      lout->Add(hShowerOpenAngleRes); 
      hPi0MomentumRes = new TH2D("i001Pi0Momentum_RES","", 20, 0, 1, 10, -0.5, 0.5);
      lout->Add(hPi0MomentumRes);
      hPi0MassRes = new TH2D("i002Pi0Mass_RES","", 20, 0, 0.5, 20, -0.5, 0.5); 
      lout->Add(hPi0MassRes);

      hTruthPi0Momentum = new TH1D("i003TruthPi0Momentum","", 20, 0, 1.5);
      lout->Add(hTruthPi0Momentum);

      hTruthldShowerTheta = new TH1D("i004TruthldShowerTheta","", 30, 0, 180);
      lout->Add(hTruthldShowerTheta);
      hTruthslShowerTheta = new TH1D("i005TruthslShowerTheta","", 30, 0, 180);
      lout->Add(hTruthslShowerTheta);

      hPi0MomentumResRaw = new TH2D("i006Pi0Momentum_RES_RAW","", 20, 0, 1, 10, -0.5, 0.5);
      lout->Add(hPi0MomentumResRaw);
      hPi0MassResRaw = new TH2D("i007Pi0Mass_RES_RAW","", 20, 0, 0.5, 20, -0.5, 0.5);
      lout->Add(hPi0MassResRaw);

      //hPi0MomentumResFit = new TH2D("i008Pi0Momentum_RES_FIT","", 20, 0, 1, 10, -0.5, 0.5);
      //lout->Add(hPi0MomentumResFit);
      //hPi0MassResFit = new TH2D("i009Pi0Mass_RES_FIT","", 20, 0, 0.5, 20, -0.5, 0.5);
      //lout->Add(hPi0MassResFit);
      //hShowerOpenAngleResFit = new TH2D("i010ShowerOpenAngle_RES_FIT","", 20, 0, 180, 20, -1.1, 1.1);
      //lout->Add(hShowerOpenAngleResFit);

      hTruthPi0Mass = new TH1D("i011hTruthPi0Mass","", 20, 0, 1);
      lout->Add(hTruthPi0Mass);

      // CVM bin histograms
      const double BinE1_1x1[] = {0,1.2};
      const double BinE2_1x1[] = {0,0.8};
      
      const double BinE1_2x2[] = {0,0.3,1.2};
      const double BinE2_2x2[] = {0,0.15,0.8};

      const double BinE1_3x3[] = {0,0.25,0.35,1.2};
      const double BinE2_3x3[] = {0,0.12,0.2,0.8};
      
      const double BinE1_4x4[] = {0,0.2,0.3,0.4,1.2};
      const double BinE2_4x4[] = {0,0.1,0.16,0.26,0.8};


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


      hShowerE1PreFitRes = new TH2D("x001hShowerE1Pre_RES"," ",20, 0, 0.8, 20, -1.1, 1.1);
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

      hShowerE1Compare = new TH1D("y001hShowerE1Compare","", 20, -1.1, 1.1);
      lout->Add(hShowerE1Compare);
      hShowerE2Compare = new TH1D("y002hShowerE2Compare","", 20, -1.1, 1.1);
      lout->Add(hShowerE2Compare);
      hShowerOACompare = new TH1D("y003hShowerOACompare","", 20, -50, 50);
      lout->Add(hShowerOACompare);
      hPi0MassCompare = new TH1D("y004hPi0MassCompare","", 20, 0, 0.3);
      lout->Add(hPi0MassCompare);

      hShowerE1ComparePost = new TH1D("y001hShowerE1ComparePost","", 20, -1.1, 1.1);
      lout->Add(hShowerE1ComparePost);
      hShowerE2ComparePost = new TH1D("y002hShowerE2ComparePost","", 20, -1.1, 1.1);
      lout->Add(hShowerE2ComparePost);
      hShowerOAComparePost = new TH1D("y003hShowerOAComparePost","", 20, -50, 50);
      lout->Add(hShowerOAComparePost);
      hPi0MassComparePost = new TH1D("y004hPi0MassComparePost","", 20, 0, 0.3);
      lout->Add(hPi0MassComparePost);

      hTrackPurityVSnHits = new TH2D("z001hTrackPurityVSnHits_REG","", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hTrackPurityVSnHits);
      hTrackCompletenessVSnHits = new TH2D("z002hTrackCompletenessVSnHits_REG","", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hTrackCompletenessVSnHits);

      hShowerPurityVSnHits = new TH2D("z003hShowerPurityVSnHits_REG","", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hShowerPurityVSnHits);
      hShowerCompletenessVSnHits = new TH2D("z004hShowerCompletenessVSnHits_REG","", 50, 0, 1000, 20, 0, 1.2);
      lout->Add(hShowerCompletenessVSnHits);

      hTrackCompleteness = new TH2D("z005hTrackCompleteness_STK","", 20, 0, 1, nparType, parTypemin, parTypemax );
      lout->Add(hTrackCompleteness);
      hShowerCompleteness = new TH2D("z006hShowerCompleteness_STK","", 20, 0, 1, nparType, parTypemin, parTypemax );
      lout->Add(hShowerCompleteness);


    }
  }// End of IniHist

} // End of namespace
#endif
