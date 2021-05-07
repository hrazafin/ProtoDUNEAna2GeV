#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "PlotUtils.h"
#include "AnaFunctions.h"

namespace AnaIO
{

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
  vector<double>  *reco_beam_calibrated_dEdX = 0x0;

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
  TH1I * hEvent = 0x0;
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

  // Declare histograms
  TH1I * hTruthBeamType = 0x0;
  TH1I * hTruthFSParticleNumber = 0x0;  
  TH1I * hTruthFSParticleType = 0x0;
  TH1I * hTruthFSPi0Number = 0x0;
  TH1I * hTruthFSMultiPi0 = 0x0;
  TH1D * hTruthLeadingProtonP = 0x0;
  TH1D * hTruthSubLeadingProtonP = 0x0;
  TH1D * hTruthGammaMaxE = 0x0;
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
    tree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);

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

    return tree;
  } // End of GetInputTree
  
  // Initialise histograms
  void IniHist(TList * lout, const bool kMC)
  {
    // Binning 
    const int nparType = 11;
    const double parTypemin = 0.5;
    const double parTypemax = 11.5;

    //====================== Reco (MC and Data)======================//
    hEvent = new TH1I("Event", "", 10000, 0, 100000); 
    lout->Add(hEvent);
    hCutBeamIDPass = new TH1I("CutBeamIDPass","", 2, -0.5, 1.5); 
    lout->Add(hCutBeamIDPass);
    hCutBeamType = new TH1I("CutBeamType","",30, -4.5, 25.5);
    lout->Add(hCutBeamType); 
    hCutBeamPosPass = new TH1I("CutBeamPosPass", "", 4, -0.5, 3.5); 
    lout->Add(hCutBeamPosPass);
    hCutBeamEndZPass = new TH1I("CutBeamEndZPass", "", 4, -0.5, 3.5);
    lout->Add(hCutBeamEndZPass); 
    hCutBeamEndZ = new TH1D("CutBeamEndZ","",50, 0, 500);
    lout->Add(hCutBeamEndZ);
    hRecBeamTheta = new TH2D("RecBeamTheta_STK","", 80 , 0, 60, 3, -0.5, 2.5); 
    lout->Add(hRecBeamTheta);
    hRecBeamMomentum = new TH2D("RecBeamMomentum_STK","", 50, 0, 2, 3, -0.5, 2.5); 
    lout->Add(hRecBeamMomentum);
    hRecProtonTheta = new TH2D("RecProtonTheta_STK","",  20, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecProtonTheta);
    hRecProtonMomentum = new TH2D("RecProtonMomentum_STK","",  20, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecProtonMomentum);
    hRecPiPlusTheta = new TH2D("RecPiPlusTheta_STK","", 15, 0, 180, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusTheta);
    hRecPiPlusMomentum = new TH2D("RecPiPlusMomentum_STK","",  20, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecPiPlusMomentum);
    hRecShowerEnergy = new TH2D("RecShowerEnergy_STK","", 15, 0, 1.2, nparType, parTypemin, parTypemax); 
    lout->Add(hRecShowerEnergy);
    hRecShowerEnergyRaw = new TH2D("RecShowerEnergy_STK_RAW","", 15, 0, 1.2, nparType, parTypemin, parTypemax);
    lout->Add(hRecShowerEnergyRaw);

    hRecShowerTheta = new TH2D("RecShowerTheta_STK","", 15, 0, 180, nparType, parTypemin, parTypemax);
    lout->Add(hRecShowerTheta);
    hRecShowerLength = new TH2D("RecShowerLength_STK","", 40, 0, 200, nparType, parTypemin, parTypemax); 
    lout->Add(hRecShowerLength);

    hRecPi0Nshower = new TH2D("RecPi0Nshower_STK","", 10, -0.5, 9.5, 3, -0.5, 2.5); 
    lout->Add(hRecPi0Nshower);
    hRecShowerOpenAngle = new TH2D("RecShowerOpenAngle_STK","", 20, 0, 180, 3, -0.5, 2.5); 
    lout->Add(hRecShowerOpenAngle);
    hRecPi0Mass = new TH2D("RecPi0Mass_STK","", 15, 0, 0.5, 3, -0.5, 2.5); 
    lout->Add(hRecPi0Mass);
    hRecPi0MassRaw = new TH2D("RecPi0Mass_STK_RAW","", 15, 0, 0.5, 3, -0.5, 2.5);
    lout->Add(hRecPi0MassRaw);
    hRecPi0MassFit = new TH2D("RecPi0Mass_STK_FIT","", 15, 0, 0.5, 3, -0.5, 2.5);
    lout->Add(hRecPi0MassFit);
    hRecPi0Momentum = new TH2D("RecPi0Momentum_STK", "", 15, 0, 1, 3, -0.5, 2.5);
    lout->Add(hRecPi0Momentum);
    hRecPi0MomentumRaw = new TH2D("RecPi0Momentum_STK_RAW", "", 15, 0, 1, 3, -0.5, 2.5);
    lout->Add(hRecPi0MomentumRaw);

    hRecPi0ShowerSep = new TH2D("RecPi0ShowerSep_STK","", 50, 0, 200, 3, -0.5, 2.5);
    lout->Add(hRecPi0ShowerSep);
    hRecLeadingShowerEnergy = new TH2D("RecLeadingShowerEnergy_STK","", 15, 0, 0.8, nparType, parTypemin, parTypemax);
    lout->Add(hRecLeadingShowerEnergy);
    hRecSubLeadingShowerEnergy = new TH2D("RecSubLeadingShowerEnergy_STK","", 15, 0, 0.5, nparType, parTypemin, parTypemax);
    lout->Add(hRecSubLeadingShowerEnergy); 
    hRecLeadingShowerEnergyRaw = new TH2D("RecLeadingShowerEnergy_STK_RAW","", 15, 0, 0.8, nparType, parTypemin, parTypemax);
    lout->Add(hRecLeadingShowerEnergyRaw);
    hRecSubLeadingShowerEnergyRaw = new TH2D("RecSubLeadingShowerEnergy_STK_RAW","", 15, 0, 0.5, nparType, parTypemin, parTypemax);
    lout->Add(hRecSubLeadingShowerEnergyRaw);

  const int nevtType = 3;
  const double evtTypemin = -0.5;
  const double evtTypemax = 2.5;

    // Commom reco TKI
    hRecdalphat = new TH2D("Recdalphat_STK","", 9, 0, 180,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdalphat);
    hRecdphit = new TH2D("Recphit_STK","", 9, 0, 180,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdphit);
    hRecdpt = new TH2D("Recdpt_STK","", 8, 0, 0.8,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecdpt);
    hRecpn = new TH2D("Recpn_STK","", 8, 0, 0.8,nevtType, evtTypemin, evtTypemax); 
    lout->Add(hRecpn);

    hCutTracknHits = new TH2D("CutTracknHits_STK","", 50, 0, 1000, nparType, parTypemin, parTypemax); 
    lout->Add(hCutTracknHits);
    hCutTrackScore = new TH2D("CutTrackScore_STK","", 50, 0, 1, nparType, parTypemin, parTypemax);
    lout->Add(hCutTrackScore);
    hCutlastTME = new TH2D("CutlastTME_STK","", 60, 0, 10, nparType, parTypemin, parTypemax); 
    lout->Add(hCutlastTME);
    hCutChi2NDF = new TH2D("CutChi2NDF_STK","", 30, 0, 500, nparType, parTypemin, parTypemax);   
    lout->Add(hCutChi2NDF);
    hCutemScore = new TH2D("CutemScore_STK","", 50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutemScore);
    hCutmichelScore = new TH2D("CutmichelScore_STK","",50, 0, 1, nparType, parTypemin, parTypemax); 
    lout->Add(hCutmichelScore);
    hCutShowerDist = new TH2D("CutShowerDist_STK","", 31, 0, 93, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerDist);
    hCutShowerIP = new TH2D("CutShowerIP_STK","", 30, 0, 30, nparType, parTypemin, parTypemax); 
    lout->Add(hCutShowerIP);

    //====================== Truth (MC only)======================//
    if(kMC){
      hTruthBeamType = new TH1I("TruthBeamType",  "", 20, -0.5, 19.5); 
      lout->Add(hTruthBeamType);
      hTruthFSParticleNumber = new TH1I("TruthFSParticleNumber", "", 160, -0.5, 159.5); 
      lout->Add(hTruthFSParticleNumber);
      hTruthFSParticleType = new TH1I("TruthFSParticleType", "", 20, -0.5, 19.5);
      lout->Add(hTruthFSParticleType);
      hTruthFSPi0Number = new TH1I("TruthFSPi0Number", "", 5, 0, 5);
      lout->Add(hTruthFSPi0Number);
      hTruthFSMultiPi0 = new TH1I("TruthFSMultiPi0", "", 5, 0, 5);
      lout->Add(hTruthFSMultiPi0);
      hTruthLeadingProtonP = new TH1D("TruthLeadingProtonP", "", 20, 0, 2);
      lout->Add(hTruthLeadingProtonP);
      hTruthSubLeadingProtonP = new TH1D("TruthSubLeadingProtonP", "", 20, 0, 2);
      lout->Add(hTruthSubLeadingProtonP);
      hTruthGammaMaxE = new TH1D("TruthGammaMaxE", "", 20, 0, 2);
      lout->Add(hTruthGammaMaxE);
      hBeamThetaRes = new TH2D("BeamTheta_RES","", 80 , 0, 60, 25, -20, 30);
      lout->Add(hBeamThetaRes);
      hBeamMomentumRes  = new TH2D("BeamMomentum_RES","", 50, 0, 2, 20, -0.5, 0.5);
      lout->Add(hBeamMomentumRes);
      hProtonThetaRes = new TH2D("ProtonTheta_RES","", 20, 0, 180, 25, -20, 30); 
      lout->Add(hProtonThetaRes);
      hProtonMomentumRes = new TH2D("ProtonMomentum_RES","",20, 0, 1.2, 20, -0.2, 0.2); 
      lout->Add(hProtonMomentumRes);
      hPiPlusThetaRes = new TH2D("PiPlusTheta_RES","", 15, 0, 180, 25, -20, 30); 
      lout->Add(hPiPlusThetaRes);
      hPiPlusMomentumRes = new TH2D("PiPlusMomentum_RES","", 20, 0, 1.2, 20, -0.2, 0.2); 
      lout->Add(hPiPlusMomentumRes);
      hShowerEnergyRes = new TH2D("ShowerEnergy_RES","", 20, 0, 1.2, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyRes);
      hShowerEnergyResRaw = new TH2D("ShowerEnergy_RES_RAW","", 20, 0, 1.2, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResRaw);
      hShowerThetaRes = new TH2D("ShowerThetaRes_RES","", 15, 0, 180, 25, -20, 30);
      lout->Add(hShowerThetaRes);
      hShowerEnergyRecVSTruth_REG = new TH2D("ShowerEnergyRecVSTruth_REG","", 20, 0, 1, 20, 0, 1);
      lout->Add(hShowerEnergyRecVSTruth_REG);
      hShowerEnergyRawRecVSTruth_REG = new TH2D("ShowerEnergyRawRecVSTruth_REG","", 20, 0, 1, 20, 0, 1);
      lout->Add(hShowerEnergyRawRecVSTruth_REG);
      hShowerEnergyResVSnHits_REG = new TH2D("ShowerEnergyResVSnHits_REG","", 50, 0, 1000, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResVSnHits_REG); 
      hShowerEnergyResVSIP_REG = new TH2D("ShowerEnergyResVSIP_REG","", 20, 0, 30, 30, -1.1, 0.7);
      lout->Add(hShowerEnergyResVSIP_REG);
      hLeadingShowerEnergyRes = new TH2D("LeadingShowerEnergy_RES","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hLeadingShowerEnergyRes);
      hSubLeadingShowerEnergyRes = new TH2D("SubLeadingShowerEnergy_RES","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hSubLeadingShowerEnergyRes);  
      hLeadingShowerEnergyResRaw = new TH2D("LeadingShowerEnergy_RES_RAW","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hLeadingShowerEnergyResRaw);
      hSubLeadingShowerEnergyResRaw = new TH2D("SubLeadingShowerEnergy_RES_RAW","", 10, 0, 1.5, 20, -1.1, 1.1);
      lout->Add(hSubLeadingShowerEnergyResRaw);
      hShowerOpenAngleRes = new TH2D("ShowerOpenAngle_RES","", 20, 0, 180, 20, -1.1, 1.1);
      lout->Add(hShowerOpenAngleRes); 
      hPi0MomentumRes = new TH2D("Pi0Momentum_RES","", 20, 0, 1, 10, -0.5, 0.5);
      lout->Add(hPi0MomentumRes);
      hPi0MassRes = new TH2D("Pi0Mass_RES","", 20, 0, 0.5, 20, -0.5, 0.5); 
      lout->Add(hPi0MassRes);

      hPi0MomentumResRaw = new TH2D("Pi0Momentum_RES_RAW","", 20, 0, 1, 10, -0.5, 0.5);
      lout->Add(hPi0MomentumResRaw);
      hPi0MassResRaw = new TH2D("Pi0Mass_RES_RAW","", 20, 0, 0.5, 20, -0.5, 0.5);
      lout->Add(hPi0MassResRaw);

      hPi0MomentumResFit = new TH2D("Pi0Momentum_RES_FIT","", 20, 0, 1, 10, -0.5, 0.5);
      lout->Add(hPi0MomentumResFit);
      hPi0MassResFit = new TH2D("Pi0Mass_RES_FIT","", 20, 0, 0.5, 20, -0.5, 0.5);
      lout->Add(hPi0MassResFit);
      hShowerOpenAngleResFit = new TH2D("ShowerOpenAngle_RES_FIT","", 20, 0, 180, 20, -1.1, 1.1);
      lout->Add(hShowerOpenAngleResFit);
    }
  }// End of IniHist

} // End of namespace
#endif
