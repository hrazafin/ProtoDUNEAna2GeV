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
  // Declare histograms
  TH1I * hEvent = 0x0;
  TH1I * hCutBeamIDPass = 0x0;
  TH1I * hCutBeamType = 0x0;
  TH1I * hCutBeamPosPass = 0x0;
  TH1I * hCutBeamEndZPass = 0x0;
  TH1D * hCutBeamEndZ = 0x0; 
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

  vector<double> *true_beam_daughter_startPx = 0x0;
  vector<double> *true_beam_daughter_startPy = 0x0;
  vector<double> *true_beam_daughter_startPz = 0x0;
  vector<double> *true_beam_daughter_startX = 0x0;
  vector<double> *true_beam_daughter_startY = 0x0;
  vector<double> *true_beam_daughter_startZ = 0x0;
  
  vector<int> *true_beam_daughter_PDG = 0x0;

  // Declare histograms
  TH1I * hTruthBeamType = 0x0;
  TH1I * hTruthFSParticleNumber = 0x0;  
  TH1I * hTruthFSParticleType = 0x0;
  TH1I * hTruthFSPi0Number = 0x0;
  TH1I * hTruthFSMultiPi0 = 0x0;
  TH1D * hTruthLeadingProtonP = 0x0;
  TH1D * hTruthSubLeadingProtonP = 0x0;
  TH1D * hTruthGammaMaxE = 0x0;

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

    tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
    tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
    tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

    tree->SetBranchAddress("true_beam_daughter_startX", &true_beam_daughter_startX);
    tree->SetBranchAddress("true_beam_daughter_startY", &true_beam_daughter_startY);
    tree->SetBranchAddress("true_beam_daughter_startZ", &true_beam_daughter_startZ);
    
    tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
    
    return tree;
  } // End of GetInputTree
  
  // Initialise histograms
  void IniHist(TList * lout, const TString tag, const bool kMC)
  {
    //====================== Reco (MC and Data)======================//
    hEvent = new TH1I("Event_"+tag, "", 10000, 0, 100000); 
    lout->Add(hEvent);
    hCutBeamIDPass = new TH1I("CutBeamIDPass_"+tag,"", 2, -0.5, 1.5); 
    lout->Add(hCutBeamIDPass);
    hCutBeamType = new TH1I("CutBeamType_"+tag,"",30, -4.5, 25.5);
    lout->Add(hCutBeamType); 
    hCutBeamPosPass = new TH1I("CutBeamPosPass_"+tag, "", 4, -0.5, 3.5); 
    lout->Add(hCutBeamPosPass);
    hCutBeamEndZPass = new TH1I("CutBeamEndZPass_"+tag, "", 4, -0.5, 3.5);
    lout->Add(hCutBeamEndZPass); 
    hCutBeamEndZ = new TH1D("CutBeamEndZ_"+tag,"",50, 0, 500);
    lout->Add(hCutBeamEndZ);
    //====================== Truth (MC only)======================//
    if(kMC){
      hTruthBeamType = new TH1I("TruthBeamType_"+tag,  "", 20, -0.5, 19.5); 
      lout->Add(hTruthBeamType);
      hTruthFSParticleNumber = new TH1I("TruthFSParticleNumber_"+tag, "", 160, -0.5, 159.5); 
      lout->Add(hTruthFSParticleNumber);
      hTruthFSParticleType = new TH1I("TruthFSParticleType_"+tag, "", 20, -0.5, 19.5);
      lout->Add(hTruthFSParticleType);
      hTruthFSPi0Number = new TH1I("TruthFSPi0Number_"+tag, "", 5, 0, 5);
      lout->Add(hTruthFSPi0Number);
      hTruthFSMultiPi0 = new TH1I("TruthFSMultiPi0_"+tag, "", 5, 0, 5);
      lout->Add(hTruthFSMultiPi0);
      hTruthLeadingProtonP = new TH1D("TruthLeadingProtonP_"+tag, "", 20, 0, 2);
      lout->Add(hTruthLeadingProtonP);
      hTruthSubLeadingProtonP = new TH1D("TruthSubLeadingProtonP_"+tag, "", 20, 0, 2);
      lout->Add(hTruthSubLeadingProtonP);
      hTruthGammaMaxE = new TH1D("TruthGammaMaxE_"+tag, "", 20, 0, 2);
      lout->Add(hTruthGammaMaxE);
    }
  }// End of IniHist

} // End of namespace
#endif
