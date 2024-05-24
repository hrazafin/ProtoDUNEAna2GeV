#include "../include/AnaCut.h"

bool AnaCut::CutBeamAllInOne(const bool kMC, bool kFill)
{ 
  // Standard procedures of protoDUNE-SP pion beam selections
   // Calculate event weight
  //kFill = true;
  //double weight = 1.0;//anaUtils.CalWeight(kMC);
  double weight = anaUtils.CalWeight(kMC);

  //Richie begin
  for (long unsigned int calo=0; calo<AnaIO::reco_beam_calo_X->size();calo++){
     if (calo==0)      AnaIO::reco_beam_startZ=AnaIO::reco_beam_calo_Z->at(0);
     if (AnaIO::reco_beam_calo_Z->at(calo)>40) break;
     if (abs(AnaIO::reco_beam_startZ-30)<abs(AnaIO::reco_beam_calo_Z->at(calo)-30)) continue;
     AnaIO::reco_beam_startX=AnaIO::reco_beam_calo_X->at(calo);
     AnaIO::reco_beam_startY=AnaIO::reco_beam_calo_Y->at(calo);
     AnaIO::reco_beam_startZ=AnaIO::reco_beam_calo_Z->at(calo);
  }

  AnaIO::reco_beam_calo_startX = AnaIO::reco_beam_startX;
  AnaIO::reco_beam_calo_startY = AnaIO::reco_beam_startY;
  AnaIO::reco_beam_calo_startZ = AnaIO::reco_beam_startZ;

  AnaIO::reco_beam_endX = AnaIO::reco_beam_calo_endX;
  AnaIO::reco_beam_endY = AnaIO::reco_beam_calo_endY;
  AnaIO::reco_beam_endZ = AnaIO::reco_beam_calo_endZ;

  double cosDirX=AnaIO::beam_inst_dirX/AnaIO::beam_inst_dirZ;
  double cosDirY=AnaIO::beam_inst_dirY/AnaIO::beam_inst_dirZ;
  double newX=(30-AnaIO::beam_inst_Z)*cosDirX+AnaIO::beam_inst_X;
  double newY=(30-AnaIO::beam_inst_Z)*cosDirY+AnaIO::beam_inst_Y;
     AnaIO::beam_inst_X=newX;
     AnaIO::beam_inst_Y=newY;
     AnaIO::beam_inst_Z=30;
  double dZCalo=AnaIO::reco_beam_endZ-AnaIO::reco_beam_startZ;
  double dXCalo=AnaIO::reco_beam_endX-AnaIO::reco_beam_startX;
  double dYCalo=AnaIO::reco_beam_endY-AnaIO::reco_beam_startY;
  double reco_beam_trkLen_30cm=TMath::Sqrt(dXCalo*dXCalo+dYCalo*dYCalo+dZCalo*dZCalo);
     AnaIO::reco_beam_trackDirX=dXCalo/reco_beam_trkLen_30cm;
     AnaIO::reco_beam_trackDirY=dYCalo/reco_beam_trkLen_30cm;
     AnaIO::reco_beam_trackDirZ=dZCalo/reco_beam_trkLen_30cm;
  //Richie end

  // Get the beam particle type/event channel/XS event using the truth level information
  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
  const int evtXStype = kMC ? anaUtils.GetFillXSEventType() : anaUtils.gkXSBmBkg;

  // Get the reconstructed beam endZ after calo correction (variable of interest to test each cut performance -- can be changed to other variable)
  const double beamEndZ = AnaIO::reco_beam_calo_endZ;

  // 1. Beam PDG cut (pion beam/other)
  const bool passBeamPDG = CutBeamPDG(kMC);
  if(kFill) AnaIO::hCutBeamPDGPass->Fill(passBeamPDG);
  if(!passBeamPDG) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutPDG, beamEndZ, parType, weight);
  if(kFill && parType == anaUtils.gkBeamPiPlus) AnaIO::hTruthMatchedBeamSample->Fill(beamEndZ, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamPDG, beamEndZ, channelType, weight);
  if(kFill && channelType == anaUtils.gkChargeExchange) AnaIO::hTruthMatchedChannelCEXSample->Fill(beamEndZ, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamPDG, beamEndZ, evtXStype, weight);
  if(kFill && evtXStype == anaUtils.gkXSSignal) AnaIO::hTruthMatchedXsEvtCEXSample->Fill(beamEndZ, weight);

  
  // Get the beam inst and front-face KE (upstream energy loss study)
  double beam_inst_KE = -999.0, true_ffKE = -999.0;
  if(kMC) anaUtils.SetBeamInstKEandFrontFaceKE(beam_inst_KE,true_ffKE,kFill);

  // Fill beam quality cut histograms (use fitted value for beam quality cut latter)
  if(kFill) anaUtils.FillBeamQualityHist();

  
  // 2. Pandora slice cut (track like/shower like)
  const bool passPandoraSlice = CutPandoraSlice();
  if(kFill) AnaIO::hCutPandoraSlicePass->Fill(passPandoraSlice);
  if(!passPandoraSlice) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutPandora, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamPandora, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamPandora, beamEndZ, evtXStype, weight);

  // 3. Calo size cut 
  const bool passCaloSize = CutCaloSize();
  if(kFill) AnaIO::hCutCaloSizePass->Fill(passCaloSize);
  if(!passCaloSize) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutCaloSize, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamCalo, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamCalo, beamEndZ, evtXStype, weight);
  

  // 4. Beam quality cut
  const bool passBeamQuality = CutBeamQuality(kMC,true,kFill);
  if(kFill) AnaIO::hCutBeamQualityPass->Fill(passBeamQuality);
  if(!passBeamQuality) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutBeamQuality, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamQuality, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamQuality, beamEndZ, evtXStype, weight);


  // 5. APA3 cut
  const bool passAPA3EndZ = CutAPA3EndZ(kMC, kFill, weight);
  if(kFill) AnaIO::hCutAPA3EndZPass->Fill(passAPA3EndZ);
  if(!passAPA3EndZ) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutAPA3, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamAPA3, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamAPA3, beamEndZ, evtXStype, weight);


  // 6. Michel Score cut (remove muons)
  const bool passMichelScore = CutMichelScore(kMC, kFill, weight);
  if(kFill) AnaIO::hCutMichelScorePass->Fill(passMichelScore);
  if(!passMichelScore) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutMichelScore, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamMichel, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamMichel, beamEndZ, evtXStype, weight);

/*
  // 7. Median dEdx cut (remove protons) replace by Chi2/DOF
  const bool passMediandEdx = CutMediandEdx(kMC);
  AnaIO::hCutMediandEdxPass->Fill(passMediandEdx);
  if(!passMediandEdx) return false;
*/

  // 7. Median dEdx cut (remove protons)
  const bool passProtonChi2 = CutProtonChi2DOF(kMC,kFill, weight);
  if(kFill) AnaIO::hCutProtonChi2Pass->Fill(passProtonChi2);
  if(!passProtonChi2) return false;
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutChi2DOF, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamChi2, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamChi2, beamEndZ, evtXStype, weight);

  // Fill upstream energy loss before BeamScraper cut
  if(kMC && kFill) anaUtils.FillUpStreamEnergyLossHistBeforeCut(beam_inst_KE, true_ffKE);

  // 8. CutBeamScraper
  /*
  const bool passBeamScraper = CutBeamScraper(kMC);
  if(kFill) AnaIO::hCutBeamScraperPass->Fill(passBeamScraper);
  if(!passBeamScraper) return false;

  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_CutBeamScraper, beamEndZ, parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_BeamScraper, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_BeamScraper, beamEndZ, evtXStype, weight);
  // Fill upstream energy loss after BeamScraper cut
  if(kMC && kFill) anaUtils.FillUpStreamEnergyLossHistAfterCut(beam_inst_KE, true_ffKE);
  */

  //========= Event Passed All Pion Beam Cuts =======//
  if(kFill) anaUtils.FillBeamVariablesAfterAllCuts(kMC, parType, channelType);

  // Selected beam rec. efficiency 
  if(kFill && parType == anaUtils.gkBeamPiPlus){
    //const TVector3 truthBeamFull = anaUtils.GetTruthBeamFull();
    AnaIO::hMatchedTruthBeamMomentum->Fill(AnaIO::true_beam_endP, weight);
  }
  
  return true;
}

bool AnaCut::CutBeamPDG(const bool kMC)
{
  if(kMC){ // MC
    // In data the beam instrumentation cannot distinguish 1GeV muon and pion, need to select both in MC
    const bool passMC_beamPDG = (AnaIO::true_beam_PDG == 211 || AnaIO::true_beam_PDG == -13);
    if(!passMC_beamPDG){
      return false;
    }
  }
  else{ // Data
    const bool passData_beamPDG = ((std::find(AnaIO::beam_inst_PDG_candidates->begin(), AnaIO::beam_inst_PDG_candidates->end(), 211)) != AnaIO::beam_inst_PDG_candidates->end());
    if(!passData_beamPDG){
      return false;
    }
    if(AnaIO::beam_inst_trigger == 8){
      return false;
    }
    if(AnaIO::beam_inst_nMomenta != 1){
      return false;
    }
    if(AnaIO::beam_inst_nTracks != 1){
      return false;
    }
    if(AnaIO::reco_reconstructable_beam_event == 0){
      return false;
    }
  }
  
  return true;
}

bool AnaCut::CutPandoraSlice()
{
  const int beam_recoType = AnaIO::reco_beam_type;
  if(beam_recoType!=13){ //13: Pandora "track like"
    return false;
  }

  return true;
}

bool AnaCut::CutCaloSize()
{
  if(AnaIO::reco_beam_calo_wire->empty()) return false;
  
  return true;
}

bool AnaCut::CutBeamQuality(bool kMC, bool DoAngleCut, bool kFill)
{
  // More info: https://docs.dunescience.org/cgi-bin/private/RetrieveFile?docid=23122&filename=ProtoDUNE_CNN_Hit_Tagging_Technote.pdf&version=1

  ///////const double cut_beamQuality_TPC_xy_low = -1.;
  const double cut_beamQuality_TPC_xy = 3.5;
  ///////const double cut_beamQuality_TPC_y = 0.;
  const double cut_beamQuality_TPC_z = 3.;
  const double cut_beamQuality_TPC_cosTheta = 0.95;
  //double weight = anaUtils.CalWeight(kMC);
  double weight = 1.; 
  //Richie parameters
  
  for (long unsigned int calo=0; calo<AnaIO::reco_beam_calo_X->size();calo++){
     if (calo==0)      AnaIO::reco_beam_startZ=AnaIO::reco_beam_calo_Z->at(0);
     if (AnaIO::reco_beam_calo_Z->at(calo)>40) break;
     if (abs(AnaIO::reco_beam_startZ-30)<abs(AnaIO::reco_beam_calo_Z->at(calo)-30)) continue;
     AnaIO::reco_beam_startX=AnaIO::reco_beam_calo_X->at(calo);
     AnaIO::reco_beam_startY=AnaIO::reco_beam_calo_Y->at(calo);
     AnaIO::reco_beam_startZ=AnaIO::reco_beam_calo_Z->at(calo);
  }
  
  AnaIO::reco_beam_calo_startX = AnaIO::reco_beam_startX;
  AnaIO::reco_beam_calo_startY = AnaIO::reco_beam_startY;
  AnaIO::reco_beam_calo_startZ = AnaIO::reco_beam_startZ;

  AnaIO::reco_beam_endX = AnaIO::reco_beam_calo_endX;
  AnaIO::reco_beam_endY = AnaIO::reco_beam_calo_endY;
  AnaIO::reco_beam_endZ = AnaIO::reco_beam_calo_endZ;
 
  double cosDirX=AnaIO::beam_inst_dirX/AnaIO::beam_inst_dirZ;
  double cosDirY=AnaIO::beam_inst_dirY/AnaIO::beam_inst_dirZ;
  double newX=(30-AnaIO::beam_inst_Z)*cosDirX+AnaIO::beam_inst_X;
  double newY=(30-AnaIO::beam_inst_Z)*cosDirY+AnaIO::beam_inst_Y;
     AnaIO::beam_inst_X=newX;
     AnaIO::beam_inst_Y=newY;
     AnaIO::beam_inst_Z=30;
  double dZCalo=AnaIO::reco_beam_endZ-AnaIO::reco_beam_startZ;
  double dXCalo=AnaIO::reco_beam_endX-AnaIO::reco_beam_startX;
  double dYCalo=AnaIO::reco_beam_endY-AnaIO::reco_beam_startY;
  double reco_beam_trkLen_30cm=TMath::Sqrt(dXCalo*dXCalo+dYCalo*dYCalo+dZCalo*dZCalo);
     AnaIO::reco_beam_trackDirX=dXCalo/reco_beam_trkLen_30cm;
     AnaIO::reco_beam_trackDirY=dYCalo/reco_beam_trkLen_30cm;
     AnaIO::reco_beam_trackDirZ=dZCalo/reco_beam_trkLen_30cm;
  
  //end parameters

  TVector3 BeamDir(AnaIO::reco_beam_calo_endX-AnaIO::reco_beam_calo_startX,
                  AnaIO::reco_beam_calo_endY-AnaIO::reco_beam_calo_startY,
                  AnaIO::reco_beam_calo_endZ-AnaIO::reco_beam_calo_startZ);
  TVector3 DetX(1,0,0);
  TVector3 DetY(0,1,0);
  TVector3 DetZ(0,0,1);

  const double thetaX = BeamDir.Angle(DetX);
  const double thetaY = BeamDir.Angle(DetY);
  const double thetaZ = BeamDir.Angle(DetZ);

  const double deltaX = kMC ? AnaIO::reco_beam_calo_startX - MCmeanStartX : AnaIO::reco_beam_calo_startX - DATAmeanStartX;
  const double deltaY = kMC ? AnaIO::reco_beam_calo_startY - MCmeanStartY : AnaIO::reco_beam_calo_startY - DATAmeanStartY;
  const double deltaZ = kMC ? AnaIO::reco_beam_calo_startZ - MCmeanStartZ : AnaIO::reco_beam_calo_startZ - DATAmeanStartZ;
  
  const double cut_xy = kMC ? sqrt(pow(deltaX/MCsigmaStartX,2)+pow(deltaY/MCsigmaStartY,2)) : sqrt(pow(deltaX/DATAsigmaStartX,2)+pow(deltaY/DATAsigmaStartY,2));

  const double cut_x = kMC ? deltaX/MCsigmaStartX : deltaX/DATAsigmaStartX;
  const double cut_y = kMC ? deltaY/MCsigmaStartY : deltaY/DATAsigmaStartY;
  const double cut_z = kMC ? deltaZ/MCsigmaStartZ : deltaZ/DATAsigmaStartZ;
  


  const double cut_theta = kMC ? cos(thetaX)*cos(MCmeanThetaX)+cos(thetaY)*cos(MCmeanThetaY)+cos(thetaZ)*cos(MCmeanThetaZ)
                               : cos(thetaX)*cos(DATAmeanThetaX)+cos(thetaY)*cos(DATAmeanThetaY)+cos(thetaZ)*cos(DATAmeanThetaZ);
  
  AnaIO::hBeamDeltaXYSigma->Fill(cut_xy);
  AnaIO::hxstartystartBeam->Fill(AnaIO::reco_beam_calo_startX, AnaIO::reco_beam_calo_startY);
  //cout << "======= cut_xy = "<< cut_xy << " cut_z = "<< cut_z << " cut_theta = "<< cut_theta << endl;
  if(kFill) AnaIO::hCutBeamDeltaXYSigma->Fill(cut_xy);
  if(kFill) AnaIO::hCutBeamDeltaZSigma->Fill(cut_z);
  if(kFill) AnaIO::hCutBeamCosTheta->Fill(cut_theta);

  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamQualityXY, cut_xy , parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamQualityX, cut_x , parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamQualityY, cut_y , parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamQualityZ, cut_z , parType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamQualityTheta, cut_theta , parType, weight);
  
  if(cut_xy > cut_beamQuality_TPC_xy) return false;

  /////if(abs(cut_y) > cut_beamQuality_TPC_y) return false;
  if(abs(cut_z) > cut_beamQuality_TPC_z) return false;

  if(DoAngleCut){
    if(cut_theta < cut_beamQuality_TPC_cosTheta) return false;
  }

  //if(CutBeamInstQuality(kMC) == false) return false;

  if(AnaIO::beam_inst_P > 2.2 || AnaIO::beam_inst_P < 1.8) return false;

  return true;
}

bool AnaCut::CutBeamInstQuality(bool kMC)
{
  const double cut_beamQuality_TPC_xy = 2.5;
  // Richie begin
  for (long unsigned int calo=0; calo<AnaIO::reco_beam_calo_X->size();calo++){
     if (calo==0)      AnaIO::reco_beam_startZ=AnaIO::reco_beam_calo_Z->at(0);
     if (AnaIO::reco_beam_calo_Z->at(calo)>40) break;
     if (abs(AnaIO::reco_beam_startZ-30)<abs(AnaIO::reco_beam_calo_Z->at(calo)-30)) continue;
     AnaIO::reco_beam_startX=AnaIO::reco_beam_calo_X->at(calo);
     AnaIO::reco_beam_startY=AnaIO::reco_beam_calo_Y->at(calo);
     AnaIO::reco_beam_startZ=AnaIO::reco_beam_calo_Z->at(calo);
  }
  AnaIO::reco_beam_calo_startX = AnaIO::reco_beam_startX;
  AnaIO::reco_beam_calo_startY = AnaIO::reco_beam_startY;
  AnaIO::reco_beam_calo_startZ = AnaIO::reco_beam_startZ;

  AnaIO::reco_beam_endX = AnaIO::reco_beam_calo_endX;
  AnaIO::reco_beam_endY = AnaIO::reco_beam_calo_endY;
  AnaIO::reco_beam_endZ = AnaIO::reco_beam_calo_endZ;

  double cosDirX=AnaIO::beam_inst_dirX/AnaIO::beam_inst_dirZ;
  double cosDirY=AnaIO::beam_inst_dirY/AnaIO::beam_inst_dirZ;
  double newX=(30-AnaIO::beam_inst_Z)*cosDirX+AnaIO::beam_inst_X;
  double newY=(30-AnaIO::beam_inst_Z)*cosDirY+AnaIO::beam_inst_Y;
     AnaIO::beam_inst_X=newX;
     AnaIO::beam_inst_Y=newY;
     AnaIO::beam_inst_Z=30;
  double dZCalo=AnaIO::reco_beam_endZ-AnaIO::reco_beam_startZ;
  double dXCalo=AnaIO::reco_beam_endX-AnaIO::reco_beam_startX;
  double dYCalo=AnaIO::reco_beam_endY-AnaIO::reco_beam_startY;
  double reco_beam_trkLen_30cm=TMath::Sqrt(dXCalo*dXCalo+dYCalo*dYCalo+dZCalo*dZCalo);
     AnaIO::reco_beam_trackDirX=dXCalo/reco_beam_trkLen_30cm;
     AnaIO::reco_beam_trackDirY=dYCalo/reco_beam_trkLen_30cm;
     AnaIO::reco_beam_trackDirZ=dZCalo/reco_beam_trkLen_30cm;
  //Richie end
  const double deltaX = kMC ? AnaIO::beam_inst_X - MCmeanInstX : AnaIO::beam_inst_X - DATAmeanInstX;
  const double deltaY = kMC ? AnaIO::beam_inst_Y - MCmeanInstY : AnaIO::beam_inst_Y - DATAmeanInstY;
  
  const double cut_xy = kMC ? sqrt(pow(deltaX/MCsigmaInstX,2)+pow(deltaY/MCsigmaInstY,2)) : sqrt(pow(deltaX/DATAsigmaInstX,2)+pow(deltaY/DATAsigmaInstY,2));
  const double cut_x = kMC ? deltaX/MCsigmaInstX : deltaX/DATAsigmaInstX;
  const double cut_y = kMC ? deltaY/MCsigmaInstY : deltaY/DATAsigmaInstY;

  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  plotUtils.FillHist(AnaIO::hCutBeamQualityInstXY, cut_xy , parType);
  plotUtils.FillHist(AnaIO::hCutBeamQualityInstX, cut_x , parType);
  plotUtils.FillHist(AnaIO::hCutBeamQualityInstY, cut_y , parType);
  
  if(cut_xy > cut_beamQuality_TPC_xy) return false;

  return true;
}


bool AnaCut::CutAPA3EndZ(const bool kMC, bool kFill, const double weight){

  //AnaIO::hCutBeamEndZ->Fill(AnaIO::reco_beam_calo_endZ);
  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  const double endZ = AnaIO::reco_beam_calo_endZ;
  if(kFill) plotUtils.FillHist(AnaIO::hCutAPA3EndZ, endZ , parType, weight);
  //if(endZ >= cut_EndZ_APA3){
  if(endZ >= cut_EndZ_APA3){
    return false;
  }

  return true;
}

bool AnaCut::CutMichelScore(const bool kMC, bool kFill, const double weight)
{
  if( AnaIO::reco_beam_vertex_michel_score == 0) return true;
  //const double MichelScorePerHits = AnaIO::reco_beam_vertex_michel_score/AnaIO::reco_beam_vertex_nHits;
  const double MichelScoreWBC = AnaIO::reco_beam_vertex_michel_score_weight_by_charge;

  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamMichelScoreNhits, MichelScoreWBC , parType, weight);
  // There is a muon candidate 
  if(MichelScoreWBC > cut_MichelScore) return false;

  return true;
}
/* // Replaced by the chi2/dof cut
bool AnaCut::CutMediandEdx(const bool kMC)
{
  auto median_it = AnaIO::reco_beam_calibrated_dEdX_SCE->begin() + AnaIO::reco_beam_calibrated_dEdX_SCE->size()/2;
  std::nth_element(AnaIO::reco_beam_calibrated_dEdX_SCE->begin(), median_it, AnaIO::reco_beam_calibrated_dEdX_SCE->end());
  const double median = (*AnaIO::reco_beam_calibrated_dEdX_SCE)[AnaIO::reco_beam_calibrated_dEdX_SCE->size()/2];

  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  plotUtils.FillHist(AnaIO::hCutBeamMediandEdx, median , parType);
  
  if(kMC){
    AnaIO::hCutBeamMedianType->Fill(parType);
  }

  if(median > 2.4) return false;

  return true;
}
*/
bool AnaCut::CutProtonChi2DOF(const bool kMC, bool kFill, const double weight)
{
  const double Chi2DOF = AnaIO::reco_beam_Chi2_proton/AnaIO::reco_beam_Chi2_ndof;
  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  if(kFill) plotUtils.FillHist(AnaIO::hCutBeamProtonChi2DOF, Chi2DOF , parType, weight);
  
  if(kMC){
    if(kFill) AnaIO::hCutBeamProtonChi2DOFType->Fill(parType, weight);
  }

  if(Chi2DOF < 80) return false;

  return true;
}


bool AnaCut::CutBeamScraper(const bool kMC){

  // Define beam plug circle in [cm] unit
  double center_x = 0.;
  double center_y = 0.;
  double radious = 0.;
  double N_sigma = 0.;
  
  // for 2 GeV :
  if(kMC){
    center_x = -3615;
    center_y = 422.;
    radious = 4.8;
    N_sigma = 1.4;
  }
  else{
    center_x = -32.5;
    center_y = 424.7;
    radious = 4.8;
    N_sigma = 1.2;
  }
  //default for 1GeV
  /*if(kMC){
    center_x = -29.6;
    center_y = 422.;
    radious = 4.8;
    N_sigma = 1.4;
  }
  else{
    center_x = -32.16;
    center_y = 422.7;
    radious = 4.8;
    N_sigma = 1.2;
  }*/

  double this_distance = sqrt( pow(AnaIO::beam_inst_X - center_x , 2.) + pow(AnaIO::beam_inst_Y - center_y, 2.) );
  if(this_distance > radious * N_sigma) return false;

  return true;
}

// --------------------------------- Legacy Cuts (not used anymore) ----------------------- //

bool AnaCut::CutMCBeamPos()
{
  return Manual_beamPos_mc(AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                           AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                           AnaIO::true_beam_startDirX, AnaIO::true_beam_startDirY, AnaIO::true_beam_startDirZ,
                           AnaIO::true_beam_startX, AnaIO::true_beam_startY, AnaIO::true_beam_startZ);
}

bool AnaCut::CutDataBeamPos()
{
  return Manual_beamPos_data(AnaIO::event, AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                             AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                             AnaIO::beam_inst_X, AnaIO::beam_inst_Y, AnaIO::beam_inst_dirX, AnaIO::beam_inst_dirY,
                             AnaIO::beam_inst_dirZ, AnaIO::beam_inst_nMomenta, AnaIO::beam_inst_nTracks);
}

bool AnaCut::Manual_beamPos_mc(const double beam_startX, const double beam_startY, const double beam_startZ, const double beam_dirX, const double beam_dirY,   const double beam_dirZ, const double true_dirX,   const double true_dirY, const double true_dirZ,   const double true_startX, const double true_startY, const double true_startZ)
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.h

  //For MC from Owen Goodwins studies
  const double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
  const double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

  const double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  const double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  const double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;

  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;

  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;

  if ( cos < coslow)
    return false;
    
  return true;
}

bool AnaCut::Manual_beamPos_data(const int event,            const double data_startX,
                         const double data_startY,   const double data_startZ,
                         const double data_dirX,     const double data_dirY,
                         const double data_dirZ,     const double beam_inst_X,
                         const double beam_inst_Y,     const double beam_inst_dirX,
                         const double beam_inst_dirY,  const double beam_inst_dirZ,
                         const int beam_inst_nMomenta, const int beam_inst_nTracks)
{

  //For Data from Owen Goodwin
  const double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
  const double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

  const double deltaX = data_startX - beam_inst_X;
  const double deltaY = data_startY - beam_inst_Y;
  const double cos = beam_inst_dirX*data_dirX + beam_inst_dirY*data_dirY +
               beam_inst_dirZ*data_dirZ;

  if(beam_inst_nMomenta != 1 || beam_inst_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

}

// --------------------------------- Legacy Cuts (not used anymore) ----------------------- //


bool AnaCut::CutTopology(const bool kMC, double & pi0KineticE, double & pi0costheta, bool kFill)
{ 
  // Calculate event weight
  //double weight = 1.0;//anaUtils.CalWeight(kMC);
  double weight = anaUtils.CalWeight(kMC);
  
  pi0KineticE = -999;
  pi0costheta = -999;
  // Count reco final state particles and determine types 
  CountPFP(kMC,kFill);
 
  /*cout << "Rec nproton: " << nproton << " Rec npi0shower: " << npi0shower << " Rec npiplus: " << npiplus << " Rec nmichel: " << nmichel << 
  " Rec npi0: " << npi0 << endl;
  cout << "anaUtils.RecProtonLTVet: " << anaUtils.RecProtonLTVet.P() << " anaUtils.TruthProtonLTVet: " << anaUtils.TruthProtonLTVet.P() << endl;
  cout << "anaUtils.RecPi0LTVet: " << anaUtils.RecPi0LTVet.P() << " anaUtils.TruthPi0LTVet: " << anaUtils.TruthPi0LTVet.P() << endl;
  cout << "anaUtils.RecPi0LTVet Mass: " << anaUtils.RecPi0LTVet.M() << " anaUtils.TruthPi0LTVet Mass: " << anaUtils.TruthPi0LTVet.M() << endl;
*/
  const double beamEndZ = AnaIO::reco_beam_calo_endZ;
  const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
  const int evtXStype = kMC ? anaUtils.GetFillXSEventType() : anaUtils.gkXSBmBkg;

  // ----------------------- Do cuts below ------------------------ // 
  // Get event type
  const int evtType = anaUtils.GetFillEventType();
  // Proton
  if(kFill) plotUtils.FillHist(AnaIO::hCutnproton, nproton, evtType);  
  //if(nproton!=1) return false;
  // Showers  
  if(kFill) plotUtils.FillHist(AnaIO::hCutnshower, npi0shower, evtType); 
  // Cex signal 
  if(npi0shower < 2) return false;
  // Pion prod. background (No pi0) 
  //if(npi0shower > 1) return false;
  // Pion prod. background (One pi0) 
  //if(npi0shower < 2) return false;
  
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_TwoShowers, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_TwoShowers, beamEndZ, evtXStype, weight);

  // Piplus
  if(kFill) plotUtils.FillHist(AnaIO::hCutnpiplus, npiplus, evtType);  
  // Cex signal 
  if(npiplus!=0) return false;
  // Pion prod. background (No pi0) 
  //if(npiplus==0) return false;
  // Pion prod. background (One pi0) 
  //if(npiplus==0) return false;

  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_NoPiPlus, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_NoPiPlus, beamEndZ, evtXStype, weight);

  //cout << "npiplus: " << npiplus << endl;
  //cout << "nmichel: " << nmichel << endl;

  // Michel electron
  if(kFill) plotUtils.FillHist(AnaIO::hCutnmichel, nmichel, evtType); 
  // Cex signal  
  if(nmichel!=0) return false;
  // Pion prod. background (No pi0) No michel restriction
  // Pion prod. background (One pi0) No michel restriction

  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_NoMichel, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_NoMichel, beamEndZ, evtXStype, weight);

  // PiZero cut
  if(kFill) plotUtils.FillHist(AnaIO::hCutnpi0, npi0, evtType);
  // Cex signal    
  if(npi0 != 1) return false;
  // Pion prod. background (No pi0)    
  //if(npi0 != 0) return false;
  // Pion prod. background (One pi0)    
  //if(npi0 < 1) return false;
  
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_OnePi0, beamEndZ, channelType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hBeamEndZ_XsEvt_OnePi0, beamEndZ, evtXStype, weight);

  if(channelType == anaUtils.gkInelastic) {
    const int size = (*AnaIO::true_beam_daughter_PDG).size();
    for (int ii = 0; ii < size; ++ii) {
      if (abs((*AnaIO::true_beam_daughter_PDG)[ii]) == 211) {
        if(kMC && kFill) AnaIO::hPiMom_InelasticChannels->Fill((*AnaIO::true_beam_daughter_startP)[ii], weight);
      }
    }
  }

  // KF Pass
  if(kFill) plotUtils.FillHist(AnaIO::hCutKFPass, anaUtils.GoodFit, evtType, weight);
  //if(anaUtils.GoodFit == false) return false;

  
  // Before KF
  if(kFill) plotUtils.FillHist(AnaIO::hRecPi0TotalEafterSel_Default, PiZeroVec.E() - AnaFunctions::PiZeroMass(), evtXStype, weight);

  // Rewrite the pi zero vector this time need to turn on KF to improve pi0 energy resolution
  PiZeroVec = anaUtils.GetRecPiZeroFromShowers(OA,kMC,false,true,truthPi0Type);
  // After KF
  if(kFill) plotUtils.FillHist(AnaIO::hRecPi0TotalEafterSel_Fitted, PiZeroVec.E() - AnaFunctions::PiZeroMass(), evtXStype, weight);

  double pi0KE = PiZeroVec.E() - AnaFunctions::PiZeroMass();
  double costheta = PiZeroVec.Pz()/PiZeroVec.P();
  pi0KineticE = pi0KE; 
  pi0costheta = costheta;
  //if(pi0costheta < 0. && pi0costheta > -0.2)  cout << "Cut pi0cos: " << pi0costheta << " PiZeroVec.P(): " << PiZeroVec.P() << " PiZeroVec.Pz(): " << PiZeroVec.Pz() << endl; 
  if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_AfterTOP_EVT, pi0KE, evtXStype, weight);

/* //fake-data control
  if(kFill) {

    bool isFakeData = anaUtils.IsFakeData();

    if(isFakeData) plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_AfterTOP_EVT, pi0KE, evtXStype);
    if(!isFakeData) plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_AfterTOP_EVT, pi0KE, evtXStype);

  }
*/
  return true;
}

void AnaCut::CountPFP(const bool kMC, const bool kFill)
{ 
  // Calculate event weight
  //double weight = 1.0;//anaUtils.CalWeight(kMC);
  double weight = anaUtils.CalWeight(kMC);

  // Get the size of reco final state particles
  const int recsize = AnaIO::reco_daughter_PFP_ID->size();
  // Initialize number of reco final state particles (beam reco daughter particles counter)
  nproton = 0;
  npiplus = 0;
  nshower = 0;
  npi0shower = 0;
  nmichel = 0;
  npi0 = 0;
  int nPFP = 0;
  // Need to clear vector for each event
  anaUtils.CleanShowerArray();

  anaUtils.RecPi0LTVet = TLorentzVector();
  anaUtils.TruthPi0LTVet = TLorentzVector();
  anaUtils.RecProtonLTVet = TLorentzVector();
  anaUtils.TruthProtonLTVet = TLorentzVector();

  vector<double> bufferPiZeroGammamom;
  vector<double> bufferProtonmom;
  vector<double> bufferPiPlusmom;

  double ProtonMom = -999;
  // Loop over all reco FS particles
  for(int ii=0; ii<recsize; ii++){
    // Get the truth particle type of this reco particle
    truthParticleType = kMC ? anaUtils.GetTruthParticleInfoFromRec(ii) : anaUtils.gkOthers;

    const int nhits_coll = (*AnaIO::reco_daughter_PFP_nHits_collection)[ii];
    if(kFill) plotUtils.FillHist(AnaIO::hNhitsCollection,nhits_coll, truthParticleType, weight);

    const double trackScore_coll = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];
    const double emScore_coll = (*AnaIO::reco_daughter_PFP_emScore_collection)[ii];
    const double michelScore_coll = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];

    if(kFill) plotUtils.FillHist(AnaIO::htrackScoreCollection,trackScore_coll, truthParticleType, weight);
    if(kFill) plotUtils.FillHist(AnaIO::hemScoreCollection,emScore_coll, truthParticleType, weight);
    if(kFill) plotUtils.FillHist(AnaIO::hmichelScoreCollection,michelScore_coll, truthParticleType, weight);

    recParticleType = -999;
    // Proton candidates selection
    if(IsProton(ii,kMC,weight,kFill)){
      nproton++;
      recParticleType = anaUtils.gkProton;
      // Fill rec and truth-matching info
      if(kFill) anaUtils.FillFSParticleKinematics(ii, kMC, truthParticleType, recParticleType);

      const TVector3 recMomProton = anaUtils.GetRecTrackVectLab(ii, true);
      TLorentzVector rec4MomProton;
      const TVector3 truthMomProton = anaUtils.GetTruthMatchedTrackVectLab(ii);
      TLorentzVector truth4MomProton;
      rec4MomProton.SetVectM(recMomProton,AnaFunctions::ProtonMass());
      truth4MomProton.SetVectM(truthMomProton,AnaFunctions::ProtonMass());
      if(rec4MomProton.P() > ProtonMom) {
        anaUtils.RecProtonLTVet = rec4MomProton;
        anaUtils.TruthProtonLTVet = truth4MomProton;
        ProtonMom = rec4MomProton.P();
      }
      if(kFill && truthParticleType == anaUtils.gkProton){
        bufferProtonmom.push_back(truth4MomProton.P());
        AnaIO::hMatchedTruthProtonMomentum->Fill(truth4MomProton.P(), weight);
      }
    }
    // Shower candidates selection
    if(IsShower(ii,kMC,weight,kFill)){
      nshower++;
      // Select good shower candidates for reconstructing pi0   
      if(IsPiZeroShower(ii,kMC,weight,kFill)){
        npi0shower++;
        recParticleType = anaUtils.gkGamma;

        // Fill rec and truth-matching info
        if(kFill) anaUtils.FillFSParticleKinematics(ii, kMC, truthParticleType, recParticleType);
        
        const double startX = (*AnaIO::reco_daughter_allShower_startX)[ii];
        const double startY = (*AnaIO::reco_daughter_allShower_startY)[ii];
        const double startZ = (*AnaIO::reco_daughter_allShower_startZ)[ii];
        const double showerLength = (*AnaIO::reco_daughter_allShower_len)[ii];  
        
        if(kFill) plotUtils.FillHist(AnaIO::hRecShowerLength, showerLength, truthParticleType, weight);
        if(kFill) plotUtils.FillHist(AnaIO::hShowerStartX, startX, truthParticleType, weight);
        if(kFill) plotUtils.FillHist(AnaIO::hShowerStartY, startY, truthParticleType, weight);
        if(kFill) plotUtils.FillHist(AnaIO::hShowerStartZ, startZ, truthParticleType, weight);

        const int nhits_coll = (*AnaIO::reco_daughter_PFP_nHits_collection)[ii];
        if(kFill) plotUtils.FillHist(AnaIO::hShowernHitsColl, nhits_coll, truthParticleType, weight);

        const double emScore = (*AnaIO::reco_daughter_PFP_emScore_collection)[ii];
        const double showerE = (*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3;
        // 2D map positionVSemScore
        if(kFill) plotUtils.FillHist(AnaIO::hCutemScoreVSStartX,startX,emScore, weight);
        if(kFill) plotUtils.FillHist(AnaIO::hCutemScoreVSStartY,startY,emScore, weight);
        if(kFill) plotUtils.FillHist(AnaIO::hCutemScoreVSStartZ,startZ,emScore, weight);

        if(kFill) plotUtils.FillHist(AnaIO::hCutemScoreVSEnergy,showerE,emScore, weight);

        if(kFill) plotUtils.FillHist(AnaIO::hCutemScore_AfterCut,emScore,truthParticleType, weight);

        if(kFill) plotUtils.FillHist(AnaIO::hCutStartXVSStartY,startX,startY, weight);
        if(kFill) plotUtils.FillHist(AnaIO::hCutStartZVSStartX,startZ,startX, weight);

        const TLorentzVector truthShowerMom = anaUtils.GetTruthMatchedShowerLTVectLab(ii);
        if(kFill && truthParticleType == anaUtils.gkGamma) bufferPiZeroGammamom.push_back(truthShowerMom.E());

      }
    }

    // Piplus candidates selection
    if(IsPiplus(ii,kMC,weight,kFill)){
      npiplus++;
      recParticleType = anaUtils.gkPiPlus;
      // Fill rec and truth-matching info
      if(kFill) anaUtils.FillFSParticleKinematics(ii, kMC, truthParticleType, recParticleType);

      const TVector3 truthMomPiPlus = anaUtils.GetTruthMatchedTrackVectLab(ii);
      TLorentzVector truth4MomPiPlus;
      truth4MomPiPlus.SetVectM(truthMomPiPlus,AnaFunctions::PionMass());
      
      if(kFill && truthParticleType == anaUtils.gkPiPlus){
        bufferPiPlusmom.push_back(truth4MomPiPlus.P());
        AnaIO::hMatchedTruthPiPlusMomentum->Fill(truth4MomPiPlus.P(), weight);
      }
      
    }

    // Michel candidates selection
    if(IsMichel(ii,kMC,weight,kFill)){
      nmichel++;
    }
    // Count total number of FS particles
    nPFP++; 
  }
  if(recsize!=nPFP) cout << "CountPFP not looping all FS particles!!" << endl;

  // Leading shower enregy (truth-matched)
  int size = bufferPiZeroGammamom.size();
  if(kFill && npi0shower > 0 && size > 0){
    double p = -999;
    int leadingID = -999;
    for(int id = 0; id < size; id++){
      if(bufferPiZeroGammamom[id] > p){
        p = bufferPiZeroGammamom[id];
        leadingID = id;
      }
    }
    AnaIO::hMatchedTruthLeadingShowerEnergy->Fill(bufferPiZeroGammamom[leadingID], weight);
  }

  // Leading proton start p (truth-matched)
  int size_np = bufferProtonmom.size();
  if(kFill && nproton > 0 && size_np > 0){
    double p = -999;
    int leadingID = -999;
    for(int id = 0; id < size_np; id++){
      if(bufferProtonmom[id] > p){
        p = bufferProtonmom[id];
        leadingID = id;
      }
    }
    AnaIO::hMatchedTruthLeadingProtonMomentum->Fill(bufferProtonmom[leadingID], weight);
  }

  // Leading piplus start p (truth-matched)
  int size_npi = bufferPiPlusmom.size();
  if(kFill && npiplus > 0 && size_npi > 0){
    double p = -999;
    int leadingID = -999;
    for(int id = 0; id < size_npi; id++){
      if(bufferPiPlusmom[id] > p){
        p = bufferPiPlusmom[id];
        leadingID = id;
      }
    }
    AnaIO::hMatchedTruthLeadingPiPlusMomentum->Fill(bufferPiPlusmom[leadingID], weight);
  }
 
  // Pi zero candidates selection (!one per event! Get pi0 from at least two pi0 showers)
  if(IsPizero(kMC,kFill,weight)){
    npi0++;
    
    if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After, PiZeroVec.E(), truthPi0Type, weight);
    if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Theta_OVERLAY_After, PiZeroVec.Theta()*TMath::RadToDeg(), truthPi0Type, weight);
    if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Phi_OVERLAY_After, PiZeroVec.Phi()*TMath::RadToDeg(), truthPi0Type, weight);

    double pi0KE = PiZeroVec.E() - AnaFunctions::PiZeroMass();
    const int evtType = anaUtils.GetFillXSEventType();

    if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_EVT, pi0KE, evtType, weight);
    if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_Kinetic, pi0KE, truthPi0Type, weight);

  } 
  //if(npi0shower > 1 && nproton > 0) printf("CountPFP PFP size %d nlooped %d nshower %d npi0shower %d nmichel %d npiplus %d nproton %d\n", recsize, nPFP, nshower, npi0shower, nmichel, npiplus, nproton);

}

bool AnaCut::IsProton(const int ii, const bool kMC, const double weight, const bool kFill)
{
  // Proton candidate must be a track like particle
  if(!IsTrack(ii,kMC,weight,kFill)) return false; 
  // Proton candidate must pass the folowing cuts 
  if(!PassProtonSubPID(ii,weight,kFill)) return false;

  return true;
}

bool AnaCut::IsTrack(const int ii,const bool kMC, const double weight, const bool kFill)
{ 
  // Get the number of hits of this reco particle
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  const int nhits_coll = (*AnaIO::reco_daughter_PFP_nHits_collection)[ii];
  // Get the track score
  const double trackScore = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];
  // Get the shower purity and completeness
  const double Purity = (*AnaIO::reco_daughter_PFP_true_byHits_purity)[ii];
  const double Completeness = (*AnaIO::reco_daughter_PFP_true_byHits_completeness)[ii];

  const double ProtonMom = (*AnaIO::reco_daughter_allTrack_momByRange_alt_proton)[ii];

  const double startX = (*AnaIO::reco_daughter_allTrack_startX)[ii];
  const double startY = (*AnaIO::reco_daughter_allTrack_startY)[ii];
  const double startZ = (*AnaIO::reco_daughter_allTrack_startZ)[ii];
  if(kFill) plotUtils.FillHist(AnaIO::hAllTrackCutStartXVSStartY,startX,startY, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hAllTrackCutStartZVSStartX,startZ,startX, weight);


  // Fill Cut histograms for those two variables
  if(kFill) plotUtils.FillHist(AnaIO::hCutTracknHits, nhits, truthParticleType, weight);  
  if(kFill) plotUtils.FillHist(AnaIO::hCutTrackScore, trackScore, truthParticleType, weight);
  if(kFill && kMC) plotUtils.FillHist(AnaIO::hTrackCompleteness, Completeness, truthParticleType, weight);

  if(kFill) plotUtils.FillHist(AnaIO::hProtonMom_NoCut,ProtonMom,truthParticleType, weight);
  if(kFill && truthParticleType == anaUtils.gkProton) AnaIO::hTruthMatchedProtonSample->Fill(ProtonMom,weight);

  // Check if this reco particle has forced track ID
  if((*AnaIO::reco_daughter_allTrack_ID)[ii]==-1) return false;
  // Cut on track score
  if(trackScore <= 0.5){
    if(kFill) AnaIO::hCutDaughterTrackScorePass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterTrackScorePass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hProtonMom_CutTrackScore,ProtonMom,truthParticleType, weight);

  if(kFill) {
    if(kMC) plotUtils.FillHist(AnaIO::hTrackPurityVSnHits, nhits, Purity, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hTrackCompletenessVSnHits, nhits, Completeness, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hTrackPurityVSnHitsColl, nhits_coll, Purity, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hTrackCompletenessVSnHitsColl, nhits_coll, Completeness, weight);
  }
  // Cut on number of hits
  if(nhits <= 20){
    if(kFill) AnaIO::hCutDaughterTracknHitsPass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterTracknHitsPass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hProtonMom_CutnHits,ProtonMom,truthParticleType, weight);

  return true;
}


bool AnaCut::IsPionTrack(const int ii,const bool kMC, const double weight, const bool kFill)
{
  // Get the number of hits of this reco particle
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  
  // Get the track score
  const double trackScore = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];
  const double PiPlusMom = (*AnaIO::reco_daughter_allTrack_momByRange_alt_muon)[ii];

  if(kFill) plotUtils.FillHist(AnaIO::hPiPlusMom_NoCut,PiPlusMom,truthParticleType, weight);
  if(kFill && truthParticleType == anaUtils.gkPiPlus) AnaIO::hTruthMatchedPiPlusSample->Fill(PiPlusMom,weight);

  // Check if this reco particle has forced track ID
  if((*AnaIO::reco_daughter_allTrack_ID)[ii]==-1) return false;

  // Cut on track score
  if(trackScore <= 0.5){
    if(kFill) AnaIO::hCutDaughterPionTrackScorePass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterPionTrackScorePass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hPiPlusMom_CutTrackScore,PiPlusMom,truthParticleType, weight);
  
  // Cut on number of hits
  if(nhits <= 20){
    if(kFill) AnaIO::hCutDaughterPionTracknHitsPass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterPionTracknHitsPass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hPiPlusMom_CutnHits,PiPlusMom,truthParticleType, weight);

  return true;
}


bool AnaCut::PassProtonSubPID(const int ii, const double weight, const bool kFill)
{
  // Get the Truncated mean dEdx backwards of this reco particle
  const double lastTME = anaUtils.GetFSParticleTME(ii,false);
  // Get the Chi2/NDF value
  const double Chi2NDF = anaUtils.GetChi2NDF(ii);

  const double ProtonMom = (*AnaIO::reco_daughter_allTrack_momByRange_alt_proton)[ii];

  // Fill Cut histograms for those two variables
  if(kFill) plotUtils.FillHist(AnaIO::hCutlastTME, lastTME, truthParticleType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hCutChi2NDF, Chi2NDF, truthParticleType, weight);
  
  // Use these two variables for proton selections
  if(Chi2NDF < 70 || lastTME > 2.8){
  //if(Chi2NDF < 60 || lastTME > 3){
    if(kFill) AnaIO::hCutDaughterProtonSubPIDPass->Fill(true);
    if(kFill) plotUtils.FillHist(AnaIO::hProtonMom_CutSubPID,ProtonMom,truthParticleType, weight);
    return true;
  }
  else {
    if(kFill) AnaIO::hCutDaughterProtonSubPIDPass->Fill(false);
  }
  return false;
}

bool AnaCut::PassPionSubPID(const int ii, const double weight, const bool kFill)
{
  // Get the Truncated mean dEdx backwards of this reco particle
  const double lastTME = anaUtils.GetFSParticleTME(ii,false);
  // Get the Chi2/NDF value
  const double Chi2NDF = anaUtils.GetChi2NDF(ii);
  const double PiPlusMom = (*AnaIO::reco_daughter_allTrack_momByRange_alt_muon)[ii];
  
  // Use these two variables for pion selections
  if((lastTME < 2.8 && lastTME > 0.5) || (lastTME < 3.4 && Chi2NDF >70)){
    if(kFill) AnaIO::hCutDaughterPionSubPIDPass->Fill(true);
    if(kFill) plotUtils.FillHist(AnaIO::hPiPlusMom_CutSubPID,PiPlusMom,truthParticleType, weight);

    return true;
  }
  else {
    if(kFill) AnaIO::hCutDaughterPionSubPIDPass->Fill(false);
  }
  return false;
}


bool AnaCut::IsPiplus(const int ii, const bool kMC, const double weight, const bool kFill)
{
  // PiPlus candidate must be a track like particle  
  //if(!IsTrack(ii,kMC)) return false;
  if(!IsPionTrack(ii,kMC,weight,kFill)) return false;
  // Assume particle not pass ProtonSubPID is piplus 
  if(!PassPionSubPID(ii,weight,kFill)) return false;

  
  return true;  
}

bool AnaCut::IsShower(const int ii, const bool kMC, const double weight, const bool kFill)
{
  
  // Get the em score and nhits of this particle
  const double emScore = (*AnaIO::reco_daughter_PFP_emScore_collection)[ii];
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  const int nhits_coll = (*AnaIO::reco_daughter_PFP_nHits_collection)[ii];
  const double showerE = (*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3;
  //const int typePandora = (*AnaIO::reco_daughter_pandora_type)[ii];

  //if(typePandora != 11){
    //if(kFill) AnaIO::hCutDaughterPandoraShowerPass->Fill(false);
    //return false;
  //}
  //else {
  //  if(kFill) AnaIO::hCutDaughterPandoraShowerPass->Fill(true);
  //}
  // Get the shower purity and completeness
  const double Purity = (*AnaIO::reco_daughter_PFP_true_byHits_purity)[ii];
  const double Completeness = (*AnaIO::reco_daughter_PFP_true_byHits_completeness)[ii];

  // Fill Cut histogram
  if(kFill) plotUtils.FillHist(AnaIO::hCutemScore, emScore, truthParticleType, anaUtils.CalWeight(kMC));
  // Check SCE affect shower score
  const double startX = (*AnaIO::reco_daughter_allShower_startX)[ii];
  const double startY = (*AnaIO::reco_daughter_allShower_startY)[ii];
  const double startZ = (*AnaIO::reco_daughter_allShower_startZ)[ii];
  if(kFill){
    if(startX < -50) plotUtils.FillHist(AnaIO::hCutemScore_R1, emScore, truthParticleType, weight);
    else plotUtils.FillHist(AnaIO::hCutemScore_R2, emScore, truthParticleType, weight);

    if(startY > 415) plotUtils.FillHist(AnaIO::hCutemScore_R3, emScore, truthParticleType, weight);
    else plotUtils.FillHist(AnaIO::hCutemScore_R4, emScore, truthParticleType, weight);
  }
  //Herilala Hist
  AnaIO::hxstartystart->Fill(startX,startY);
  // 2D map positionVSemScore
  if(kFill) plotUtils.FillHist(AnaIO::hAllCutemScoreVSStartX,startX,emScore, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hAllCutemScoreVSStartY,startY,emScore, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hAllCutemScoreVSStartZ,startZ,emScore, weight);

  if(kFill) plotUtils.FillHist(AnaIO::hAllCutStartXVSStartY,startX,startY, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hAllCutStartZVSStartX,startZ,startX, weight);

  if(kFill) plotUtils.FillHist(AnaIO::hShowernHitsVSEnergy, nhits, showerE, weight);
  //if(kFill) plotUtils.FillHist(AnaIO::hShowernHitsColl, nhits_coll, truthParticleType);
  
  if(kMC && kFill) plotUtils.FillHist(AnaIO::hShowerCompleteness, Completeness, truthParticleType, weight);

  
  if(kFill) plotUtils.FillHist(AnaIO::hShowerEnergy_NoCut, showerE, truthParticleType, weight);
  if(kFill && truthParticleType == anaUtils.gkGamma) AnaIO::hTruthMatchedShowerSample->Fill(showerE, weight);
  
  if((*AnaIO::reco_daughter_allShower_ID)[ii]==-1) return false;

  // Cut on em score
  if(emScore <= 0.5){
    if(kFill) AnaIO::hCutDaughterShowerScorePass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterShowerScorePass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hShowerEnergy_CutEMscore, showerE, truthParticleType, weight);

  if(kFill) {
    if(kMC) plotUtils.FillHist(AnaIO::hShowerPurityVSnHits, nhits, Purity, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hShowerCompletenessVSnHits, nhits, Completeness, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hShowerPurityVSEnergy, showerE, Purity, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hShowerCompletenessVSEnergy, showerE, Completeness, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hShowerPurityVSnHitsColl, nhits_coll, Purity, weight);
    if(kMC) plotUtils.FillHist(AnaIO::hShowerCompletenessVSnHitsColl, nhits_coll, Completeness, weight);
  }
  // Cut on number of hits
  if(nhits <= 80){
    if(kFill) AnaIO::hCutDaughterShowernHitsPass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterShowernHitsPass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hShowerEnergy_CutnHits, showerE, truthParticleType, weight);

  // Cut on shower stratZ (in cm)
  if(kFill) plotUtils.FillHist(AnaIO::hCutShowerStartZ, startZ, truthParticleType, weight);
  if(startZ <= -999){
    if(kFill) AnaIO::hCutDaughterShowerStartZPass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterShowerStartZPass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hShowerEnergy_CutStartZ, showerE, truthParticleType, weight);

  return true; 
}

bool AnaCut::IsMichel(const int ii, const bool kMC, const double weight, const bool kFill)
{
  // Get Michel score of this particle
  const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];
  if(kFill) plotUtils.FillHist(AnaIO::hCutmichelScore, michelScore, truthParticleType, weight);
  // Cut on Michel score
  if(michelScore<=0.5) return false;

  return true;
}

bool AnaCut::IsPiZeroShower(const int ii, const bool kMC, const double weight, const bool kFill)
{
  const double showerE = (*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3;
  if(showerE == -999*1E-3){
    if(kFill) AnaIO::hCutDaughterShowerNonEmptyEPass->Fill(false);
    return false;
  }
  else{
    if(kFill) AnaIO::hCutDaughterShowerNonEmptyEPass->Fill(true);
  }

  const double Purity = (*AnaIO::reco_daughter_PFP_true_byHits_purity)[ii];
  const double Completeness = (*AnaIO::reco_daughter_PFP_true_byHits_completeness)[ii];
  
  // Get the position vector point from vertex to shower start position
  const TVector3 dist = anaUtils.GetRecShowerDistVector(ii);
  // Get reco and truth shower momentum vector in both lab frame
  const TLorentzVector recShowerMom = anaUtils.GetRecShowerLTVectLab(ii);//GetRecShowerRefBeam(false,ii);//GetRecShowerLTVectLab(ii);
  const TLorentzVector recShowerMomRaw = anaUtils.GetRecShowerLTVectLab(ii,false);//GetRecShowerRefBeam(false,ii,false);//GetRecShowerLTVectLab(ii,false);
  const TLorentzVector truthShowerMom = anaUtils.GetTruthMatchedShowerLTVectLab(ii);//GetRecShowerRefBeam(true,ii);//GetTruthMatchedShowerLTVectLab(ii);
  // Calculate shower impact parameter
  const double IP = dist.Mag()*TMath::Sin((recShowerMom.Angle(dist)));
  // Get the shower length
  //const double showerLength = (*AnaIO::reco_daughter_allShower_len)[ii];  
  // Get the shower position vector
  const TVector3 showerPosition((*AnaIO::reco_daughter_allShower_startX)[ii],(*AnaIO::reco_daughter_allShower_startY)[ii],(*AnaIO::reco_daughter_allShower_startZ)[ii]);
  // Fill Cut histogram
  if(kFill) plotUtils.FillHist(AnaIO::hCutShowerDist, dist.Mag(), truthParticleType, weight);
  if(kFill) plotUtils.FillHist(AnaIO::hCutShowerIP, IP, truthParticleType, weight);
  //plotUtils.FillHist(AnaIO::hRecShowerLength, showerLength, truthParticleType);
  
  if(kMC && kFill) plotUtils.FillHist(AnaIO::hShowerPurityVSIP, IP, Purity, weight);
  if(kMC && kFill) plotUtils.FillHist(AnaIO::hShowerCompletenessVSIP, IP, Completeness, weight);

  // In unit of cm
  if( dist.Mag() < 3 || dist.Mag() > 90 ) {
    if(kFill) AnaIO::hCutDaughterShowerDistPass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterShowerDistPass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hShowerEnergy_CutDistance, showerE, truthParticleType, weight);

  // Impact Parameter Cut
  if( IP > 20 ){
    if(kFill) AnaIO::hCutDaughterShowerIPPass->Fill(false);
    return false;
  }
  else {
    if(kFill) AnaIO::hCutDaughterShowerIPPass->Fill(true);
  }
  if(kFill) plotUtils.FillHist(AnaIO::hShowerEnergy_CutIP, showerE, truthParticleType, weight);
  
  // Calculate the angle between dist vector and Pandora shower direction
  const double angleOffset = recShowerMom.Angle(dist)*TMath::RadToDeg();
  if(kFill) plotUtils.FillHist(AnaIO::hShowerAngleOffset, angleOffset, truthParticleType, weight);


  // Get Michel score of this particle
  const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];
  if(kFill) plotUtils.FillHist(AnaIO::hShowerMichelScore, michelScore, truthParticleType, weight);

  //if(michelScore > 0.5) {cout << "michel cut false" << endl; return false;}

  // Need to save all pizero shower candidates to reconstruct pizero
  anaUtils.SavePiZeroShower(recShowerMom, recShowerMomRaw, truthShowerMom, recShowerMom.E(), truthShowerMom.E(), showerPosition, truthParticleType);
  return true;
}


bool AnaCut::IsPizero(const bool kMC, const bool kFill, const double weight)
{
  // Pi zero selection
  OA = -999;
  truthPi0Type = -999;
  
  // Turn this on when need to KF sideband tune
  PiZeroVec = anaUtils.GetRecPiZeroFromShowers(OA,kMC,kFill,true,truthPi0Type);

  // Get pi0 vector no need to do KF in this step
  PiZeroVec = anaUtils.GetRecPiZeroFromShowers(OA,kMC,kFill,false,truthPi0Type);

  // Output tree for kinematic fitting
  if(kFill && kMC) anaUtils.SavePi0ShowersForKF();

  //if(kMC) cout << "Cut Purity num: " << anaUtils.selected << endl;
  //if(kMC) cout << "Cut Purity demo: " << anaUtils.total << endl;

  double mass = PiZeroVec.M();

  if(npi0shower < 2) {
    AnaIO::hCutDaughterPi0NOCutsPass->Fill(false);
    return false;
  }
  else AnaIO::hCutDaughterPi0NOCutsPass->Fill(true);
  if(kFill) plotUtils.FillHist(AnaIO::hPi0Energy_NoCut,PiZeroVec.E(),truthPi0Type, weight);
  if(kFill && truthPi0Type == anaUtils.gkTwoGammasSamePi0) AnaIO::hTruthMatchedPi0Sample->Fill(PiZeroVec.E(), weight);

  // Invariant mass cut
  if(mass < 0.05 || mass > 0.25){
    AnaIO::hCutDaughterPi0MassPass->Fill(false);
    return false;
  }
  else AnaIO::hCutDaughterPi0MassPass->Fill(true);
  if(kFill) plotUtils.FillHist(AnaIO::hPi0Energy_CutMass,PiZeroVec.E(),truthPi0Type, weight);

  // Opening angle cut
  if(OA < 10 || OA > 80){
    AnaIO::hCutDaughterPi0OAPass->Fill(false);
    return false;
  }
  else AnaIO::hCutDaughterPi0OAPass->Fill(true);
  if(kFill) plotUtils.FillHist(AnaIO::hPi0Energy_CutOA,PiZeroVec.E(),truthPi0Type, weight);

  return true;
}




