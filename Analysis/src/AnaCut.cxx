#include "../include/AnaCut.h"

bool AnaCut::CutBeamAllInOne(const bool kMC)
{
  // Standard procedure from pion analyses (updates on 16 Nov.)

  // 1. Beam PDG cut (pion beam/other)
  const bool passBeamPDG = CutBeamPDG(kMC);
  AnaIO::hCutBeamPDGPass->Fill(passBeamPDG);
  if(!passBeamPDG) return false;

  // Fill beam quality cut histograms (use fitted value for beam quality cut latter)
  FillBeamQualityHist();

  // 2. Pandora slice cut (track like/shower like)
  const bool passPandoraSlice = CutPandoraSlice();
  AnaIO::hCutPandoraSlicePass->Fill(passPandoraSlice);
  if(!passPandoraSlice) return false;

  // 3. Calo size cut 
  const bool passCaloSize = CutCaloSize();
  AnaIO::hCutCaloSizePass->Fill(passCaloSize);
  if(!passCaloSize) return false;

  // 4. Beam quality cut
  const bool passBeamQuality = CutBeamQuality(kMC);
  AnaIO::hCutBeamQualityPass->Fill(passBeamQuality);
  if(!passBeamQuality) return false;

  // 5. APA3 cut
  const bool passAPA3EndZ = CutAPA3EndZ();
  AnaIO::hCutAPA3EndZPass->Fill(passAPA3EndZ);
  if(!passAPA3EndZ) return false;

  // 6. Michel Score cut (remove muons)
  const bool passMichelScore = CutMichelScore(kMC);
  AnaIO::hCutMichelScorePass->Fill(passMichelScore);
  if(!passMichelScore) return false;

  // 7. Median dEdx cut (remove protons)
  const bool passMediandEdx = CutMediandEdx(kMC);
  AnaIO::hCutMediandEdxPass->Fill(passMediandEdx);
  if(!passMediandEdx) return false;
  
  return true;
}

void AnaCut::FillBeamQualityHist()
{
  // Beam direction vector
  TVector3 BeamDir(AnaIO::reco_beam_calo_endX-AnaIO::reco_beam_calo_startX,
                  AnaIO::reco_beam_calo_endY-AnaIO::reco_beam_calo_startY,AnaIO::reco_beam_calo_endZ-AnaIO::reco_beam_calo_startZ);
  // Detector axis
  TVector3 DetX(1,0,0);
  TVector3 DetY(0,1,0);
  TVector3 DetZ(0,0,1);
  // Fill histograms
  AnaIO::hRecBeamStartX->Fill(AnaIO::reco_beam_calo_startX);
  AnaIO::hRecBeamStartY->Fill(AnaIO::reco_beam_calo_startY);
  AnaIO::hRecBeamStartZ->Fill(AnaIO::reco_beam_calo_startZ);

  AnaIO::hRecBeamThetaX->Fill(BeamDir.Angle(DetX)*TMath::RadToDeg());
  AnaIO::hRecBeamThetaY->Fill(BeamDir.Angle(DetY)*TMath::RadToDeg());
  AnaIO::hRecBeamThetaZ->Fill(BeamDir.Angle(DetZ)*TMath::RadToDeg());
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

bool AnaCut::CutBeamQuality(bool kMC, bool DoAngleCut)
{
  // More info: https://docs.dunescience.org/cgi-bin/private/RetrieveFile?docid=23122&filename=ProtoDUNE_CNN_Hit_Tagging_Technote.pdf&version=1

  const double cut_beamQuality_TPC_xy = 3.;
  const double cut_beamQuality_TPC_z = 3.;
  const double cut_beamQuality_TPC_cosTheta = 0.95;

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
  
  AnaIO::hCutBeamDeltaXYSigma->Fill(cut_xy);
  AnaIO::hCutBeamDeltaZSigma->Fill(cut_z);
  AnaIO::hCutBeamCosTheta->Fill(cut_theta);

  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  plotUtils.FillHist(AnaIO::hCutBeamQualityXY, cut_xy , parType);
  plotUtils.FillHist(AnaIO::hCutBeamQualityX, cut_x , parType);
  plotUtils.FillHist(AnaIO::hCutBeamQualityY, cut_y , parType);
  plotUtils.FillHist(AnaIO::hCutBeamQualityZ, cut_z , parType);
  plotUtils.FillHist(AnaIO::hCutBeamQualityTheta, cut_theta , parType);
  
  


  if(cut_xy > cut_beamQuality_TPC_xy) return false;

  if(cut_z > cut_beamQuality_TPC_z) return false;

  if(DoAngleCut){
    if(cut_theta < cut_beamQuality_TPC_cosTheta) return false;
  }

  return true;
}

bool AnaCut::CutAPA3EndZ(){

  //AnaIO::hCutBeamEndZ->Fill(AnaIO::reco_beam_calo_endZ);
  if(AnaIO::reco_beam_calo_endZ >= cut_EndZ_APA3){
    return false;
  }

  return true;
}

bool AnaCut::CutMichelScore(const bool kMC)
{
  if( AnaIO::reco_beam_vertex_michel_score == 0) return true;
  const double MichelScorePerHits = AnaIO::reco_beam_vertex_michel_score/AnaIO::reco_beam_vertex_nHits;
  
  const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
  plotUtils.FillHist(AnaIO::hCutBeamMichelScoreNhits, MichelScorePerHits , parType);
  // There is a muon candidate 
  if(MichelScorePerHits > cut_MichelScore) return false;

  return true;
}

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

int AnaCut::GetTruthPDGFromID(const int inID, const vector<int> * idarray, const vector<int> * pdgarray)
{ 
  int outpdg = -999;
  for(unsigned int ii = 0; ii<idarray->size(); ii++){
    if((*idarray)[ii] == inID ){
      outpdg = (*pdgarray)[ii];
    }
  }
  return outpdg;
}

// Aim to find secondary daughter particles
int AnaCut::GetTruthParticleInfoFromRec(const int recidx)
{
  // Get true PDG of this reco particle at index recidx
  const int directPDG = (*AnaIO::reco_daughter_PFP_true_byHits_PDG)[recidx];
  // Get true ID array of this event
  const vector<int> *trueIDarray = AnaIO::reco_daughter_PFP_true_byHits_ID;
  // Get true ID of this reco particle
  const int truthID = (*trueIDarray)[recidx];
  // Initialize variables
  bool isPrimary = false;  
  int pdg = -999;
  
  //-------------- First search direct daughter ----------------// 
  //(Loop over true beam daughther ID and PDG by function GetTruthPDGFromID) 
  pdg = GetTruthPDGFromID(truthID, AnaIO::true_beam_daughter_ID, AnaIO::true_beam_daughter_PDG);

  if(pdg!=-999){ //1. is direct daughter
    // Not a proton, pion, electron, muon, photon or kaon
    if(pdg!=2212 && TMath::Abs(pdg)!=211 && TMath::Abs(pdg)!=11 && TMath::Abs(pdg)!=13 && TMath::Abs(pdg)!=22 && TMath::Abs(pdg)!=321 && pdg<1000000000){
      printf("GetTruthParticleInfoFromRec reconstructed truth daughter not proton or pion! %d %d\n", pdg, directPDG); exit(1);
    } 
    isPrimary = true;

  }// End of is direct daughter

  else{ // 1. not direct daughter
    //-------------- Then search pi0 daughter -----------------// 
    pdg = GetTruthPDGFromID(truthID, AnaIO::true_beam_Pi0_decay_ID, AnaIO::true_beam_Pi0_decay_PDG);
    
    if(pdg!=-999){//2. is pi0 daughter
      // Not a photon or electron
      if(pdg!=22 && TMath::Abs(pdg)!=11){
      printf("GetTruthParticleInfoFromRec Pi0 decay not to gamma! %d %d\n", pdg, directPDG); exit(1);
      }
      //pi0 direct daughter also primary
      isPrimary = true;
    }
    else{// 2. not pi0 daughter

      //----------- Then search grand daugher (This will be secondary particle candidates)
      pdg = GetTruthPDGFromID(truthID, AnaIO::true_beam_grand_daughter_ID, AnaIO::true_beam_grand_daughter_PDG);
      if(pdg!=-999){
      
      }
      else{
        //--- lump great grand daughter here
        if(TMath::Abs(directPDG)==11 || TMath::Abs(directPDG)==13 || TMath::Abs(directPDG)==22 || TMath::Abs(directPDG)==211 || TMath::Abs(directPDG)==321 || directPDG==2212 || directPDG==2112 || directPDG==-1 || directPDG>1000000000){//when no true match found pdg = -1
          pdg = directPDG;
        }
        else{
          printf("AnaCut::GetTruthFromRec search not done! %d %d\n", recidx, directPDG); exit(1);
        }
      }
    } // End of not pi0 daughter
  } // End of not direct daughter
  // The directPDG should be equal to the pdg we found
  if(directPDG!=pdg){
    printf("GetTruthFromRec inconsistent PDG %d %d\n", pdg, directPDG); exit(1);
  }
  // Default truthParticleType is others
  int truthParticleType = anaUtils.gkOthers;
  if(isPrimary){
    if(pdg==2212){//proton
      truthParticleType = anaUtils.gkProton;
    }
    else if(pdg==211){//pi+
      truthParticleType = anaUtils.gkPiPlus;
    }
    else if(pdg==-211){//pi-
      truthParticleType = anaUtils.gkPiMinus;
    }
    else if(pdg==22){//gamma
      truthParticleType = anaUtils.gkGamma;
    }
  }
  else{
    if(pdg==2212){//proton
      truthParticleType = anaUtils.gkSecondaryProton;
    }
    else if(pdg==211){//pi+
      truthParticleType = anaUtils.gkSecondaryPiPlus;
    }
    else if(pdg==-211){//pi-
      truthParticleType = anaUtils.gkSecondaryPiMinus;
    }
    else if(pdg==22){//gamma
      truthParticleType = anaUtils.gkSecondaryGamma;
    }
    else if(TMath::Abs(pdg)==11){//e+/-
      truthParticleType = anaUtils.gkSecondaryEplusEminus;
    }
    else if(TMath::Abs(pdg)==13){//mu+/-
      truthParticleType = anaUtils.gkSecondaryMuon;
    }
  }

  return truthParticleType;
}



bool AnaCut::CutTopology(const bool kMC)
{
  // Count reco final state particles and determine types 
  CountPFP(kMC,true);
  // ----------------------- Do cuts below ------------------------ // 
  // Get event type
  const int evtType = anaUtils.GetFillEventType();
  // Proton
  plotUtils.FillHist(AnaIO::hCutnproton, nproton, evtType);  
  if(nproton<1) return false;
  // Showers
  plotUtils.FillHist(AnaIO::hCutnshower, npi0shower, evtType);  
  if(npi0shower<2) return false;
  // Piplus
  plotUtils.FillHist(AnaIO::hCutnpiplus, npiplus, evtType);  
  if(npiplus!=0) return false;
  // Michel electron
  plotUtils.FillHist(AnaIO::hCutnmichel, nmichel, evtType);  
  if(nmichel!=0) return false;

  return true;
}

void AnaCut::CountPFP(const bool kMC, const bool kFill)
{
  // Get the size of reco final state particles
  const int recsize = AnaIO::reco_daughter_PFP_ID->size();
  // Initialize number of reco final state particles (beam reco daughter particles counter)
  nproton = 0;
  npiplus = 0;
  nshower = 0;
  npi0shower = 0;
  nmichel = 0;
  int nPFP = 0;
  // Need to clear vector for each event
  anaUtils.CleanShowerArray();

  double ProtonMom = -999;
  // Loop over all reco FS particles
  for(int ii=0; ii<recsize; ii++){
    // Get the truth particle type of this reco particle
    truthParticleType = kMC? GetTruthParticleInfoFromRec(ii) : anaUtils.gkOthers;
    recParticleType = -999;
    // Proton candidates selection
    if(IsProton(ii,kMC)){
      nproton++;
      recParticleType = anaUtils.gkProton;
      // Fill rec and truth-matching info
      if(kFill) anaUtils.FillFSParticleKinematics(ii, truthParticleType, recParticleType);

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
    }
    // Piplus candidates selection
    if(IsPiplus(ii,kMC)){
      npiplus++;
      recParticleType = anaUtils.gkPiPlus;
      // Fill rec and truth-matching info
      if(kFill) anaUtils.FillFSParticleKinematics(ii, truthParticleType, recParticleType);
    }
    // Shower candidates selection
    if(IsShower(ii,kMC)){
      nshower++;
      // Select good shower candidates for reconstructing pi0   
      if(IsPiZeroShower(ii,kMC)){
        npi0shower++;
        recParticleType = anaUtils.gkGamma;
        // Fill rec and truth-matching info
        if(kFill) anaUtils.FillFSParticleKinematics(ii, truthParticleType, recParticleType);
      }
    }
    // Michel candidates selection
    if(IsMichel(ii,kMC)){
      nmichel++;
    }
    // Count total number of FS particles
    nPFP++; 
  }
  if(recsize!=nPFP) cout << "CountPFP not looping all FS particles!!" << endl;
 
  if(kFill) anaUtils.GetRecPiZeroFromShowers();
  // Output tree for kinematic fitting
  if(kFill && kMC) anaUtils.SavePi0ShowersForKF();
  
  //if(npi0shower > 1 && nproton > 0) printf("CountPFP PFP size %d nlooped %d nshower %d npi0shower %d nmichel %d npiplus %d nproton %d\n", recsize, nPFP, nshower, npi0shower, nmichel, npiplus, nproton);
}

bool AnaCut::IsProton(const int ii, const bool kMC)
{
  // Proton candidate must be a track like particle
  if(!IsTrack(ii,kMC)) return false; 
  // Proton candidate must pass the folowing cuts 
  if(!PassProtonSubPID(ii)) return false;

  return true;
}

bool AnaCut::IsTrack(const int ii,const bool kMC)
{
  // Get the number of hits of this reco particle
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  // Get the track score
  const double trackScore = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];
  // Get the shower purity and completeness
  const double Purity = (*AnaIO::reco_daughter_PFP_true_byHits_purity)[ii];
  const double Completeness = (*AnaIO::reco_daughter_PFP_true_byHits_completeness)[ii];

  // Fill Cut histograms for those two variables
  plotUtils.FillHist(AnaIO::hCutTracknHits, nhits, truthParticleType);  
  plotUtils.FillHist(AnaIO::hCutTrackScore, trackScore, truthParticleType);
  if(kMC) plotUtils.FillHist(AnaIO::hTrackCompleteness, Completeness, truthParticleType);
  // Check if this reco particle has forced track ID
  if((*AnaIO::reco_daughter_allTrack_ID)[ii]==-1) return false;
  // Cut on track score
  if(trackScore <= 0.5) return false;

  if(kMC) plotUtils.FillHist(AnaIO::hTrackPurityVSnHits, nhits, Purity);
  if(kMC) plotUtils.FillHist(AnaIO::hTrackCompletenessVSnHits, nhits, Completeness);

  // Cut on number of hits
  if(nhits <= 40) return false;

  return true;
}

bool AnaCut::PassProtonSubPID(const int ii)
{
  // Get the Truncated mean dEdx backwards of this reco particle
  const double lastTME = anaUtils.GetFSParticleTME(ii,false);
  // Get the Chi2/NDF value
  const double Chi2NDF = anaUtils.GetChi2NDF(ii);
  // Fill Cut histograms for those two variables
  plotUtils.FillHist(AnaIO::hCutlastTME, lastTME, truthParticleType);
  plotUtils.FillHist(AnaIO::hCutChi2NDF, Chi2NDF, truthParticleType);
  // Use these two variables for proton selections
  if(Chi2NDF < 50 || lastTME > 3.5){
    return true;
  }
  return false;
}

bool AnaCut::IsPiplus(const int ii, const bool kMC)
{
  // PiPlus candidate must be a track like particle  
  if(!IsTrack(ii,kMC)) return false;
  // Assume particle not pass ProtonSubPID is piplus 
  if(PassProtonSubPID(ii)) return false;
  
  return true;  
}

bool AnaCut::IsShower(const int ii, const bool kMC)
{
  // Get the em score and nhits of this particle
  const double emScore = (*AnaIO::reco_daughter_PFP_emScore_collection)[ii];
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  // Get the shower purity and completeness
  const double Purity = (*AnaIO::reco_daughter_PFP_true_byHits_purity)[ii];
  const double Completeness = (*AnaIO::reco_daughter_PFP_true_byHits_completeness)[ii];

  // Fill Cut histogram
  plotUtils.FillHist(AnaIO::hCutemScore, emScore, truthParticleType);
  if(kMC) plotUtils.FillHist(AnaIO::hShowerCompleteness, Completeness, truthParticleType);

  if((*AnaIO::reco_daughter_allShower_ID)[ii]==-1) return false;
  // Cut on em score
  if(emScore <= 0.5) return false;

  if(kMC) plotUtils.FillHist(AnaIO::hShowerPurityVSnHits, nhits, Purity);
  if(kMC) plotUtils.FillHist(AnaIO::hShowerCompletenessVSnHits, nhits, Completeness);
  
  // Cut on number of hits
  if(nhits <= 80) return false;

  return true; 
}

bool AnaCut::IsMichel(const int ii, const bool kMC)
{
  // Get Michel score of this particle
  const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];
  plotUtils.FillHist(AnaIO::hCutmichelScore, michelScore, truthParticleType);
  // Cut on Michel score
  if(michelScore<=0.5) return false;

  return true;
}

bool AnaCut::IsPiZeroShower(const int ii, const bool kMC)
{
  if((*AnaIO::reco_daughter_allShower_energy)[ii] == -999) return false;

  // Get the position vector point from vertex to shower start position
  const TVector3 dist = anaUtils.GetRecShowerDistVector(ii);

  // Get reco and truth shower momentum vector in both lab frame
  const TLorentzVector recShowerMom = anaUtils.GetRecShowerLTVectLab(ii);//GetRecShowerRefBeam(false,ii);//GetRecShowerLTVectLab(ii);
  const TLorentzVector recShowerMomRaw = anaUtils.GetRecShowerLTVectLab(ii,false);//GetRecShowerRefBeam(false,ii,false);//GetRecShowerLTVectLab(ii,false);
  const TLorentzVector truthShowerMom = anaUtils.GetTruthMatchedShowerLTVectLab(ii);//GetRecShowerRefBeam(true,ii);//GetTruthMatchedShowerLTVectLab(ii);
  // Calculate shower impact parameter
  const double IP = dist.Mag()*TMath::Sin((recShowerMom.Angle(dist)));
  // Get the shower length
  const double showerLength = (*AnaIO::reco_daughter_allShower_len)[ii];  
  // Get the shower position vector
  const TVector3 showerPosition((*AnaIO::reco_daughter_allShower_startX)[ii],(*AnaIO::reco_daughter_allShower_startY)[ii],(*AnaIO::reco_daughter_allShower_startZ)[ii]);
  // Fill Cut histogram
  plotUtils.FillHist(AnaIO::hCutShowerDist, dist.Mag(), truthParticleType);
  plotUtils.FillHist(AnaIO::hCutShowerIP, IP, truthParticleType);
  plotUtils.FillHist(AnaIO::hRecShowerLength, showerLength, truthParticleType);
  // In unit of cm
  if( dist.Mag() < 2 || dist.Mag() > 90 ) return false;
  // Impact Parameter Cut
  if( IP > 15 ) return false;
  // Need to save all pizero shower candidates to reconstruct pizero
  anaUtils.SavePiZeroShower(recShowerMom, recShowerMomRaw, truthShowerMom, recShowerMom.E(), truthShowerMom.E(), showerPosition, truthParticleType);
  return true;
}






