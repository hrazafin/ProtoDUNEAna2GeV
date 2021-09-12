#include "../include/AnaCut.h"

bool AnaCut::CutBeamAllInOne(const bool kMC)
{
  // Standard procedure from pion analyses
  
  // 1. Beam ID cut (pion beam/other)
  if(kMC){ // MC
    // In data the beam instrumentation is able to filter for these events but is not inside the MC
    const bool mc_beamID =(AnaIO::true_beam_PDG == 211 || AnaIO::true_beam_PDG == -13);
    AnaIO::hCutBeamIDPass->Fill(mc_beamID);
    if(!mc_beamID){
      return false;
    }
  }
  else{ // Data
    const bool data_beamID = CutBeamID((*AnaIO::beam_inst_PDG_candidates));
    AnaIO::hCutBeamIDPass->Fill(data_beamID);
    if(!data_beamID){
      return false;
    }
  }
  
  //2. Primary beam type cut (track like/shower like)
  const int beam_recoType = AnaIO::reco_beam_type;
  AnaIO::hCutBeamType->Fill(beam_recoType);
  if(beam_recoType!=13){//13: Pandora "track like"
    return false;
  }

  //3. Beam position cut
  const bool kBeamPosPass = kMC ? CutMCBeamPos() : CutDataBeamPos();
  AnaIO::hCutBeamPosPass->Fill(kBeamPosPass);
  if(!kBeamPosPass){
    return false;
  }

  //4. APA3 cut
  const double endzcut = 226;
  AnaIO::hCutBeamEndZ->Fill(AnaIO::reco_beam_endZ);
  AnaIO::hCutBeamEndZPass->Fill(!(AnaIO::reco_beam_endZ>=endzcut));
  if(AnaIO::reco_beam_endZ >= endzcut){
    return false;
  }
  
  return true;
}


bool AnaCut::CutBeamID(const std::vector<int> & pidCandidates)
{
  return ( (std::find(pidCandidates.begin(), pidCandidates.end(), 211)) != pidCandidates.end());
}

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
  // Loop over all reco FS particles
  for(int ii=0; ii<recsize; ii++){
    // Get the truth particle type of this reco particle
    truthParticleType = kMC? GetTruthParticleInfoFromRec(ii) : anaUtils.gkOthers;
    recParticleType = -999;
    // Proton candidates selection
    if(IsProton(ii)){
      nproton++;
      recParticleType = anaUtils.gkProton;
      // Fill rec and truth-matching info
      if(kFill) anaUtils.FillFSParticleKinematics(ii, truthParticleType, recParticleType);
    }
    // Piplus candidates selection
    if(IsPiplus(ii)){
      npiplus++;
      recParticleType = anaUtils.gkPiPlus;
      // Fill rec and truth-matching info
      if(kFill) anaUtils.FillFSParticleKinematics(ii, truthParticleType, recParticleType);
    }
    // Shower candidates selection
    if(IsShower(ii)){
      nshower++;
      // Select good shower candidates for reconstructing pi0   
      if(IsPiZeroShower(ii)){
        npi0shower++;
        recParticleType = anaUtils.gkGamma;
        // Fill rec and truth-matching info
        if(kFill) anaUtils.FillFSParticleKinematics(ii, truthParticleType, recParticleType);
      }
    }
    // Michel candidates selection
    if(IsMichel(ii)){
      nmichel++;
    }
    // Count total number of FS particles
    nPFP++; 
  }
  if(recsize!=nPFP) cout << "CountPFP not looping all FS particles!!" << endl;
  if(kFill && npi0shower >= 2) anaUtils.GetPiZero();
  // Output tree for kinematic fitting
  if(kFill && npi0shower >= 2 && kMC) anaUtils.GetPi0Showers();
  
  //if(npi0shower > 1 && nproton > 0) printf("CountPFP PFP size %d nlooped %d nshower %d npi0shower %d nmichel %d npiplus %d nproton %d\n", recsize, nPFP, nshower, npi0shower, nmichel, npiplus, nproton);
}

bool AnaCut::IsProton(const int ii)
{
  // Proton candidate must be a track like particle
  if(!IsTrack(ii)) return false; 
  // Proton candidate must pass the folowing cuts 
  if(!PassProtonSubPID(ii)) return false;

  return true;
}

bool AnaCut::IsTrack(const int ii)
{
  // Get the number of hits of this reco particle
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  // Get the track score
  const double trackScore = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];
  // Fill Cut histograms for those two variables
  plotUtils.FillHist(AnaIO::hCutTracknHits, nhits, truthParticleType);  
  plotUtils.FillHist(AnaIO::hCutTrackScore, trackScore, truthParticleType);

  // Check if this reco particle has forced track ID
  if((*AnaIO::reco_daughter_allTrack_ID)[ii]==-1) return false;
  // Cut on number of hits
  if(nhits <= 40) return false;
  // Cut on track score
  if(trackScore <= 0.5) return false;
  
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

bool AnaCut::IsPiplus(const int ii)
{
  // PiPlus candidate must be a track like particle  
  if(!IsTrack(ii)) return false;
  // Assume particle not pass ProtonSubPID is piplus 
  if(PassProtonSubPID(ii)) return false;
  
  return true;  
}

bool AnaCut::IsShower(const int ii)
{
  // Get the em score and nhits of this particle
  const double emScore = (*AnaIO::reco_daughter_PFP_emScore_collection)[ii];
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  // Fill Cut histogram
  plotUtils.FillHist(AnaIO::hCutemScore, emScore, truthParticleType);

  if((*AnaIO::reco_daughter_allShower_ID)[ii]==-1) return false;
  // Cut on number of hits
  if(nhits <= 80) return false;
  // Cut on em score
  if(emScore <= 0.5) return false;
  
  return true; 
}

bool AnaCut::IsMichel(const int ii)
{
  // Get Michel score of this particle
  const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];
  plotUtils.FillHist(AnaIO::hCutmichelScore, michelScore, truthParticleType);
  // Cut on Michel score
  if(michelScore<=0.5) return false;

  return true;
}

bool AnaCut::IsPiZeroShower(const int ii)
{
  // Get the position vector point from vertex to shower start position
  const TVector3 dist = anaUtils.GetRecShowerDistVector(ii);
  // Get reco and truth shower momentum vector in both lab frame
  const TLorentzVector recShowerMom = anaUtils.GetRecShowerLTVectLab(ii);
  const TLorentzVector recShowerMomRaw = anaUtils.GetRecShowerLTVectLab(ii,false);
  const TLorentzVector truthShowerMom = anaUtils.GetTruthMatchedShowerLTVectLab(ii);
  // Calculate shower impact parameter
  const double IP = dist.Mag()*TMath::Sin((recShowerMom.Angle(dist)));
  // Get the shower length
  const double showerLength = (*AnaIO::reco_daughter_allShower_len)[ii];  
  // Get the shower position vector
  const TVector3 showerPosition((*AnaIO::reco_daughter_allShower_startX)[ii],(*AnaIO::reco_daughter_allShower_startY)[ii],(*AnaIO::reco_daughter_allShower_startZ)[ii]);
 
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







