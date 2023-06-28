#include "../include/AnaCut.h"

bool AnaCut::CutBeamAllInOne(const bool kMC, bool kStopMuon, bool kFill)
{ 
  // Standard procedures of protoDUNE-SP pion beam selections
  double weight = 1.0;//anaUtils.CalWeight(kMC);

  // 1. Beam PDG cut (pion beam/other)
  if(!kStopMuon){
    const bool passBeamPDG = CutBeamPDG(kMC);
    if(!passBeamPDG) return false;
  }
  else{
    const bool passBeamPDG = CutBeamPDG_StopMuon(kMC);
    if(!passBeamPDG) return false;
  }
  // 2. Pandora slice cut (track like/shower like)
  const bool passPandoraSlice = CutPandoraSlice();
  if(!passPandoraSlice) return false;

  // 3. Calo size cut 
  const bool passCaloSize = CutCaloSize();
  if(!passCaloSize) return false;
  

  // 4. Beam quality cut
  const bool passBeamQuality = CutBeamQuality(kMC,true,kFill);
  if(!passBeamQuality) return false;

/*
  // 5. APA3 cut
  if(!kStopMuon){
    const bool passAPA3EndZ = CutAPA3EndZ(kMC, kFill, weight);
    if(!passAPA3EndZ) return false;
  }
  else{
    const bool passAPA3EndZ = CutAPA3EndZ(kMC, kFill, weight);
    if(!passAPA3EndZ) return false;
  }
*/
  // 6. Michel Score cut (remove muons)
  if(!kStopMuon){
    const bool passMichelScore = CutMichelScore(kMC, kFill, weight);
    if(!passMichelScore) return false;
  }
  else{ 
    const bool passMichelScore = CutMichelScore_StopMuon(kMC, kFill, weight);
    if(!passMichelScore) return false;
  }

  // 7. Median dEdx cut (remove protons)
  const bool passProtonChi2 = CutProtonChi2DOF(kMC,kFill, weight);
  if(!passProtonChi2) return false;


  // 8. CutBeamScraper
  const bool passBeamScraper = CutBeamScraper(kMC);
  if(!passBeamScraper) return false;


  //========= Event Passed All Pion Beam Cuts =======//
  
  return true;
}


bool AnaCut::CutBeamLongTrack(const bool kMC, bool kFill)
{ 
  // Standard procedures of protoDUNE-SP pion beam selections (for long track)
  //double weight = 1.0;//anaUtils.CalWeight(kMC);

  // 1. Beam PDG cut (pion beam/other)
  const bool passBeamPDG = CutBeamPDG(kMC);
  if(!passBeamPDG) return false;
  // 2. Pandora slice cut (track like/shower like)
  const bool passPandoraSlice = CutPandoraSlice();
  if(!passPandoraSlice) return false;
  // 3. Calo size cut 
  const bool passCaloSize = CutCaloSize();
  if(!passCaloSize) return false;
  // 4. Beam quality cut
  const bool passBeamQuality = CutBeamQuality(kMC,true,kFill);
  if(!passBeamQuality) return false;

  //========= Event Passed All Pion Beam Cuts =======//
  
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

bool AnaCut::CutBeamPDG_StopMuon(const bool kMC)
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

  //const double cut_beamQuality_TPC_xy_low = -1.;
  const double cut_beamQuality_TPC_xy = 3.;
  //const double cut_beamQuality_TPC_y = 3.;
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

  const double cut_z = kMC ? deltaZ/MCsigmaStartZ : deltaZ/DATAsigmaStartZ;
  


  const double cut_theta = kMC ? cos(thetaX)*cos(MCmeanThetaX)+cos(thetaY)*cos(MCmeanThetaY)+cos(thetaZ)*cos(MCmeanThetaZ)
                               : cos(thetaX)*cos(DATAmeanThetaX)+cos(thetaY)*cos(DATAmeanThetaY)+cos(thetaZ)*cos(DATAmeanThetaZ);
  

  if(cut_xy > cut_beamQuality_TPC_xy) return false;

  //if(abs(cut_y) > cut_beamQuality_TPC_y) return false;
  if(abs(cut_z) > cut_beamQuality_TPC_z) return false;
  //if(cut_z > cut_beamQuality_TPC_z) return false;

  if(DoAngleCut){
    if(cut_theta < cut_beamQuality_TPC_cosTheta) return false;
  }

  //if(CutBeamInstQuality(kMC) == false) return false;

  if(AnaIO::beam_inst_P > 1.2 || AnaIO::beam_inst_P < 0.8) return false;

  return true;
}

bool AnaCut::CutBeamInstQuality(bool kMC)
{
  const double cut_beamQuality_TPC_xy = 2.5;

  const double deltaX = kMC ? AnaIO::beam_inst_X - MCmeanInstX : AnaIO::beam_inst_X - DATAmeanInstX;
  const double deltaY = kMC ? AnaIO::beam_inst_Y - MCmeanInstY : AnaIO::beam_inst_Y - DATAmeanInstY;
  
  const double cut_xy = kMC ? sqrt(pow(deltaX/MCsigmaInstX,2)+pow(deltaY/MCsigmaInstY,2)) : sqrt(pow(deltaX/DATAsigmaInstX,2)+pow(deltaY/DATAsigmaInstY,2));

  if(cut_xy > cut_beamQuality_TPC_xy) return false;

  return true;
}


bool AnaCut::CutAPA3EndZ(const bool kMC, bool kFill, const double weight){

  const double endZ = AnaIO::reco_beam_calo_endZ;
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
  // There is a muon candidate 
  if(MichelScoreWBC > cut_MichelScore) return false;

  return true;
}

bool AnaCut::CutMichelScore_StopMuon(const bool kMC, bool kFill, const double weight)
{
  if( AnaIO::reco_beam_vertex_michel_score == 0) return true;
  const double MichelScoreWBC = AnaIO::reco_beam_vertex_michel_score_weight_by_charge;
  // There is a muon candidate 
  //if(MichelScoreWBC > cut_MichelScore) return true;
  if(MichelScoreWBC > 0.60) return true;

  return false;
}

bool AnaCut::CutProtonChi2DOF(const bool kMC, bool kFill, const double weight)
{
  const double Chi2DOF = AnaIO::reco_beam_Chi2_proton/AnaIO::reco_beam_Chi2_ndof;
  if(Chi2DOF < 80) return false;

  return true;
}


bool AnaCut::CutBeamScraper(const bool kMC){

  // Define beam plug circle in [cm] unit
  double center_x = 0.;
  double center_y = 0.;
  double radious = 0.;
  double N_sigma = 0.;
  if(kMC){
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
  }

  double this_distance = sqrt( pow(AnaIO::beam_inst_X - center_x , 2.) + pow(AnaIO::beam_inst_Y - center_y, 2.) );
  if(this_distance > radious * N_sigma) return false;

  return true;
}

bool AnaCut::CutTopology(const bool kMC, bool kFill)
{ 
  pi0KineticE = -999;
  pi0costheta = -999;
  // Count reco final state particles and determine types 
  CountPFP(kMC,kFill);
 
  // ----------------------- Do cuts below ------------------------ // 
  // Proton
  //if(nproton!=1) return false;
  // Showers  
  // Cex signal 
  if(npi0shower < 2) return false;
  // Pion prod. background (No pi0) 
  //if(npi0shower > 1) return false;
  // Pion prod. background (One pi0) 
  //if(npi0shower < 2) return false;
  
  // Piplus
  // Cex signal 
  if(npiplus!=0) return false;
  // Pion prod. background (No pi0) 
  //if(npiplus==0) return false;
  // Pion prod. background (One pi0) 
  //if(npiplus==0) return false;


  //cout << "npiplus: " << npiplus << endl;
  //cout << "nmichel: " << nmichel << endl;

  // Michel electron
  // Cex signal  
  if(nmichel!=0) return false;
  // Pion prod. background (No pi0) No michel restriction
  // Pion prod. background (One pi0) No michel restriction

  // PiZero cut
  // Cex signal    
  if(npi0 != 1) return false;
  // Pion prod. background (No pi0)    
  //if(npi0 != 0) return false;
  // Pion prod. background (One pi0)    
  //if(npi0 < 1) return false;


  // Rewrite the pi zero vector this time need to turn on KF to improve pi0 energy resolution
  PiZeroVec = anaUtils.GetRecPiZeroFromShowers(OA,kMC,false,true,truthPi0Type);
  // After KF

  double pi0KE = PiZeroVec.E() - AnaFunctions::PiZeroMass();
  double costheta = PiZeroVec.Pz()/PiZeroVec.P();
  pi0KineticE = pi0KE; 
  pi0costheta = costheta;
  //if(pi0costheta < 0. && pi0costheta > -0.2)  cout << "Cut pi0cos: " << pi0costheta << " PiZeroVec.P(): " << PiZeroVec.P() << " PiZeroVec.Pz(): " << PiZeroVec.Pz() << endl; 

  return true;
}

void AnaCut::CountPFP(const bool kMC, const bool kFill)
{ 
  // Calculate event weight
  double weight = 1.0;//anaUtils.CalWeight(kMC);

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

  // Loop over all reco FS particles
  for(int ii=0; ii<recsize; ii++){
    // Get the truth particle type of this reco particle
    truthParticleType = kMC ? anaUtils.GetTruthParticleInfoFromRec(ii) : anaUtils.gkOthers;

    // Proton candidates selection
    if(IsProton(ii,kMC,weight,kFill)){
      nproton++;
    }
    // Shower candidates selection
    if(IsShower(ii,kMC,weight,kFill)){
      nshower++;
      // Select good shower candidates for reconstructing pi0   
      if(IsPiZeroShower(ii,kMC,weight,kFill)){
        npi0shower++;
      }
    }

    // Piplus candidates selection
    if(IsPiplus(ii,kMC,weight,kFill)){
      npiplus++;
    }

    // Michel candidates selection
    if(IsMichel(ii,kMC,weight,kFill)){
      nmichel++;
    }
    // Count total number of FS particles
    nPFP++; 
  }
  if(recsize!=nPFP) cout << "CountPFP not looping all FS particles!!" << endl;

  // Pi zero candidates selection (!one per event! Get pi0 from at least two pi0 showers)
  if(IsPizero(kMC,kFill,weight)){
    npi0++;
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
  // Get the track score
  const double trackScore = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];

  // Check if this reco particle has forced track ID
  if((*AnaIO::reco_daughter_allTrack_ID)[ii]==-1) return false;
  // Cut on track score
  if(trackScore <= 0.5){
    return false;
  }

  // Cut on number of hits
  if(nhits <= 20){
    return false;
  }
  
  return true;
}


bool AnaCut::IsPionTrack(const int ii,const bool kMC, const double weight, const bool kFill)
{
  // Get the number of hits of this reco particle
  const int nhits = (*AnaIO::reco_daughter_PFP_nHits)[ii];
  
  // Get the track score
  const double trackScore = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];

  // Check if this reco particle has forced track ID
  if((*AnaIO::reco_daughter_allTrack_ID)[ii]==-1) return false;

  // Cut on track score
  if(trackScore <= 0.5){
    return false;
  }
  // Cut on number of hits
  if(nhits <= 20){
    return false;
  }

  return true;
}


bool AnaCut::PassProtonSubPID(const int ii, const double weight, const bool kFill)
{
  // Get the Truncated mean dEdx backwards of this reco particle
  const double lastTME = anaUtils.GetFSParticleTME(ii,false);
  // Get the Chi2/NDF value
  const double Chi2NDF = anaUtils.GetChi2NDF(ii);

  // Use these two variables for proton selections
  if(Chi2NDF < 70 || lastTME > 2.8){
    return true;
  }
  
  return false;
}

bool AnaCut::PassPionSubPID(const int ii, const double weight, const bool kFill)
{
  // Get the Truncated mean dEdx backwards of this reco particle
  const double lastTME = anaUtils.GetFSParticleTME(ii,false);
  // Get the Chi2/NDF value
  const double Chi2NDF = anaUtils.GetChi2NDF(ii);
  
  // Use these two variables for pion selections
  if((lastTME < 2.8 && lastTME > 0.5) || (lastTME < 3.4 && Chi2NDF >70)){
    return true;
  }
  
  return false;
}


bool AnaCut::IsPiplus(const int ii, const bool kMC, const double weight, const bool kFill)
{
  // PiPlus candidate must be a track like particle  
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
  const double startZ = (*AnaIO::reco_daughter_allShower_startZ)[ii];
  
  if((*AnaIO::reco_daughter_allShower_ID)[ii]==-1) return false;

  // Cut on em score
  if(emScore <= 0.5){
    return false;
  }
  // Cut on number of hits
  if(nhits <= 80){
    return false;
  }
  // Cut on shower stratZ (in cm)
  if(startZ <= -999){
    return false;
  }
  
  return true; 
}

bool AnaCut::IsMichel(const int ii, const bool kMC, const double weight, const bool kFill)
{
  // Get Michel score of this particle
  const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];
  // Cut on Michel score
  if(michelScore<=0.5) return false;
  return true;
}

bool AnaCut::IsPiZeroShower(const int ii, const bool kMC, const double weight, const bool kFill)
{
  const double showerE = (*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3;
  if(showerE == -999*1E-3){
    return false;
  }
  
  // Get the position vector point from vertex to shower start position
  const TVector3 dist = anaUtils.GetRecShowerDistVector(ii);
  // Get reco and truth shower momentum vector in both lab frame
  const TLorentzVector recShowerMom = anaUtils.GetRecShowerLTVectLab(ii);//GetRecShowerRefBeam(false,ii);//GetRecShowerLTVectLab(ii);
  const TLorentzVector recShowerMomRaw = anaUtils.GetRecShowerLTVectLab(ii,false);//GetRecShowerRefBeam(false,ii,false);//GetRecShowerLTVectLab(ii,false);
  const TLorentzVector truthShowerMom = anaUtils.GetTruthMatchedShowerLTVectLab(ii);//GetRecShowerRefBeam(true,ii);//GetTruthMatchedShowerLTVectLab(ii);
  // Calculate shower impact parameter
  const double IP = dist.Mag()*TMath::Sin((recShowerMom.Angle(dist)));
  // Get the shower position vector
  const TVector3 showerPosition((*AnaIO::reco_daughter_allShower_startX)[ii],(*AnaIO::reco_daughter_allShower_startY)[ii],(*AnaIO::reco_daughter_allShower_startZ)[ii]);
  // Fill Cut histogram

  // In unit of cm
  if( dist.Mag() < 3 || dist.Mag() > 90 ) {
    return false;
  }
  // Impact Parameter Cut
  if( IP > 20 ){
    return false;
  }
  // Need to save all pizero shower candidates to reconstruct pizero
  anaUtils.SavePiZeroShower(recShowerMom, recShowerMomRaw, truthShowerMom, recShowerMom.E(), truthShowerMom.E(), showerPosition, truthParticleType);
  return true;
}


bool AnaCut::IsPizero(const bool kMC, const bool kFill, const double weight)
{
  // Pi zero selection
  OA = -999;
  truthPi0Type = -999;
  
  // Get pi0 vector no need to do KF in this step
  PiZeroVec = anaUtils.GetRecPiZeroFromShowers(OA,kMC,kFill,false,truthPi0Type);

  double mass = PiZeroVec.M();

  if(npi0shower < 2) {
    return false;
  }
  // Invariant mass cut
  if(mass < 0.05 || mass > 0.25){
    return false;
  }

  // Opening angle cut
  if(OA < 10 || OA > 80){
    return false;
  }

  return true;
}




