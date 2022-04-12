#include "../include/AnaUtils.h"
#include "AnaFit.cxx"

int AnaUtils::GetBeamParticleType(const int pdg)
{
  int type = -999;
  if(pdg==2212){
    type = gkBeamProton;
  }
  else if(pdg==211){
    type = gkBeamPiPlus;
  }
  else if(pdg==-13){
    type = gkBeamMuPlus;
  }
  else if(pdg==11 || pdg==-11 || pdg==22){
    type = gkBeamElectronGamma;
  }
  else if(pdg==-211){
    type = gkBeamPiMinus;
  }
  else{
    type = gkBeamOthers;
  }
  return type;
}

int AnaUtils::GetParticleType(const int pdg)
{
  int type = -999;
  if(pdg==2212){
    type = gkProton;
  }
  else if(pdg==211){
    type = gkPiPlus;
  }
  else if(pdg==-211){
    type = gkPiMinus;
  }
  else if(pdg==111){
    type = gkPiZero;
  }
  else if(pdg==-11){
    type = gkPositron;
  }
  else if(pdg==11){
    type = gkElectron;
  }
  else if(pdg==-13){
    type = gkMuPlus;
  }
  else if(pdg==13){
    type = gkMuMinus;
  }
  else if(pdg==321||pdg==-321||pdg==310||pdg==130||pdg==331){
    type = gkKaon;
  }
 else if(pdg==2112){
    type = gkNeutron;
  }
  else if(pdg==22){
    type = gkGamma;
  }
  else if(pdg==14 || pdg==-14 || pdg==12 || pdg==-12){
    type = gkNeutrino;
  }
  else if(pdg==3112||pdg==3122||pdg==3212||pdg==3222||pdg==3312||pdg==3322){
    type = gkHyperon;
  }
  else if(pdg==221){
    type = gkMesons;
  }
  else if(pdg>9999){
    type = gkNucleus;
  }
  else if(pdg==-1){
    type = gkNoTruth;
  }
  else{
    cout<<"AnaUtils::GetParticleType unknown pdg "<<pdg<<endl; exit(1);
  }

  return type;
}

int AnaUtils::GetFillEventType()
{
  int filleventtype = -999;
  if(AnaIO::Signal){
    filleventtype = gkSignal;
  }
  else if(AnaIO::true_beam_PDG==211){
    filleventtype = gkEvtBkg;
  }
  else{
    filleventtype = gkBmBkg;
  }

  return filleventtype;
}

void AnaUtils::SetFullSignal()
{
  // Phase space cut on protons 
  // 0.45 GeV/c is the reconstruction threshold, 1 GeV/c is limit where momentum by range is reliable
  // Leading proton momentum 0.45 - 1 GeV/c
  // Subleading proton momentum < 0.45 GeV/c
  // No cuts on pions

  // Make sure we start from false
  AnaIO::Signal = false;
  // Get final state particles vector in this event
  vector<TLorentzVector> vecFSParticle = GetFSParticlesTruth();
  // Check the pi0 daughters 
  if(nPiZero > 0) GetFSPiZeroDecayDaughterTruth();
  // Get the FS particles momentum
  double LeadingPiZeroP = vecFSParticle[0].P();
  double LeadingProtonP = vecFSParticle[1].P();
  double SubLeadingProtonP = vecFSParticle[2].P();

  // Check event topology (TKI event)
  // Initial pion beam and at least one proton and one pizero, no other mesons in final state (but not consider number of neutrons)
  
  if( AnaIO::true_beam_PDG==211 && IsSignal(nProton,nPiZero,nPiPlus,nParticleBkg) == true){
    // Proton momentum selection (below 0.45 GeV/c is not detectbale)
    if(LeadingProtonP < 1 && LeadingProtonP > 0.45 && SubLeadingProtonP < 0.45){
      // No restrictions on pizero momentum
      if(LeadingPiZeroP > 0){} //cout << "This is a good TKI event" << endl; //AnaIO::Signal = true; 
    }
  }
  
  // Check Pi0 inclusive events
  if( AnaIO::true_beam_PDG == 211 && IsSignal(nProton,nPiZero,nPiPlus,nParticleBkg) == true){
    AnaIO::Signal = true; 
  }
  
}

vector<TLorentzVector> AnaUtils::GetFSParticlesTruth()
{
  // Get true beam daughter information
  const vector<int> * pdg = AnaIO::true_beam_daughter_PDG;
  const vector<double> * px = AnaIO::true_beam_daughter_startPx; // GeV/c
  const vector<double> * py = AnaIO::true_beam_daughter_startPy; // GeV/c
  const vector<double> * pz = AnaIO::true_beam_daughter_startPz; // GeV/c

  // Class member variables (beam truth daughter particles counter)
  nProton = 0;
  nNeutron = 0;
  nPiPlus = 0;
  nPiZero = 0;
  nGamma = 0;
  nParticleBkg = 0;
  // Vectors to save particle info
  TLorentzVector pPiZero, pProton, pSecondaryProton;

  // Get the size of final state particles
  const int np = pdg->size();
  // Fill the hTruthFSParticleNumber hitogram
  AnaIO::hTruthFSParticleNumber->Fill(np);

  double Protonmom[np];
  double PiZeromom[np];
  double Gammamom[np];
  vector<TVector3> bufferProtonmom;
  vector<TVector3> bufferPiZeromom;
  bufferType.clear();
  // Now loop over FS particles
  for(int ii=0; ii<np; ii++){
    // Get the FS particle type
    const int itype = GetParticleType((*pdg)[ii]);
    // Get the FS particle 3-momentum
    const TVector3 tmpp( (*px)[ii], (*py)[ii], (*pz)[ii] );
    // Fill the FS particle type
    AnaIO::hTruthFSParticleType->Fill(itype); 

    // Check each FS type and save info
    // Proton
    if(itype == gkProton){
      // Save particle's momentum magnitude
      Protonmom[nProton] = tmpp.Mag();
      // Save momentum vector
      bufferProtonmom.push_back(tmpp);
      // Increase proton number 
      nProton++; 
    } 
    // PiZero
    else if(itype == gkPiZero){
      // Save particle's momentum magnitude
      PiZeromom[nPiZero] = tmpp.Mag();
      // Save momentum vector
      bufferPiZeromom.push_back(tmpp);      
      // Increase proton number 
      nPiZero++;
    }
    // Gamma
    else if(itype == gkGamma){
      Gammamom[nGamma] = tmpp.Mag();
      nGamma++;
    }
    // PiPlus
    else if(itype == gkPiPlus){
      nPiPlus++;
    }
    // Neutron
    else if(itype == gkNeutron){
      nNeutron++;
    }
    // PiMinus and Kaon
    else if(itype==gkPiMinus||itype==gkKaon){
      nParticleBkg++;
    }
    // Store FS particle type info
    bufferType.push_back(itype);

  } // End of loop over FS particles

  //======================= Proton =======================
  int leadingProtonID = 0, subldProtonID = -999;
  // Sort protons when more than one protons are found
  if(nProton>1){
    int Protonsortid[nProton];
    // Sort index according to it's momentum
    TMath::Sort(nProton, Protonmom, Protonsortid);
    // Save sorted index
    leadingProtonID = Protonsortid[0];
    subldProtonID = Protonsortid[1];
  }
  // At least one proton is found
  if(nProton>0){
    // Save info to leading proton TLorentzVector
    pProton.SetVectM(bufferProtonmom[leadingProtonID], AnaFunctions::ProtonMass());
    // Fill leading proton momentum
    AnaIO::hTruthLeadingProtonP->Fill(Protonmom[leadingProtonID]);
  }
  // At least two proton is found
  if(nProton>1){
    // Save info to subleading proton TLorentzVector
    pSecondaryProton.SetVectM(bufferProtonmom[subldProtonID], AnaFunctions::ProtonMass());
    // Fill subleading proton momentum
    AnaIO::hTruthSubLeadingProtonP->Fill(Protonmom[subldProtonID]);
  }

  //======================== PiZero ========================
  int leadingPiZeroID = 0, subldPiZeroID = -999;
  if(nPiZero>1){
    // Fill the FS pi0 number (at least two pi0)
    AnaIO::hTruthFSMultiPi0->Fill(nPiZero);
    int PiZerosortid[nPiZero];
    // Sort index according to it's momentum
    TMath::Sort(nPiZero, PiZeromom, PiZerosortid);
    // Save sorted index
    leadingPiZeroID = PiZerosortid[0];
    subldPiZeroID = PiZerosortid[1];
  }
  if(nPiZero>0){
    // Fill histogram for FS pi0 number
    AnaIO::hTruthFSPi0Number->Fill(nPiZero);
    // Save info to pi0 TLorentzVector
    pPiZero.SetVectM(bufferPiZeromom[leadingPiZeroID], AnaFunctions::PiZeroMass());
    AnaIO::hTruthLeadingPiZeroP->Fill(PiZeromom[leadingPiZeroID]);
    AnaIO::hTruthLeadingPiZeroE->Fill(pPiZero.E());
  }
  if(nPiZero>1){
    AnaIO::hTruthSubLeadingPiZeroP->Fill(PiZeromom[subldPiZeroID]);
  }

  //======================== Gamma ========================
  int leadingGammaID = 0;
  if(nGamma>1){
    int Gammasortid[nGamma];
    TMath::Sort(nGamma, Gammamom, Gammasortid);
    leadingGammaID = Gammasortid[0];
  }
  if(nGamma>0) AnaIO::hTruthGammaMaxE->Fill(Gammamom[leadingGammaID]);

  // Fill vector of FS particles 
  vector<TLorentzVector> vec;
  vec.push_back(pPiZero);
  vec.push_back(pProton);
  vec.push_back(pSecondaryProton);
  return vec;
}

vector<TLorentzVector> AnaUtils::GetFSPiZeroDecayDaughterTruth()
{
  // Get beam vertex information
  const TVector3 vtx(AnaIO::true_beam_endX, AnaIO::true_beam_endY, AnaIO::true_beam_endZ);

  // Get true pi0 daughter information
  const vector<int> * pdg = AnaIO::true_beam_Pi0_decay_PDG;
  const vector<double> * px = AnaIO::true_beam_Pi0_decay_startPx; // GeV/c
  const vector<double> * py = AnaIO::true_beam_Pi0_decay_startPy; // GeV/c
  const vector<double> * pz = AnaIO::true_beam_Pi0_decay_startPz; // GeV/c

  const vector<double> * startX = AnaIO::true_beam_Pi0_decay_startX;
  const vector<double> * startY = AnaIO::true_beam_Pi0_decay_startY;
  const vector<double> * startZ = AnaIO::true_beam_Pi0_decay_startZ;

   // Counter
  nPiZeroGamma = 0;
  nPiZeroElectron = 0;
  nPiZeroPositron = 0;
  // Pi0 decay gamma momentum
  TLorentzVector LeadingGamma, SubldGamma, Electron, Positron;
  // Get the size of pi0 daughter particles
  const int np = pdg->size();
  AnaIO::hTruthPi0DecayParticleNumber->Fill(np);

  double PiZeroGammamom[np];
  double PiZeroElectronmom[np];
  double PiZeroPositronmom[np];
  double ShowerDist[np];
  vector<TVector3> bufferPiZeroGammamom;
  vector<TVector3> bufferPiZeroElectronmom;
  vector<TVector3> bufferPiZeroPositronmom;

  // Now loop over pi0 daughter particles
  for(int ii=0; ii<np; ii++){
    // Get the pi0 daughter particle type
    const int itype = GetParticleType((*pdg)[ii]);
    // Get the pi0 daughter particle 3-momentum
    const TVector3 tmpp( (*px)[ii], (*py)[ii], (*pz)[ii] );
    // Start position of beam daughter
    const TVector3 startPos((*startX)[ii], (*startY)[ii], (*startZ)[ii]);
    const TVector3 dist = vtx - startPos;
    // Gamma
    if(itype == gkGamma){
      // Save particle's momentum magnitude
      PiZeroGammamom[nPiZeroGamma] = tmpp.Mag();
      // Save momentum vector
      bufferPiZeroGammamom.push_back(tmpp);
      ShowerDist[nPiZeroGamma] = dist.Mag();
      // Increase proton number 
      nPiZeroGamma++; 
    }
    
    else if(itype == gkElectron){
      PiZeroElectronmom[nPiZeroElectron] = tmpp.Mag();
      bufferPiZeroElectronmom.push_back(tmpp);
      nPiZeroElectron++;
    }

    else if(itype == gkPositron){
      PiZeroPositronmom[nPiZeroElectron] = tmpp.Mag();
      bufferPiZeroPositronmom.push_back(tmpp);
      nPiZeroPositron++;
    }
    
  } // end of loop

  //======================== Pi0 Gammas ========================
  int leadingGammaID = 0, subldGammaID = -999;
  // Decays into two gammas
  if(nPiZeroGamma>1){
    int PiZeroGammasortid[nPiZeroGamma];
    // Sort index according to it's momentum
    TMath::Sort(nPiZeroGamma, PiZeroGammamom, PiZeroGammasortid);
    // Save sorted index
    leadingGammaID = PiZeroGammasortid[0];
    subldGammaID = PiZeroGammasortid[1];
    // Save info to gammas TLorentzVector
    LeadingGamma.SetVectM(bufferPiZeroGammamom[leadingGammaID], 0);
    SubldGamma.SetVectM(bufferPiZeroGammamom[subldGammaID], 0);
    AnaIO::hTruthLeadingPi0GammaP->Fill(PiZeroGammamom[leadingGammaID]);
    AnaIO::hTruthSubLeadingPi0GammaP->Fill(PiZeroGammamom[subldGammaID]);
    // Opening angle
    double OA = (bufferPiZeroGammamom[0].Angle(bufferPiZeroGammamom[1]))*TMath::RadToDeg();
    AnaIO::hTruthPi0OA->Fill(OA);
    AnaIO::hTruthLeadingPi0GammaOA->Fill(OA);
    AnaIO::hTruthSubLeadingPi0GammaOA->Fill(OA);
    plotUtils.FillHist(AnaIO::hTruthPi0GammaEnergy,LeadingGamma.E(),0);
    plotUtils.FillHist(AnaIO::hTruthPi0GammaEnergy,SubldGamma.E(),1);

    AnaIO::hTruthLeadingPiZeroGammaDist->Fill(ShowerDist[leadingGammaID]);
    AnaIO::hTruthSubLeadingPiZeroGammaDist->Fill(ShowerDist[subldGammaID]);

  }
  // Rare decay - one gamma one electron and one positron
  else {
    LeadingGamma.SetVectM(bufferPiZeroGammamom[leadingGammaID], 0);
    Electron.SetVectM(bufferPiZeroElectronmom[0], AnaFunctions::ElectronMass());
    Positron.SetVectM(bufferPiZeroPositronmom[0], AnaFunctions::ElectronMass());
    AnaIO::hTruthRarePi0GammaP->Fill(PiZeroGammamom[leadingGammaID]);
    AnaIO::hTruthRarePi0ElectronP->Fill(PiZeroElectronmom[0]);
    AnaIO::hTruthRarePi0PositronP->Fill(PiZeroPositronmom[0]);
  }

  // Fill vector of FS particles 
  vector<TLorentzVector> vec;
  vec.push_back(LeadingGamma);
  vec.push_back(SubldGamma);
  return vec; 

}

bool AnaUtils::IsSignal(const int nProton, const int nPiZero, const int nPiPlus, const int nParticleBkg)
{
  bool tmpSig = false;
  // TKI event selection 
  //if(nProton > 0 && nPiZero > 0 && nPiPlus == 0 && nParticleBkg == 0) tmpSig = true;
  // Pi0 inclusive event selection
  if(nPiZero > 0) tmpSig = true;
  return tmpSig;
} 

TVector3 AnaUtils::GetRecBeamFull(){

  TVector3 beamdir;
  // Set beam end direction
  // index 1 uses direction from line projected between last 2 points;
  beamdir.SetXYZ((*AnaIO::reco_beam_calo_endDirX)[1], 
                 (*AnaIO::reco_beam_calo_endDirY)[1], 
                 (*AnaIO::reco_beam_calo_endDirZ)[1] );
  // Get the beam end kinetic energy
  double ke = AnaIO::reco_beam_interactingEnergy/1E3;
  if(ke<0) ke = 1E-10;
  // Get pion mass since signal is pion beam
  const double mpi = AnaFunctions::PionMass();
  // Calculate pion beam end momentum
  const double pionEndP = TMath::Sqrt(ke*ke+2*ke*mpi);
  // Get momentum 3 vector
  const TVector3 fullbeam = beamdir.Unit()*pionEndP;
  return fullbeam;
}

TVector3 AnaUtils::GetTruthBeamFull()
{
  //const TVector3 tmpbeam(AnaIO::true_beam_endPx,
  //                       AnaIO::true_beam_endPy,
  //                       AnaIO::true_beam_endPz );
  const TVector3 tmpbeam(AnaIO::reco_beam_true_byHits_endPx,
                         AnaIO::reco_beam_true_byHits_endPy,
                         AnaIO::reco_beam_true_byHits_endPz );
   
  return tmpbeam;
}

void AnaUtils::FillBeamKinematics(const int kMC)
{
  // Get event type
  const int evtType = GetFillEventType();
  
  const TVector3 recBeamFull = GetRecBeamFull();
  if(kMC){
    const TVector3 truthBeamFull = GetTruthBeamFull();

    const double beamthetaRes    = (recBeamFull.Theta()-truthBeamFull.Theta())*TMath::RadToDeg();//use absolute difference 
    const double beammomentumRes = recBeamFull.Mag()/truthBeamFull.Mag()-1;

    plotUtils.FillHist(AnaIO::hBeamThetaRes,    truthBeamFull.Theta()*TMath::RadToDeg(), beamthetaRes);
    plotUtils.FillHist(AnaIO::hBeamMomentumRes, truthBeamFull.Mag(),                     beammomentumRes);
  }
  // This evtType only works for MC, data will not have this info but fill it anyway
  plotUtils.FillHist(AnaIO::hRecBeamTheta,    recBeamFull.Theta()*TMath::RadToDeg(), evtType);
  plotUtils.FillHist(AnaIO::hRecBeamMomentum, recBeamFull.Mag(),                     evtType);

}

vector<double> AnaUtils::GetdEdxVector(const vector<double> &arraydEdx, const bool kForward)
{
  // Get the size of this raw dEdx vector
  const unsigned int ncls = arraydEdx.size();

  vector<double> dEdxVec;
  // Start from [2] because [0] and [1] in both start and last are weird
  for(unsigned int kk=2; kk<ncls; kk++){
    // Fill the vector start from beginning of arraydEdx
    if(kForward) dEdxVec.push_back(arraydEdx[kk]);
    // Fill it backwards
    else{
     const double endpe = arraydEdx[ncls-1-kk];
     dEdxVec.push_back(endpe);
    }
  }
  return dEdxVec;
}

double AnaUtils::GetTruncatedMean(const vector<double> & dEdxVec, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac)
{
  // Require nsample0 < nsample1 < dEdxVec size
  if(nsample1>=dEdxVec.size() || nsample0>=nsample1) return -999;
  
  vector<double> dEdxVecTruncated;
  // Get the truncated dEdx vector specified by nsample0 and nsample1
  for(unsigned int ii=nsample0; ii<=nsample1; ii++){
    dEdxVecTruncated.push_back(dEdxVec[ii]);
  }
  // Get the upper and lower bound  
  const int iter0 = dEdxVecTruncated.size()*lowerFrac;
  const int iter1 = dEdxVecTruncated.size()*upperFrac;
  const int nterm = iter1-iter0;
  if( nterm<=0 ) return -999;
  // Sort this truncated dEdx vector 
  std::sort(dEdxVecTruncated.begin(), dEdxVecTruncated.end());
  double sum = 0.0;
  // Calcuate the mean value dEdx specified by lowerFrac and upperFrac
  for(int ii=iter0; ii< iter1; ii++){
    sum += dEdxVecTruncated[ii];
  }
  return sum / (nterm+1E-10);
}

double AnaUtils::GetFSParticleTME(const unsigned int ii, const bool kForward)
{
  // Get the raw dEdx vector
  const vector<double> recodEdxarray = (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
  // Get the meaningful dEdx vector
  vector<double> dEdxVec = GetdEdxVector(recodEdxarray,kForward);
  // Get the size of dEdxVec 
  const int ndEdx = dEdxVec.size();
  double TME = -999;
  
  if(kForward) TME = GetTruncatedMean(dEdxVec, 2, 6, 0.4,  0.95); 
  else TME  = GetTruncatedMean(dEdxVec, 2, ndEdx-8, 0.05, 0.6);

  return TME;
}

double AnaUtils::GetChi2NDF(const int ii)
{
  const double chi2       = (*AnaIO::reco_daughter_allTrack_Chi2_proton)[ii];
  const double ndof       = (*AnaIO::reco_daughter_allTrack_Chi2_ndof)[ii];

  const double chnf = chi2/(ndof+1E-10);

  return chnf;
}

TVector3 AnaUtils::GetTruthMatchedTrackVectLab(const int ii)
{
  const TVector3 trackVectLab((*AnaIO::reco_daughter_PFP_true_byHits_startPx)[ii],
                              (*AnaIO::reco_daughter_PFP_true_byHits_startPy)[ii],
                              (*AnaIO::reco_daughter_PFP_true_byHits_startPz)[ii] );
  return trackVectLab;
}

TVector3 AnaUtils::GetRecTrackVectLab(const int ii, const bool kProton, bool DoCorrection)
{ 
  // MBR: momentum by range
  double trackMBR = -999;
  if (DoCorrection && (*AnaIO::reco_daughter_allTrack_momByRange_alt_proton)[ii] > 0.45) {
    trackMBR = kProton? GetProtonCorrectedMom((*AnaIO::reco_daughter_allTrack_momByRange_alt_proton)[ii]) : (*AnaIO::reco_daughter_allTrack_momByRange_alt_muon)[ii]; 
  }
  else trackMBR = kProton? (*AnaIO::reco_daughter_allTrack_momByRange_alt_proton)[ii] : (*AnaIO::reco_daughter_allTrack_momByRange_alt_muon)[ii]; 

  TVector3 trackVectLab;
  // Get this reco particle momentum vector in lab frame
  
  //if(kProton && DoCorrection)  trackVectLab.SetMagThetaPhi(trackMBR, GetProtonCorrectedTheta((*AnaIO::reco_daughter_allTrack_Theta)[ii]), GetProtonCorrectedPhi((*AnaIO::reco_daughter_allTrack_Phi)[ii]));
  //else 
  trackVectLab.SetMagThetaPhi(trackMBR, (*AnaIO::reco_daughter_allTrack_Theta)[ii], (*AnaIO::reco_daughter_allTrack_Phi)[ii]);


  return trackVectLab;
}

TLorentzVector AnaUtils::GetMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton, bool DoProtonMomCorrection)
{
  // Get true/reco beam momentum vector
  const TVector3 tmpBeam = isTruth ? GetTruthBeamFull() : GetRecBeamFull();
  // Get true/reco FS particle momentum vector
  const TVector3 vecLab = isTruth ? GetTruthMatchedTrackVectLab(recIndex) : GetRecTrackVectLab(recIndex, kProton, DoProtonMomCorrection);
  // Get theta angle reletive to the beam
  const double thetaRefBeam = AnaFunctions::GetThetaRef(vecLab, tmpBeam.Unit());

  TVector3 vectRefBeam;
  vectRefBeam.SetMagThetaPhi(vecLab.Mag(), thetaRefBeam, 0);

  TLorentzVector momentumRefBeam;
  momentumRefBeam.SetVectM(vectRefBeam, kProton? AnaFunctions::ProtonMass() : AnaFunctions::PionMass() );

  return momentumRefBeam;
}

double AnaUtils::GetTransverseMomentumRefBeam(const bool isTruth, const int recIndex, const bool kProton, bool DoProtonMomCorrection)
{
  // Get true/reco beam momentum vector
  const TVector3 tmpBeam = isTruth ? GetTruthBeamFull() : GetRecBeamFull();
  // Get true/reco FS particle momentum vector
  const TVector3 vecLab = isTruth ? GetTruthMatchedTrackVectLab(recIndex) : GetRecTrackVectLab(recIndex, kProton, DoProtonMomCorrection);
  // Get theta angle reletive to the beam
  const double thetaRefBeam = AnaFunctions::GetThetaRef(vecLab, tmpBeam.Unit());

  double pT = vecLab.Mag()*sin(thetaRefBeam);

  return pT;
}


TLorentzVector AnaUtils::GetPi0MomentumRefBeam(const TLorentzVector RecPi0, const TLorentzVector TruthPi0, const bool isTruth)
{
  // Get reco beam momentum vector
  const TVector3 tmpBeam = isTruth? GetTruthBeamFull() : GetRecBeamFull();
  // Get reco Pi0 particle momentum vector
  const TVector3 vecLab = isTruth? TruthPi0.Vect() : RecPi0.Vect();
  // Get theta angle reletive to the beam
  const double thetaRefBeam = AnaFunctions::GetThetaRef(vecLab, tmpBeam.Unit());

  TVector3 vectRefBeam;
  vectRefBeam.SetMagThetaPhi(vecLab.Mag(), thetaRefBeam, 0);

  TLorentzVector momentumRefBeam;
  momentumRefBeam.SetVectM(vectRefBeam, AnaFunctions::PiZeroMass() );

  return momentumRefBeam;
}


void AnaUtils::FillFSParticleKinematics(const int recIndex, const int truthParticleType, const int recParticleType)
{
  // ---------------- Fill proton kinematics ----------------//
  if(recParticleType == gkProton){
    // Get the uncorrected shower momentum vector for comparison 
    const TLorentzVector recMomRawRefBeam = GetMomentumRefBeam(false /*=>reco*/, recIndex, true /*=>proton*/,false /*=>no correction*/);
    // Get this reco particle momentum vector relative to beam (corrected)
    const TLorentzVector recMomRefBeam = GetMomentumRefBeam(false /*=>reco*/, recIndex, true /*=>proton*/);
    // Get this reco particle momentum vector in lab frame
    const TVector3 recMomLab = GetRecTrackVectLab(recIndex, true, false);
    plotUtils.FillHist(AnaIO::hRecProtonMomentum,recMomRawRefBeam.P(), truthParticleType);
    plotUtils.FillHist(AnaIO::hRecProtonTheta, recMomRawRefBeam.Theta()*TMath::RadToDeg(), truthParticleType);
    // Truth-matching primary proton
    if(truthParticleType == gkProton){ // MC only (Data loop won't pass this)
      // Get the truth-matched particle truth momentum vector relative to beam
      const TLorentzVector truthMomRefBeam = GetMomentumRefBeam(true /*=>truth-matched*/, recIndex, true /*=>proton*/);
      // Get the truth-matched particle truth momentum vector in lab frame
      const TVector3 truthMomLab = GetTruthMatchedTrackVectLab(recIndex);
      // Relative to beam
      const double momentumRes = recMomRefBeam.P()/truthMomRefBeam.P()-1;
      const double momentumRawRes = recMomRawRefBeam.P()/truthMomRefBeam.P()-1;
      const double thetaRes    = (recMomRefBeam.Theta()-truthMomRefBeam.Theta())*TMath::RadToDeg();
      //const double thetaResRaw    = (recMomRawRefBeam.Theta()-truthMomRefBeam.Theta())*TMath::RadToDeg();
      // Lab frame
      const double thetaLabRes = (recMomLab.Theta()-truthMomLab.Theta())*TMath::RadToDeg();
      const double phiLabRes = (recMomLab.Phi()-truthMomLab.Phi())*TMath::RadToDeg();
      plotUtils.FillHist(AnaIO::hProtonMomentumRes, truthMomRefBeam.P(), momentumRes);
      plotUtils.FillHist(AnaIO::hProtonMomentumRawRes, truthMomRefBeam.P(), momentumRawRes);
      plotUtils.FillHist(AnaIO::hProtonThetaRes, truthMomRefBeam.Theta()*TMath::RadToDeg(), thetaRes);
      plotUtils.FillHist(AnaIO::hProtonThetaLabRes, truthMomLab.Theta()*TMath::RadToDeg(), thetaLabRes);
      plotUtils.FillHist(AnaIO::hProtonPhiLabRes, truthMomLab.Phi()*TMath::RadToDeg(), phiLabRes);
      // Rec VS truth
      plotUtils.FillHist(AnaIO::hProtonMomentumRawRecVSTruth_REG,recMomRawRefBeam.P(),truthMomRefBeam.P());
      plotUtils.FillHist(AnaIO::hProtonMomentumRecVSTruth_REG, recMomRefBeam.P(),truthMomRefBeam.P());
      plotUtils.FillHist(AnaIO::hProtonTransverseMomentumRawRecVSTruth_REG,GetTransverseMomentumRefBeam(false,recIndex,true,false),
                         GetTransverseMomentumRefBeam(true, recIndex, true));
      plotUtils.FillHist(AnaIO::hProtonTransverseMomentumRecVSTruth_REG,GetTransverseMomentumRefBeam(false,recIndex,true),
                         GetTransverseMomentumRefBeam(true, recIndex, true));
        
      plotUtils.FillHist(AnaIO::hProtonMomentumRecVSTruth_REG_Correction,recMomRawRefBeam.P(),momentumRawRes);
      plotUtils.FillHist(AnaIO::hProtonMomentumRecVSTruth_REG_AfterCor,recMomRefBeam.P(),truthMomRefBeam.P());

      plotUtils.FillHist(AnaIO::hProtonThetaRecVSTruth_REG_Correction,recMomLab.Theta()*TMath::RadToDeg(),thetaLabRes);
      plotUtils.FillHist(AnaIO::hProtonThetaRecVSTruth_REG_AfterCor,recMomLab.Theta()*TMath::RadToDeg(),truthMomLab.Theta()*TMath::RadToDeg());

      plotUtils.FillHist(AnaIO::hProtonPhiRecVSTruth_REG_Correction,recMomLab.Phi()*TMath::RadToDeg(),phiLabRes);
      plotUtils.FillHist(AnaIO::hProtonPhiRecVSTruth_REG_AfterCor,recMomLab.Phi()*TMath::RadToDeg(),truthMomLab.Phi()*TMath::RadToDeg());

      // Transverse momentum
      const double pTRaw = GetTransverseMomentumRefBeam(false,recIndex,true,false);
      const double pTTruth = GetTransverseMomentumRefBeam(true, recIndex, true);
      plotUtils.FillHist(AnaIO::hProtonTransverseMomentumRecVSTruth_REG_Correction, pTRaw, pTRaw/pTTruth - 1 );
      // Save proton truth and raw momentum
      //ProtonMomTruth.push_back(GetTransverseMomentumRefBeam(true,recIndex,true));
      //ProtonMomRaw.push_back(GetTransverseMomentumRefBeam(false,recIndex,true));
      ProtonMomTruth.push_back(truthMomRefBeam.P());
      ProtonMomRaw.push_back(recMomRefBeam.P());
    }
  }
  // ---------------- Fill piplus kinematics ----------------//
  if(recParticleType == gkPiPlus){
    // Get this reco particle momentum vector relative to beam
    const TLorentzVector recMomRefBeam = GetMomentumRefBeam(false /*=>reco*/, recIndex, false/*=>piplus*/);
    plotUtils.FillHist(AnaIO::hRecPiPlusMomentum,recMomRefBeam.P(), truthParticleType);
    plotUtils.FillHist(AnaIO::hRecPiPlusTheta, recMomRefBeam.Theta()*TMath::RadToDeg(), truthParticleType);
    // Truth-matching primary piplus
    if(truthParticleType == gkPiPlus){ // MC only (Data loop won't pass this)
      // Get the truth-matched particle truth momentum vector relative to beam
      const TLorentzVector truthMomRefBeam = GetMomentumRefBeam(true /*=>truth-matched*/, recIndex, false /*=>piplus*/);
      const double momentumRes = recMomRefBeam.P()/truthMomRefBeam.P()-1;
      const double thetaRes    = (recMomRefBeam.Theta()-truthMomRefBeam.Theta())*TMath::RadToDeg();
      plotUtils.FillHist(AnaIO::hPiPlusMomentumRes, truthMomRefBeam.P(), momentumRes);
      plotUtils.FillHist(AnaIO::hPiPlusThetaRes, truthMomRefBeam.Theta()*TMath::RadToDeg(), thetaRes);
    }  
  }
  // ---------------- Fill shower kinematics ----------------//
  if(recParticleType == gkGamma){
    // Get this reco particle momentum vector in lab frame 
    const TLorentzVector recShowerMom = GetRecShowerLTVectLab(recIndex);//GetRecShowerRefBeam(false,recIndex); //GetRecShowerLTVectLab(recIndex);
    // Get the uncorrected shower momentum vector for comparison 
    const TLorentzVector recShowerMomRaw = GetRecShowerLTVectLab(recIndex,false);//GetRecShowerRefBeam(false,recIndex,false); //GetRecShowerLTVectLab(recIndex,false);
    plotUtils.FillHist(AnaIO::hRecShowerEnergy, recShowerMom.E(), truthParticleType);
    plotUtils.FillHist(AnaIO::hRecShowerTheta, recShowerMom.Theta()*TMath::RadToDeg(), truthParticleType);
    plotUtils.FillHist(AnaIO::hRecShowerEnergyRaw, recShowerMomRaw.E(), truthParticleType);
    // Truth-matching
    if(truthParticleType == gkGamma){ // MC only (Data loop won't pass this)
      // Get the truth-matched particle truth momentum vector in lab frame
      const TLorentzVector truthShowerMom = GetTruthMatchedShowerLTVectLab(recIndex); //GetRecShowerRefBeam(true,recIndex); //GetTruthMatchedShowerLTVectLab(recIndex); 
      const double momentumRes = recShowerMom.E()/truthShowerMom.E()-1;
      const double thetaRes    = (recShowerMom.Theta()-truthShowerMom.Theta())*TMath::RadToDeg();
      const double momentumResRaw = recShowerMomRaw.E()/truthShowerMom.E()-1;
      const double thetaResRaw    = (recShowerMomRaw.Theta()-truthShowerMom.Theta())*TMath::RadToDeg();
      const double phiRes    = (recShowerMom.Phi()-truthShowerMom.Phi())*TMath::RadToDeg();
      const double phiResRaw    = (recShowerMomRaw.Phi()-truthShowerMom.Phi())*TMath::RadToDeg();

      plotUtils.FillHist(AnaIO::hShowerEnergyRes, truthShowerMom.E(), momentumRes);
      plotUtils.FillHist(AnaIO::hShowerEnergyResRaw, truthShowerMom.E(), momentumResRaw);
      plotUtils.FillHist(AnaIO::hShowerThetaRes, truthShowerMom.Theta()*TMath::RadToDeg(), thetaRes);
      plotUtils.FillHist(AnaIO::hShowerThetaResRaw, truthShowerMom.Theta()*TMath::RadToDeg(), thetaResRaw);
      plotUtils.FillHist(AnaIO::hShowerPhiRes, truthShowerMom.Phi()*TMath::RadToDeg(), phiRes);
      plotUtils.FillHist(AnaIO::hShowerPhiResRaw, truthShowerMom.Phi()*TMath::RadToDeg(), phiResRaw);
      // Rec VS truth
      plotUtils.FillHist(AnaIO::hShowerEnergyRecVSTruth_REG, truthShowerMom.E(), recShowerMom.E());
      plotUtils.FillHist(AnaIO::hShowerEnergyRawRecVSTruth_REG, truthShowerMom.E(), recShowerMomRaw.E());
      const int nhits = (*AnaIO::reco_daughter_PFP_nHits_collection)[recIndex];
      plotUtils.FillHist(AnaIO::hShowerEnergyResVSnHits_REG, nhits, momentumRes);

      plotUtils.FillHist(AnaIO::hShowerEnergyRecVSTruth_REG_Correction, recShowerMomRaw.E(), momentumResRaw);
      plotUtils.FillHist(AnaIO::hShowerEnergyRecVSTruth_REG_AfterCor, recShowerMom.E(), truthShowerMom.E());

      plotUtils.FillHist(AnaIO::hShowerThetaRecVSTruth_REG_Correction, recShowerMomRaw.Theta()*TMath::RadToDeg(), thetaResRaw);
      plotUtils.FillHist(AnaIO::hShowerThetaRecVSTruth_REG_AfterCor, recShowerMom.Theta()*TMath::RadToDeg(), truthShowerMom.Theta()*TMath::RadToDeg());

      plotUtils.FillHist(AnaIO::hShowerPhiRecVSTruth_REG_Correction, recShowerMomRaw.Phi()*TMath::RadToDeg(), phiResRaw);
      plotUtils.FillHist(AnaIO::hShowerPhiRecVSTruth_REG_AfterCor, recShowerMom.Phi()*TMath::RadToDeg(), truthShowerMom.Phi()*TMath::RadToDeg());

      const TVector3 dist = GetRecShowerDistVector(recIndex);
      // Calculate shower impact parameter
      const double IP = dist.Mag()*TMath::Sin((recShowerMom.Angle(dist)));
      plotUtils.FillHist(AnaIO::hShowerEnergyResVSIP_REG, IP, momentumRes);
      /*
      cout << "event: " << AnaIO::event << endl;
      cout << "Truth Shower energy: " << truthShowerMom.E() << endl;
      cout << "Reco Raw Shower energy: " << recShowerMomRaw.E() << endl;
      cout << "Reco Shower energy: " << recShowerMom.E() << endl;
      */
    }
  }
}

TVector3 AnaUtils::GetRecShowerDistVector(const int ii)
{
  // Get reco beam end point
  const TVector3 vtx(AnaIO::reco_beam_endX, AnaIO::reco_beam_endY, AnaIO::reco_beam_endZ);
  // Get reco shower start position
  const TVector3 shw((*AnaIO::reco_daughter_allShower_startX)[ii], (*AnaIO::reco_daughter_allShower_startY)[ii], (*AnaIO::reco_daughter_allShower_startZ)[ii]);
  // Calculate distance vector
  const TVector3 dist=shw-vtx;
  return dist;
}

TVector3 AnaUtils::GetTruthMatchedShower3VectLab(const int ii)
{
  const TVector3 showerVectLab((*AnaIO::reco_daughter_PFP_true_byHits_startPx)[ii],
                              (*AnaIO::reco_daughter_PFP_true_byHits_startPy)[ii],
                              (*AnaIO::reco_daughter_PFP_true_byHits_startPz)[ii] );

  return showerVectLab;
}

TLorentzVector AnaUtils::GetTruthMatchedShowerLTVectLab(const int ii)
{
  const TVector3 showerVectLab((*AnaIO::reco_daughter_PFP_true_byHits_startPx)[ii],
                              (*AnaIO::reco_daughter_PFP_true_byHits_startPy)[ii],
                              (*AnaIO::reco_daughter_PFP_true_byHits_startPz)[ii] );

  TLorentzVector showerLTVectLab;
  showerLTVectLab.SetVectM(showerVectLab, 0 );

  return showerLTVectLab;
}

TVector3 AnaUtils::GetRecShower3VectLab(const int ii, bool DoCorrection)
{
  // Get reco shower direction
  const TVector3 showerDir((*AnaIO::reco_daughter_allShower_dirX)[ii],(*AnaIO::reco_daughter_allShower_dirY)[ii],(*AnaIO::reco_daughter_allShower_dirZ)[ii] );
  // Get reco shower energy in GeV
  const double showerE = (*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3;
  double showerE_corrected = -999;
  if (DoCorrection) showerE_corrected = GetShowerCorrectedE(showerE); //(showerE - 0.010848)/0.6846; //(showerE - 0.038008)/0.71237;
  else showerE_corrected = showerE;//GetShowerCorrectedE(showerE); //showerE;
  // Combine diretion and energy to get shower 3 momemtum
  const TVector3 showerVectLab = showerDir*showerE_corrected;
  
  return showerVectLab;
}

TLorentzVector AnaUtils::GetRecShowerLTVectLab(const int ii, bool DoCorrection)
{
  // Get reco shower direction
  TVector3 showerDir((*AnaIO::reco_daughter_allShower_dirX)[ii],(*AnaIO::reco_daughter_allShower_dirY)[ii],(*AnaIO::reco_daughter_allShower_dirZ)[ii] );
  // Get reco shower energy in GeV
  const double showerE = (*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3;
  double showerE_corrected = -999;
  if (DoCorrection /*&& showerE > 0.05*/){
    double Theta_corrected = GetShowerCorrectedTheta(showerDir.Theta()*TMath::RadToDeg())*TMath::DegToRad();
    double Phi_corrected = GetShowerCorrectedPhi(showerDir.Phi()*TMath::RadToDeg())*TMath::DegToRad();
    showerDir.SetTheta(Theta_corrected);
    showerDir.SetPhi(Phi_corrected);
    showerE_corrected = GetShowerCorrectedE(showerE); //(showerE - 0.010848)/0.6846; //(showerE - 0.038008)/0.71237;
  }
  else showerE_corrected = showerE;
  // Combine diretion and energy to get shower 3 momemtum
  const TVector3 showerVectLab = showerDir*showerE_corrected;
  
  TLorentzVector showerLTVectLab;
  showerLTVectLab.SetVectM(showerVectLab, 0 );

  return showerLTVectLab;
}

TLorentzVector AnaUtils::GetRecShowerRefBeam(const bool isTruth, const int ii, bool DoCorrection)
{
  // Get true/reco beam momentum vector
  const TVector3 tmpBeam = isTruth ? GetTruthBeamFull() : GetRecBeamFull();
  // Get true/reco FS particle momentum vector
  const TVector3 vecLab = isTruth ? GetTruthMatchedShower3VectLab(ii) : GetRecShower3VectLab(ii, DoCorrection);
  //const TVector3 vecLab = LTvecLab.Vect();
  // Get theta angle reletive to the beam
  const double thetaRefBeam = AnaFunctions::GetThetaRef(vecLab, tmpBeam.Unit());

  TVector3 vectRefBeam;
  vectRefBeam.SetMagThetaPhi(vecLab.Mag(), thetaRefBeam, 0);

  TLorentzVector momentumRefBeam;
  momentumRefBeam.SetVectM(vectRefBeam, 0);

  return momentumRefBeam;
}

TLorentzVector AnaUtils::GetRecPiZeroFromShowers(double &OA, bool kMC, bool kFill, bool kDoKF, int &truthPi0Type)
{
  // Get the true event type (signal,background or beam background) only works for MC
  const int truthEventType = GetFillEventType();
  // Get the size of shower array
  const int showerSize = showerArray.size();
  plotUtils.FillHist(AnaIO::hRecPi0Nshower, showerSize, truthEventType); 
  // Declare PiZero vector
  TLorentzVector PiZeroVec, PiZeroVec_MassCal;
  // Raw PiZero vector without energy correction
  TLorentzVector PiZeroVecRaw;
  // Truth-Matched PiZero vector
  TLorentzVector PiZeroTruthVec;

  // Need to have at least two showers to reconstruct pi0
  if(showerSize>=2){
    total++;
    double separation = -999;
    // Get two pi0 showers with leading and subleading energy
    vector<TLorentzVector> RecPi0Showers = GetTwoPi0Showers(separation,kMC,kFill,kDoKF);

    // Combine leading shower and subleading shower to get pi0 vector
    PiZeroVec = RecPi0Showers[0] + RecPi0Showers[1];
    // Without energy correction
    PiZeroVecRaw = RecPi0Showers[2] + RecPi0Showers[3];
    // Without KF correction
    PiZeroVec_MassCal = RecPi0Showers[4] + RecPi0Showers[5];

    if(kFill) AnaIO::hPi0Total->Fill(PiZeroVec.P());

    // TKI analysis pi0 vector
    RecPi0LTVet = PiZeroVec;
    // Save the truth type of two reco showers
    truthPi0Type = -999;
    // Get two truth-matched particles (maynot be photons indicated by the truthPi0Type)
    vector<TLorentzVector> TruthPi0Showers = GetTwoTruthMatchedPi0Showers(truthPi0Type,kFill);
    // Fill pi0 info
    if(kFill){
      const double openingAngle = RecPi0Showers[0].Angle(RecPi0Showers[1].Vect())*TMath::RadToDeg();
      const double openingAngleRaw = RecPi0Showers[2].Angle(RecPi0Showers[3].Vect())*TMath::RadToDeg();
      OA = openingAngle;
      // Fill pi0 mass and momentum
      plotUtils.FillHist(AnaIO::hRecPi0Mass, PiZeroVec_MassCal.M(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0MassRaw, PiZeroVecRaw.M(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0Momentum, PiZeroVec.P(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0MomentumRaw, PiZeroVecRaw.P(), truthPi0Type);

      plotUtils.FillHist(AnaIO::hRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0MassRaw_OVERLAY, PiZeroVecRaw.M(), truthPi0Type);

      plotUtils.FillHist(AnaIO::hRecPi0Momentum_OVERLAY, PiZeroVec.P(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0MomentumRaw_OVERLAY, PiZeroVecRaw.P(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0OA_OVERLAY, openingAngle, truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0OARaw_OVERLAY, openingAngleRaw, truthPi0Type);

      plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY, PiZeroVec.E(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0ShowerSep_OVERLAY, separation, truthPi0Type);

      plotUtils.FillHist(AnaIO::hRecPi0Theta_OVERLAY, PiZeroVec.Theta()*TMath::RadToDeg(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0Phi_OVERLAY, PiZeroVec.Phi()*TMath::RadToDeg(), truthPi0Type);

      plotUtils.FillHist(AnaIO::hRecPi0ThetaRaw_OVERLAY, PiZeroVecRaw.Theta()*TMath::RadToDeg(), truthPi0Type);
      plotUtils.FillHist(AnaIO::hRecPi0PhiRaw_OVERLAY, PiZeroVecRaw.Phi()*TMath::RadToDeg(), truthPi0Type);

      if(truthPi0Type == gkTwoGammasSamePi0){

        const double ldPhotonAngle = RecPi0Showers[0].Angle(TruthPi0Showers[0].Vect())*TMath::RadToDeg();
        const double slPhotonAngle = RecPi0Showers[1].Angle(TruthPi0Showers[1].Vect())*TMath::RadToDeg();

        plotUtils.FillHist(AnaIO::hLeadingPhotonAngleRes, TruthPi0Showers[0].E(), ldPhotonAngle);
        plotUtils.FillHist(AnaIO::hSubLeadingPhotonAngleRes, TruthPi0Showers[1].E(), slPhotonAngle);

        PiZeroTruthVec = TruthPi0Showers[0] + TruthPi0Showers[1];

        const double openingAngleTruth = TruthPi0Showers[0].Angle(TruthPi0Showers[1].Vect())*TMath::RadToDeg();

        plotUtils.FillHist(AnaIO::hOpeningAngleRes, PiZeroTruthVec.P(), openingAngle - openingAngleTruth);

        TruthPi0LTVet = PiZeroTruthVec;

        // Calculate pi0 mass resolution
        const double mpi0Res = PiZeroVec_MassCal.M()/PiZeroTruthVec.M() -1;
        const double mpi0ResRaw = PiZeroVecRaw.M()/PiZeroTruthVec.M() -1;
        plotUtils.FillHist(AnaIO::hPi0MassRes, PiZeroTruthVec.M(), mpi0Res);
        plotUtils.FillHist(AnaIO::hPi0MassResRaw, PiZeroTruthVec.M(), mpi0ResRaw);

        // Calculate pi0 momentum resolution
        const double pi0momRes = PiZeroVec.E()/PiZeroTruthVec.E()-1;
        const double pi0momResRaw = PiZeroVecRaw.E()/PiZeroTruthVec.E()-1;
        const double pi0ThetaRes = (PiZeroVec.Theta() - PiZeroTruthVec.Theta())*TMath::RadToDeg();
        const double pi0PhiRes = (PiZeroVec.Phi() - PiZeroTruthVec.Phi())*TMath::RadToDeg();
        const double pi0ThetaResRaw = (PiZeroVecRaw.Theta() - PiZeroTruthVec.Theta())*TMath::RadToDeg();
        const double pi0PhiResRaw = (PiZeroVecRaw.Phi() - PiZeroTruthVec.Phi())*TMath::RadToDeg();

        plotUtils.FillHist(AnaIO::hPi0MomentumRes, PiZeroTruthVec.E(), pi0momRes);
        plotUtils.FillHist(AnaIO::hPi0MomentumResRaw, PiZeroTruthVec.E(), pi0momResRaw);
        plotUtils.FillHist(AnaIO::hPi0ThetaRes, PiZeroTruthVec.Theta()*TMath::RadToDeg(), pi0ThetaRes); 
        plotUtils.FillHist(AnaIO::hPi0PhiRes, PiZeroTruthVec.Phi()*TMath::RadToDeg(), pi0PhiRes); 
        plotUtils.FillHist(AnaIO::hPi0ThetaResRaw, PiZeroTruthVec.Theta()*TMath::RadToDeg(), pi0ThetaResRaw); 
        plotUtils.FillHist(AnaIO::hPi0PhiResRaw, PiZeroTruthVec.Phi()*TMath::RadToDeg(), pi0PhiResRaw); 
    
        if(PiZeroTruthVec.M() > 0.134 && PiZeroTruthVec.M() < 0.135){
          selected++;
          if(kFill) AnaIO::hPi0Selected->Fill(PiZeroVec.P());
          AnaIO::hMatchedTruthPi0Momentum->Fill(PiZeroTruthVec.E());
          
          AnaIO::hMatchedTruthPi0OA->Fill(openingAngleTruth);
          AnaIO::hMatchedTruthPi0Mass->Fill(PiZeroTruthVec.M());
        }
      }
    }
    // Turn this on to save event number if you want to do event display
    //if( nProton > 1 ) cout << "event : " << AnaIO::event << endl;
  }

  return PiZeroVec;
}

vector<TLorentzVector> AnaUtils::GetTwoPi0Showers(double &separation, bool kMC, bool kFill, bool kDoKF)
{
  vector<TLorentzVector> pi0Showers;

  const int showerSize = showerArray.size();
  // Get the true event type (signal,background or beam background) only works for MC
  const int truthEventType = GetFillEventType();
  // Get [0] element of shower energy vector
  const double* shE = &(showerEarr[0]);
  int *nindex = new int[showerSize];
  // Sort shower energy
  TMath::Sort(showerSize, shE, nindex, true);
  // Set reco leading and subleading shower vector
  TLorentzVector ldShower = showerArray[nindex[0]];
  TLorentzVector slShower = showerArray[nindex[1]];

  // Set reco leading and subleading shower vector (these vectors will not be changed by KF - pi0 mass calculation)
  TLorentzVector ldShower_MassCal = showerArray[nindex[0]];
  TLorentzVector slShower_MassCal = showerArray[nindex[1]];

  // Set raw reco leading and subleading shower vector
  const TLorentzVector ldShowerRaw = showerArrayRaw[nindex[0]];
  const TLorentzVector slShowerRaw = showerArrayRaw[nindex[1]];
  
  if(kFill){
    plotUtils.FillHist(AnaIO::hRecLeadingShowerEnergy, ldShower.E(), showerTypeArray[nindex[0]]); 
    plotUtils.FillHist(AnaIO::hRecSubLeadingShowerEnergy, slShower.E(), showerTypeArray[nindex[1]]);
    plotUtils.FillHist(AnaIO::hRecLeadingShowerEnergyRaw, ldShowerRaw.E(), showerTypeArray[nindex[0]]);
    plotUtils.FillHist(AnaIO::hRecSubLeadingShowerEnergyRaw, slShowerRaw.E(), showerTypeArray[nindex[1]]);
    // Calculate reco opening angle
    const double openingAngle = ldShower.Angle(slShower.Vect())*TMath::RadToDeg();
    plotUtils.FillHist(AnaIO::hRecShowerOpenAngle, openingAngle, truthEventType);

    // Get position vector of two reco showers 
    const TVector3 ldShowerPos = showerPos[nindex[0]];
    const TVector3 slShowerPos = showerPos[nindex[1]];
    const TVector3 distVect = ldShowerPos - slShowerPos;
    // Calculate two reco showers start separation   
    separation = distVect.Mag();
  }

  if(kDoKF){
    vector<double> FittedVars;
    KF(ldShower, slShower, FittedVars);
    // Set the shower energy after KF
    ldShower.SetE(FittedVars[0]);
    slShower.SetE(FittedVars[1]);
  }

  // Output two showers
  pi0Showers.push_back(ldShower);
  pi0Showers.push_back(slShower);
  pi0Showers.push_back(ldShowerRaw);
  pi0Showers.push_back(slShowerRaw);
  pi0Showers.push_back(ldShower_MassCal);
  pi0Showers.push_back(slShower_MassCal);

  return pi0Showers;
}

vector<TLorentzVector> AnaUtils::GetTwoTruthMatchedPi0Showers(int &truthPi0Type, const bool &kFill)
{ 
  vector<TLorentzVector> pi0ShowersTruthMatched;
  const int showerSize = showerArray.size();
  // Get [0] element of shower energy vector
  const double* shE = &(showerEarr[0]);
  int *nindex = new int[showerSize];
  // Sort shower energy
  TMath::Sort(showerSize, shE, nindex, true);

  // Set reco leading and subleading shower vector
   TLorentzVector ldShower = showerArray[nindex[0]];
   TLorentzVector slShower = showerArray[nindex[1]];
  // Set raw reco leading and subleading shower vector
  const TLorentzVector ldShowerRaw = showerArrayRaw[nindex[0]];
  const TLorentzVector slShowerRaw = showerArrayRaw[nindex[1]];

  // Truth-matching
  // MC only (Data loop won't pass this)
  if(showerTypeArray[nindex[0]] == gkGamma && showerTypeArray[nindex[1]] == gkGamma){
    const TLorentzVector ld = showerTruthArray[nindex[0]];
    const TLorentzVector sl = showerTruthArray[nindex[1]];
    TLorentzVector pi0 = ld + sl;
    if(pi0.M() > 0.134 && pi0.M() < 0.135) truthPi0Type = gkTwoGammasSamePi0;
    else truthPi0Type = gkTwoGammasDiffPi0;
  }
  else if (showerTypeArray[nindex[0]] == gkGamma || showerTypeArray[nindex[1]] == gkGamma){
    truthPi0Type = gkOneGamma;
  }
  else truthPi0Type = gkNoGammas;

  // Get truth leading and subleading shower energy
  const TLorentzVector ldShowerTruth = showerTruthArray[nindex[0]];
  const TLorentzVector slShowerTruth = showerTruthArray[nindex[1]];

  // Truth-matched particles are photons and their invariant mass is pi0 mass
  if(truthPi0Type == gkTwoGammasSamePi0 && kFill){
    
    // Get the theta angle (relative to z axis)
    TVector3 unitZ(0,0,1);
    double TruthldTheta = (ldShowerTruth.Vect()).Angle(unitZ);
    double TruthslTheta = (slShowerTruth.Vect()).Angle(unitZ);
    AnaIO::hMatchedTruthldShowerTheta->Fill(TruthldTheta*TMath::RadToDeg());
    AnaIO::hMatchedTruthslShowerTheta->Fill(TruthslTheta*TMath::RadToDeg());
    const TLorentzVector PiZeroTruthVec = ldShowerTruth + slShowerTruth;

    // Calculate leading and subleading shower energy resolution
    const double ldShowerRes = (ldShower.E()/ldShowerTruth.E())-1;
    const double slShowerRes = (slShower.E()/slShowerTruth.E())-1;
    const double ldShowerResRaw = (ldShowerRaw.E()/ldShowerTruth.E())-1;
    const double slShowerResRaw = (slShowerRaw.E()/slShowerTruth.E())-1;
    
    
    plotUtils.FillHist(AnaIO::hLeadingShowerEnergyRes, ldShowerTruth.E(), ldShowerRes);
    plotUtils.FillHist(AnaIO::hSubLeadingShowerEnergyRes, slShowerTruth.E(), slShowerRes);
    plotUtils.FillHist(AnaIO::hLeadingShowerEnergyResRaw, ldShowerTruth.E(), ldShowerResRaw);
    plotUtils.FillHist(AnaIO::hSubLeadingShowerEnergyResRaw, slShowerTruth.E(), slShowerResRaw);
    
    // Calculate opening angle resolution
    const double openingAngleTruth = ldShowerTruth.Angle(slShowerTruth.Vect())*TMath::RadToDeg();
    const double openingAngle = ldShower.Angle(slShower.Vect())*TMath::RadToDeg();
    const double openingAngleRes = openingAngle - openingAngleTruth;
    plotUtils.FillHist(AnaIO::hShowerOpenAngleRes, openingAngleTruth, openingAngleRes); 

    // Correction? (no need)
    plotUtils.FillHist(AnaIO::hLeadingShowerEnergyRecVSTruth_REG_Correction, ldShowerRaw.E(), ldShowerResRaw);
    plotUtils.FillHist(AnaIO::hSubLeadingShowerEnergyRecVSTruth_REG_Correction, slShowerRaw.E(), slShowerResRaw);
    plotUtils.FillHist(AnaIO::hOpeningAngleRecVSTruth_REG_Correction, openingAngle, openingAngleRes);
    
    AnaIO::hMatchedTruthldShowerEnergy->Fill(ldShowerTruth.E());
    AnaIO::hMatchedTruthslShowerEnergy->Fill(slShowerTruth.E());
    AnaIO::hMatchedTruthldShowerOA->Fill(openingAngleTruth);
    AnaIO::hMatchedTruthslShowerOA->Fill(openingAngleTruth);
  
  }

  pi0ShowersTruthMatched.push_back(ldShowerTruth);
  pi0ShowersTruthMatched.push_back(slShowerTruth);

  return pi0ShowersTruthMatched;
}

void AnaUtils::SavePi0ShowersForKF()
{

  AnaIO::LeadingShowerEnergy = -1;
  AnaIO::SubLeadingShowerEnergy = -1;
  AnaIO::LeadingShowerEnergyRaw = -1;
  AnaIO::SubLeadingShowerEnergyRaw = -1;
  AnaIO::OpeningAngle = -1;
  AnaIO::LeadingShowerEnergyTruth = -1;
  AnaIO::SubLeadingShowerEnergyTruth = -1;
  AnaIO::OpeningAngleTruth = -1;

  // Get the size of truth shower array
  const int showerSize = showerArray.size();
  // Need to have at least two showers to reconstruct pi0
  if(showerSize>=2){
    // Get [0] element of shower energy vector
    const double* shE = &(showerEarr[0]);
    int *nindex = new int[showerSize];
    // Sort truth shower energy
    TMath::Sort(showerSize, shE, nindex, true);
    if( showerTypeArray[nindex[0]] == gkGamma && showerTypeArray[nindex[1]] == gkGamma ){
      // Get truth leading ans subleading shower energy
      const TLorentzVector ldShowerTruth = showerTruthArray[nindex[0]];
      const TLorentzVector slShowerTruth = showerTruthArray[nindex[1]];

      // Get truth pi0 vector
      TLorentzVector PiZeroTruthVec = ldShowerTruth + slShowerTruth;

      // Get the truth-mathched photon type
      int truthPi0Type = -999;
      if(showerTypeArray[nindex[0]] == gkGamma && showerTypeArray[nindex[1]] == gkGamma){
        const TLorentzVector ld = showerTruthArray[nindex[0]];
        const TLorentzVector sl = showerTruthArray[nindex[1]];
        TLorentzVector pi0 = ld + sl;
        if(pi0.M() > 0.134 && pi0.M() < 0.135) truthPi0Type = gkTwoGammasSamePi0;
        else truthPi0Type = gkTwoGammasDiffPi0;
      }
      else if (showerTypeArray[nindex[0]] == gkGamma || showerTypeArray[nindex[1]] == gkGamma){
        truthPi0Type = gkOneGamma;
      }
      else {
        truthPi0Type = gkNoGammas;
        cout << "gkNoGammas!" << endl;
      }
      

      cout << "truthPi0Type: " << truthPi0Type << endl;
      //if(truthPi0Type == gkTwoGammasSamePi0 && PiZeroTruthVec.M() < 0.1350 && PiZeroTruthVec.M() > 0.1349) {
      if(truthPi0Type == gkTwoGammasSamePi0) {
        // Set reco leading and subleading shower vector
        const TLorentzVector ldShower = showerArray[nindex[0]];
        const TLorentzVector slShower = showerArray[nindex[1]];
        // Set raw reco leading and subleading shower vector
        const TLorentzVector ldShowerRaw = showerArrayRaw[nindex[0]];
        const TLorentzVector slShowerRaw = showerArrayRaw[nindex[1]];
        // Calculate reco opening angle
        const double openingAngle = ldShower.Angle(slShower.Vect())*TMath::RadToDeg();
        // Calculate truth opening angle
        const double openingAngleTruth = ldShowerTruth.Angle(slShowerTruth.Vect())*TMath::RadToDeg();

        // Save KF shower info
        /*if(ldShowerTruth.E() < slShowerTruth.E()){
          LdShowerEnergyTruth.push_back(slShowerTruth.E());
          SlShowerEnergyTruth.push_back(ldShowerTruth.E());
          LdShowerDirTruth.push_back(slShowerTruth.Vect());
          SlShowerDirTruth.push_back(ldShowerTruth.Vect());
        }*/
        //else{
          LdShowerEnergyTruth.push_back(ldShowerTruth.E());
          SlShowerEnergyTruth.push_back(slShowerTruth.E());
          LdShowerDirTruth.push_back(ldShowerTruth.Vect());
          SlShowerDirTruth.push_back(slShowerTruth.Vect());
        //}

        OpenAngleTruth.push_back(openingAngleTruth*TMath::DegToRad());
        LdShowerEnergyRaw.push_back(ldShower.E());
        SlShowerEnergyRaw.push_back(slShower.E());
        OpenAngle.push_back(openingAngle*TMath::DegToRad());

        LdShowerDir.push_back(ldShower.Vect());
        SlShowerDir.push_back(slShower.Vect());
        
        // Set output tree variables
        if(ldShowerTruth.E() < slShowerTruth.E()){
          AnaIO::LeadingShowerEnergyTruth = slShowerTruth.E();
          AnaIO::SubLeadingShowerEnergyTruth = ldShowerTruth.E();
          AnaIO::LeadingShowerEnergyUnitDirTruth = slShowerTruth.Vect();
          AnaIO::SubLeadingShowerEnergyUnitDirTruth = ldShowerTruth.Vect();
        }
        else{
          AnaIO::LeadingShowerEnergyTruth = ldShowerTruth.E();
          AnaIO::SubLeadingShowerEnergyTruth = slShowerTruth.E();
          AnaIO::LeadingShowerEnergyUnitDirTruth = ldShowerTruth.Vect();
          AnaIO::SubLeadingShowerEnergyUnitDirTruth = slShowerTruth.Vect();
        }

        AnaIO::LeadingShowerEnergy = ldShower.E();
        AnaIO::SubLeadingShowerEnergy = slShower.E();
        AnaIO::LeadingShowerEnergyRaw = ldShowerRaw.E();
        AnaIO::SubLeadingShowerEnergyRaw = slShowerRaw.E();
        AnaIO::OpeningAngle = openingAngle*TMath::DegToRad();
        AnaIO::LeadingShowerEnergyUnitDir = ldShower.Vect();
        AnaIO::SubLeadingShowerEnergyUnitDir = slShower.Vect();
        AnaIO::OpeningAngleTruth = openingAngleTruth*TMath::DegToRad();

      }
    }
  }

}
void AnaUtils::TruthMatchingTKI(TLorentzVector dummypi0, TLorentzVector dummyproton, TLorentzVector dummypi0Truth, TLorentzVector dummyprotonTruth, const bool kMC)
{
  const int truthEventType = GetFillEventType();
  const TLorentzVector beamFullP(GetRecBeamFull(), AnaFunctions::PionMass());
  const TLorentzVector beamFullPTruth(GetTruthBeamFull(), AnaFunctions::PionMass());

  const int targetA = 40;
  const int targetZ = 18;
  double dalphat,dphit,dpt,pn,finPitheta,finProtontheta;
  double dalphat_truth,dphit_truth,dpt_truth,pn_truth,finPitheta_truth,finProtontheta_truth;

  AnaFunctions::getCommonTKI(targetA, targetZ, &beamFullP, &(dummypi0), &(dummyproton), dalphat, dphit, dpt, pn, finPitheta, finProtontheta);
  
  //cout << "dalphat,dphit,dpt,pn,finPitheta,finProtontheta"  << dalphat << " " << dphit << " " << dpt << " " << pn << " " << finPitheta << " " << finProtontheta<< endl;
  plotUtils.FillHist(AnaIO::hRecdalphat,dalphat,truthEventType);
  plotUtils.FillHist(AnaIO::hRecdphit,dphit,truthEventType);
  plotUtils.FillHist(AnaIO::hRecdpt,dpt,truthEventType);
  plotUtils.FillHist(AnaIO::hRecpn,pn,truthEventType);
  //cout << "dalphat,dphit,dpt,pn,finPitheta,finProtontheta"  << dalphat_truth << " " << dphit_truth << " " << dpt_truth << " " << pn_truth << " " << finPitheta_truth << " " << finProtontheta_truth<< endl;
  if(kMC){
    AnaFunctions::getCommonTKI(targetA, targetZ, &beamFullPTruth, &(dummypi0Truth), &(dummyprotonTruth), dalphat_truth, dphit_truth, dpt_truth, pn_truth, finPitheta_truth, finProtontheta_truth);
    plotUtils.FillHist(AnaIO::hRecdalphat_truth,dalphat_truth,1);
    plotUtils.FillHist(AnaIO::hRecdphit_truth,dphit_truth,1);
    plotUtils.FillHist(AnaIO::hRecdpt_truth,dpt_truth,1);
    plotUtils.FillHist(AnaIO::hRecpn_truth,pn_truth,1);

    plotUtils.FillHist(AnaIO::hdalphat_REG,dalphat,dalphat_truth);
    plotUtils.FillHist(AnaIO::hdphit_REG,dphit,dphit_truth);
    plotUtils.FillHist(AnaIO::hdpt_REG,dpt,dpt_truth);
    plotUtils.FillHist(AnaIO::hpn_REG,pn,pn_truth);

    double dalphatRES = dalphat - dalphat_truth;
    double dphitRES = dphit - dphit_truth;
    double dptRES = dpt/dpt_truth - 1;
    double pnRES = pn/pn_truth - 1;

    plotUtils.FillHist(AnaIO::hdalphat_RES,dalphat,dalphatRES);
    plotUtils.FillHist(AnaIO::hdphit_RES,dphit,dphitRES);
    plotUtils.FillHist(AnaIO::hdpt_RES,dpt,dptRES);
    plotUtils.FillHist(AnaIO::hpn_RES,pn,pnRES);
    
  } 
}


void AnaUtils::DoTruthTKICalculation(){
  // Get number of truth FS particles
  vector<double> NParList = GetNParticles();
  AnaIO::hTruthNproton->Fill(NParList[0]);
  AnaIO::hTruthNneutron->Fill(NParList[1]);
  AnaIO::hTruthNPiZero->Fill(NParList[2]);

  //======================== TKI calculation =====================//
  const TLorentzVector beamFullP(AnaIO::true_beam_endPx, AnaIO::true_beam_endPy, AnaIO::true_beam_endPz, AnaFunctions::PionMass());
  //const TVector3 tmpbeam(AnaIO::reco_beam_true_byHits_endPx,
  //                       AnaIO::reco_beam_true_byHits_endPy,
  //                       AnaIO::reco_beam_true_byHits_endPz );
  //const TLorentzVector beamFullP(tmpbeam, AnaFunctions::PionMass());                  
  AnaIO::iniPimomentum = beamFullP.P();
  AnaIO::iniPitheta = beamFullP.Theta()*TMath::RadToDeg();

  const int targetA = 40;
  const int targetZ = 18;

  vector<TLorentzVector> vecPiP = GetFSParticlesTruth();
  AnaIO::finPimomentum = vecPiP[0].P();
  AnaIO::finProtonmomentum = vecPiP[1].P();
  AnaIO::fin2Pmom = vecPiP[2].P();

  AnaFunctions::getCommonTKI(targetA, targetZ, &beamFullP, &(vecPiP[0]), &(vecPiP[1]), 
                              AnaIO::dalphat, AnaIO::dphit, AnaIO::dpt, AnaIO::pn, AnaIO::finPitheta, AnaIO::finProtontheta);

  AnaIO::hTruthMomIniPi->Fill(AnaIO::iniPimomentum);
  AnaIO::hTruthThetaIniPi->Fill(AnaIO::iniPitheta);
  AnaIO::hTruthMomFinPi->Fill(AnaIO::finPimomentum);
  AnaIO::hTruthThetaFinPi->Fill(AnaIO::finPitheta);
  AnaIO::hTruthMomFinProton->Fill(AnaIO::finProtonmomentum);
  AnaIO::hTruthThetaFinProton->Fill(AnaIO::finProtontheta);
  AnaIO::hTruthMomFin2Proton->Fill(AnaIO::fin2Pmom);
  AnaIO::hTruthDalphat->Fill(AnaIO::dalphat);
  AnaIO::hTruthDphit->Fill(AnaIO::dphit);
  AnaIO::hTruthDpt->Fill(AnaIO::dpt);
  AnaIO::hTruthPn->Fill(AnaIO::pn);

  // 1p0n
  if(NParList[0] == 1 && NParList[1] == 0){
    AnaIO::hTruthDalphat1p0n->Fill(AnaIO::dalphat);
    AnaIO::hTruthDphit1p0n->Fill(AnaIO::dphit);
    AnaIO::hTruthDpt1p0n->Fill(AnaIO::dpt);
    AnaIO::hTruthPn1p0n->Fill(AnaIO::pn);
    AnaIO::hTruthMomIniPi1p0n->Fill(AnaIO::iniPimomentum);
    AnaIO::hTruthThetaIniPi1p0n->Fill(AnaIO::iniPitheta);
    AnaIO::hTruthMomFinPi1p0n->Fill(AnaIO::finPimomentum);
    AnaIO::hTruthThetaFinPi1p0n->Fill(AnaIO::finPitheta);
    AnaIO::hTruthMomFinProton1p0n->Fill(AnaIO::finProtonmomentum);
    AnaIO::hTruthThetaFinProton1p0n->Fill(AnaIO::finProtontheta);
    AnaIO::hTruthMomFin2Proton1p0n->Fill(AnaIO::fin2Pmom);
  }
  // Np0n
  if(NParList[0] != 1 && NParList[1] == 0){
    AnaIO::hTruthDalphatNp0n->Fill(AnaIO::dalphat);
    AnaIO::hTruthDphitNp0n->Fill(AnaIO::dphit);
    AnaIO::hTruthDptNp0n->Fill(AnaIO::dpt);
    AnaIO::hTruthPnNp0n->Fill(AnaIO::pn);
    AnaIO::hTruthMomIniPiNp0n->Fill(AnaIO::iniPimomentum);
    AnaIO::hTruthThetaIniPiNp0n->Fill(AnaIO::iniPitheta);
    AnaIO::hTruthMomFinPiNp0n->Fill(AnaIO::finPimomentum);
    AnaIO::hTruthThetaFinPiNp0n->Fill(AnaIO::finPitheta);
    AnaIO::hTruthMomFinProtonNp0n->Fill(AnaIO::finProtonmomentum);
    AnaIO::hTruthThetaFinProtonNp0n->Fill(AnaIO::finProtontheta);
    AnaIO::hTruthMomFin2ProtonNp0n->Fill(AnaIO::fin2Pmom);
  }
  // 1pMn
  if(NParList[0] == 1 && NParList[1] != 0){
    AnaIO::hTruthDalphat1pMn->Fill(AnaIO::dalphat);
    AnaIO::hTruthDphit1pMn->Fill(AnaIO::dphit);
    AnaIO::hTruthDpt1pMn->Fill(AnaIO::dpt);
    AnaIO::hTruthPn1pMn->Fill(AnaIO::pn);
    AnaIO::hTruthMomIniPi1pMn->Fill(AnaIO::iniPimomentum);
    AnaIO::hTruthThetaIniPi1pMn->Fill(AnaIO::iniPitheta);
    AnaIO::hTruthMomFinPi1pMn->Fill(AnaIO::finPimomentum);
    AnaIO::hTruthThetaFinPi1pMn->Fill(AnaIO::finPitheta);
    AnaIO::hTruthMomFinProton1pMn->Fill(AnaIO::finProtonmomentum);
    AnaIO::hTruthThetaFinProton1pMn->Fill(AnaIO::finProtontheta);
    AnaIO::hTruthMomFin2Proton1pMn->Fill(AnaIO::fin2Pmom);
  } 
  // NpMn
  if(NParList[0] != 1 && NParList[1] != 0){
    AnaIO::hTruthDalphatNpMn->Fill(AnaIO::dalphat);
    AnaIO::hTruthDphitNpMn->Fill(AnaIO::dphit);
    AnaIO::hTruthDptNpMn->Fill(AnaIO::dpt);
    AnaIO::hTruthPnNpMn->Fill(AnaIO::pn);
    AnaIO::hTruthMomIniPiNpMn->Fill(AnaIO::iniPimomentum);
    AnaIO::hTruthThetaIniPiNpMn->Fill(AnaIO::iniPitheta);
    AnaIO::hTruthMomFinPiNpMn->Fill(AnaIO::finPimomentum);
    AnaIO::hTruthThetaFinPiNpMn->Fill(AnaIO::finPitheta);
    AnaIO::hTruthMomFinProtonNpMn->Fill(AnaIO::finProtonmomentum);
    AnaIO::hTruthThetaFinProtonNpMn->Fill(AnaIO::finProtontheta);
    AnaIO::hTruthMomFin2ProtonNpMn->Fill(AnaIO::fin2Pmom);
  }
}

void AnaUtils::DoKinematicFitting(){

  cout << "DoKinematicFitting!!" << endl;
  
  // Method A (energy independent)
  //vector<double> tmp_Bin = AnaFit::DoBinCVM(AnaUtils::LdShowerEnergyTruth,AnaUtils::SlShowerEnergyTruth,AnaUtils::OpenAngleTruth,
  //              AnaUtils::LdShowerEnergyRaw,AnaUtils::SlShowerEnergyRaw,AnaUtils::OpenAngle);
  // Get the sample size
  int sampleSize = AnaUtils::LdShowerEnergyTruth.size();
  //cout << "sampleSize: " << sampleSize << endl;
  SetCVM();
  // Counters
  //int BadFitVar = 0;
  int emptyCVM = 0;
  // Method A (energy independent but now with 2 bins)
    //vector<double> CVM_Dir = AnaFit::GetCVM(AnaUtils::LdShowerEnergyTruth,AnaUtils::SlShowerEnergyTruth,AnaUtils::OpenAngleTruth,
    //            AnaUtils::LdShowerEnergyRaw,AnaUtils::SlShowerEnergyRaw,AnaUtils::OpenAngle,
    //            AnaUtils::LdShowerEnergyTruth[0],AnaUtils::SlShowerEnergyTruth[0],AnaUtils::OpenAngleTruth[0],
    //            AnaUtils::LdShowerEnergyRaw[0],AnaUtils::SlShowerEnergyRaw[0],AnaUtils::OpenAngle[0]);
  // Loop over shower vetor
  for(int ii = 0; ii < sampleSize; ii++){

    // Get the CVM matrix
    // Method A (energy independent but now with 2 bins)
    vector<double> CVM_Dir = AnaFit::GetCVM(AnaUtils::LdShowerEnergyTruth,AnaUtils::SlShowerEnergyTruth,AnaUtils::OpenAngleTruth,
                AnaUtils::LdShowerEnergyRaw,AnaUtils::SlShowerEnergyRaw,AnaUtils::OpenAngle,
                AnaUtils::LdShowerEnergyTruth[ii],AnaUtils::SlShowerEnergyTruth[ii],AnaUtils::OpenAngleTruth[ii],
                AnaUtils::LdShowerEnergyRaw[ii],AnaUtils::SlShowerEnergyRaw[ii],AnaUtils::OpenAngle[ii]);

    // Method B (energy dependent)
    vector<double> CVM_Bin = AnaFit::GetBinCVM(AnaUtils::LdShowerEnergyTruth[ii],AnaUtils::SlShowerEnergyTruth[ii],AnaUtils::OpenAngleTruth[ii],
                AnaUtils::LdShowerEnergyRaw[ii],AnaUtils::SlShowerEnergyRaw[ii],AnaUtils::OpenAngle[ii],
                AnaUtils::LdShowerEnergyTruth,AnaUtils::SlShowerEnergyTruth,AnaUtils::OpenAngleTruth,
                AnaUtils::LdShowerEnergyRaw,AnaUtils::SlShowerEnergyRaw,AnaUtils::OpenAngle);
    if(CVM_Dir.size()==0) emptyCVM++;

    //if(CVM_Bin.size()!=0){
    cout << "CVM_Dir: " << endl;
      for(auto i : CVM_Dir){
        cout << "ele: " << i << endl;
      }
      
      // Get the fitted variables
      bool GoodFit = false;
      vector<double> FittedVars = DoKF(AnaUtils::LdShowerEnergyRaw[ii],AnaUtils::SlShowerEnergyRaw[ii],AnaUtils::OpenAngle[ii], CVM_Dir, GoodFit);
      AnaIO::hKFPassRate->Fill(GoodFit);

      if(GoodFit){
/*
  vector<double> BinE1={0,0.3,2};
  vector<double> BinE2={0,0.15,1};
std::pair <int,int> bin;
for(unsigned int xx = 1; xx <= BinE1.size() - 1; xx++){
    for(unsigned int yy = 1; yy <= BinE2.size() - 1; yy++){
      if(AnaUtils::LdShowerEnergyRaw[ii] > BinE1[xx - 1]  && AnaUtils::LdShowerEnergyRaw[ii] < BinE1[xx] 
        && AnaUtils::SlShowerEnergyRaw[ii] > BinE2[yy - 1] && AnaUtils::SlShowerEnergyRaw[ii] < BinE2[yy]){
        bin = std::make_pair (xx,yy);
      }
    }
  }

if(AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii]-1 < -0.1){
  AnaIO::hPi0MassLowE1->Fill(sqrt(2*AnaUtils::LdShowerEnergyRaw[ii]*AnaUtils::SlShowerEnergyRaw[ii]*(1-cos(AnaUtils::OpenAngle[ii]))));
}
if(AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii]-1 > -0.1){
  AnaIO::hPi0MassHighE1->Fill(sqrt(2*AnaUtils::LdShowerEnergyRaw[ii]*AnaUtils::SlShowerEnergyRaw[ii]*(1-cos(AnaUtils::OpenAngle[ii]))));
}
//if(AnaUtils::CVM[bin].size()!=0){
vector<double> FittedVars = DoKF(AnaUtils::LdShowerEnergyRaw[ii],AnaUtils::SlShowerEnergyRaw[ii],AnaUtils::OpenAngle[ii],AnaUtils::CVM[bin]);
*/     
      // Fill histograms
      AnaIO::hShowerE1PreFitRes->Fill(AnaUtils::LdShowerEnergyTruth[ii],AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii]-1);
      AnaIO::hShowerE1PostFitRes->Fill(AnaUtils::LdShowerEnergyTruth[ii],FittedVars[0]/AnaUtils::LdShowerEnergyTruth[ii]-1);

      AnaIO::hShowerE2PreFitRes->Fill(AnaUtils::SlShowerEnergyTruth[ii],AnaUtils::SlShowerEnergyRaw[ii]/AnaUtils::SlShowerEnergyTruth[ii]-1);
      AnaIO::hShowerE2PostFitRes->Fill(AnaUtils::SlShowerEnergyTruth[ii],FittedVars[1]/AnaUtils::SlShowerEnergyTruth[ii]-1);

      AnaIO::hShowerOAPreFitRes->Fill(AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg(),AnaUtils::OpenAngle[ii]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg());
      AnaIO::hShowerOAPostFitRes->Fill(AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg(),FittedVars[2]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg());

      AnaIO::hShowerE1Compare->Fill(AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii]-1);
      AnaIO::hShowerE1ComparePost->Fill(FittedVars[0]/AnaUtils::LdShowerEnergyTruth[ii]-1);

      AnaIO::hShowerE2Compare->Fill(AnaUtils::SlShowerEnergyRaw[ii]/AnaUtils::SlShowerEnergyTruth[ii]-1);
      AnaIO::hShowerE2ComparePost->Fill(FittedVars[1]/AnaUtils::SlShowerEnergyTruth[ii]-1);

      AnaIO::hShowerOACompare->Fill(AnaUtils::OpenAngle[ii]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg());
      AnaIO::hShowerOAComparePost->Fill(FittedVars[2]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg());

      AnaIO::hPi0MassCompare->Fill(sqrt(2*AnaUtils::LdShowerEnergyRaw[ii]*AnaUtils::SlShowerEnergyRaw[ii]*(1-cos(AnaUtils::OpenAngle[ii]))));
      AnaIO::hPi0MassComparePost->Fill(sqrt(2*FittedVars[0]*FittedVars[1]*(1-cos(FittedVars[2]))));

      AnaIO::hShowerE1Compare_REG->Fill(AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii]-1,FittedVars[0]/AnaUtils::LdShowerEnergyTruth[ii]-1);
      AnaIO::hShowerE2Compare_REG->Fill(AnaUtils::SlShowerEnergyRaw[ii]/AnaUtils::SlShowerEnergyTruth[ii]-1,FittedVars[1]/AnaUtils::SlShowerEnergyTruth[ii]-1);
      AnaIO::hShowerOACompare_REG->Fill(AnaUtils::OpenAngle[ii]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg(), FittedVars[2]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg());
      // Print out some info
      cout << "LD shower E (truth): " << AnaUtils::LdShowerEnergyTruth[ii] << endl;
      cout << "LD shower E: " << AnaUtils::LdShowerEnergyRaw[ii] << endl;
      cout << "LD shower E (fitted): " << FittedVars[0] << endl;

      cout << "\nSL shower E (truth): " << AnaUtils::SlShowerEnergyTruth[ii] << endl;
      cout << "SL shower E: " << AnaUtils::SlShowerEnergyRaw[ii] << endl;
      cout << "SL shower E (fitted): " << FittedVars[1] << endl;

      cout << "\nOpenAngle (truth): " << AnaUtils::OpenAngleTruth[ii] << endl;
      cout << "OpenAngle: " << AnaUtils::OpenAngle[ii] << endl;
      cout << "OpenAngle (fitted): " << FittedVars[2] << endl;

      cout << "\nMass (truth): " << sqrt(2*AnaUtils::LdShowerEnergyTruth[ii]*AnaUtils::SlShowerEnergyTruth[ii]*(1-cos(AnaUtils::OpenAngleTruth[ii]))) << endl;
      cout << "Mass: " << sqrt(2*AnaUtils::LdShowerEnergyRaw[ii]*AnaUtils::SlShowerEnergyRaw[ii]*(1-cos(AnaUtils::OpenAngle[ii]))) << endl;
      cout << "Mass (fitted): " << sqrt(2*FittedVars[0]*FittedVars[1]*(1-cos(FittedVars[2]))) << endl;

      const TVector3 ldShower3Vect = AnaUtils::LdShowerDir[ii].Unit()*FittedVars[0];
      const TVector3 slShower3Vect = AnaUtils::SlShowerDir[ii].Unit()*FittedVars[1];

      const TVector3 ldShower3VectPre = AnaUtils::LdShowerDir[ii].Unit()*AnaUtils::LdShowerEnergyRaw[ii];
      const TVector3 slShower3VectPre = AnaUtils::SlShowerDir[ii].Unit()*AnaUtils::SlShowerEnergyRaw[ii];

      const TVector3 ldShower3VectTruth = AnaUtils::LdShowerDirTruth[ii].Unit()*AnaUtils::LdShowerEnergyTruth[ii];
      const TVector3 slShower3VectTruth = AnaUtils::SlShowerDirTruth[ii].Unit()*AnaUtils::SlShowerEnergyTruth[ii];

      TLorentzVector ldShowerLTVectLab;
      ldShowerLTVectLab.SetVectM(ldShower3Vect, 0 );

      TLorentzVector slShowerLTVectLab;
      slShowerLTVectLab.SetVectM(slShower3Vect, 0 );

      TLorentzVector ldShowerLTVectLabPre;
      ldShowerLTVectLabPre.SetVectM(ldShower3VectPre, 0 );

      TLorentzVector slShowerLTVectLabPre;
      slShowerLTVectLabPre.SetVectM(slShower3VectPre, 0 );
      
      TLorentzVector ldShowerLTVectLab_Truth;
      ldShowerLTVectLab_Truth.SetVectM(ldShower3VectTruth, 0 );

      TLorentzVector slShowerLTVectLab_Truth;
      slShowerLTVectLab_Truth.SetVectM(slShower3VectTruth, 0 );

      TLorentzVector RecPi0 = ldShowerLTVectLab + slShowerLTVectLab;

      TLorentzVector RecPi0Pre = ldShowerLTVectLabPre + slShowerLTVectLabPre;

      TLorentzVector TruthPi0 = ldShowerLTVectLab_Truth + slShowerLTVectLab_Truth;

      
      double LDRes = AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii] - 1;
      double SLRes = AnaUtils::SlShowerEnergyRaw[ii]/AnaUtils::SlShowerEnergyTruth[ii] - 1;
      double OARes = (AnaUtils::OpenAngle[ii]/AnaUtils::OpenAngleTruth[ii]) -1;

      AnaIO::hLDShowerE->Fill(LDRes);
      AnaIO::hSLShowerE->Fill(SLRes);
      AnaIO::hOAShower->Fill(OARes);

      cout << "LD: " << AnaUtils::LdShowerEnergyRaw[ii] - AnaUtils::LdShowerEnergyTruth[ii] << endl;
      cout << "SL: " << AnaUtils::SlShowerEnergyRaw[ii] - AnaUtils::SlShowerEnergyTruth[ii] << endl;
      cout << "OA: " << (AnaUtils::OpenAngle[ii] - AnaUtils::OpenAngleTruth[ii])*TMath::RadToDeg() << endl;

      // Get Variables
      double E1 = AnaUtils::LdShowerEnergyRaw[ii];
      double E2 = AnaUtils::SlShowerEnergyRaw[ii];
      double OA = AnaUtils::OpenAngle[ii];

      double E1_Fit = FittedVars[0];
      double E2_Fit = FittedVars[1];
      double OA_Fit = FittedVars[2];

      double E1_Truth = AnaUtils::LdShowerEnergyTruth[ii];
      double E2_Truth = AnaUtils::SlShowerEnergyTruth[ii];
      double OA_Truth = AnaUtils::OpenAngleTruth[ii];

      // Method A (Pre Fit)
      double pion_E_pre_m1 = E1 + E2;
      // Method A (Post Fit)
      double pion_E_post_m1 = E1_Fit + E2_Fit;

      // Method B (Pre Fit)
      double pion_E_pre_m2 = E1 + 0.134977*0.134977/(2*E1*(1-cos(OA)));
      // Method B (Post Fit)
      double pion_E_post_m2 = E1_Fit + 0.134977*0.134977/(2*E1_Fit*(1-cos(OA_Fit)));

      // Method C (Pre Fit)
      double alpha = abs(E1-E2)/(E1+E2);
      double p_pi0 = 0.134977*sqrt(2/((1-alpha*alpha)*(1-cos(OA))));
      double e_pi0 = p_pi0;//sqrt(0.134977*0.134977+p_pi0*p_pi0);
      double pion_E_pre_m3 = e_pi0;
      
      // Method C (Post Fit)
      double alpha_post = abs(E1_Fit-E2_Fit)/(E1_Fit+E2_Fit);
      double p_pi0_post = 0.134977*sqrt(2/((1-alpha_post*alpha_post)*(1-cos(OA_Fit))));
      double e_pi0_post = p_pi0_post;//sqrt(0.134977*0.134977+p_pi0_post*p_pi0_post);
      double pion_E_post_m3 = e_pi0_post;

      double pion_E_Truth = E1_Truth + 0.134977*0.134977/(2*E1_Truth*(1-cos(OA_Truth)));
      double pion_E_Truth_dir = E1_Truth + E2_Truth;

      double alpha_Truth = abs(E1_Truth-E2_Truth)/(E1_Truth+E2_Truth);
      double p_pi0_Truth = 0.134977*sqrt(2/((1-alpha_Truth*alpha_Truth)*(1-cos(OA_Truth))));
      double e_pi0_Truth = p_pi0_Truth;//sqrt(0.134977*0.134977+p_pi0_Truth*p_pi0_Truth);
      cout << "pion_E_Truth: " << pion_E_Truth << " pion_E_Truth_dir: " << pion_E_Truth_dir << " e_pi0_Truth: " << e_pi0_Truth << endl;

/*
      double pion_E = E1 + 0.134977*0.134977/(2*E1*(1-cos(OA)));
      double pion_P = sqrt(pion_E*pion_E - 0.134977*0.134977);

      double pion_E_Pre = E1 + 0.134977*0.134977/(2*E1*(1-cos(AnaUtils::OpenAngle[ii])));
      double pion_P_Pre = sqrt(pion_E_Pre*pion_E_Pre - 0.134977*0.134977);

      double pion_E_Truth = E1_Truth + 0.134977*0.134977/(2*E1_Truth*(1-cos(OA_truth)));

      double pion_E_dir_Truth = E1_Truth + E2_Truth;

      double pion_P_Truth = sqrt(pion_E_Truth*pion_E_Truth - 0.134977*0.134977);

      // Paper Equation
      double alpha = abs(E1-E2)/(E1+E2);
      double p_pi0 = 0.134977*sqrt(2/((1-alpha*alpha)*(1-cos(AnaUtils::OpenAngle[ii]))));
      double e_pi0 = sqrt(0.134977*0.134977+p_pi0*p_pi0);

      double alpha_post = abs(FittedVars[0]-FittedVars[1])/(FittedVars[0]+FittedVars[1]);
      double p_pi0_post = 0.134977*sqrt(2/((1-alpha_post*alpha_post)*(1-cos(FittedVars[2]))));
      double e_pi0_post = sqrt(0.134977*0.134977+p_pi0_post*p_pi0_post);
*/

      AnaIO::hPrePi0MomNorm->Fill(pion_E_pre_m1/pion_E_Truth - 1);
      AnaIO::hPrePi0MomEOA->Fill(pion_E_pre_m2/pion_E_Truth - 1);
      AnaIO::hPrePi0MomAsym->Fill(pion_E_pre_m3/pion_E_Truth - 1);

      AnaIO::hPi0MomNorm->Fill(pion_E_post_m1/pion_E_Truth - 1);
      AnaIO::hPi0MomEOA->Fill(pion_E_post_m2/pion_E_Truth - 1);
      AnaIO::hPi0MomAsym->Fill(pion_E_post_m3/pion_E_Truth - 1);

/*
      cout << "RecPi0Pre.P(): " << RecPi0Pre.E() << "RecPi0.E() : " << RecPi0.E() << "TruthPi0.E() : " << TruthPi0.E() << endl;
      cout << "pion_E_Truth: " << pion_P_Truth << endl;
      cout << "pion_E_Truth.E(): " << TruthPi0.E() << endl;
      cout << "pion_E_dir_Truth: " << pion_E_dir_Truth << endl;
      cout << "pion_E: " << pion_E << endl;
      cout << "E2: " << E2_Truth << " " << 0.134977*0.134977/(2*E1_Truth*(1-cos(OA_truth))) << endl;
      cout << "pion_E_Pre: " << pion_P_Pre << endl;
      cout << "pion_P: " << pion_P << endl;
      cout << "TruthPi0Pre.M(): " << TruthPi0.M() << endl;
*/

      AnaIO::hPi0MomCompareStandard->Fill(pion_E_pre_m1/pion_E_Truth - 1);
      AnaIO::hPi0MomComparePreStandard->Fill(pion_E_pre_m1/pion_E_Truth - 1);
      AnaIO::hPi0MomComparePostStandard->Fill(pion_E_post_m1/pion_E_Truth - 1);

      AnaIO::hPi0MomCompare->Fill(pion_E_pre_m2/pion_E_Truth - 1);
      AnaIO::hPi0MomComparePre->Fill(pion_E_pre_m2/pion_E_Truth - 1);
      AnaIO::hPi0MomComparePost->Fill(pion_E_post_m2/pion_E_Truth - 1);

      AnaIO::hPi0MomPreFitRes->Fill(pion_E_pre_m2,pion_E_Truth);
      AnaIO::hPi0MomPostFitRes->Fill(pion_E_post_m2,pion_E_Truth);

      AnaIO::hPi0MomCompare_REG->Fill(pion_E_pre_m2/pion_E_Truth - 1,pion_E_post_m2/pion_E_Truth - 1);

      AnaIO::hPi0MomComparePaper->Fill(pion_E_pre_m3/pion_E_Truth - 1);
      AnaIO::hPi0MomComparePrePaper->Fill(pion_E_pre_m3/pion_E_Truth - 1);
      AnaIO::hPi0MomComparePostPaper->Fill(pion_E_post_m3/pion_E_Truth - 1);

      


      AnaIO::hPi0MomPreFitResPaper->Fill(pion_E_pre_m3,pion_E_Truth);
      AnaIO::hPi0MomPostFitResPaper->Fill(pion_E_post_m3,pion_E_Truth);

      AnaIO::hPi0MomCompare_REGPaper->Fill(pion_E_pre_m3/pion_E_Truth - 1,pion_E_post_m3/pion_E_Truth - 1);
      AnaIO::hPi0MomCompare_REGStand->Fill(pion_E_pre_m1/pion_E_Truth - 1,pion_E_post_m1/pion_E_Truth - 1);
      

      const TLorentzVector recMomRefBeam = GetPi0MomentumRefBeam(RecPi0,TruthPi0,false);
      const TLorentzVector truthMomRefBeam = GetPi0MomentumRefBeam(RecPi0,TruthPi0,true);

      AnaIO::hPi0ThetaFitRes->Fill(truthMomRefBeam.Theta()*TMath::RadToDeg(),recMomRefBeam.Theta()*TMath::RadToDeg() - truthMomRefBeam.Theta()*TMath::RadToDeg());
      AnaIO::hLDShowerThetaFitRes->Fill(ldShowerLTVectLab_Truth.Theta()*TMath::RadToDeg(), ldShowerLTVectLab.Theta()*TMath::RadToDeg() - ldShowerLTVectLab_Truth.Theta()*TMath::RadToDeg());
      AnaIO::hPi0ThetaResFit->Fill(TruthPi0.Theta()*TMath::RadToDeg(), RecPi0.Theta()*TMath::RadToDeg() - TruthPi0.Theta()*TMath::RadToDeg());
      AnaIO::hPi0PhiResFit->Fill(TruthPi0.Phi()*TMath::RadToDeg(), RecPi0.Phi()*TMath::RadToDeg() - TruthPi0.Phi()*TMath::RadToDeg());
      AnaIO::hPi0MomRes->Fill(pion_E_Truth,pion_E_post_m2/pion_E_Truth - 1);
      AnaIO::hPi0MomPreRes->Fill(pion_E_Truth,pion_E_pre_m2/pion_E_Truth - 1);
    //}
      }

  }
  cout << "emptyCVM: " << emptyCVM << endl;
  
}


void AnaUtils::SetCVM(){
  // Dimension of the CVM
  int dim = 2;
  // Set empty vector to CVM
  for(int xx = 1; xx <= dim; xx++){
    for(int yy = 1; yy <= dim; yy++){
      // Create an empty vector
      vector<double> Vect;
      // Make the key
      std::pair <int,int> bin = std::make_pair (xx,yy);
      if(xx==1 && yy==1){
        Vect.push_back(0.0081142);
        Vect.push_back(-8.38321e-06);
        Vect.push_back(-0.00331411);
        Vect.push_back(-8.38321e-06);
        Vect.push_back(0.00906629);
        Vect.push_back(-0.00466705);
        Vect.push_back(-0.00331411);
        Vect.push_back(-0.00466705);
        Vect.push_back(0.118713);
      }
      else if(xx==1 && yy==2){
        Vect.push_back(0.00810685);
        Vect.push_back(-7.04183e-05);
        Vect.push_back(-0.00585613);
        Vect.push_back(-7.04183e-05);
        Vect.push_back(0.00684988);
        Vect.push_back(-0.000956773);
        Vect.push_back(-0.00585613);
        Vect.push_back(-0.000956773);
        Vect.push_back(0.110837);
      }
      else if(xx==2 && yy==1){
        Vect.push_back(0.00734139);
        Vect.push_back(-0.00293472);
        Vect.push_back(0.0015594);
        Vect.push_back(-0.00293472);
        Vect.push_back(0.00851379);
        Vect.push_back(-0.00856854);
        Vect.push_back(0.0015594);
        Vect.push_back(-0.00856854);
        Vect.push_back(0.147113);
      }
      else if(xx==2 && yy==2){
        Vect.push_back(0.00605665);
        Vect.push_back(-0.00139198);
        Vect.push_back(-0.00107769);
        Vect.push_back(-0.00139198);
        Vect.push_back(0.00523178);
        Vect.push_back(0.000339517);
        Vect.push_back(-0.00107769);
        Vect.push_back(0.000339517);
        Vect.push_back(0.127052);
      }
      CVM[bin] = Vect;
    }
  }
}

/*
void AnaUtils::SetCVM(){
  // Dimension of the CVM
  int dim = 3;
  // Set empty vector to CVM
  for(int xx = 1; xx <= dim; xx++){
    for(int yy = 1; yy <= dim; yy++){
      // Create an empty vector
      vector<double> Vect;
      // Make the key
      std::pair <int,int> bin = std::make_pair (xx,yy);
      if(xx==1 && yy==1){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00206809);
        Vect.push_back(0.00244881);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00244881);
        Vect.push_back(0.24885);
      }
      else if(xx==1 && yy==2){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00110825);
        Vect.push_back(-0.000985235);
        Vect.push_back(0.00000001);
        Vect.push_back(-0.000985235);
        Vect.push_back(0.200011);
      }
      else if(xx==1 && yy==3){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00366147);
        Vect.push_back(-0.0201989);
        Vect.push_back(0.00000001);
        Vect.push_back(-0.0201989);
        Vect.push_back(0.345658);
      }
      else if(xx==2 && yy==1){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00245149);
        Vect.push_back(0.00237991);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00237991);
        Vect.push_back(0.165784);
      }
      else if(xx==2 && yy==2){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00135179);
        Vect.push_back(0.000211054);
        Vect.push_back(0.00000001);
        Vect.push_back(0.000211054);
        Vect.push_back(0.220865);
      }
      else if(xx==2 && yy==3){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.0062439);
        Vect.push_back(-0.0107125);
        Vect.push_back(0.00000001);
        Vect.push_back(-0.0107125);
        Vect.push_back(0.17041);
      }
      else if(xx==3 && yy==1){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00426994);
        Vect.push_back(0.0122522);
        Vect.push_back(0.00000001);
        Vect.push_back(0.0122522);
        Vect.push_back(0.23304);
      }
      else if(xx==3 && yy==2){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.0018513);
        Vect.push_back(0.00174401);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00174401);
        Vect.push_back(0.199405);
      }
      else if(xx==3 && yy==3){
        Vect.push_back(1.0);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00000001);
        Vect.push_back(0.00994019);
        Vect.push_back(-0.0197883);
        Vect.push_back(0.00000001);
        Vect.push_back(-0.0197883);
        Vect.push_back(0.291611);
      }
      CVM[bin] = Vect;
    }
  }
}
*/

void AnaUtils::KF(const TLorentzVector &ldShower, const TLorentzVector &slShower, vector<double> &FittedVars){

  const double openingAngleRad = ldShower.Angle(slShower.Vect());
  vector<double> BinE1={0,0.3,5};
  vector<double> BinE2={0,0.15,5};
  std::pair <int,int> bin;
  
  for(unsigned int xx = 1; xx <= BinE1.size() - 1; xx++){
    for(unsigned int yy = 1; yy <= BinE2.size() - 1; yy++){
      if(ldShower.E() > BinE1[xx - 1]  && ldShower.E() < BinE1[xx] 
        && slShower.E() > BinE2[yy - 1] && slShower.E() < BinE2[yy]){
        bin = std::make_pair (xx,yy);
      }
    }
  }
  //cout << "ldShower.E(): " << ldShower.E() << endl;
  //cout << "slShower.E(): " << slShower.E() << endl;

  cout << "bin CVM: "<< endl;
  for(auto i : AnaUtils::CVM[bin]){
    cout << "ele: " << i << endl;
  }

  SetCVM();
  bool GoodFit = false;
  FittedVars = DoKF(ldShower.E(),slShower.E(),openingAngleRad,AnaUtils::CVM[bin],GoodFit);

  AnaIO::hKFPassRate->Fill(GoodFit);

}