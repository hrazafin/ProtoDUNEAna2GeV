#include "../include/AnaUtils.h"
#include "AnaFit.cxx"

int AnaUtils::GetBeamParticleType(const int pdg)
{
  int type = -999;
  // Info about this beam particle (truth matched)
  const int origin = AnaIO::reco_beam_true_byHits_origin; // 0: kUnknown 1: kBeamNeutrino 2: kCosmicRay 3: kSuperNovaNeutrino 4: kSingleParticle (single particles thrown at the detector) 
  const bool matched = AnaIO::reco_beam_true_byHits_matched; 
  const TString process = (*AnaIO::reco_beam_true_byHits_process); 

  if (process == "primary" && matched && origin == 4 && pdg == 211) {
    type = gkBeamPiPlus;
  }
  else if (process == "primary" && matched && origin == 4 && pdg == -13) {
    type = gkBeamMuon;
  }
  // Secondary particles misidentified as beam particle
  else if(pdg==2212){
    type = gkMisIDProton;
  }
  else if(abs(pdg)==211){
    type = gkMisIDPion;
  }
  else if(pdg==-13){
    type = gkMisIDMuon;
  }
  else if(pdg==11 || pdg==-11 || pdg==22){
    type = gkMisIDElectronGamma;
  }
  else if (origin == 2) {
    type = gkCosmicMuon;
  }
  else{
    type = gkBeamOthers;
  }
  return type;
}

int AnaUtils::GetChannelType(const TString beamEndprocess, const int beamPDG, const int npi0Daughters, const int npiplusDaughters, const int npiminusDaughters)
{
  int type = -999;
  // Pi+ beam event
  if(beamPDG == 211){
    // End process is inelastic 
    if(beamEndprocess == "pi+Inelastic"){ // This is the primary pion beam event without elastic interaction 
      // Inelastic + Charge Exchange + Absorption
/*
      // Find daughter pion with momentum above detection threshold
      bool has_daughter_pion_above_threshold = false;
      const int size = (*AnaIO::true_beam_daughter_PDG).size();
      for (int ii = 0; ii < size; ++ii) {
        if (abs((*AnaIO::true_beam_daughter_PDG)[ii]) == 211 && (*AnaIO::true_beam_daughter_startP)[ii] > 0.15) {
          has_daughter_pion_above_threshold = true;
          break;
        }
      }
*/
      // Inelastic channel (have daughter charged pions above threshold)
      if(npiplusDaughters + npiminusDaughters > 0 /*&& has_daughter_pion_above_threshold == true*/){
        type = gkInelastic;
      }
      // Charge Exchange + Absorption channels (No daughter charged pions)
      else{
        /*bool has_daughter_pi0_above_threshold = false;
        const int size = (*AnaIO::true_beam_daughter_PDG).size();
        for (int ii = 0; ii < size; ++ii) {
          double P = (*AnaIO::true_beam_daughter_startP)[ii];
          double E = sqrt(P*P+AnaFunctions::PiZeroMass()*AnaFunctions::PiZeroMass());
          if (abs((*AnaIO::true_beam_daughter_PDG)[ii]) == 111 && E > 0.2) {
            has_daughter_pi0_above_threshold = true;
            break;
          }
        }*/
        // Charge Exchange channel (have daughter neutral pions)
        if(npi0Daughters > 0){
          type = gkChargeExchange;
        }
        // Absorption channel (no pions)
        else{
          type = gkAbsorption;
        }
      }
    }
    // End process is decay 
    else if (beamEndprocess == "Decay"){
      type = gkPionDecays;
    }
    // End process is others 
    else{
      type = gkOtherChannels;
    }
  }
  else if(abs(beamPDG) == 13){ 
    type = gkBeamMuons;
  }
  // Non pi+ beam event
  else {
    type = gkBackground;
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

int AnaUtils::GetFillXSEventType()
{
  int filleventtype = -999;
  if(AnaIO::Signal){
    filleventtype = gkXSSignal;
  }
  else if(AnaIO::true_beam_PDG==211){
    if(AnaIO::true_daughter_nPiPlus == 0 && AnaIO::true_daughter_nPiMinus == 0 && AnaIO::true_daughter_nPi0 == 0) filleventtype = gkXSEvtBkgAbs;
    else if((AnaIO::true_daughter_nPiPlus != 0 || AnaIO::true_daughter_nPiMinus != 0) && AnaIO::true_daughter_nPi0 == 0) filleventtype = gkXSEvtBkgInel;
    else if(AnaIO::true_daughter_nPi0 == 1) filleventtype = gkXSEvtBkgSinglePi0;
    else /*if(AnaIO::true_daughter_nPi0 > 1)*/ filleventtype = gkXSEvtBkgMultiPi0;
    //else filleventtype = gkXSEvtBkgOthers;
  }
  else{
    filleventtype = gkXSBmBkg;
  }

  return filleventtype;
}

int AnaUtils::GetFillTKIEventType()
{
  cout << "TKI nProton: " << nProton << " nNeutron: " << nNeutron << " nPiPlus: " << nPiPlus 
  << " nPiZero: " << nPiZero << " nGamma: " << nGamma << " nParticleBkg: " << nParticleBkg << endl;
  cout << endl;

  int filleventTKItype = -999;

  if(AnaIO::Signal){
    if(nProton == 1 && nNeutron == 0){
      filleventTKItype = gk1p0n;
    }
    else if(nProton == 1 && nNeutron > 0 ){
      filleventTKItype = gk1pMn;
    }
    //else if(nProton > 0 && nNeutron == 0){
    //  filleventTKItype = gkNp0n;
    //}
    else{
      cout << "gkNpMn Cat: nProton " << nProton << " nNeutron: " << nNeutron << endl;
      filleventTKItype = gkNpMn;
    }
  }
  else if(AnaIO::true_beam_PDG==211){
    filleventTKItype = gkEvtTKIBkg;
  }
  else{
    filleventTKItype = gkBmTKIBkg;
  }

  return filleventTKItype;
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

  // We only want need beam pion daughters
  if(AnaIO::true_beam_PDG == 211){
    // Get final state particles vector in this event
    vector<TLorentzVector> vecFSParticle = GetFSParticlesTruth();
    
    // Get the FS particles momentum
    double LeadingPiZeroP = vecFSParticle[0].P();
    double LeadingProtonP = vecFSParticle[1].P();
    double SubLeadingProtonP = vecFSParticle[2].P();

    // Check event topology (TKI event)
    // Initial pion beam and at least one proton and one pizero, no other mesons in final state (but not consider number of neutrons)
    /*cout << "Sig nProton: " << nProton << " nNeutron: " << nNeutron << " nPiPlus: " << nPiPlus 
    << " nPiZero: " << nPiZero << " nGamma: " << nGamma << " nParticleBkg: " << nParticleBkg << endl;
    cout << endl;
    cout << "ldP: " << LeadingProtonP << " slP: " << SubLeadingProtonP << " ldpi0P: " << LeadingPiZeroP << endl;
    cout << "AnaIO::true_beam_PDG: " << AnaIO::true_beam_PDG << endl;
    */
    // Check the pi0 daughters 
    if(nPiZero > 0){
      vector<TLorentzVector> vecPi0FSParticle = GetFSPiZeroDecayDaughterTruth();
      //double LeadingShowerE = vecPi0FSParticle[0].E();
      //double SubLeadingShowerE = vecPi0FSParticle[1].E();
      //cout << "Sig Showers: " << nPiZeroGamma  << endl;
      //cout << "ld Shower E: " << LeadingShowerE << " sl Shower E: " << SubLeadingShowerE << endl;
    }

    if( AnaIO::true_beam_PDG==211 && IsSignal(nProton,nPiZero,nPiPlus,nPiMinus,nParticleBkg,nPionAboveThreshold) == true){
      // "LeadingProtonP: " << LeadingProtonP << " SubLeadingProtonP: " << SubLeadingProtonP << " LeadingPiZeroP: " << LeadingPiZeroP << endl;
      // Proton momentum selection (below 0.45 GeV/c is not detectbale)
      if(LeadingProtonP < 1 && LeadingProtonP > 0.45 && SubLeadingProtonP < 0.45){
        // No restrictions on pizero momentum
        if(LeadingPiZeroP > 0){ /*AnaIO::Signal = true; cout << "This is a good TKI event" << endl;*/} //cout << "This is a good TKI event" << endl; //AnaIO::Signal = true; 
      }
    }
    
    // Check Pi0 inclusive events
    if(AnaIO::true_beam_PDG == 211 && (*AnaIO::true_beam_endProcess) == "pi+Inelastic" && IsSignal(nProton,nPiZero,nPiPlus,nPiMinus,nParticleBkg,nPionAboveThreshold) == true){
      //double LeadingPiZeroE = sqrt(LeadingPiZeroP*LeadingPiZeroP+AnaFunctions::PiZeroMass()*AnaFunctions::PiZeroMass());
      if(LeadingPiZeroP > 0) AnaIO::Signal = true; 
    }
    // Check TKI event 
    if(nProton > 0 && nPiZero > 0 && nPiPlus == 0 && nParticleBkg == 0){
      if(LeadingProtonP < 1.0 && LeadingProtonP > 0.45 && SubLeadingProtonP < 0.45) DoTruthTKICalculation();
    }
  }
}

vector<TLorentzVector> AnaUtils::GetFSParticlesTruth(bool kFill)
{
  // Get true beam daughter information
  const vector<int> * pdg = AnaIO::true_beam_daughter_PDG;
  const vector<double> * px = AnaIO::true_beam_daughter_startPx; // GeV/c
  const vector<double> * py = AnaIO::true_beam_daughter_startPy; // GeV/c
  const vector<double> * pz = AnaIO::true_beam_daughter_startPz; // GeV/c
  //const vector<string> * process = AnaIO::true_beam_daughter_Process;
  //const vector<string> * EndProcess = AnaIO::true_beam_daughter_endProcess;

  // Class member variables (beam truth daughter particles counter)
  nProton = 0;
  nNeutron = 0;
  nPiPlus = 0;
  nPiMinus = 0;
  nPiZero = 0;
  nGamma = 0;
  nParticleBkg = 0;
  nPionAboveThreshold = 0;
  // Vectors to save particle info
  TLorentzVector pPiZero, pProton, pSecondaryProton, pPiPlus;

  // Get the size of final state particles
  const int np = pdg->size();
  // Fill the hTruthFSParticleNumber hitogram
  if(kFill) AnaIO::hTruthFSParticleNumber->Fill(np);

  double Protonmom[np];
  double PiZeromom[np];
  double PiPlusmom[np];

  double Gammamom[np];
  vector<TVector3> bufferProtonmom;
  vector<TVector3> bufferPiZeromom;
  vector<TVector3> bufferPiPlusmom;

  bufferType.clear();
  // Now loop over FS particles
  for(int ii=0; ii<np; ii++){
    // Get the FS particle type
    const int itype = GetParticleType((*pdg)[ii]);
    // Get the FS particle 3-momentum
    const TVector3 tmpp( (*px)[ii], (*py)[ii], (*pz)[ii] );
    // Fill the FS particle type
    if(kFill) AnaIO::hTruthFSParticleType->Fill(itype); 

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
      // Save particle's momentum magnitude
      PiPlusmom[nPiPlus] = tmpp.Mag();
      // Save momentum vector
      bufferPiPlusmom.push_back(tmpp);

      nPiPlus++;
      if(tmpp.Mag() > 0.10) nPionAboveThreshold++;
    }
    // PiMinus
    else if(itype == gkPiMinus){
      if(kFill) AnaIO::hTruthPiMinusP->Fill(tmpp.Mag());
      nPiMinus++;
      if(tmpp.Mag() > 0.10) nPionAboveThreshold++;
    }
    // Neutron
    else if(itype == gkNeutron){
      nNeutron++;
    }
    // PiMinus and Kaon
    else if(itype==gkKaon){
      nParticleBkg++;
    }
    // Others
    else{
      //cout << "itype: " << itype << endl;
      //cout << "pdg: " << (*pdg)[ii] << endl;
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
    if(kFill) AnaIO::hTruthLeadingProtonP->Fill(Protonmom[leadingProtonID]);
    // Fill all proton momentum
    for(unsigned int ii = 0; ii < bufferProtonmom.size(); ii++){
      if(kFill) AnaIO::hTruthProtonP->Fill(bufferProtonmom[ii].Mag());
    }
  }
  // At least two proton is found
  if(nProton>1){
    // Save info to subleading proton TLorentzVector
    pSecondaryProton.SetVectM(bufferProtonmom[subldProtonID], AnaFunctions::ProtonMass());
    // Fill subleading proton momentum
    if(kFill) AnaIO::hTruthSubLeadingProtonP->Fill(Protonmom[subldProtonID]);
  }

  //======================= PiPlus =======================
  int leadingPiPlusID = 0, subldPiPlusID = -999;
  // Sort piplus when more than one piplus are found
  if(nPiPlus>1){
    int PiPlussortid[nPiPlus];
    // Sort index according to it's momentum
    TMath::Sort(nPiPlus, PiPlusmom, PiPlussortid);
    // Save sorted index
    leadingPiPlusID = PiPlussortid[0];
    subldPiPlusID = PiPlussortid[1];
  }
  // At least one piplus is found
  if(nPiPlus>0){
    // Save info to leading piplus TLorentzVector
    pPiPlus.SetVectM(bufferPiPlusmom[leadingPiPlusID], AnaFunctions::PionMass());
    // Fill leading piplus momentum
    if(kFill) AnaIO::hTruthLeadingPiPlusP->Fill(PiPlusmom[leadingPiPlusID]);
    // Fill all piplus momentum
    for(unsigned int ii = 0; ii < bufferPiPlusmom.size(); ii++){
      if(kFill) AnaIO::hTruthPiPlusP->Fill(bufferPiPlusmom[ii].Mag());
    }
  }
  // At least two piplus is found
  if(nPiPlus>1){
    // Fill subleading piplus momentum
    if(kFill) AnaIO::hTruthSubLeadingPiPlusP->Fill(PiPlusmom[subldPiPlusID]);
  }

  //======================== PiZero ========================
  int leadingPiZeroID = 0, subldPiZeroID = -999;
  if(nPiZero>1){
    // Fill the FS pi0 number (at least two pi0)
    if(kFill) AnaIO::hTruthFSMultiPi0->Fill(nPiZero);
    int PiZerosortid[nPiZero];
    // Sort index according to it's momentum
    TMath::Sort(nPiZero, PiZeromom, PiZerosortid);
    // Save sorted index
    leadingPiZeroID = PiZerosortid[0];
    subldPiZeroID = PiZerosortid[1];
  }
  if(nPiZero>0){
    // Fill histogram for FS pi0 number
    if(kFill) AnaIO::hTruthFSPi0Number->Fill(nPiZero);
    // Save info to pi0 TLorentzVector
    pPiZero.SetVectM(bufferPiZeromom[leadingPiZeroID], AnaFunctions::PiZeroMass());
    if(kFill) AnaIO::hTruthLeadingPiZeroP->Fill(PiZeromom[leadingPiZeroID]);
    //if(kFill && pPiZero.E() > 0.2) AnaIO::hTruthLeadingPiZeroE->Fill(pPiZero.E());
    if(kFill) AnaIO::hTruthLeadingPiZeroE->Fill(pPiZero.E());

  }
  if(nPiZero>1){
    if(kFill) AnaIO::hTruthSubLeadingPiZeroP->Fill(PiZeromom[subldPiZeroID]);
  }

  //======================== Gamma ========================
  int leadingGammaID = 0;
  if(nGamma>1){
    int Gammasortid[nGamma];
    TMath::Sort(nGamma, Gammamom, Gammasortid);
    leadingGammaID = Gammasortid[0];
  }
  // Fill in MeV!
  if(nGamma > 0 && kFill) AnaIO::hTruthGammaMaxE->Fill(Gammamom[leadingGammaID]*1000);
  //if(nGamma > 0 kFill) AnaIO::hTruthLeadingPi0GammaP->Fill(Gammamom[leadingGammaID]);
  // Fill vector of FS particles 
  vector<TLorentzVector> vec;
  vec.push_back(pPiZero);
  vec.push_back(pProton);
  vec.push_back(pSecondaryProton);
  return vec;
}

vector<TLorentzVector> AnaUtils::GetFSPiZeroDecayDaughterTruth(const bool kFill)
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
  if(kFill) AnaIO::hTruthPi0DecayParticleNumber->Fill(np);

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
    if(kFill)AnaIO::hTruthLeadingPi0GammaP->Fill(PiZeroGammamom[leadingGammaID]);
    if(kFill)AnaIO::hTruthLeadingPi0GammaPBin1->Fill(PiZeroGammamom[leadingGammaID]);
    if(kFill)AnaIO::hTruthSubLeadingPi0GammaP->Fill(PiZeroGammamom[subldGammaID]);
    // Opening angle
    double OA = (bufferPiZeroGammamom[0].Angle(bufferPiZeroGammamom[1]))*TMath::RadToDeg();
    if(kFill)AnaIO::hTruthPi0OA->Fill(OA);
    if(kFill)AnaIO::hTruthLeadingPi0GammaOA->Fill(OA);
    if(kFill)AnaIO::hTruthSubLeadingPi0GammaOA->Fill(OA);
    if(kFill)plotUtils.FillHist(AnaIO::hTruthPi0GammaEnergy,LeadingGamma.E(),0);
    if(kFill)plotUtils.FillHist(AnaIO::hTruthPi0GammaEnergy,SubldGamma.E(),1);

    if(kFill)AnaIO::hTruthLeadingPiZeroGammaDist->Fill(ShowerDist[leadingGammaID]);
    if(kFill)AnaIO::hTruthSubLeadingPiZeroGammaDist->Fill(ShowerDist[subldGammaID]);

    if(kFill)plotUtils.FillHist(AnaIO::hTruthPiZeroGammaE1E2,PiZeroGammamom[leadingGammaID],PiZeroGammamom[subldGammaID]);
    if(kFill)plotUtils.FillHist(AnaIO::hTruthPiZeroGammaE1OA,PiZeroGammamom[leadingGammaID],OA);


  }
  // Rare decay - one gamma one electron and one positron
  else {
    LeadingGamma.SetVectM(bufferPiZeroGammamom[leadingGammaID], 0);
    Electron.SetVectM(bufferPiZeroElectronmom[0], AnaFunctions::ElectronMass());
    Positron.SetVectM(bufferPiZeroPositronmom[0], AnaFunctions::ElectronMass());
    if(kFill)AnaIO::hTruthRarePi0GammaP->Fill(PiZeroGammamom[leadingGammaID]);
    if(kFill)AnaIO::hTruthRarePi0ElectronP->Fill(PiZeroElectronmom[0]);
    if(kFill)AnaIO::hTruthRarePi0PositronP->Fill(PiZeroPositronmom[0]);
  }

  // Fill vector of FS particles 
  vector<TLorentzVector> vec;
  vec.push_back(LeadingGamma);
  vec.push_back(SubldGamma);
  return vec; 

}

int AnaUtils::GetTruthPDGFromID(const int inID, const vector<int> * idarray, const vector<int> * pdgarray)
{ 
  int outpdg = -999;
  for(unsigned int ii = 0; ii<idarray->size(); ii++){
    if((*idarray)[ii] == inID ){
      outpdg = (*pdgarray)[ii];
    }
  }
  return outpdg;
}

string AnaUtils::GetTruthProcessFromID(const int inID, const vector<int> * idarray, const vector<string> * processarray)
{ 
  string outprocess = "Null";
  for(unsigned int ii = 0; ii<idarray->size(); ii++){
    if((*idarray)[ii] == inID ){
      outprocess = (*processarray)[ii];
    }
  }
  return outprocess;
}

string AnaUtils::GetTruthEndProcessFromID(const int inID, const vector<int> * idarray, const vector<string> * Endprocessarray)
{ 
  string outEndprocess = "Null";
  for(unsigned int ii = 0; ii<idarray->size(); ii++){
    if((*idarray)[ii] == inID ){
      outEndprocess = (*Endprocessarray)[ii];
    }
  }
  return outEndprocess;
}

// Aim to find secondary daughter particles
int AnaUtils::GetTruthParticleInfoFromRec(const int recidx)
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

  //const string process = (*AnaIO::true_beam_daughter_Process)[truthID];
  //const string EndProcess = (*AnaIO::true_beam_daughter_endProcess)[truthID];
  
  //-------------- First search direct daughter ----------------// 
  //(Loop over true beam daughther ID and PDG by function GetTruthPDGFromID) 
  pdg = GetTruthPDGFromID(truthID, AnaIO::true_beam_daughter_ID, AnaIO::true_beam_daughter_PDG);

  if(pdg!=-999){ //1. is direct daughter
    // Not a proton, pion, electron, muon, photon or kaon
    if(pdg!=2212 && TMath::Abs(pdg)!=211 && TMath::Abs(pdg)!=11 && TMath::Abs(pdg)!=13 && TMath::Abs(pdg)!=22 && TMath::Abs(pdg)!=321 && pdg<1000000000){
      // Not a strange baryon
      if(pdg!=3112 && pdg!=3222 && pdg!=3212 && pdg!=3122){
        printf("GetTruthParticleInfoFromRec reconstructed truth daughter not proton or pion! %d %d\n", pdg, directPDG); exit(1);
      }
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
      //string process = GetTruthProcessFromID(truthID, AnaIO::true_beam_grand_daughter_ID, AnaIO::true_beam_grand_daughter_Process);
      //string EndProcess = GetTruthEndProcessFromID(truthID, AnaIO::true_beam_grand_daughter_ID, AnaIO::true_beam_grand_daughter_endProcess);
      if(pdg!=-999){
        /*
        if(TMath::Abs(pdg)==22) {
          cout << "2nd gamma" << endl;
          cout << "Process: " << process << endl;
          cout << "Endprocess: " << EndProcess << endl;
        }
        if(TMath::Abs(pdg)==13) {
          cout << "2nd muon" << endl;
          cout << "Process: " << process << endl;
          cout << "Endprocess: " << EndProcess << endl;
        }
        if(TMath::Abs(pdg)==11) {
          cout << "2nd ele" << endl;
          cout << "Process: " << process << endl;
          cout << "Endprocess: " << EndProcess << endl;
        }
        if(TMath::Abs(pdg)==2212) {
          cout << "2nd proton" << endl;
          cout << "Process: " << process << endl;
          cout << "Endprocess: " << EndProcess << endl;
        }
        if(pdg==211) {
          cout << "2nd pi+" << endl;
          cout << "Process: " << process << endl;
          cout << "Endprocess: " << EndProcess << endl;
        }
        if(pdg==-211) {
          cout << "2nd pi-" << endl;
          cout << "Process: " << process << endl;
          cout << "Endprocess: " << EndProcess << endl;
        }
        */
      }
      else{
        //--- lump great grand daughter here
        if(TMath::Abs(directPDG)==11 || TMath::Abs(directPDG)==13 || TMath::Abs(directPDG)==22 || TMath::Abs(directPDG)==211 || TMath::Abs(directPDG)==321 || directPDG==2212 || directPDG==2112 || directPDG==-1 || directPDG==3112 || directPDG==3222 || directPDG==3212 || directPDG==3122 || directPDG>1000000000){//when no true match found pdg = -1
          pdg = directPDG;
        }
        else{
          printf("AnaUtils::GetTruthFromRec search not done! %d %d\n", recidx, directPDG); exit(1);
        }
      }
    } // End of not pi0 daughter
  } // End of not direct daughter
  // The directPDG should be equal to the pdg we found
  if(directPDG!=pdg){
    printf("GetTruthFromRec inconsistent PDG %d %d\n", pdg, directPDG); exit(1);
  }
  // Default truthParticleType is others
  int truthParticleType = gkOthers;
  if(isPrimary){
    if(pdg==2212){//proton
      truthParticleType = gkProton;
    }
    else if(pdg==211){//pi+
      truthParticleType = gkPiPlus;
    }
    else if(pdg==-211){//pi-
      truthParticleType = gkPiMinus;
    }
    else if(pdg==22){//gamma
      truthParticleType = gkGamma;
    }
  }
  else{
    if(pdg==2212){//proton
      truthParticleType = gkSecondaryProton;
    }
    else if(pdg==211){//pi+
      truthParticleType = gkSecondaryPiPlus;
    }
    else if(pdg==-211){//pi-
      truthParticleType = gkSecondaryPiMinus;
    }
    else if(pdg==22){//gamma
      truthParticleType = gkSecondaryGamma;
    }
    else if(TMath::Abs(pdg)==11){//e+/-
      truthParticleType = gkSecondaryEplusEminus;
    }
    else if(TMath::Abs(pdg)==13){//mu+/-
      truthParticleType = gkSecondaryMuon;
    }
  }

  return truthParticleType;
}



double AnaUtils::MakeTrueIncidentEnergies(vector<double> *true_beam_traj_Z, vector<double> *true_beam_traj_KE, vector<double> *true_beam_new_incidentEnergies)
{
  // Only include trajectory points starting in the active volume
  double fTrajZStart = -0.49375;
  // Only include trajectory points less than slice 464 (the end of APA3)
  int fSliceCut = 464;
  // ProtoDUNE TPC wire pitch [cm]
  double fPitch = 0.4794; 

  double true_beam_new_interactingEnergy = -999;
  double next_slice_z = fTrajZStart;
  int next_slice_num = 0;
  
  for (size_t j = 1; j < true_beam_traj_Z->size() - 1; ++j) {
    double z = true_beam_traj_Z->at(j);
    double ke = true_beam_traj_KE->at(j);

    if (z < fTrajZStart) continue;

    if (z >= next_slice_z) {
      double temp_z = true_beam_traj_Z->at(j-1);
      double temp_e = true_beam_traj_KE->at(j-1);

      while (next_slice_z < z && next_slice_num < fSliceCut) {
        double sub_z = next_slice_z - temp_z;
        double delta_e = true_beam_traj_KE->at(j-1) - ke;
        double delta_z = z - true_beam_traj_Z->at(j-1);
        temp_e -= (sub_z/delta_z)*delta_e;

        true_beam_new_incidentEnergies->push_back(temp_e);
        temp_z = next_slice_z;
        next_slice_z += fPitch;
        ++next_slice_num;
      }
    }
  }
  // If the trajectory does not reach the end of the fiducial slices it must have interacted.
  // The interacting energy will be the last incident energy.
  if( next_slice_num < fSliceCut && true_beam_new_incidentEnergies->size() > 0 ) {
    true_beam_new_interactingEnergy = true_beam_new_incidentEnergies->back();
  } // if not true_beam_new_interactingEnergy will = -999
  return true_beam_new_interactingEnergy;
}


double AnaUtils::MakeRecoIncidentEnergies(vector<double> *reco_beam_traj_Z, vector<double> *reco_beam_traj_KE, vector<double> *reco_beam_new_incidentEnergies)
{
  // Only include trajectory points starting in the active volume
  double fTrajZStart = -0.49375;
  // Only include trajectory points less than slice 464 (the end of APA3)
  int fSliceCut = 464;
  //int fSliceCut = 1000000000;
  // ProtoDUNE TPC wire pitch [cm]
  double fPitch = 0.4794; 
  //double fPitch = 20; 


  double reco_beam_new_interactingEnergy = -999;
  double next_slice_z = fTrajZStart;
  int next_slice_num = 0;
  //int size = reco_beam_traj_Z->size();
  //cout << "size: " << size << endl;
  //cout << "reco_beam_traj_KE->at(size): " << reco_beam_traj_KE->at(size-1) << endl;
  //cout << "reco_beam_traj_KE->at(size-2): " << reco_beam_traj_KE->at(size-2) << endl;
  //cout << "traj points size: " << reco_beam_traj_Z->size() << endl;

  for (size_t j = 1; j < reco_beam_traj_Z->size() - 1; ++j) {
    double z = reco_beam_traj_Z->at(j);
    double ke = reco_beam_traj_KE->at(j);

    if (z < fTrajZStart) continue;

    if (z >= next_slice_z) {
      double temp_z = reco_beam_traj_Z->at(j-1);
      double temp_e = reco_beam_traj_KE->at(j-1);
      //cout << "new while loop true: " << endl;
      //cout << "tmp_e: " << temp_e << endl;
      //cout << "tmp_z: " << temp_z << endl;

      //cout << "next_slice_z: " << next_slice_z << endl;
      //cout << "next_slice_num: " << next_slice_num << endl;

      while (next_slice_z < z && next_slice_num < fSliceCut) {
        double sub_z = next_slice_z - temp_z;
        //cout << "sub_z true: " << sub_z << endl;
        double delta_e = reco_beam_traj_KE->at(j-1) - ke;
        //cout << "delta_e: " << delta_e << endl;
        double delta_z = z - reco_beam_traj_Z->at(j-1);
        //cout << "delta_z: " << delta_z << endl;
        temp_e -= (sub_z/delta_z)*delta_e;
        //cout << "(sub_z/delta_z)*delta_e: " << (sub_z/delta_z)*delta_e << endl;
        //cout << "AnaIO::reco_beam_slices_deltaE : " << (*AnaIO::reco_beam_slices_deltaE)[j] << endl;
        reco_beam_new_incidentEnergies->push_back(temp_e);
        //cout << "temp_e loop: " << temp_e << endl;
        temp_z = next_slice_z;
        next_slice_z += fPitch;
        ++next_slice_num;
        //cout << "next_slice_num: " << next_slice_num << endl;
      }
    }
  }
  //cout << "new size: " << reco_beam_new_incidentEnergies->size() << endl;
  // If the trajectory does not reach the end of the fiducial slices it must have interacted.
  // The interacting energy will be the last incident energy.
  if( next_slice_num < fSliceCut && reco_beam_new_incidentEnergies->size() > 0 ) {
    reco_beam_new_interactingEnergy = reco_beam_new_incidentEnergies->back();//FIXME upper Syst
    //cout << "reco_beam_new_interactingEnergy: " << reco_beam_new_interactingEnergy << endl;
  }

  return reco_beam_new_interactingEnergy;
}

double AnaUtils::GetInteractingE_CaloBased(const double & ffe){

  const vector<double> * trackPitch_SCE = AnaIO::reco_beam_TrkPitch_SCE;
  const vector<double> * dEdx_SCE = AnaIO::reco_beam_calibrated_dEdX_SCE;
  const vector<double> * caloZ = AnaIO::reco_beam_calo_Z;

  double intE = ffe;
  for(unsigned int ii = 0; ii < trackPitch_SCE->size()-1; ii++){ //-1 to not count the last slice
    double trkpit = (*trackPitch_SCE)[ii];
    double dEdx = (*dEdx_SCE)[ii];
    double z = (*caloZ)[ii];
    if(dEdx < 0) continue;
    double DeltaE = dEdx*trkpit;
    
    /*cout << "\nii: " << ii << endl;
    cout << "trkpit: " << trkpit << endl;
    cout << "dEdx: " << dEdx << endl;
    cout << "DeltaE: " << DeltaE << endl;*/
    if(z < 220.0) intE -= DeltaE;
  }
  
  return intE;

}


double AnaUtils::GetRecoTrackLength()
{
  double reco_trklen = -999;
  reco_trklen_accum.clear();
  int size = AnaIO::reco_beam_calo_Z->size();
  //cout << "size: " << size << endl;
  //cout << "size1: " << AnaIO::reco_beam_incidentEnergies->size() << endl;
  //cout << "size2: " << AnaIO::reco_beam_calo_wire->size() << endl;
  for (int i=1; i<size; i++){
    //cout << "reco_trklen: " << reco_trklen << endl;
    if (i == 1) reco_trklen = 0;
    reco_trklen += sqrt( pow( (*AnaIO::reco_beam_calo_X)[i]-(*AnaIO::reco_beam_calo_X)[i-1], 2)
                        + pow( (*AnaIO::reco_beam_calo_Y)[i]-(*AnaIO::reco_beam_calo_Y)[i-1], 2)
                        + pow( (*AnaIO::reco_beam_calo_Z)[i]-(*AnaIO::reco_beam_calo_Z)[i-1], 2)
                        );
    reco_trklen_accum.push_back(reco_trklen);

  }
  return reco_trklen;
}

double AnaUtils::GetTrueTrackLength()
{
  int start_idx = -1;
  true_trklen_accum.clear();
  int size = AnaIO::true_beam_traj_Z->size();
  for (int i=0; i<size; i++){
    if ((*AnaIO::true_beam_traj_Z)[i] >= 0){
      start_idx = i-1; // the trajectory point before entering the TPC
      if (start_idx < 0) start_idx = -1;
      break;
    }
  }
  double true_trklen = -999; // initialize
  if (start_idx >= 0){
    for (int i=start_idx+1; i<size; i++){
      if (i == start_idx+1) {
        true_trklen = sqrt( pow( (*AnaIO::true_beam_traj_X)[i]-(*AnaIO::true_beam_traj_X)[i-1], 2)
                            + pow( (*AnaIO::true_beam_traj_Y)[i]-(*AnaIO::true_beam_traj_Y)[i-1], 2)
                            + pow( (*AnaIO::true_beam_traj_Z)[i]-(*AnaIO::true_beam_traj_Z)[i-1], 2)
                            ) * (*AnaIO::true_beam_traj_Z)[i]/((*AnaIO::true_beam_traj_Z)[i]-(*AnaIO::true_beam_traj_Z)[i-1]);
      }
      else{
        true_trklen += sqrt( pow( (*AnaIO::true_beam_traj_X)[i]-(*AnaIO::true_beam_traj_X)[i-1], 2)
                            + pow( (*AnaIO::true_beam_traj_Y)[i]-(*AnaIO::true_beam_traj_Y)[i-1], 2)
                            + pow( (*AnaIO::true_beam_traj_Z)[i]-(*AnaIO::true_beam_traj_Z)[i-1], 2)
                            );
      }
      true_trklen_accum.push_back(true_trklen);
    }

    true_ffKE = -999;
    //cout << "(*AnaIO::true_beam_traj_KE)[start_idx+1]: " << (*AnaIO::true_beam_traj_KE)[start_idx+1] << endl;
    //cout << "(true_trklen_accum)[start_idx+1]: " << (true_trklen_accum)[start_idx+1] << endl;

    true_ffKE = (*AnaIO::true_beam_traj_KE)[start_idx];// + 2.18*(true_trklen_accum)[start_idx+1];
    //cout << "true_ffKE in loop: " << true_ffKE << endl;

    int_energy_true = -999;
    int traj_max = AnaIO::true_beam_traj_Z->size()-1;
    if ((*AnaIO::true_beam_traj_KE)[traj_max] != 0) {
      int_energy_true = (*AnaIO::true_beam_traj_KE)[traj_max];
    }
    else {
      int temp = traj_max-1;
      while ((*AnaIO::true_beam_traj_KE)[temp] == 0) temp--;
      int_energy_true = (*AnaIO::true_beam_traj_KE)[temp] - 2.1*((true_trklen_accum)[traj_max]-(true_trklen_accum)[temp]); // 2.1 MeV/cm
    }
    // Sunbin's method
    /*for(unsigned int i = 0; i < AnaIO::true_beam_traj_KE->size(); i++){
      if((*AnaIO::true_beam_traj_KE)[i] > 0.1){
        int_energy_true = (*AnaIO::true_beam_traj_KE)[i];
      }
      //cout << "[true_beam_traj_KE.at(" << i << ") : " << (*AnaIO::true_beam_traj_KE)->at(i) << endl;
    }*/
  }

  return true_trklen;
}

int AnaUtils::GetInitialSliceID(){
  // true initial sliceID
  int true_ini_sliceID = int(ceil( (Eklim - true_ffKE)/Eslicewidth )); // ignore incomplete slices
  if (true_ini_sliceID < 0) true_ini_sliceID = -1; // both physical and unphysical underflow
  if (true_ini_sliceID >= nthinslices) true_ini_sliceID = nthinslices; // overflow (Eff<pi::Eslicewidth)
  cout << "true_ini_sliceID: " << true_ini_sliceID << endl;
  return true_ini_sliceID;
}

int AnaUtils::GetInteractionSliceID(){

  // true interaction sliceID
  int true_sliceID = int(floor( (Eklim-int_energy_true)/Eslicewidth ));
  //if (true_sliceID <= -99) true_sliceID = -99;
  if (true_sliceID < 0) true_sliceID = -1; // unphysical underflow
  if (true_sliceID >= nthinslices) true_sliceID = nthinslices; // overflow (int_energy_true <= 0)
  // ignore incomplete slices
  double true_ini_sliceID = GetInitialSliceID();
  if (true_sliceID < true_ini_sliceID) {
    true_ini_sliceID = -1;
    true_sliceID = -1;
  } // if true_sliceID==-1, this event should not be used when calculating true XS (but should it be used in unfolding???)

  return true_sliceID;
}

bool AnaUtils::IsSignal(const int nProton, const int nPiZero, const int nPiPlus, const int nPiMinus, const int nParticleBkg, const int nPionAboveThreshold)
{
  bool tmpSig = false;
  // TKI event selection 
  //if(nProton > 0 && nPiZero > 0 && nPiPlus == 0 && nParticleBkg == 0) tmpSig = true;
  // Pi0 inclusive event selection
  if(nPiZero == 1 && nPionAboveThreshold == 0 /*&& nParticleBkg == 0*/) tmpSig = true;

  //if(nPiZero == 1 && nPiPlus == 0 && nPiMinus == 0 /*&& nParticleBkg == 0*/) tmpSig = true;
  return tmpSig;
} 

TVector3 AnaUtils::GetRecBeamFull(){

  TVector3 beamdir;
  // Set beam end direction
  // index 1 uses direction from line projected between last 2 points;
  // index 0 uses direction from line projected between first and last point;
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
/*
  double P_cor = GetBeamCorrectedE(pionEndP);
  double Theta_cor = GetBeamCorrectedTheta(beamdir.Theta()*TMath::RadToDeg())*TMath::DegToRad();
  double Phi_cor = GetBeamCorrectedPhi(beamdir.Phi()*TMath::RadToDeg())*TMath::DegToRad();

  beamdir.SetTheta(Theta_cor);
  beamdir.SetPhi(Phi_cor);
  // Get momentum 3 vector
  const TVector3 fullbeam = beamdir.Unit()*P_cor;
*/
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
  const int evtXStype = kMC ? GetFillXSEventType() : gkXSBmBkg;

  const TVector3 recBeamFull = GetRecBeamFull();
  if(kMC){
    const TVector3 truthBeamFull = GetTruthBeamFull();

    const double beamphiRes    = (recBeamFull.Phi()-truthBeamFull.Phi())*TMath::RadToDeg();//use absolute difference 
    const double beamthetaRes    = (recBeamFull.Theta()-truthBeamFull.Theta())*TMath::RadToDeg();//use absolute difference 
    const double beammomentumRes = recBeamFull.Mag()/truthBeamFull.Mag()-1;

    plotUtils.FillHist(AnaIO::hBeamPhiRes,    truthBeamFull.Phi()*TMath::RadToDeg(), beamphiRes);
    plotUtils.FillHist(AnaIO::hBeamThetaRes,    truthBeamFull.Theta()*TMath::RadToDeg(), beamthetaRes);
    plotUtils.FillHist(AnaIO::hBeamMomentumRes, truthBeamFull.Mag(),                     beammomentumRes);

    // Correction
    plotUtils.FillHist(AnaIO::hBeamMomentumRecVSTruth_REG_Correction,recBeamFull.Mag(),beammomentumRes);
    plotUtils.FillHist(AnaIO::hBeamPhiRecVSTruth_REG_Correction,recBeamFull.Phi()*TMath::RadToDeg(),beamphiRes);
    plotUtils.FillHist(AnaIO::hBeamThetaRecVSTruth_REG_Correction,recBeamFull.Theta()*TMath::RadToDeg(),beamthetaRes);

    const double recoLen = GetRecoTrackLength();
    const double trueLen = GetTrueTrackLength();

    const double beamtrackLenRes = (recoLen-trueLen)/recoLen;
    plotUtils.FillHist(AnaIO::hBeamTrackLengthRes,recoLen,beamtrackLenRes);

  }
  // This evtType only works for MC, data will not have this info but fill it anyway
  plotUtils.FillHist(AnaIO::hRecBeamPhi,    recBeamFull.Phi()*TMath::RadToDeg(), evtType);
  plotUtils.FillHist(AnaIO::hRecBeamTheta,    recBeamFull.Theta()*TMath::RadToDeg(), evtType);
  plotUtils.FillHist(AnaIO::hRecBeamMomentum, recBeamFull.Mag(),                     evtType);

  plotUtils.FillHist(AnaIO::hRecBeamPhi_XsEvt,    recBeamFull.Phi()*TMath::RadToDeg(), evtXStype);
  plotUtils.FillHist(AnaIO::hRecBeamTheta_XsEvt,    recBeamFull.Theta()*TMath::RadToDeg(), evtXStype);
  plotUtils.FillHist(AnaIO::hRecBeamMomentum_XsEvt, recBeamFull.Mag(),                     evtXStype);

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
  
  if(kProton && DoCorrection) {
    double Theta_corrected = GetProtonCorrectedTheta((*AnaIO::reco_daughter_allTrack_Theta)[ii]*TMath::RadToDeg())*TMath::DegToRad();
    double Phi_corrected = GetProtonCorrectedPhi((*AnaIO::reco_daughter_allTrack_Phi)[ii]*TMath::RadToDeg())*TMath::DegToRad();
    trackVectLab.SetMagThetaPhi(trackMBR, Theta_corrected, Phi_corrected);
  }
  else trackVectLab.SetMagThetaPhi(trackMBR, (*AnaIO::reco_daughter_allTrack_Theta)[ii], (*AnaIO::reco_daughter_allTrack_Phi)[ii]);


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


void AnaUtils::FillFSParticleKinematics(const int recIndex, const bool kMC, const int truthParticleType, const int recParticleType)
{
  double weight = 1.0;//CalWeight(kMC);
  // ---------------- Fill proton kinematics ----------------//
  if(recParticleType == gkProton){
    // Get the uncorrected proton momentum vector for comparison 
    const TLorentzVector recMomRawRefBeam = GetMomentumRefBeam(false /*=>reco*/, recIndex, true /*=>proton*/,false /*=>no correction*/);
    // Get this reco particle momentum vector relative to beam (corrected)
    const TLorentzVector recMomRefBeam = GetMomentumRefBeam(false /*=>reco*/, recIndex, true /*=>proton*/);
    // Get this reco particle momentum vector in lab frame (uncorrected)
    const TVector3 recMomLabRaw = GetRecTrackVectLab(recIndex, true, false);
    // Get this reco particle momentum vector in lab frame (corrected)
    const TVector3 recMomLab = GetRecTrackVectLab(recIndex, true, true);

    plotUtils.FillHist(AnaIO::hRecProtonMomentum,recMomRawRefBeam.P(), truthParticleType, weight);
    plotUtils.FillHist(AnaIO::hRecProtonTheta, recMomRawRefBeam.Theta()*TMath::RadToDeg(), truthParticleType, weight);
    // Truth-matching primary proton
    if(truthParticleType == gkProton){ // MC only (Data loop won't pass this)
      // Get the truth-matched particle truth momentum vector relative to beam
      const TLorentzVector truthMomRefBeam = GetMomentumRefBeam(true /*=>truth-matched*/, recIndex, true /*=>proton*/);
      // Get the truth-matched particle truth momentum vector in lab frame
      const TVector3 truthMomLab = GetTruthMatchedTrackVectLab(recIndex);
      // Relative to beam
      const double momentumRes = recMomRefBeam.P()/truthMomRefBeam.P()-1;
      const double momentumRawRes = recMomRawRefBeam.P()/truthMomRefBeam.P()-1;
      //const double thetaRes    = (recMomRefBeam.Theta()-truthMomRefBeam.Theta())*TMath::RadToDeg();
      const double thetaResRaw    = (recMomRawRefBeam.Theta()-truthMomRefBeam.Theta())*TMath::RadToDeg();
      // Lab frame
      const double thetaLabResRaw = (recMomLabRaw.Theta()-truthMomLab.Theta())*TMath::RadToDeg();
      const double phiLabResRaw = (recMomLabRaw.Phi()-truthMomLab.Phi())*TMath::RadToDeg();

      const double thetaLabRes = (recMomLab.Theta()-truthMomLab.Theta())*TMath::RadToDeg();
      const double phiLabRes = (recMomLab.Phi()-truthMomLab.Phi())*TMath::RadToDeg();

      plotUtils.FillHist(AnaIO::hProtonMomentumRes, truthMomRefBeam.P(), momentumRes);
      plotUtils.FillHist(AnaIO::hProtonMomentumRawRes, truthMomRefBeam.P(), momentumRawRes);
      plotUtils.FillHist(AnaIO::hProtonThetaRes, truthMomRefBeam.Theta()*TMath::RadToDeg(), thetaResRaw);
      plotUtils.FillHist(AnaIO::hProtonThetaLabRawRes, truthMomLab.Theta()*TMath::RadToDeg(), thetaLabResRaw);
      plotUtils.FillHist(AnaIO::hProtonPhiLabRawRes, truthMomLab.Phi()*TMath::RadToDeg(), phiLabResRaw);
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

      plotUtils.FillHist(AnaIO::hProtonThetaRecVSTruth_REG_Correction,recMomLabRaw.Theta()*TMath::RadToDeg(),thetaLabResRaw);
      plotUtils.FillHist(AnaIO::hProtonThetaRecVSTruth_REG_AfterCor,recMomLab.Theta()*TMath::RadToDeg(),truthMomLab.Theta()*TMath::RadToDeg());

      plotUtils.FillHist(AnaIO::hProtonPhiRecVSTruth_REG_Correction,recMomLabRaw.Phi()*TMath::RadToDeg(),phiLabResRaw);
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
    const TLorentzVector recMomRefBeam = GetMomentumRefBeam(false /*=>reco*/, recIndex, false/*=>piplus*/, false/*=>no correction*/);
    plotUtils.FillHist(AnaIO::hRecPiPlusMomentum,recMomRefBeam.P(), truthParticleType, weight);
    plotUtils.FillHist(AnaIO::hRecPiPlusTheta, recMomRefBeam.Theta()*TMath::RadToDeg(), truthParticleType, weight);
    // Truth-matching primary piplus
    if(truthParticleType == gkPiPlus){ // MC only (Data loop won't pass this)
      // Get the truth-matched particle truth momentum vector relative to beam
      const TLorentzVector truthMomRefBeam = GetMomentumRefBeam(true /*=>truth-matched*/, recIndex, false /*=>piplus*/, false/*=>no correction*/);
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
    plotUtils.FillHist(AnaIO::hRecShowerEnergy, recShowerMom.E(), truthParticleType, weight);
    plotUtils.FillHist(AnaIO::hRecShowerTheta, recShowerMom.Theta()*TMath::RadToDeg(), truthParticleType, weight);
    plotUtils.FillHist(AnaIO::hRecShowerEnergyRaw, recShowerMomRaw.E(), truthParticleType, weight);
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
  if(kFill) plotUtils.FillHist(AnaIO::hRecPi0Nshower, showerSize, truthEventType); 
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
/*
    cout << "rec shower ld: " << RecPi0Showers[0].E() << "rec shower sl: " << RecPi0Showers[1].E() << endl;
    cout << "truth shower ld: " << TruthPi0Showers[0].E() << "truth shower sl: " << TruthPi0Showers[1].E() << endl;

    cout << "rec shower ld raw: " << RecPi0Showers[2].E() << "rec shower sl raw: " << RecPi0Showers[3].E() << endl;
*/
    const double openingAngle = RecPi0Showers[0].Angle(RecPi0Showers[1].Vect())*TMath::RadToDeg();
    const double openingAngleRaw = RecPi0Showers[2].Angle(RecPi0Showers[3].Vect())*TMath::RadToDeg();
    OA = openingAngle;

    if(kDoKF && kFill) {
      const int evtXStype = kMC ? GetFillXSEventType() : gkXSBmBkg;
      double weight = 1.0;//CalWeight(kMC);

      double Pi0KE = PiZeroVec.E() - AnaFunctions::PiZeroMass();

      if(Pi0KE < 0.30){
        plotUtils.FillHist(AnaIO::testhRecPi0OA_OVERLAY, openingAngle, truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::testhRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::testhRecPi0Mass_OVERLAY_EVT, PiZeroVec_MassCal.M(), evtXStype, weight);

      }
      if(Pi0KE > 0.30){
        plotUtils.FillHist(AnaIO::test1hRecPi0OA_OVERLAY, openingAngle, truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::test1hRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::test1hRecPi0Mass_OVERLAY_EVT, PiZeroVec_MassCal.M(), evtXStype, weight);
      }

      plotUtils.FillHist(AnaIO::hRecPi0KineticEnergyVSPi0Mass, Pi0KE, PiZeroVec_MassCal.M(), weight);

    }
    // KF fit result
    if(kDoKF && kFill){
      if(truthPi0Type == gkTwoGammasSamePi0){
        PiZeroTruthVec = TruthPi0Showers[0] + TruthPi0Showers[1];
        const double pi0EResKF = PiZeroVec.E()/PiZeroTruthVec.E()-1;

        plotUtils.FillHist(AnaIO::hPi0MomentumResFitKF, PiZeroTruthVec.E(), pi0EResKF);
      }
    }
    // Fill pi0 info
    if(kFill){
      const int evtXStype = kMC ? GetFillXSEventType() : gkXSBmBkg;
      

      double weight = 1.0;//CalWeight(kMC);
      Pi0weight = 1.0;//CalPi0OAWeight(kMC,OA);

      weight = weight*Pi0weight;

      if(kFill) AnaIO::hRecPi0TotalE_Fitted->Fill(PiZeroVec.E(),truthPi0Type, weight);
      if(kFill) AnaIO::hRecPi0TotalE_Default->Fill(PiZeroVec_MassCal.E(), truthPi0Type, weight);
      
      //plotUtils.FillHist(AnaIO::hPi0Energy_NoCut,PiZeroVec.E(),truthPi0Type);

      // Fill pi0 mass and momentum
      plotUtils.FillHist(AnaIO::hRecPi0Mass, PiZeroVec_MassCal.M(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0MassRaw, PiZeroVecRaw.M(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0Momentum, PiZeroVec.P(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0MomentumRaw, PiZeroVecRaw.P(), truthPi0Type, weight);
      /*//fake-data control
      bool isFakeData = IsFakeData();

      if(isFakeData) plotUtils.FillHist(AnaIO::hRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type);
      if(!isFakeData) plotUtils.FillHist(AnaIO::hRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type);
      */
      plotUtils.FillHist(AnaIO::hRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0Mass_OVERLAY_EVT, PiZeroVec_MassCal.M(), evtXStype, weight);

      plotUtils.FillHist(AnaIO::hRecPi0MassRaw_OVERLAY, PiZeroVecRaw.M(), truthPi0Type, weight);

      plotUtils.FillHist(AnaIO::hRecPi0Momentum_OVERLAY, PiZeroVec.P(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0MomentumRaw_OVERLAY, PiZeroVecRaw.P(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0OA_OVERLAY, openingAngle, truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0OARaw_OVERLAY, openingAngleRaw, truthPi0Type, weight);

      plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY, PiZeroVec.E(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0ShowerSep_OVERLAY, separation, truthPi0Type, weight);

      plotUtils.FillHist(AnaIO::hRecPi0Theta_OVERLAY, PiZeroVec.Theta()*TMath::RadToDeg(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0Phi_OVERLAY, PiZeroVec.Phi()*TMath::RadToDeg(), truthPi0Type, weight);

      plotUtils.FillHist(AnaIO::hRecPi0ThetaRaw_OVERLAY, PiZeroVecRaw.Theta()*TMath::RadToDeg(), truthPi0Type, weight);
      plotUtils.FillHist(AnaIO::hRecPi0PhiRaw_OVERLAY, PiZeroVecRaw.Phi()*TMath::RadToDeg(), truthPi0Type, weight);
/*
      if(PiZeroVec.E() - AnaFunctions::PiZeroMass() < 0.35){
        plotUtils.FillHist(AnaIO::testhRecPi0OA_OVERLAY, openingAngle, truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::testhRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::testhRecPi0Mass_OVERLAY_EVT, PiZeroVec_MassCal.M(), evtXStype, weight);

      }
      if(PiZeroVec.E() - AnaFunctions::PiZeroMass() > 0.35){
        plotUtils.FillHist(AnaIO::test1hRecPi0OA_OVERLAY, openingAngle, truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::test1hRecPi0Mass_OVERLAY, PiZeroVec_MassCal.M(), truthPi0Type, weight);
        plotUtils.FillHist(AnaIO::test1hRecPi0Mass_OVERLAY_EVT, PiZeroVec_MassCal.M(), evtXStype, weight);

      }
*/
      GoodTruthMatch = false;
      if(truthPi0Type == gkTwoGammasSamePi0){

        const double ldPhotonAngle = RecPi0Showers[0].Angle(TruthPi0Showers[0].Vect())*TMath::RadToDeg();
        const double slPhotonAngle = RecPi0Showers[1].Angle(TruthPi0Showers[1].Vect())*TMath::RadToDeg();

        plotUtils.FillHist(AnaIO::hLeadingPhotonAngleRes, TruthPi0Showers[0].E(), ldPhotonAngle);
        plotUtils.FillHist(AnaIO::hSubLeadingPhotonAngleRes, TruthPi0Showers[1].E(), slPhotonAngle);

        PiZeroTruthVec = TruthPi0Showers[0] + TruthPi0Showers[1];

        const double openingAngleTruth = TruthPi0Showers[0].Angle(TruthPi0Showers[1].Vect())*TMath::RadToDeg();

        plotUtils.FillHist(AnaIO::hOpeningAngleRes, PiZeroTruthVec.P(), openingAngle - openingAngleTruth);

        GoodTruthMatch = true;
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
          AnaIO::hMatchedTruthPi0Energy->Fill(PiZeroTruthVec.E());
          
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
    double weight = 1.0;//CalWeight(kMC);

    plotUtils.FillHist(AnaIO::hRecLeadingShowerEnergy, ldShower.E(), showerTypeArray[nindex[0]], weight); 
    plotUtils.FillHist(AnaIO::hRecSubLeadingShowerEnergy, slShower.E(), showerTypeArray[nindex[1]], weight);
    plotUtils.FillHist(AnaIO::hRecLeadingShowerEnergyRaw, ldShowerRaw.E(), showerTypeArray[nindex[0]], weight);
    plotUtils.FillHist(AnaIO::hRecSubLeadingShowerEnergyRaw, slShowerRaw.E(), showerTypeArray[nindex[1]], weight);
    // Calculate reco opening angle
    const double openingAngle = ldShower.Angle(slShower.Vect())*TMath::RadToDeg();
    plotUtils.FillHist(AnaIO::hRecShowerOpenAngle, openingAngle, truthEventType, weight);

    // Get position vector of two reco showers 
    const TVector3 ldShowerPos = showerPos[nindex[0]];
    const TVector3 slShowerPos = showerPos[nindex[1]];
    const TVector3 distVect = ldShowerPos - slShowerPos;
    // Calculate two reco showers start separation   
    separation = distVect.Mag();
  }

  //if(false){ //FIXME -- speed up 
  if(kDoKF){
    vector<double> FittedVars;
    KF(ldShower, slShower, FittedVars);
    //cout << "default energy: " << ldShower.E() + slShower.E() << endl;
    //const double openingAngle = ldShower.Angle(slShower.Vect());
    //cout << "default energy M2: " << ldShower.E() + 0.134977*0.134977/(2*slShower.E()*(1-cos(openingAngle))) << endl;
    //cout << "default OA: " << openingAngle*TMath::RadToDeg() << endl;
    if(kMC && kFill){

      int truthPi0Type_tmp = -999;
      vector<TLorentzVector> TruthPi0Showers = GetTwoTruthMatchedPi0Showers(truthPi0Type_tmp,kFill);
      double pi0Energy_Truth = TruthPi0Showers[0].E() + TruthPi0Showers[1].E();

      const double openingAngleRad = ldShower.Angle(slShower.Vect());
      const double openingAngleRad_Truth = TruthPi0Showers[0].Angle(TruthPi0Showers[1].Vect());

      // Method A (Pre Fit)
      double pion_E_pre_m1 = ldShower.E() + slShower.E();
      // Method A (Post Fit)
      double pion_E_post_m1 = FittedVars[0] + FittedVars[1];

      // Method B (Pre Fit)
      double pion_E_pre_m2 = ldShower.E() + 0.134977*0.134977/(2*ldShower.E()*(1-cos(openingAngleRad)));
      // Method B (Post Fit)
      double pion_E_post_m2 = FittedVars[0] + 0.134977*0.134977/(2*FittedVars[0]*(1-cos(FittedVars[2])));

      // Method C (Pre Fit)
      double alpha = abs(ldShower.E()-slShower.E())/(ldShower.E()+slShower.E());
      double p_pi0 = 0.134977*sqrt(2/((1-alpha*alpha)*(1-cos(openingAngleRad))));
      double e_pi0 = p_pi0;//sqrt(0.134977*0.134977+p_pi0*p_pi0);
      double pion_E_pre_m3 = e_pi0;
      
      // Method C (Post Fit)
      double alpha_post = abs(FittedVars[0]-FittedVars[1])/(FittedVars[0]+FittedVars[1]);
      double p_pi0_post = 0.134977*sqrt(2/((1-alpha_post*alpha_post)*(1-cos(FittedVars[2]))));
      double e_pi0_post = p_pi0_post;//sqrt(0.134977*0.134977+p_pi0_post*p_pi0_post);
      double pion_E_post_m3 = e_pi0_post;

      if(truthPi0Type_tmp == gkTwoGammasSamePi0){

        AnaIO::hLDShower_PreFit->Fill(ldShower.E()/TruthPi0Showers[0].E()-1);
        AnaIO::hLDShower_PostFit->Fill(FittedVars[0]/TruthPi0Showers[0].E()-1);

        AnaIO::hSLShower_PreFit->Fill(slShower.E()/TruthPi0Showers[1].E()-1);
        AnaIO::hSLShower_PostFit->Fill(FittedVars[1]/TruthPi0Showers[1].E()-1);

        AnaIO::hOAShower_PreFit->Fill(openingAngleRad/openingAngleRad_Truth-1);
        AnaIO::hOAShower_PostFit->Fill(FittedVars[2]/openingAngleRad_Truth-1);
        
        AnaIO::hPi0EnergyE1E2_PreFit->Fill(pion_E_pre_m1/pi0Energy_Truth-1);
        AnaIO::hPi0EnergyE1E2_PostFit->Fill(pion_E_post_m1/pi0Energy_Truth-1);

        AnaIO::hPi0EnergyE1OA_PreFit->Fill(pion_E_pre_m2/pi0Energy_Truth-1);
        AnaIO::hPi0EnergyE1OA_PostFit->Fill(pion_E_post_m2/pi0Energy_Truth-1);

        AnaIO::hPi0EnergyAsym_PreFit->Fill(pion_E_pre_m3/pi0Energy_Truth-1);
        AnaIO::hPi0EnergyAsym_PostFit->Fill(pion_E_post_m3/pi0Energy_Truth-1);


        AnaIO::hShowerE1Compare_REG->Fill(ldShower.E()/TruthPi0Showers[0].E()-1,FittedVars[0]/TruthPi0Showers[0].E()-1);
        AnaIO::hShowerE2Compare_REG->Fill(slShower.E()/TruthPi0Showers[1].E()-1,FittedVars[1]/TruthPi0Showers[1].E()-1);
        AnaIO::hShowerOACompare_REG->Fill(openingAngleRad*TMath::RadToDeg()-openingAngleRad_Truth*TMath::RadToDeg(), FittedVars[2]*TMath::RadToDeg()-openingAngleRad_Truth*TMath::RadToDeg());

        AnaIO::hPi0EnergyComparem1_REG->Fill(pion_E_pre_m1/pi0Energy_Truth-1,pion_E_post_m1/pi0Energy_Truth-1);
        AnaIO::hPi0EnergyComparem2_REG->Fill(pion_E_pre_m2/pi0Energy_Truth-1,pion_E_post_m2/pi0Energy_Truth-1);


        const double pi0EResKF = pion_E_post_m2/pi0Energy_Truth-1;
        plotUtils.FillHist(AnaIO::hPi0MomentumResFitKF_M2, pi0Energy_Truth, pi0EResKF);


      }
    }
    
    // Set the shower energy after KF
    ldShower.SetE(FittedVars[0]);
    slShower.SetE(FittedVars[1]);
    //cout << "KF energy: " << FittedVars[0] + 0.134977*0.134977/(2*FittedVars[0]*(1-cos(FittedVars[2])))  << endl;
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
  TLorentzVector ldShowerTruth;// = showerTruthArray[nindex[0]];
  TLorentzVector slShowerTruth;// = showerTruthArray[nindex[1]];

  if(showerTruthArray[nindex[0]].E() < showerTruthArray[nindex[1]].E()){
    ldShowerTruth = showerTruthArray[nindex[1]];
    slShowerTruth = showerTruthArray[nindex[0]];
  }

  else{
    ldShowerTruth = showerTruthArray[nindex[0]];
    slShowerTruth = showerTruthArray[nindex[1]];
  }
  truthPi0TypeWeight = truthPi0Type;
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
    // Sort reco shower energy
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
      

      //cout << "truthPi0Type: " << truthPi0Type << endl;
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
        if(ldShowerTruth.E() < slShowerTruth.E()){
          LdShowerEnergyTruth.push_back(slShowerTruth.E());
          SlShowerEnergyTruth.push_back(ldShowerTruth.E());
          LdShowerDirTruth.push_back(slShowerTruth.Vect());
          SlShowerDirTruth.push_back(ldShowerTruth.Vect());
        }
        else{
          LdShowerEnergyTruth.push_back(ldShowerTruth.E());
          SlShowerEnergyTruth.push_back(slShowerTruth.E());
          LdShowerDirTruth.push_back(ldShowerTruth.Vect());
          SlShowerDirTruth.push_back(slShowerTruth.Vect());
        }

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
void AnaUtils::TruthMatchingTKI(TLorentzVector dummypi0, TLorentzVector dummyproton, TLorentzVector dummypi0Truth, TLorentzVector dummyprotonTruth, const bool kMC, const bool GoodTruthMatch)
{
  //const int truthEventType = GetFillEventType();
  const int truthEventType = GetFillTKIEventType();
  /*cout << "truthEventType: " << truthEventType << endl;
  cout << "dummyproton: " << dummyproton.P() << endl;
  cout << "dummyprotonTruth: " << dummyprotonTruth.P() << endl;
  cout << "dummypi0: " << dummypi0.P() << endl;
  cout << "dummypi0Truth: " << dummypi0Truth.P() << endl;
*/
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
  if(kMC && GoodTruthMatch){
    AnaFunctions::getCommonTKI(targetA, targetZ, &beamFullPTruth, &(dummypi0Truth), &(dummyprotonTruth), dalphat_truth, dphit_truth, dpt_truth, pn_truth, finPitheta_truth, finProtontheta_truth);
    plotUtils.FillHist(AnaIO::hRecdalphat_truth,dalphat_truth,truthEventType);
    plotUtils.FillHist(AnaIO::hRecdphit_truth,dphit_truth,truthEventType);
    plotUtils.FillHist(AnaIO::hRecdpt_truth,dpt_truth,truthEventType);
    plotUtils.FillHist(AnaIO::hRecpn_truth,pn_truth,truthEventType);

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

      AnaIO::hPi0ShowerE1Compare_REG->Fill(AnaUtils::LdShowerEnergyRaw[ii]/AnaUtils::LdShowerEnergyTruth[ii]-1,FittedVars[0]/AnaUtils::LdShowerEnergyTruth[ii]-1);
      AnaIO::hPi0ShowerE2Compare_REG->Fill(AnaUtils::SlShowerEnergyRaw[ii]/AnaUtils::SlShowerEnergyTruth[ii]-1,FittedVars[1]/AnaUtils::SlShowerEnergyTruth[ii]-1);
      AnaIO::hPi0ShowerOACompare_REG->Fill(AnaUtils::OpenAngle[ii]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg(), FittedVars[2]*TMath::RadToDeg()-AnaUtils::OpenAngleTruth[ii]*TMath::RadToDeg());
      // Print out some info
      /*cout << "LD shower E (truth): " << AnaUtils::LdShowerEnergyTruth[ii] << endl;
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
*/
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

      //cout << "LD: " << AnaUtils::LdShowerEnergyRaw[ii] - AnaUtils::LdShowerEnergyTruth[ii] << endl;
      //cout << "SL: " << AnaUtils::SlShowerEnergyRaw[ii] - AnaUtils::SlShowerEnergyTruth[ii] << endl;
      //cout << "OA: " << (AnaUtils::OpenAngle[ii] - AnaUtils::OpenAngleTruth[ii])*TMath::RadToDeg() << endl;

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
      //cout << "pion_E_Truth: " << pion_E_Truth << " pion_E_Truth_dir: " << pion_E_Truth_dir << " e_pi0_Truth: " << e_pi0_Truth << endl;
      if(pion_E_Truth_dir > 0 && e_pi0_Truth > 0)
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
        Vect.push_back(0.0110251);
        Vect.push_back(0.00148956);
        Vect.push_back(-0.00596679);
        Vect.push_back(0.00148956);
        Vect.push_back(0.00315956);
        Vect.push_back(-0.00201437);
        Vect.push_back(-0.00596679);
        Vect.push_back(-0.00201437);
        Vect.push_back(0.118713);
      }
      else if(xx==1 && yy==2){
        Vect.push_back(0.0106343);
        Vect.push_back(0.00094108);
        Vect.push_back(-0.00810455);
        Vect.push_back(0.00094108);
        Vect.push_back(0.00229943);
        Vect.push_back(0.00129165);
        Vect.push_back(-0.00810455);
        Vect.push_back(0.00129165);
        Vect.push_back(0.110837);
      }
      else if(xx==2 && yy==1){
        Vect.push_back(0.00622274);
        Vect.push_back(-0.000985182);
        Vect.push_back(-0.000145423);
        Vect.push_back(-0.000985182);
        Vect.push_back(0.00573335);
        Vect.push_back(-0.00686372);
        Vect.push_back(-0.000145423);
        Vect.push_back(-0.00686372);
        Vect.push_back(0.147113);
      }
      else if(xx==2 && yy==2){
        Vect.push_back(0.00550302);
        Vect.push_back(-0.000419183);
        Vect.push_back(-0.0019177);
        Vect.push_back(-0.000419183);
        Vect.push_back(0.00383982);
        Vect.push_back(0.00117953);
        Vect.push_back(-0.0019177);
        Vect.push_back(0.00117953);
        Vect.push_back(0.127052);
      }
      CVM[bin] = Vect;
    }
  }
}

void AnaUtils::SetCVMind(){

  CVMind.clear();

  CVMind.push_back(0.00825215);
  CVMind.push_back(0.000305574);
  CVMind.push_back(-0.00439822);
  CVMind.push_back(0.000305574);
  CVMind.push_back(0.00371694);
  CVMind.push_back(-0.00167137);
  CVMind.push_back(-0.00439822);
  CVMind.push_back(-0.00167137);
  CVMind.push_back(0.118882);

}


/* Without swap truth showers
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
*/
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
  /*vector<double> BinE1={0,0.3,5};
  vector<double> BinE2={0,0.15,5};
  std::pair <int,int> bin;
  
  for(unsigned int xx = 1; xx <= BinE1.size() - 1; xx++){
    for(unsigned int yy = 1; yy <= BinE2.size() - 1; yy++){
      if(ldShower.E() > BinE1[xx - 1]  && ldShower.E() < BinE1[xx] 
        && slShower.E() > BinE2[yy - 1] && slShower.E() < BinE2[yy]){
        bin = std::make_pair (xx,yy);
      }
    }
  }*/
  //cout << "ldShower.E(): " << ldShower.E() << endl;
  //cout << "slShower.E(): " << slShower.E() << endl;
/*
  cout << "bin CVM: "<< endl;
  for(auto i : AnaUtils::CVM[bin]){
    cout << "ele: " << i << endl;
  }

  cout << "bin CVM: "<< endl;
  for(auto i : CVMind){
    cout << "ele: " << i << endl;
  }
*/
  //SetCVM();
  //SetCVMind();

 

  GoodFit = false;
  FittedVars = DoKF(ldShower.E(),slShower.E(),openingAngleRad,CVMind,GoodFit);

  AnaIO::hKFPassRate->Fill(GoodFit);

}

void AnaUtils::FillXSTrueHistograms(int &true_avaPionBeam, int &true_avaPionCEXevt, int &true_avaDiffCEXevt){
  // Option to control fake data sample
  bool isFakeData = IsFakeData();
  // Get the MC weight for each event
  double weight = CalWeight(true);
  double g4rw = CalG4RW();
  weight = weight*g4rw;

  // Select true pion beam
  if(AnaIO::true_beam_PDG == 211 && !isFakeData){
    // Clear the incident beam energy vector
    AnaIO::true_beam_incidentEnergies->clear();
    // Fill the beam traj Z position with/without SCE correction
    for(unsigned int ii = 0; ii < AnaIO::true_beam_traj_Z->size(); ii ++){
      AnaIO::hTruthIncidentZ->Fill((*AnaIO::true_beam_traj_Z)[ii]);
      AnaIO::hTruthIncidentZ_SCE->Fill((*AnaIO::true_beam_traj_Z_SCE)[ii]);
    }
    // Calculate the true track length at each space point and the front-face energy and interacting energy (Method 1)
    double trackLen = GetTrueTrackLength();
    // Get the true track length vector at each point
    vector<double> trackLenAccum = GetTrueTrackLengthAccumVect();
    // Calculate the new beam incident energy vector
    // Last element is the interacting energy
    double interactingE = MakeTrueIncidentEnergies(AnaIO::true_beam_traj_Z, AnaIO::true_beam_traj_KE, AnaIO::true_beam_incidentEnergies);

    if(AnaIO::true_beam_incidentEnergies->size() != 0){
      true_avaPionBeam++;

      // First element is the initial energy
      // One method
      double initialE = (*AnaIO::true_beam_incidentEnergies)[0];
      // Another method
      //double initialE = GetTrueFrontFaceEnergy();
      //interactingE = GetTrueIntEnergy();
      //if(initialE > 1050.0) return;
      
      double ffe = -999, intE = -999;
      if(trackLen!= -999){
        ffe = GetTrueFrontFaceEnergy();
        intE = GetTrueIntEnergy();
        AnaIO::hTruthTestFFEnergyM1->Fill(ffe, weight);
        AnaIO::hTruthTestFFEnergyM2->Fill(initialE, weight);
        AnaIO::hTruthTestIntEnergyM1->Fill(intE, weight);
        AnaIO::hTruthTestIntEnergyM2->Fill(interactingE, weight);
        AnaIO::hTruthTestFFEnergyDiff->Fill(ffe-initialE, weight);
        AnaIO::hTruthTestIntEnergyDiff->Fill(intE-interactingE, weight);
      } 

      //vector<vector<double>> incE_Info = ComputeTrueIncidentHist(initialE, interactingE, ffe, intE, AnaIO::true_beam_incidentEnergies);
      vector<vector<double>> incE_Info = ComputeTrueIncidentHist(initialE, interactingE, initialE, interactingE, AnaIO::true_beam_incidentEnergies);
      
      // == Binnings for cross sections
      //int N_binning_100MeV = 20;
      //vector<double> binning_100MeV = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000.};
      //vector<double> binning_100MeV = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};
      //double interactingE_tmp = interactingE;
      //if(interactingE == -999) interactingE_tmp = AnaIO::true_beam_incidentEnergies->back();

      //if(interactingE != -999) FillEsliceHistograms(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthCEXInteractingHist, initialE, interactingE_tmp, -1, weight, binning_100MeV, N_binning_100MeV, false);

      if(incE_Info.size() != 0){
        int incE_size = incE_Info[0].size();
        for(int kk = 0; kk < incE_size; kk++){
          double incE = incE_Info[0][kk];
          double wt = incE_Info[1][kk];
          //AnaIO::hTruthIncidentHist->Fill(incE,wt);
          AnaIO::hTruthBeamIncidentHist->Fill(incE,wt*weight);
        }
      }

      for(unsigned int ii = 0; ii < AnaIO::true_beam_incidentEnergies->size(); ii ++){
        double incidentE = (*AnaIO::true_beam_incidentEnergies)[ii];
        AnaIO::hTruthIncidentHist->Fill(incidentE, weight);
        AnaIO::hTruthBeamIncidentHistOldM->Fill(incidentE, weight);
        AnaIO::hTruthSingleIncidentHist->Fill(incidentE, weight);       
      }
      // Select all true inelastic events
      if((*AnaIO::true_beam_endProcess) == "pi+Inelastic"){
        // Fill the interacting histogram for inelastic events
        if(interactingE != -999) AnaIO::hTruthInteractingHist->Fill(interactingE, weight);
        //(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthInteractingHist, initialE, interactingE, interactingE, weight, binning_100MeV, N_binning_100MeV, true);
        //if(interactingE != -999) FillEsliceHistograms(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthInteractingHist, initialE, interactingE, interactingE, weight, binning_100MeV, N_binning_100MeV, true);

        AnaIO::hTruthSingleInteractingHist->Fill(interactingE, weight);
      }
      // Select all true charge exchange (CEX) events i.e. 1-pi0 and 0-pi+/-
      if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" &&  AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0) {
        true_avaPionCEXevt++;
        // Fill the interacting histogram for CEX events
        if(interactingE != -999) AnaIO::hTruthCEXInteractingHist->Fill(interactingE, weight);
        //if(interactingE != -999) FillEsliceHistograms(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthCEXInteractingHist, initialE, interactingE, interactingE, weight, binning_100MeV, N_binning_100MeV, true);

        // Get the true pi0 4-vectors (0 element of this vector is pion)
        vector<TLorentzVector> vecFSParticle = GetFSParticlesTruth();
        // Get the FS pi0 momentum
        double LeadingPiZeroP = vecFSParticle[0].P();
        double LeadingPiZeroPz = vecFSParticle[0].Pz();
        // Calculate the kinetic energy
        double LeadingPiZeroKE = sqrt(LeadingPiZeroP*LeadingPiZeroP + pow(AnaFunctions::PiZeroMass(),2)) - AnaFunctions::PiZeroMass(); 
        // Fill the pi0 KE spectrum --Todo change the name of this histogram?
        AnaIO::hTruthDiffCEXInteractingHist->Fill(LeadingPiZeroKE*1000, weight);
        // Compute the cos theta
        double costheta = LeadingPiZeroPz/LeadingPiZeroP;
        // Fill the pi0 KE VS cos theta (all sample)
        AnaIO::hTruthOutgoingKEcosTheta->Fill(LeadingPiZeroKE*1000, costheta, weight);
        // Get the theta angle between pz and p in degree
        double theta = TMath::RadToDeg() * TMath::ACos(costheta);

        // Get true beam momentum vector
        //const TVector3 truthBeamFull = GetTruthBeamFull();
        // Get the angle between pi0 vector and pi+ beam
        //double angle = truthBeamFull.Angle(vecFSParticle[0].Vect());
        //double cos_angle = TMath::Cos(angle);
        
        // Select different pion beam interacting slice for Diff. xsec calculation
        if(interactingE > 650 && interactingE < 800){
          true_avaDiffCEXevt++;
          AnaIO::hTruthDiffCEXInteractingHist_650to800MeV->Fill(LeadingPiZeroKE*1000, weight);
          AnaIO::hTruthDiffCEXInteractingHistTheta_650to800MeV->Fill(theta, weight);
          AnaIO::hTruthDiffCEXInteractingHistCosTheta_650to800MeV->Fill(costheta, weight);

          // Select pi0 KE for three regions
          if(LeadingPiZeroKE*1000 > 650){ // High pi0 KE region

            AnaIO::hTruthCEXDaughters_High->Fill(AnaIO::true_daughter_nProton,0);
            AnaIO::hTruthCEXDaughters_High->Fill(AnaIO::true_daughter_nNeutron,1);
            AnaIO::hTruthCEXDaughters_High->Fill(AnaIO::true_daughter_nNucleus,2);
            AnaIO::hTruthPi0DaughtersCosTheta_High->Fill(LeadingPiZeroPz/LeadingPiZeroP);

            /*cout << "\nhigh region" << endl;
            cout << AnaIO::true_daughter_nNucleus << endl;
            for(auto pdg : (*AnaIO::true_beam_daughter_PDG)){
              cout << pdg << endl;
            }
            cout << "grand daughter" << endl;
            for(auto pdg : (*AnaIO::true_beam_grand_daughter_PDG)){
              cout << pdg << endl;
            }*/
          }
          
          if(LeadingPiZeroKE*1000 > 200 && LeadingPiZeroKE*1000 < 400){ // Middle pi0 KE region
          
            AnaIO::hTruthCEXDaughters_Middle->Fill(AnaIO::true_daughter_nProton,0);
            AnaIO::hTruthCEXDaughters_Middle->Fill(AnaIO::true_daughter_nNeutron,1);
            AnaIO::hTruthCEXDaughters_Middle->Fill(AnaIO::true_daughter_nNucleus,2);
            AnaIO::hTruthPi0DaughtersCosTheta_Middle->Fill(LeadingPiZeroPz/LeadingPiZeroP);
            /*cout << "\nmid region" << endl;
            cout << AnaIO::true_daughter_nNucleus << endl;

            for(auto pdg : (*AnaIO::true_beam_daughter_PDG)){
              cout << pdg << endl;
            }
            cout << "grand daughter" << endl;
            for(auto pdg : (*AnaIO::true_beam_grand_daughter_PDG)){
              cout << pdg << endl;
            }*/
          }

          if(LeadingPiZeroKE*1000 < 150){ // Low pi0 KE Region
            
            AnaIO::hTruthCEXDaughters_Low->Fill(AnaIO::true_daughter_nProton,0);
            AnaIO::hTruthCEXDaughters_Low->Fill(AnaIO::true_daughter_nNeutron,1);
            AnaIO::hTruthCEXDaughters_Low->Fill(AnaIO::true_daughter_nNucleus,2);
            AnaIO::hTruthPi0DaughtersCosTheta_Low->Fill(LeadingPiZeroPz/LeadingPiZeroP);
            /*cout << "\nlow region" << endl;
            cout << AnaIO::true_daughter_nNucleus << endl;

            for(auto pdg : (*AnaIO::true_beam_daughter_PDG)){
              cout << pdg << endl;
            }
            cout << "grand daughter" << endl;
            for(auto pdg : (*AnaIO::true_beam_grand_daughter_PDG)){
              cout << pdg << endl;
            }*/
          }

        }
        if(interactingE > 650 && interactingE < 700) {
          // Fill the pi0 KE VS cos theta (intE = 650-700 MeV)
          AnaIO::hTruthOutgoingKEcosTheta_Mid->Fill(LeadingPiZeroKE*1000, LeadingPiZeroPz/LeadingPiZeroP, weight);
          // Fill the pi0 KE spectrum and angular info
          AnaIO::hTruthDiffCEXInteractingHist_700MeV->Fill(LeadingPiZeroKE*1000, weight);
          AnaIO::hTruthDiffCEXInteractingHistTheta_700MeV->Fill(theta, weight);
          AnaIO::hTruthDiffCEXInteractingHistCosTheta_700MeV->Fill(costheta, weight);
        }
        if(interactingE > 750 && interactingE < 800) {
          // Fill the pi0 KE spectrum and angular info
          AnaIO::hTruthDiffCEXInteractingHist_800MeV->Fill(LeadingPiZeroKE*1000, weight);
          AnaIO::hTruthSingleDiffCEXInteractingHist_800MeV->Fill(LeadingPiZeroKE*1000, weight);
          AnaIO::hTruthDiffCEXInteractingHistTheta_800MeV->Fill(theta, weight);
          AnaIO::hTruthDiffCEXInteractingHistCosTheta_800MeV->Fill(costheta, weight);
        }
        if(interactingE > 850 && interactingE < 900) {
          // Fill the pi0 KE VS cos theta (intE = 850-900 MeV)
          AnaIO::hTruthOutgoingKEcosTheta_High->Fill(LeadingPiZeroKE*1000, LeadingPiZeroPz/LeadingPiZeroP, weight);
          // Fill the pi0 KE spectrum and angular info
          AnaIO::hTruthDiffCEXInteractingHist_900MeV->Fill(LeadingPiZeroKE*1000, weight);
          AnaIO::hTruthDiffCEXInteractingHistTheta_900MeV->Fill(theta, weight);
          AnaIO::hTruthDiffCEXInteractingHistCosTheta_900MeV->Fill(costheta, weight);
          /*
          // Select pi0 KE for three regions
          if(LeadingPiZeroKE*1000 > 800){ // High pi0 KE region

            AnaIO::hTruthCEXDaughters_High->Fill(AnaIO::true_daughter_nProton,0);
            AnaIO::hTruthCEXDaughters_High->Fill(AnaIO::true_daughter_nNeutron,1);
            AnaIO::hTruthCEXDaughters_High->Fill(AnaIO::true_daughter_nNucleus,2);
            AnaIO::hTruthPi0DaughtersCosTheta_High->Fill(LeadingPiZeroPz/LeadingPiZeroP);
          }
          
          if(LeadingPiZeroKE*1000 > 200 && LeadingPiZeroKE*1000 < 400){ // Middle pi0 KE region
          
            AnaIO::hTruthCEXDaughters_Middle->Fill(AnaIO::true_daughter_nProton,0);
            AnaIO::hTruthCEXDaughters_Middle->Fill(AnaIO::true_daughter_nNeutron,1);
            AnaIO::hTruthCEXDaughters_Middle->Fill(AnaIO::true_daughter_nNucleus,2);
            AnaIO::hTruthPi0DaughtersCosTheta_Middle->Fill(LeadingPiZeroPz/LeadingPiZeroP);
          }

          if(LeadingPiZeroKE*1000 < 150){ // Low pi0 KE Region
            
            AnaIO::hTruthCEXDaughters_Low->Fill(AnaIO::true_daughter_nProton,0);
            AnaIO::hTruthCEXDaughters_Low->Fill(AnaIO::true_daughter_nNeutron,1);
            AnaIO::hTruthCEXDaughters_Low->Fill(AnaIO::true_daughter_nNucleus,2);
            AnaIO::hTruthPi0DaughtersCosTheta_Low->Fill(LeadingPiZeroPz/LeadingPiZeroP);
          }*/
        }

        if(interactingE > 450 && interactingE < 500) {
          AnaIO::hTruthOutgoingKEcosTheta_Low->Fill(LeadingPiZeroKE*1000, LeadingPiZeroPz/LeadingPiZeroP);
          AnaIO::hTruthDiffCEXInteractingHist_500MeV->Fill(LeadingPiZeroKE*1000);
        }
        if(interactingE > 250 && interactingE < 300) {
          AnaIO::hTruthDiffCEXInteractingHist_300MeV->Fill(LeadingPiZeroKE*1000);
        }
        if(interactingE > 150 && interactingE < 200) {
          AnaIO::hTruthDiffCEXInteractingHist_200MeV->Fill(LeadingPiZeroKE*1000);
        }
        if(interactingE > 50 && interactingE < 100) {
          AnaIO::hTruthDiffCEXInteractingHist_100MeV->Fill(LeadingPiZeroKE*1000);
        } 
      }
    }
  } // End of truth xsec measurements
}


void AnaUtils::FillXSNewTrueHistograms(){

  // Option to control fake data sample
  bool isFakeData = IsFakeData();
  // Get the MC weight for each event
  double weight = CalWeight(true);
  double g4rw = CalG4RW();
  weight = weight*g4rw;

  // Select true pion beam
  if(AnaIO::true_beam_PDG == 211 && !isFakeData){

    //AnaIO::true_beam_incidentEnergies->clear();
    // Calculate the true track length at each space point and the front-face energy and interacting energy 
    GetTrueTrackLength();
    // Get the true track length vector at each point
    vector<double> trackLenAccum = GetTrueTrackLengthAccumVect();

    double initialE = GetTrueFrontFaceEnergy();
    double interactingE = GetTrueIntEnergy();

    //double interactingE = MakeTrueIncidentEnergies(AnaIO::true_beam_traj_Z, AnaIO::true_beam_traj_KE, AnaIO::true_beam_incidentEnergies);
    //double initialE = (*AnaIO::true_beam_incidentEnergies)[0];
    //double interactingE_tmp = interactingE;
    //if(interactingE==-999) interactingE_tmp = AnaIO::true_beam_incidentEnergies->back();

    // Fill initial and beam int histograms
    FillEsliceHistograms(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthCEXInteractingHist, initialE, interactingE, interactingE, weight, binning_100MeV, N_binning_100MeV, false);

    // Select true pion inelastic events 
    if((*AnaIO::true_beam_endProcess) == "pi+Inelastic"){
      // Fill total inelastic interacting histogram
      FillEsliceHistograms(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthInteractingHist, initialE, interactingE, interactingE, weight, binning_100MeV, N_binning_100MeV, true);
      //AnaIO::hTruthSingleInteractingHist->Fill(interactingE, weight);
    }

    // Select all true charge exchange (CEX) events i.e. 1-pi0 and 0-pi+/-
    if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" &&  AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0) {
      // Fill CEX interacting histogram
      FillEsliceHistograms(AnaIO::hNewTruthBeamInitialHist, AnaIO::hNewTruthBeamInteractingHist, AnaIO::hNewTruthBeamIncidentHist, AnaIO::hNewTruthCEXInteractingHist, initialE, interactingE, interactingE, weight, binning_100MeV, N_binning_100MeV, true);
    }
  }

}


vector<vector<double>> AnaUtils::ComputeTrueIncidentHist(const double & initialE, const double & interactingE, const double & ffe, const double & intE, vector<double> *true_beam_incidentEnergies)
{
  vector<vector<double>> incE_Info;
  vector<double> incEvect;
  vector<double> weightvect;
  // Get the MC weight for each event
  double weight = CalWeight(true);
  double g4rw = CalG4RW();
  weight = weight*g4rw;

  double interactingE_tmp = interactingE;
  double intE_tmp = intE;
  double size = true_beam_incidentEnergies->size();
  //cout << size << endl;
  //double size = AnaIO::true_beam_traj_KE->size();
  // If the beam interacts after the APA3 i.e. it's value is -999. Then set the last element as the interactingE to calculate incident Hist
  if(interactingE == -999) interactingE_tmp = true_beam_incidentEnergies->back();
  if(intE == -999) intE_tmp = true_beam_incidentEnergies->back();
  
  // Calculate the incident entry
  //double slice_width = slice_width;
  //double interval = initialE - interactingE_tmp; // 10MeV bin
  //double nbin = interval/slice_width; // Number of bins 
  double inibin = initialE/slice_width;
  double intbin = interactingE_tmp/slice_width;

  //int y = (int)nbin + 1;
  //double intbin = interactingE/10.;
  int y0 = (int)intbin + 1;
  //int y1 = y0 + y;
  int y1 = (int)inibin + 1;

  double wt = 1.0, wt1 = 1.0;

  if(size != 0){

    wt = intbin - floor(intbin);
    //cout << "interactingE_tmp: " << interactingE_tmp << "intbin: " << intbin << "wt: " << wt << endl;
    wt1 = inibin - floor(inibin);
    //cout << "initialE: " << initialE << "inibin: " << inibin << "wt1: " << wt1 << endl;

    Int_t initialE_bin = AnaIO::hTruthInitialHist->GetXaxis()->FindBin(initialE);
    Int_t interactingE_bin = AnaIO::hTruthBeamInteractingHist->GetXaxis()->FindBin(interactingE_tmp);
    if(initialE_bin==interactingE_bin) AnaIO::hTruthTestSameBin->Fill(initialE_bin);

    AnaIO::hTruthInitialHist->Fill(ffe, weight);
    AnaIO::hTruthBeamInitialHist->Fill(ffe, weight);
    AnaIO::hTruthBeamInteractingHist->Fill(intE_tmp, weight);

    AnaIO::hTruthBeamInitialHist_50MeVbin->Fill(ffe, weight);
    AnaIO::hTruthBeamInteractingHist_50MeVbin->Fill(intE_tmp, weight);
    
    for(int idy = y0; idy <= y1; idy++){
      double incE = idy * slice_width - slice_width/2.0;
      if(idy == y0 && y0 != y1) {incEvect.push_back(incE); weightvect.push_back(wt);}// AnaIO::hTruthSingleIncidentHist->Fill(incE,wt);
      else if(idy == y0 && y0 == y1) {incEvect.push_back(incE); weightvect.push_back(wt1-wt);}//AnaIO::hTruthSingleIncidentHist->Fill(incE,wt1-wt);
      else if(idy != y1) {incEvect.push_back(incE); weightvect.push_back(1.);}//AnaIO::hTruthSingleIncidentHist->Fill(incE);
      else {incEvect.push_back(incE); weightvect.push_back(wt1);}//AnaIO::hTruthSingleIncidentHist->Fill(incE,wt1);
    }
    incE_Info.push_back(incEvect);
    incE_Info.push_back(weightvect);

  }
  return incE_Info;
}

double AnaUtils::CalWeight(const bool & kMC){
  double weight = 1.;
  //return weight;
  
  double mufrac = 1.53;
  
  double mom_mu0 = 1.00362;//1.0033;
  double mom_sigma0 = 0.0595377;//0.0609;
  double mom_mu = 1.01604;//1.01818;
  double mom_sigma = 0.0708354;//0.07192;
  double wlimit = 3.5;//3; 
  
  if(kMC){
    // momentum reweight (outlier weights are set to 1e-5)
    double deno = exp(-pow((AnaIO::true_beam_startP-mom_mu0)/mom_sigma0,2)/2);
    //if (deno < wlimit) deno = wlimit;
    double numo = exp(-pow((AnaIO::true_beam_startP-mom_mu)/mom_sigma,2)/2);
    //if (numo < wlimit) numo = wlimit;
    weight *= numo;
    weight /= deno;
    if (weight>wlimit) weight=wlimit;
    if (weight<1./wlimit) weight=1./wlimit;
    if (AnaIO::true_beam_PDG == -13) weight = 1;
    // muon reweight
    if (AnaIO::true_beam_PDG == -13 && AnaIO::reco_beam_true_byE_matched) { // kMuon
      weight *= mufrac;
    } 
  }
  //double g4rw = CalG4RW();

  //return weight*g4rw;
  return weight;

} 

double AnaUtils::CalBckWeight(const bool & kMC){

  double weight = 1.;

  if (kMC){

    const int origin = AnaIO::reco_beam_true_byHits_origin; // 0: kUnknown 1: kBeamNeutrino 2: kCosmicRay 3: kSuperNovaNeutrino 4: kSingleParticle (single particles thrown at the detector) 
    const bool matched = AnaIO::reco_beam_true_byHits_matched; 
    const TString process = (*AnaIO::reco_beam_true_byHits_process); 
    const int pdg = AnaIO::reco_beam_true_byHits_PDG;

    // Background weight
    double mu_weight = 0.65, p_weight = 1.65, pi_weight = 1.47;
    
    if(process == "primary" && matched && origin == 4 && pdg == -13) { // beam muon bkg  
      weight = mu_weight;     
    }
    if(process == "primary" && matched && origin == 4 && pdg == 211) { // beam pion bkg
      weight = 1.0; 
    }
    else {
      if(pdg == -13) { // secondary muon bkg
        weight = mu_weight;
      }
      else if(pdg == 2212) { // secondary proton bkg
        weight = p_weight;
      }
      else if(abs(pdg) == 211) { // secondary pion bkg
        weight = pi_weight;
      }
    }
  }
  // No need to include bck weight for now
  weight = 1.;
  return weight;
}


double AnaUtils::CalXSEvtWeight(const bool & kMC, const double & intE, const int & evtXStype){

  double weight = 1.;

  // With 12MeV E Loss
  /*if (kMC){
    // 0pi0
    if(evtXStype == gkXSEvtBkgInel){
      if(intE < 650.0) weight = 1.01219;
      else weight = 0.966249;
    }
    // 1pi0
    else if(evtXStype == gkXSEvtBkgSinglePi0 || evtXStype == gkXSEvtBkgMultiPi0){
      if(intE < 700.0) weight = 2.13924;
      else weight = 0.688961;
    }
    // No weight for other bck
    else weight = 1.;
  }*/

  // With 2MeV E Loss
  /*if (kMC){
    // 0pi0
    if(evtXStype == gkXSEvtBkgInel){
      if(intE < 650.0) weight = 1.0101;
      else weight = 0.968285;
    }
    // 1pi0
    else if(evtXStype == gkXSEvtBkgSinglePi0 || evtXStype == gkXSEvtBkgMultiPi0){
      if(intE < 700.0) weight = 2.54864;
      else weight = 0.622455;
    }
    // No weight for other bck
    else weight = 1.;
  }*/

  // With dep E Loss Feb.
  /*if (kMC){
    // 0pi0
    if(evtXStype == gkXSEvtBkgInel){
      if(intE < 700.0) weight = 0.997976;
      else weight = 0.969012;
    }
    // 1pi0
    else if(evtXStype == gkXSEvtBkgSinglePi0 || evtXStype == gkXSEvtBkgMultiPi0){
      if(intE < 700.0) weight = 2.39891;
      else weight = 0.706501;
    }
    // No weight for other bck
    else weight = 1.;
  }*/

  // With dep E Loss Apr. (new)
  if (kMC){
    // 0pi0
    if(evtXStype == gkXSEvtBkgInel){
      if(intE < 700.0) weight = 1.01585;
      else weight = 0.965132;
    }
    // 1pi0
    else if(evtXStype == gkXSEvtBkgSinglePi0 || evtXStype == gkXSEvtBkgMultiPi0){
      if(intE < 700.0) weight = 2.20256;
      else weight = 0.713314;
    }
    // No weight for other bck
    else weight = 1.;
  }

  return weight;
}

double AnaUtils::CalPi0OAWeight(const bool & kMC, const double & OA){
  double weight = 1.;
  //return weight;
  if(kMC){
    //double IniWeight[] = {27.0/59.0638, 83.0/100.622, 135.0/194.633, 178.0/242.524, 174.0/201.365, 142.0/137.606, 125.0/97.6361, 87.0/87.0227, 79.0/70.1378, 63.0/55.6535, 61.0/46.6166, 53.0/35.1211, 
    //                 42.0/22.1165, 34.0/17.5821, 30.0/19.3663, 19.0/8.79664, 18.0/12.3029, 18.0/11.2505, 9.0/10.8019, 5.0/4.03585};
    double IniWeight[] = {27.0/51.0259, 83.0/86.9285, 135.0/168.146, 178.0/209.519, 174.0/173.961, 142.0/118.879, 125.0/84.3489, 87.0/75.1799, 79.0/60.5928, 63.0/48.0797, 61.0/40.2726, 53.0/30.3415, 
                     42.0/19.1067, 34.0/15.1894, 30.0/16.7307, 19.0/7.59951, 18.0/10.6286, 18.0/9.7194, 9.0/9.33184, 5.0/3.48662};
    
    Int_t binx = AnaIO::hRecPi0OA_OVERLAY->GetXaxis()->FindBin(OA);
    weight *= IniWeight[binx-1];
  }
  
  return weight;
}

double AnaUtils::CalBeamIniWeight(const double & iniE){

  double weight = 1.;
  //double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.30802, 0.416522, 
  //                 0.549213, 0.813086, 0.849306, 0.851946, 0.872305, 0.873036, 0.875985, 0.870969};
  //double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.262599, 0.449623, 
  //                 0.657675, 0.844744, 0.840194, 0.857892, 0.867017, 0.875814, 0.877931, 0.865566};
  // ==== new bck comp
  //double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  //                   1.0, 0.857099, 0.868815, 0.869798, 0.885231, 0.890228, 0.887002, 0.88836};
  //double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  //                   1.0, 0.839761, 0.865995, 0.870306, 0.883158, 0.887925, 0.888798, 0.888868};

  // Feb 2023
  //double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  //                   0.0, 0.801123, 0.859177, 0.86919, 0.880718, 0.886977, 0.889779, 0.887032};
  // April 2023 new
  //double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  //                   0.0, 0.858606, 0.870282, 0.882274, 0.887159, 0.892445, 0.894791, 0.889779};

  // April 2023 new latest
  double IniWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     0.0, 0.841146, 0.870546, 0.879581, 0.88791, 0.893029, 0.894007, 0.891124};

  Int_t binx = AnaIO::hRecPiPlusInitialEnergy->GetXaxis()->FindBin(iniE);
  weight *= IniWeight[binx-1];

  return weight;

}

double AnaUtils::CalBeamIntWeight(const double & intE){

  double weight = 1.;
  //double IntWeight[] = {0.32537, 0.419454, 0.146989, 0.183387, 0.382133, 0.521113, 0.526616, 0.597888, 0.666668, 0.74798, 0.80511, 0.84424, 
  //                 0.881101, 0.899525, 0.915137, 0.926683, 0.93251, 0.929668, 0.934508, 0.884038};
  //double IntWeight[] = {0.38796, 0.154826, 0.119607, 0.289206, 0.377158, 0.514024, 0.522487, 0.597579, 0.67405, 0.744799, 0.804759, 0.84912, 
  //                 0.877425, 0.900835, 0.912441, 0.927018, 0.934754, 0.934351, 0.918808, 0.8892};
  // ==== new bck comp
  //double IntWeight[] = {0.465024, 0, 0.547152, 0.205098, 0.57463, 0.529922, 0.550447, 0.622745, 0.682822, 0.785272, 0.827535, 0.869196, 
  //                   0.889707, 0.907593, 0.922458, 0.935966, 0.940945, 0.946625, 0.93507, 0.94387};
  //double IntWeight[] = {0.439768, 0, 1.0, 0.288137, 0.524855, 0.523935, 0.559962, 0.595805, 0.67353, 0.766597, 0.820003, 0.863974, 
  //                    0.884379, 0.906827, 0.920416, 0.930346, 0.941464, 0.943562, 0.944133, 0.944218};
  
  // Feb 2023
  //double IntWeight[] = {0.439768, 1.0, 0.0, 0.240104, 0.325475, 0.547874, 0.525402, 0.569011, 0.663258, 0.755661, 0.820706, 0.858978, 
  //                   0.886342, 0.901809, 0.920025, 0.931599, 0.945724, 0.945179, 0.943406, 0.938913};
  // April 2023 new
  double IntWeight[] = {0, 1.0, 0.0, 0.200557, 0.372566, 0.586148, 0.553071, 0.6047, 0.67587, 0.770301, 0.831756, 0.867062, 
                     0.891131, 0.907382, 0.923663, 0.934314, 0.946617, 0.947189, 0.947067, 0.938439};

  Int_t binx = AnaIO::hRecPiPlusInteractingEnergy->GetXaxis()->FindBin(intE);
  weight *= IntWeight[binx-1];

  return weight;
}

double AnaUtils::CalCEXIntWeight(const double & intE){

  double weight = 1.;
  //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.52362, 0.742799, 0.835278, 0.871476, 0.733481, 0.645124, 
  //                 0.728537, 0.679767, 0.597596, 0.665684, 0.65989, 0.60755, 0.728686, 0.829107};
  //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.675337, 0.701481, 0.932205, 0.816573, 0.746547, 0.705951/1.33, 
  //                 0.742581/1.33, 0.65943, 0.630317*1.33, 0.644456, 0.638316, 0.620395/1.33, 0.73646, 0.608416};
  //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.675337, 0.701481, 0.932205, 0.816573, 0.746547, 0.705951, 
  //                 0.742581, 0.65943, 0.630317, 0.644456, 0.638316, 0.620395, 0.73646, 0.608416};

  // New bck portion
  //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.52362, 0.744688, 0.800389, 0.831086, 0.643261, 0.52957, 
  //                 0.616279, 0.691406, 0.610101, 0.677686, 0.671867, 0.62, 0.737002, 0.836766};
  // ==== new bck comp
  //double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.748344, 0.602063, 0.841422, 0.815263, 0.67177, 0.553415, 
  //                  0.596665, 0.544498, 0.703856, 0.690155, 0.731224, 0.720228, 0.704154, 1};

  double IntWeight[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.69794, 0.818919, 0.775484, 0.800257, 0.53928, 
                    0.719739, 0.530114, 0.753455, 0.752041, 0.76991, 0.74918, 0.743382, 1.0};
    
  Int_t binx = AnaIO::hRecPiPlusInteractingEnergyEvt->GetXaxis()->FindBin(intE);
  weight *= IntWeight[binx-1];

  return weight;
}


double AnaUtils::CalCEXPi0KEWeight(const double & intE){

  double weight = 1.;

  double IntWeight[] = {0.0, 0.513447, 0.703031, 0.568209, 0.520239, 0.66259, 0.573659, 0.590971, 0.593988, 0.543884, 0.66232, 0.694243, 
                    0.884691, 1.0, 0.771999, 1.0, 1.0, 0.0, 0.0, 0.0};

  //double IntWeight[] = {0.0, 0.536816, 0.690962, 0.649954, 0.596227, 0.710346, 0.613647, 0.618428, 0.49974, 0.573666, 0.816561, 0.664146, 
  //                  0.946883, 0.948054, 0.77866, 1.0, 1.0, 0.0, 0.0, 0.0};

  Int_t binx = AnaIO::hRecPiZeroSliceKineticEnergyEvt->GetXaxis()->FindBin(intE);
  weight *= IntWeight[binx-1];

  return weight;
}

double AnaUtils::CalG4RW(){
  /*vector<vector<double>> cex;
  
  double g4rw = 1.0;
  if (AnaIO::true_beam_PDG == 211) {
    for (int i=0; i<11; ++i) {
      double g4rw_cex = 0;
      cex.push_back((*AnaIO::g4rw_full_grid_piplus_coeffs)[i]);
      if (cex[i].size() > 0) {
        double sum = 0; 
        for (size_t j = 0; j < cex[i].size(); ++j) {
          sum += cex[i][j];
          g4rw_cex += cex[i][j] * pow(weight[i], j);
        }
        
        if (abs(sum-1)>0.1) {
          cout<<"out"<<i<<"\t"<<sum<<endl;
          cout<<(*AnaIO::g4rw_full_grid_piplus_weights)[i][9]<<endl;
          cout<<AnaIO::run<<"\t"<<AnaIO::subrun<<"\t"<<AnaIO::event<<endl;
          g4rw_cex = 1;
        }
        g4rw *= g4rw_cex;
      }
    }
  }*/
  vector<double> cex_0_400 = (*AnaIO::g4rw_full_grid_piplus_coeffs)[4];
  vector<double> cex_400_800 = (*AnaIO::g4rw_full_grid_piplus_coeffs)[5];
  vector<double> cex_800_2000 = (*AnaIO::g4rw_full_grid_piplus_coeffs)[6];

  double w_cex_0_400 = 1.0;//w1; // KE 0 ~ 476.44931 MeV
  double w_cex_400_800 = 1.0;//w2; // KE 476.44931 ~ 1865.2940 MeV
  double w_cex_800_2000 = 1.0;//w2; // KE 476.44931 ~ 1865.2940 MeV
  
  double g4rw = 1;
  if (AnaIO::true_beam_PDG == 211 && (cex_0_400.size() > 0) && (cex_400_800.size() > 0) && (cex_800_2000.size() > 0)&& (w_cex_0_400 != 1 || w_cex_400_800 != 1 || w_cex_800_2000 != 1)){
    if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" && AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0){
    //if((*AnaIO::true_beam_endProcess) == "pi+Inelastic"){
    double g4rw_cex_0_400 = 0;
    for (size_t i = 0; i < cex_0_400.size(); ++i) {
      g4rw_cex_0_400 += cex_0_400[i] * pow(w_cex_0_400, i);
    }
    double g4rw_cex_400_800 = 0;
    for (size_t i = 0; i < cex_400_800.size(); ++i) {
      g4rw_cex_400_800 += cex_400_800[i] * pow(w_cex_400_800, i);
    }
    double g4rw_cex_800_2000 = 0;
    for (size_t i = 0; i < cex_800_2000.size(); ++i) {
      g4rw_cex_800_2000 += cex_800_2000[i] * pow(w_cex_800_2000, i);
    }
    g4rw *= g4rw_cex_0_400;
    g4rw *= g4rw_cex_400_800;
    g4rw *= g4rw_cex_800_2000;
    }
  }

  return g4rw;
}


double AnaUtils::GetUpStreamEnergyLoss(const bool & kMC, const double & InstE){
  
  double delta_E = 0.0;
  //if(kMC && InstE > 0){
  if(InstE > 0){
    /* --- tag one
    double E = InstE*1000.0;
    // Outliers 
    if(InstE > 0.9) E = 900.0;
    if(InstE < 0.7) E = 700.0;
    //double p0 = 126.0; double p1 = -0.401; double p2 = 0.000309;
    delta_E = GetLoss(E);

    if(InstE < 0.7) delta_E = 20.0;*/

    double E = InstE*1000.0;
    // Outliers 
    //if(InstE > 0.9) E = 975.0;
    if(InstE > 1.05) E = 1050;
    if(InstE < 0.7) E = 700;

    //double p0 = 126.0; double p1 = -0.401; double p2 = 0.000309;
    delta_E = GetLoss(E);

    //if(InstE < 0.7) delta_E = 13.0;

  }

  else{
    cout << "Waring InstE less than 0 !!" << endl; exit(1);
  }
  return delta_E;
}


void AnaUtils::FillEsliceHistograms(TH1D* hinit, TH1D* hend, TH1D* hinc, TH1D* hint, double KE_init, double KE_end, double KE_int, double weight, const vector<double> binning, int N_bin, bool fill_int)
{
  // == Fill KE_init, KE_inc, and KE_end distributions

  int i_init = -1;
  int i_end = -1;
  int i_int = -1;
  for(unsigned int i_bin = 0; i_bin < binning.size() - 1; i_bin++){
    if(KE_init > binning.at(i_bin) && KE_init < binning.at(i_bin + 1)) i_init = i_bin;
    if(KE_end > binning.at(i_bin) && KE_end < binning.at(i_bin + 1)) i_end = i_bin;
    if(KE_int > binning.at(i_bin) && KE_int < binning.at(i_bin + 1)) i_int = i_bin;
  }
  
  if(i_init == -1 || i_end == -1 || i_init == i_end) return;

  double bin_width = 50.0;
  double end_weight = (KE_end-binning.at(i_end))/bin_width;
  double ini_weight = (KE_init-binning.at(i_init))/bin_width;

  //if(end_weight < 0.8) return;
  
  if(fill_int){
    // == Fill KE_int
    //if(i_int > 0) AnaIO::hNewRecoInteractingHist->Fill(binning.at(i_int) + 0.1, weight);
    if(i_int > 0) hint->Fill(binning.at(i_int) + 0.1, weight);
    /*else {
      cout <<"fill_int but i_int < 0" << endl; 
      cout << "KE_init: " << KE_init << endl;
      cout << "KE_end: " << KE_end << endl;
      cout << "KE_int: " << KE_int << endl;

      //exit(1);
    }*/
  }
  else{
    // == Fill KE_init
    //AnaIO::hNewRecoInitialHist->Fill(binning.at(i_init - 1) + 0.1, weight);
    hinit->Fill(binning.at(i_init) + 0.1, weight);

    // == Fill KE_end
    //AnaIO::hNewRecoBeamInteractingHist->Fill(binning.at(i_end) + 0.1, weight);
    hend->Fill(binning.at(i_end) + 0.1, weight);

    if((*AnaIO::true_beam_endProcess) == "pi+Inelastic") AnaIO::hNewTruthInteractingHistTest->Fill(end_weight);
    if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" &&  AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0) AnaIO::hNewTruthCEXInteractingHistTest->Fill(end_weight);

    // == Fill KE_inc
    for(int i_bin = i_end; i_bin <= i_init; i_bin++){
      //AnaIO::hNewRecoIncidentHist->Fill(binning.at(i_bin) + 0.1, weight);
      if(i_bin == i_end) weight = weight*end_weight;
      if(i_bin == i_init) weight = weight*ini_weight;

      hinc->Fill(binning.at(i_bin) + 0.1, weight);

    }
  }

  return;
}

void AnaUtils::FillBeamQualityHist()
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

  AnaIO::hRecBeamInstX->Fill(AnaIO::beam_inst_X);
  AnaIO::hRecBeamInstY->Fill(AnaIO::beam_inst_Y);

}

void AnaUtils::SetBeamInstKEandFrontFaceKE(double &beam_inst_KE, double &true_ffKE, bool kFill)
{  
  int start_idx = -1;
  // Get reco beam inst KE
  beam_inst_KE = sqrt(pow(AnaIO::beam_inst_P,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
  // Get true front face KE (need to find the index when Z > 0)
  for (unsigned int i=0; i<AnaIO::true_beam_traj_Z->size(); i++){
    if ((*AnaIO::true_beam_traj_Z)[i] >= 0){
      start_idx = i-1; // the trajectory point before entering the TPC
      if (start_idx < 0) start_idx = -1;
      break;
    }
  }
  if(start_idx >= 0) true_ffKE = (*AnaIO::true_beam_traj_KE)[start_idx];
  if(AnaIO::true_beam_PDG == 211 && kFill) plotUtils.FillHist(AnaIO::hBeamInstVSTruthKEffNoCuts,beam_inst_KE*1000,true_ffKE);

}
void AnaUtils::FillUpStreamEnergyLossHistBeforeCut(double beam_inst_KE, double true_ffKE)
{
  if(AnaIO::true_beam_PDG == 211) plotUtils.FillHist(AnaIO::hBeamInstVSTruthKEffAfterCuts,beam_inst_KE*1000,true_ffKE);
    
  double UpStreamELoss = beam_inst_KE*1000-true_ffKE; 
  if(beam_inst_KE > 0.7 && beam_inst_KE < 0.8) AnaIO::hUpStreamELoss700MeV->Fill(UpStreamELoss);
  if(beam_inst_KE > 0.8 && beam_inst_KE < 0.9) AnaIO::hUpStreamELoss800MeV->Fill(UpStreamELoss);
  if(beam_inst_KE > 0.9 && beam_inst_KE < 1.0) AnaIO::hUpStreamELoss900MeV->Fill(UpStreamELoss);
  if(beam_inst_KE > 1.0 && beam_inst_KE < 1.1) AnaIO::hUpStreamELoss1000MeV->Fill(UpStreamELoss);

  double beam_inst_X = -999.0, beam_inst_Y = -999.0;
  int start_idx = -1;
  // Get reco front face KE (need to find the index when Z > 0)
  for (unsigned int i=0; i<AnaIO::true_beam_traj_Z->size(); i++){
    if ((*AnaIO::true_beam_traj_Z)[i] >= 0){
      start_idx = i-1; // the trajectory point before entering the TPC
      if (start_idx < 0) start_idx = -1;
      break;
    }
  }
  // Front face X and Y
  if(start_idx != -1){
    beam_inst_X = AnaIO::beam_inst_X;
    beam_inst_Y = AnaIO::beam_inst_Y;

    if(beam_inst_KE > 0.7 && beam_inst_KE < 0.8){
      if(UpStreamELoss < 56.1) plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYBeam700MeV,beam_inst_X,beam_inst_Y);
      else plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYScraper700MeV,beam_inst_X,beam_inst_Y);
    }
    if(beam_inst_KE > 0.8 && beam_inst_KE < 0.9){
      if(UpStreamELoss < 72.5) plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYBeam800MeV,beam_inst_X,beam_inst_Y);
      else plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYScraper800MeV,beam_inst_X,beam_inst_Y);
    }
    if(beam_inst_KE > 0.9 && beam_inst_KE < 1.0){
      if(UpStreamELoss < 88.0) plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYBeam900MeV,beam_inst_X,beam_inst_Y);
      else plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYScraper900MeV,beam_inst_X,beam_inst_Y);
    }
    if(beam_inst_KE > 1.0 && beam_inst_KE < 1.1){
      if(UpStreamELoss < 112.6) plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYBeam1000MeV,beam_inst_X,beam_inst_Y);
      else plotUtils.FillHist(AnaIO::hBeamInstXVSBeamInstYScraper1000MeV,beam_inst_X,beam_inst_Y);
    }
  }
}

void AnaUtils::FillUpStreamEnergyLossHistAfterCut(double beam_inst_KE, double true_ffKE)
{
  if(AnaIO::true_beam_PDG == 211) plotUtils.FillHist(AnaIO::hBeamInstVSTruthKEffAfterScraperCuts,beam_inst_KE*1000,true_ffKE);
  double UpStreamELoss = beam_inst_KE*1000-true_ffKE; 
  if(beam_inst_KE > 0.7 && beam_inst_KE < 0.8) AnaIO::hUpStreamELoss700MeVAfterScraperCuts->Fill(UpStreamELoss);
  if(beam_inst_KE > 0.8 && beam_inst_KE < 0.9) AnaIO::hUpStreamELoss800MeVAfterScraperCuts->Fill(UpStreamELoss);
  if(beam_inst_KE > 0.9 && beam_inst_KE < 1.0) AnaIO::hUpStreamELoss900MeVAfterScraperCuts->Fill(UpStreamELoss);
  if(beam_inst_KE > 1.0 && beam_inst_KE < 1.1) AnaIO::hUpStreamELoss1000MeVAfterScraperCuts->Fill(UpStreamELoss);
}

void AnaUtils::FillBeamVariablesAfterAllCuts(const bool kMC, int parType, int channelType)
{
  const double beamEndZ = AnaIO::reco_beam_calo_endZ;

  plotUtils.FillHist(AnaIO::hBeamInstXY, AnaIO::beam_inst_X, AnaIO::beam_inst_Y);
  if(parType == gkMisIDProton){
    plotUtils.FillHist(AnaIO::hBeamInstXY_misIDproton, AnaIO::beam_inst_X, AnaIO::beam_inst_Y);
  }
  if(parType == gkMisIDPion){
    plotUtils.FillHist(AnaIO::hBeamInstXY_misIDpion, AnaIO::beam_inst_X, AnaIO::beam_inst_Y);
  }

  plotUtils.FillHist(AnaIO::hBeamEndZ_Channels, beamEndZ, channelType);


  const double emScore = AnaIO::reco_beam_PFP_emScore_collection;
  const double emScore_wbc = AnaIO::reco_beam_PFP_emScore_collection_weight_by_charge;

  plotUtils.FillHist(AnaIO::hBeamemScore, emScore, parType);
  plotUtils.FillHist(AnaIO::hBeamemScore_wbc, emScore_wbc, parType);

  const double trackScore = AnaIO::reco_beam_PFP_trackScore_collection;
  const double trackScore_wbc = AnaIO::reco_beam_PFP_trackScore_collection_weight_by_charge;

  plotUtils.FillHist(AnaIO::hBeamtrackScore, trackScore, parType);
  plotUtils.FillHist(AnaIO::hBeamtrackScore_wbc, trackScore_wbc, parType);

  const int DaughterNumber = (*AnaIO::reco_daughter_PFP_ID).size();
  plotUtils.FillHist(AnaIO::hBeamDaughterNumbers, DaughterNumber, parType);

  // Beam start and end position study
  double startX = AnaIO::reco_beam_startX;
  double startCaloX = AnaIO::reco_beam_calo_startX;
  double startY = AnaIO::reco_beam_startY;
  double startCaloY = AnaIO::reco_beam_calo_startY;
  double startZ = AnaIO::reco_beam_startZ;
  double startCaloZ = AnaIO::reco_beam_calo_startZ;

  double endX = AnaIO::reco_beam_endX;
  double endCaloX = AnaIO::reco_beam_calo_endX;
  double endY = AnaIO::reco_beam_endY;
  double endCaloY = AnaIO::reco_beam_calo_endY;
  double endZ = AnaIO::reco_beam_endZ;
  double endCaloZ = AnaIO::reco_beam_calo_endZ;

  plotUtils.FillHist(AnaIO::hBeamStartX,startX,parType);
  plotUtils.FillHist(AnaIO::hBeamStartY,startY,parType);
  plotUtils.FillHist(AnaIO::hBeamStartZ,startZ,parType);
  plotUtils.FillHist(AnaIO::hBeamStartCaloX,startCaloX,parType);
  plotUtils.FillHist(AnaIO::hBeamStartCaloY,startCaloY,parType);
  plotUtils.FillHist(AnaIO::hBeamStartCaloZ,startCaloZ,parType);

  plotUtils.FillHist(AnaIO::hBeamEndX,endX,parType);
  plotUtils.FillHist(AnaIO::hBeamEndY,endY,parType);
  plotUtils.FillHist(AnaIO::hBeamEndZ,endZ,parType);
  plotUtils.FillHist(AnaIO::hBeamEndCaloX,endCaloX,parType);
  plotUtils.FillHist(AnaIO::hBeamEndCaloY,endCaloY,parType);
  plotUtils.FillHist(AnaIO::hBeamEndCaloZ,endCaloZ,parType);

  TVector3 BeamDir(AnaIO::reco_beam_calo_endX-AnaIO::reco_beam_calo_startX,
                  AnaIO::reco_beam_calo_endY-AnaIO::reco_beam_calo_startY,AnaIO::reco_beam_calo_endZ-AnaIO::reco_beam_calo_startZ);

  TVector3 BeamDir_Dir((*AnaIO::reco_beam_calo_endDirX)[1],(*AnaIO::reco_beam_calo_endDirY)[1],(*AnaIO::reco_beam_calo_endDirZ)[1]);

  double angleDiff = BeamDir.Angle(BeamDir_Dir)*TMath::RadToDeg();
  //cout << "angleDiff: " << angleDiff << endl;
  plotUtils.FillHist(AnaIO::hBeamDirDiff,angleDiff,parType);
  
  //if((*AnaIO::reco_beam_true_byHits_endProcess) == "Decay") plotUtils.FillHist(AnaIO::hBeamDirDiff,angleDiff,parType);
  
  // Track pitch 
  const vector<double> * trackPitch_SCE = AnaIO::reco_beam_TrkPitch_SCE;
  const vector<double> * trackPitch_NoSCE = AnaIO::reco_beam_TrkPitch_NoSCE;

  for(unsigned int ii = 0; ii < trackPitch_SCE->size(); ii++){
    double trkpit = (*trackPitch_SCE)[ii];
    AnaIO::hRecBeamTrackPitch_SCE->Fill(trkpit);

    double trkpit_noSCE = (*trackPitch_NoSCE)[ii];
    AnaIO::hRecBeamTrackPitch_NoSCE->Fill(trkpit_noSCE);
  }

  // dEdx 
  const vector<double> * dEdx_SCE = AnaIO::reco_beam_calibrated_dEdX_SCE;
  const vector<double> * dEdx_NoSCE = AnaIO::reco_beam_dEdX_NoSCE;

  for(unsigned int ii = 0; ii < dEdx_SCE->size(); ii++){

    double dEdx = (*dEdx_SCE)[ii];
    AnaIO::hRecBeamdEdx_SCE->Fill(dEdx);

    double dEdx_noSCE = (*dEdx_NoSCE)[ii];
    AnaIO::hRecBeamdEdx_NoSCE->Fill(dEdx_noSCE);

    double DeltaE = dEdx*(*trackPitch_SCE)[ii];
    double DeltaE_noSCE = dEdx_noSCE*(*trackPitch_NoSCE)[ii];

    AnaIO::hRecBeamDeltaE_SCE->Fill(DeltaE);
    AnaIO::hRecBeamDeltaE_NoSCE->Fill(DeltaE_noSCE);

  }

}

