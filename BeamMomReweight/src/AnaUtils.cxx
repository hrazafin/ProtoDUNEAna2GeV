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
      if(tmpp.Mag() > 0.15) nPionAboveThreshold++;
    }
    // PiMinus
    else if(itype == gkPiMinus){
      nPiMinus++;
      if(tmpp.Mag() > 0.15) nPionAboveThreshold++;
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
  }
  // At least two proton is found
  if(nProton>1){
    // Save info to subleading proton TLorentzVector
    pSecondaryProton.SetVectM(bufferProtonmom[subldProtonID], AnaFunctions::ProtonMass());
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
    // Fill all piplus momentum
    for(unsigned int ii = 0; ii < bufferPiPlusmom.size(); ii++){
      bufferPiPlusmom[subldPiPlusID];
    }
  }
  // At least two piplus is found
  if(nPiPlus>1){
    // Fill subleading piplus momentum
  }

  //======================== PiZero ========================
  int leadingPiZeroID = 0, subldPiZeroID = -999;
  if(nPiZero>1){
    // Fill the FS pi0 number (at least two pi0)
    int PiZerosortid[nPiZero];
    // Sort index according to it's momentum
    TMath::Sort(nPiZero, PiZeromom, PiZerosortid);
    // Save sorted index
    leadingPiZeroID = PiZerosortid[0];
    subldPiZeroID = PiZerosortid[1];
  }
  if(nPiZero>0){
    // Fill histogram for FS pi0 number
    // Save info to pi0 TLorentzVector
    pPiZero.SetVectM(bufferPiZeromom[leadingPiZeroID], AnaFunctions::PiZeroMass());
    bufferPiZeromom[subldPiZeroID];
  }
  if(nPiZero>1){
  }

  //======================== Gamma ========================
  if(nGamma>1){
    int Gammasortid[nGamma];
    TMath::Sort(nGamma, Gammamom, Gammasortid);
  }
  // Fill in MeV!
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

  double PiZeroGammamom[np];
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
      // Increase proton number 
      nPiZeroGamma++; 
    }
    
    else if(itype == gkElectron){
      bufferPiZeroElectronmom.push_back(tmpp);
      nPiZeroElectron++;
    }

    else if(itype == gkPositron){
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
    
  }
  // Rare decay - one gamma one electron and one positron
  else {
    LeadingGamma.SetVectM(bufferPiZeroGammamom[leadingGammaID], 0);
    Electron.SetVectM(bufferPiZeroElectronmom[0], AnaFunctions::ElectronMass());
    Positron.SetVectM(bufferPiZeroPositronmom[0], AnaFunctions::ElectronMass());
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
  // Get the size of shower array
  const int showerSize = showerArray.size();
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

    
    // TKI analysis pi0 vector
    RecPi0LTVet = PiZeroVec;
    // Save the truth type of two reco showers
    truthPi0Type = -999;
    // Get two truth-matched particles (maynot be photons indicated by the truthPi0Type)
    vector<TLorentzVector> TruthPi0Showers = GetTwoTruthMatchedPi0Showers(truthPi0Type,kFill);

    const double openingAngle = RecPi0Showers[0].Angle(RecPi0Showers[1].Vect())*TMath::RadToDeg();
    OA = openingAngle;

  }

  return PiZeroVec;
}

vector<TLorentzVector> AnaUtils::GetTwoPi0Showers(double &separation, bool kMC, bool kFill, bool kDoKF)
{
  vector<TLorentzVector> pi0Showers;

  const int showerSize = showerArray.size();
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

  return pi0ShowersTruthMatched;
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




void AnaUtils::KF(const TLorentzVector &ldShower, const TLorentzVector &slShower, vector<double> &FittedVars){

  const double openingAngleRad = ldShower.Angle(slShower.Vect());
 
  GoodFit = false;
  FittedVars = DoKF(ldShower.E(),slShower.E(),openingAngleRad,CVMind,GoodFit);
}

double AnaUtils::CalWeight(const bool & kMC){
  double weight = 1.;
  //return weight;
  
  double mufrac = 1.71;
  
  double mom_mu0 = 1.0033;
  double mom_sigma0 = 0.0609;
  double mom_mu = 1.01818;
  double mom_sigma = 0.07192;
  double wlimit = 3; 
  
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

  // With dep E Loss
  if (kMC){
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
  }

  return weight;
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




