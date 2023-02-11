#include "include/AnaIO.h"
#include "src/PlotUtils.cxx"
#include "src/AnaUtils.cxx"
#include "src/AnaCut.cxx"
#include "src/Unfold.cxx"

TRandom3 *grdm = new TRandom3(1); // fixed seed

int anaRec(const TString finName, TList *lout, const TString tag, const int nEntryToStop, BetheBloch & bb, Unfold & uf)
{
  //======================================== Basic Settings ======================================== 
  cout << "Input file:" << endl;
  gSystem->Exec(Form("readlink -f %s", finName.Data()));
  cout << endl;
  // Control MC or data sample
  bool kMC = true;
  if(!finName.Contains("_mc_"))  kMC = false;

  // Open the input file
  TFile * fin = new TFile(finName);
  if(!fin->IsOpen()){
    cout << "fin not open!" << endl;
    exit(1);
  }
  else cout << "fin is open!" << endl;

  // Get the TTree from input file
  TTree  * tree = AnaIO::GetInputTree(fin, "pduneana/beamana", tag);
  TTree  * tout = AnaIO::GetOutputTree(lout, tag);

  // Initialise reco histograms
  AnaIO::IniHist(lout, kMC);

  // Initialise class objects
  AnaUtils anaUtils;
  AnaCut anaCut;
  PlotUtils plotUtils;
 
  // Initialise entry and counters
  double ientry = 0;
  double BeamCount = 0;
  double CEXEventCount = 0;

  // Pion beam truth
  int true_avaPionBeam = 0;
  // Pion beam reco
  int recoMatch_avaPionBeam = 0, reco_selectedPionBeam = 0, reco_notPionBeam = 0, reco_missedBeam = 0;
  // Cex event truth
  int true_avaPionCEXevt = 0;
  // Cex event reco
  int recoMatch_avaPionCEXevt = 0, reco_selectedPionCEXevt = 0, reco_notPionCEXevt = 0, reco_missedevt = 0;
  // Differential cex truth
  int true_avaDiffCEXevt = 0;
  // Differential cex reco
  int recoMatch_avaDiffCEXevt = 0, reco_selectedDiffCEXevt = 0, reco_notDiffCEXevt = 0, reco_missedDiffevt = 0;

  // Control XS measurement
  bool doXS = true;

  // Loop over TTree
  while(tree->GetEntry(ientry)){

    // Break 
    if(nEntryToStop > 0 && ientry>=nEntryToStop){
      cout << "Break the loop after " << nEntryToStop << " entries!" << endl;
      break;
    }
    // update ientry after each loop
    ientry++;

    // Gaus smearing 
    double radGaus = 0.0;

    if(kMC && doXS){
      // Select true pion beam for the unfolding process
      if(anaCut.CutBeamAllInOne(kMC)){
        double trackLen = anaUtils.GetRecoTrackLength();
        // Only use the event where reco pion beam track length is physical
        if(trackLen != -999) {
          radGaus = grdm->Gaus(-0.00955164,0.0176993);
          //radGausffe = grdm->Gaus(0.00541435,0.0367387);
        }
      }
    }
    
    //====================== Extract truth information (MC only)======================//
    if(kMC){
      // Get true beam particle type from it's pdg code
      int TruthBeamType = anaUtils.GetParticleType(AnaIO::true_beam_PDG); 
      // Fill the TruthBeamType histogram
      AnaIO::hTruthBeamType->Fill(TruthBeamType);

      if(AnaIO::true_beam_PDG == 211 || AnaIO::true_beam_PDG == -13){
        const int pdg = AnaIO::reco_beam_true_byHits_PDG; 
        const int origin = AnaIO::reco_beam_true_byHits_origin; // 0: kUnknown 1: kBeamNeutrino 2: kCosmicRay 3: kSuperNovaNeutrino 4: kSingleParticle (single particles thrown at the detector) 
        const bool matched = AnaIO::reco_beam_true_byHits_matched; 
        const TString process = (*AnaIO::reco_beam_true_byHits_process); 
        if (process == "primary" && matched && origin == 4 && pdg == 211)  AnaIO::hBeamEndZ_TrueAvailable->Fill(AnaIO::true_beam_endZ);
      }
      
      // Fill XS hisotgrams
      anaUtils.FillXSTrueHistograms(true_avaPionBeam, true_avaPionCEXevt, true_avaDiffCEXevt);

      // Set signal using truth info
      anaUtils.SetFullSignal();
      AnaIO::hTruthSignal->Fill(AnaIO::Signal);
      // If this MC event is a signal 
      if(AnaIO::Signal){
        // Get FS particle type vector
        vector<int> SignalFSParticleType = anaUtils.GetBufferType();
        // Fill histograms
        for(int type : SignalFSParticleType){
          AnaIO::hTruthSignalFSParticleType->Fill(type);
        }
        AnaIO::hTruthSignalFSParticleNumber->Fill(SignalFSParticleType.size());

        const double beamEndZ = AnaIO::true_beam_endZ;
        //const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
        AnaIO::hBeamEndZ_ChannelsTrueSignal->Fill(beamEndZ);
        AnaIO::hBeamEndZ_ChargeExcChannel->Fill(beamEndZ);
        // Calcualte the TKI in truth level
        //anaUtils.DoTruthTKICalculation();
      }
    }
    //====================== End truth information (MC only)======================//

    //====================== MC reco information (MC cex and train response matrix)======================//

    // Thin-slice method, incident and interaction histograms calcualtion and filling
    if(kMC && doXS){
      // Divide the MC sample into half according to the event number
      bool isFakeData = anaUtils.IsFakeData();
      // Get the MC weight for each event
      double weight = anaUtils.CalWeight(kMC);
      //double bckweight = anaUtils.CalBckWeight(kMC);

      // Test truth matched beam is true pion beam (and yes it is)
      if(AnaIO::reco_beam_true_byHits_PDG == 211){
        AnaIO::hTruthPDGCode_PionMatched->Fill(AnaIO::true_beam_PDG);
      }

      // == Binnings for cross sections (SungBin's method)
      int N_binning_100MeV = 20;
      vector<double> binning_100MeV = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000.};
    
      // Select true pion beam for the unfolding process
      if(AnaIO::true_beam_PDG == 211){     
        // True beam incident energy (vector) and interaction energy calculation (use true_beam_traj_Z and true_beam_traj_KE)
        AnaIO::true_beam_incidentEnergies->clear();
        double interactingE = anaUtils.MakeTrueIncidentEnergies(AnaIO::true_beam_traj_Z, AnaIO::true_beam_traj_KE, AnaIO::true_beam_incidentEnergies);
        double LeadingPiZeroCosTheta = -999, LeadingPiZeroTheta = -999, LeadingPiZeroKE = -999, LeadingPiZeroP = -999;
        if(AnaIO::true_beam_incidentEnergies->size() != 0){
          // Get truth-matched beam particle types
          const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
          //const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
          // Get XS event types
          const int evtXStype = kMC ? anaUtils.GetFillXSEventType() : anaUtils.gkXSBmBkg;

          // Get the true FS particles (need pi0 for differential XS calcualtion)
          vector<TLorentzVector> vecFSParticle = anaUtils.GetFSParticlesTruth();
          // Get the FS pi0 particle's momentum
          LeadingPiZeroP = vecFSParticle[0].P();
          LeadingPiZeroKE = sqrt(LeadingPiZeroP*LeadingPiZeroP + pow(AnaFunctions::PiZeroMass(),2)) - AnaFunctions::PiZeroMass(); 
          LeadingPiZeroCosTheta = vecFSParticle[0].Pz()/vecFSParticle[0].P();
          LeadingPiZeroTheta = TMath::RadToDeg() * TMath::ACos(LeadingPiZeroCosTheta);

          // Get the reco track length from reco_beam_calo_X,Y,Z
          double trackLen = anaUtils.GetRecoTrackLength();
          // Get the reco track length at each point
          vector<double> trackLenAccum = anaUtils.GetRecoTrackLengthAccumVect();
          // Only use the event where reco pion beam track length is physical
          if(trackLen != -999) {
            // Count the fake data sample size
            if(isFakeData) recoMatch_avaPionBeam++;
            // Calculate the reco beam energy at the entry point
            //radGaus = grdm->Gaus(-0.00985,0.017756843);
            double beam_inst_P_smearing = AnaIO::beam_inst_P + radGaus; //grdm->Gaus(-0.00985,0.017756843);

            double beam_inst_KE = sqrt(pow(beam_inst_P_smearing,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
            // Taking into account the front face energy lost 12.74 MeV is used here
            //double Eloss = anaUtils.GetUpStreamEnergyLoss(kMC,beam_inst_KE);
            //double ff_energy_reco = beam_inst_KE*1000 - 12.74; 2.65229
            double ff_energy_reco = beam_inst_KE*1000 - 2.65229; //FIXME upper Syst

            double initialE_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[0]); 

            // Based on the reco beam accumulated length get the reco beam incident energy
            // Using BetheBloch formular to refill the AnaIO::reco_beam_incidentEnergies (analog to AnaIO::true_beam_traj_KE in true version)
            AnaIO::reco_beam_incidentEnergies->clear();
            for(unsigned int idx = 0; idx < trackLenAccum.size(); idx++){
              // Need first use BetheBloch formular to get reco beam incident energy using track length at each point
              double energy_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[idx]);
              // Fill the incident histogram
              AnaIO::reco_beam_incidentEnergies->push_back(energy_reco);
            }
            // Reco beam incident energy (vector) and interaction energy calculation (use reco_beam_calo_Z and reco_beam_incidentEnergies) 
            AnaIO::reco_beam_new_incidentEnergies->clear();
            // Slicing the beam into thiner slice using wire pitch spacing
            double interactingE_reco = anaUtils.MakeRecoIncidentEnergies(AnaIO::reco_beam_calo_Z, AnaIO::reco_beam_incidentEnergies, AnaIO::reco_beam_new_incidentEnergies);
            
            //double XSevtweight = 1.0;//anaUtils.CalXSEvtWeight(kMC, interactingE_reco, evtXStype);

            if(AnaIO::reco_beam_new_incidentEnergies->size() != 0){

              // Get the new reco beam incident energy vector size
              int size_reco = AnaIO::reco_beam_new_incidentEnergies->size();
              int size_true = AnaIO::true_beam_incidentEnergies->size();
              //int size_reco = AnaIO::reco_beam_new_incidentEnergies->size();
              //bool whichsize = size_true > size_reco;
              //int size = whichsize ? size_true : size_reco;

              int size_diff = size_reco - size_true;
              double interval_reco = (*AnaIO::reco_beam_new_incidentEnergies)[0] - interactingE_reco;
              double interval_true = (*AnaIO::true_beam_incidentEnergies)[1] - interactingE;

              double interval_diff = interval_reco - interval_true;

              AnaIO::hOldIncSizeDiff->Fill(size_diff);
              AnaIO::hNewIncIntervalDiff->Fill(interval_diff);

              // If it is not fake data
              if(!isFakeData){
                for(int i = 0; i < size_reco; i++){
                  // Fill the efficiency for incident energy vector
                  uf.eff_den_Inc->Fill((*AnaIO::true_beam_incidentEnergies)[i]);
                }
                uf.eff_den_Ini->Fill((*AnaIO::true_beam_incidentEnergies)[1]);
                if(interactingE!=-999) uf.eff_den_BeamInt->Fill(interactingE);
                else uf.eff_den_BeamInt->Fill(AnaIO::true_beam_incidentEnergies->back());
              }
                
              // Beam and event topology cuts and fill reco incident and interaction histograms
              // Do the beam cut and fill incident histogram
              if(anaCut.CutBeamAllInOne(kMC) && parType == anaUtils.gkBeamPiPlus){ // Pass the beam cuts
                // Fake data used to fill incident histogram
                if(isFakeData){
                  reco_selectedPionBeam++;
                  // Fill initial histogram
                  //AnaIO::hRecoInitialHist->Fill((*AnaIO::reco_beam_new_incidentEnergies)[0]);
                  AnaIO::hRecoInitialHist->Fill(initialE_reco, weight);
                  //cout << "(*AnaIO::reco_beam_new_incidentEnergies)[0]: " << (*AnaIO::reco_beam_new_incidentEnergies)[0] << " initialE_reco: " << initialE_reco << endl;
                  if(interactingE_reco!=-999) AnaIO::hRecoBeamInteractingHist->Fill(interactingE_reco, weight);
                  else AnaIO::hRecoBeamInteractingHist->Fill(AnaIO::reco_beam_new_incidentEnergies->back(), weight);

                  if(interactingE_reco!=-999) interactingE_reco = AnaIO::reco_beam_new_incidentEnergies->back();
                  // Sungbin's new method
                  anaUtils.FillEsliceHistograms(AnaIO::hNewRecoInitialHist, AnaIO::hNewRecoBeamInteractingHist, AnaIO::hNewRecoIncidentHist, AnaIO::hNewRecoInteractingHist, initialE_reco, interactingE_reco, -1, weight, binning_100MeV, N_binning_100MeV, false);
                  // Fill the reco incident histogram
                  for(int i = 0; i < size_reco; i++){
                    double energy_reco = (*AnaIO::reco_beam_new_incidentEnergies)[i];
                    //>!
                    AnaIO::hRecoIncidentHist->Fill(energy_reco);
                    if(parType == anaUtils.gkBeamPiPlus) AnaIO::hTruthMatchedIncidentHist->Fill(energy_reco);
                    if(parType != anaUtils.gkBeamPiPlus) AnaIO::hNonPionBeamIncidentHist->Fill(energy_reco);
                    AnaIO::hRecoBckSubIncidentHist->Fill(energy_reco);
                    //uf.response_SliceID_Inc.Fill((*AnaIO::reco_beam_new_incidentEnergies)[i],(*AnaIO::true_beam_incidentEnergies)[i]);
                  }
                }

                // Non-fake data used to train unfolding matrix
                else{
                  // Initial hist unfolding
                  //uf.response_SliceID_Ini.Fill((*AnaIO::reco_beam_new_incidentEnergies)[0],(*AnaIO::true_beam_incidentEnergies)[1]);
                  uf.response_SliceID_Ini.Fill(initialE_reco,(*AnaIO::true_beam_incidentEnergies)[1], weight);

                  if(interactingE_reco!=-999 && interactingE!=-999) uf.response_SliceID_BeamInt.Fill(interactingE_reco,interactingE, weight);
                  else if(interactingE_reco==-999 && interactingE!=-999) uf.response_SliceID_BeamInt.Fill(AnaIO::reco_beam_new_incidentEnergies->back(),interactingE, weight);
                  else if(interactingE_reco!=-999 && interactingE==-999) uf.response_SliceID_BeamInt.Fill(interactingE_reco,AnaIO::true_beam_incidentEnergies->back(), weight);
                  else uf.response_SliceID_BeamInt.Fill(AnaIO::reco_beam_new_incidentEnergies->back(),AnaIO::true_beam_incidentEnergies->back(), weight);

                  uf.eff_num_Ini->Fill((*AnaIO::true_beam_incidentEnergies)[1]);
                  if(interactingE!=-999) uf.eff_num_BeamInt->Fill(interactingE);
                  else uf.eff_num_BeamInt->Fill(AnaIO::true_beam_incidentEnergies->back());

                  /*for(int i = 0; i < size_reco; i++){
                    uf.response_SliceID_Inc.Fill((*AnaIO::reco_beam_new_incidentEnergies)[i],(*AnaIO::true_beam_incidentEnergies)[i]);
                    uf.eff_num_Inc->Fill((*AnaIO::true_beam_incidentEnergies)[i]);
                  }*/

                  if(size_diff > 0){ // Rec > Truth
                    for(int i = 0; i < size_reco; i++){
                      if(i < size_true){
                        uf.response_SliceID_Inc.Fill((*AnaIO::reco_beam_new_incidentEnergies)[i],(*AnaIO::true_beam_incidentEnergies)[i]);
                        uf.eff_num_Inc->Fill((*AnaIO::true_beam_incidentEnergies)[i]);
                        //cout << "(*AnaIO::true_beam_incidentEnergies)[i] inloop: " << (*AnaIO::true_beam_incidentEnergies)[i] << endl;
                        //cout << "(*AnaIO::reco_beam_new_incidentEnergies)[i] inloop: " << (*AnaIO::reco_beam_new_incidentEnergies)[i] << endl;

                      }
                      else{
                        //uf.response_SliceID_Inc.Fill((*AnaIO::reco_beam_new_incidentEnergies)[i],(*AnaIO::true_beam_incidentEnergies)[i]);
                        uf.response_SliceID_Inc.Fill((*AnaIO::reco_beam_new_incidentEnergies)[i],-999);

                        //cout << "size_reco: " << size_reco << " size_true: " << size_true << endl;
                        //cout << "i: " << i << endl;
                        //cout << "(*AnaIO::true_beam_incidentEnergies)[i]: " << (*AnaIO::true_beam_incidentEnergies)[i] << endl;
                        //cout << "(*AnaIO::reco_beam_new_incidentEnergies)[i]: " << (*AnaIO::reco_beam_new_incidentEnergies)[i] << endl;
                        
                      }
                    }
                  }
                  else{ // Truth > Rec
                    for(int i = 0; i < size_true; i++){
                      if(i < size_reco){
                        uf.response_SliceID_Inc.Fill((*AnaIO::reco_beam_new_incidentEnergies)[i],(*AnaIO::true_beam_incidentEnergies)[i]);
                        uf.eff_num_Inc->Fill((*AnaIO::true_beam_incidentEnergies)[i]);
                      }
                      else{
                        //uf.response_SliceID_Inc.Fill(-999,(*AnaIO::true_beam_incidentEnergies)[i]);
                        //uf.response_SliceID_Inc.Miss((*AnaIO::true_beam_incidentEnergies)[i]);
                      }
                    }
                  }
                  
                  /*
                  // Fill the reco incident histograms
                  for(int i = 0; i < size_reco; i++){
                    double energy_reco = (*AnaIO::reco_beam_new_incidentEnergies)[i];
                    AnaIO::hRecoIncidentHist->Fill(energy_reco);
                    if(parType == anaUtils.gkBeamPiPlus) AnaIO::hTruthMatchedIncidentHist->Fill(energy_reco);
                    if(parType != anaUtils.gkBeamPiPlus) AnaIO::hNonPionBeamIncidentHist->Fill(energy_reco);
                    AnaIO::hRecoBckSubIncidentHist->Fill(energy_reco);
                  }
                  */
                }

              }
              else{ // Not pass the beam cuts
                if(isFakeData) reco_notPionBeam++;
                if(!isFakeData){
                  uf.response_SliceID_Ini.Miss((*AnaIO::true_beam_incidentEnergies)[1], weight);
                  if(interactingE!=-999) uf.response_SliceID_BeamInt.Miss(interactingE, weight);
                  else uf.response_SliceID_BeamInt.Miss(AnaIO::true_beam_incidentEnergies->back(), weight);

                  for(unsigned int i = 0; i < AnaIO::true_beam_incidentEnergies->size(); i++){
                    uf.response_SliceID_Inc.Miss((*AnaIO::true_beam_incidentEnergies)[i]);
                  }
                }     
              }
              
              // Select pion inelastic events with signal topology on FS particles
              if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" &&  AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0){
                
                // Count fake data size for signal
                if(isFakeData) recoMatch_avaPionCEXevt++;
                // Non-fake data to fill efficiency for interacting energy
                if(!isFakeData) uf.eff_den_Int->Fill(interactingE);
                
                //if(interactingE > 750 && interactingE < 800) {
                if(interactingE > 650 && interactingE < 800) {
                  uf.eff_den_Pi0KE->Fill(LeadingPiZeroKE*1000);
                  uf.eff_den_Pi0CosTheta->Fill(LeadingPiZeroCosTheta);
                  uf.eff_den_Pi0Theta->Fill(LeadingPiZeroTheta);
                }

                // Declear pi0 KE
                double pi0KineticE, pi0costheta;
                // Do event topology cut and fill interacting histogram
                if(anaCut.CutBeamAllInOne(kMC) && anaCut.CutTopology(kMC, pi0KineticE, pi0costheta) && parType == anaUtils.gkBeamPiPlus && evtXStype == anaUtils.gkXSSignal) {  // Pass the event topology
                  //anaCut.CutTopology(kMC, pi0KineticE);
                  plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyBckSub, interactingE_reco, evtXStype, weight);
                  // Fake data used to fill interacting histogram
                  if(isFakeData){
                    reco_selectedPionCEXevt++;
                    AnaIO::hRecoInteractingHist->Fill(interactingE_reco, weight);
                    // Sungbin's new method
                    anaUtils.FillEsliceHistograms(AnaIO::hNewRecoInitialHist, AnaIO::hNewRecoBeamInteractingHist, AnaIO::hNewRecoIncidentHist, AnaIO::hNewRecoInteractingHist, initialE_reco, interactingE_reco, interactingE_reco, weight, binning_100MeV, N_binning_100MeV, true);

                    if(parType == anaUtils.gkBeamPiPlus) AnaIO::hTruthMatchedInteractingHist->Fill(interactingE_reco);
                    if(parType != anaUtils.gkBeamPiPlus) AnaIO::hNonPionBeamInteractingHist->Fill(interactingE_reco);
                    AnaIO::hRecoBckSubInteractingHist->Fill(interactingE_reco);
                    uf.response_SliceID_Int.Fill(interactingE_reco, interactingE, weight);
                  }
                  // Non-fake data used to train unfolding matrix
                  else {
                    uf.response_SliceID_Int.Fill(interactingE_reco, interactingE, weight);   
                    uf.eff_num_Int->Fill(interactingE);
                    //>!
                    AnaIO::hRecoInteractingHist->Fill(interactingE_reco, weight);
                    // Sungbin's new method
                    anaUtils.FillEsliceHistograms(AnaIO::hNewRecoInitialHist, AnaIO::hNewRecoBeamInteractingHist, AnaIO::hNewRecoIncidentHist, AnaIO::hNewRecoInteractingHist, initialE_reco, interactingE_reco, interactingE_reco, weight, binning_100MeV, N_binning_100MeV, true);

                  }

                  /*/ ??
                  if(!isFakeData){
                    //AnaIO::hRecoInteractingHist->Fill(interactingE_reco);
                    if(parType == anaUtils.gkBeamPiPlus) AnaIO::hTruthMatchedInteractingHist->Fill(interactingE_reco);
                    if(parType != anaUtils.gkBeamPiPlus) AnaIO::hNonPionBeamInteractingHist->Fill(interactingE_reco);
                    AnaIO::hRecoBckSubInteractingHist->Fill(interactingE_reco);
                  }
                  */
                  //if(anaCut.CutTopology(kMC, pi0KineticE) && evtXStype == anaUtils.gkXSSignal){
                    // ?? Indenpent of fake data use full MC to fill diff. XS unfolding matrix
                    //if(/*!isFakeData && */interactingE > 750 && interactingE < 800) {
                    if(/*!isFakeData && */interactingE > 650 && interactingE < 800) {
                      if(isFakeData) {recoMatch_avaDiffCEXevt++; reco_selectedDiffCEXevt++;}
                      uf.response_SliceID_Pi0KE.Fill(pi0KineticE*1000, LeadingPiZeroKE*1000, weight);
                      uf.response_SliceID_Pi0CosTheta.Fill(pi0costheta, LeadingPiZeroCosTheta, weight);
                      double pi0theta = TMath::RadToDeg() * TMath::ACos(pi0costheta);
                      uf.response_SliceID_Pi0Theta.Fill(pi0theta, LeadingPiZeroTheta, weight);

                      uf.eff_num_Pi0KE->Fill(LeadingPiZeroKE*1000);
                      uf.eff_num_Pi0CosTheta->Fill(LeadingPiZeroCosTheta);
                      uf.eff_num_Pi0Theta->Fill(LeadingPiZeroTheta);

                    //}
                    // ?? Fill diff. XS pi0 KE histogram
                      //if(interactingE_reco > 750 && interactingE_reco < 800) {
                        
                        plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_EVT_High, pi0KineticE*1000, evtXStype, weight);
                        AnaIO::hRecoPi0KEHist->Fill(pi0KineticE*1000, weight);
                        AnaIO::hRecoPi0CosThetaHist->Fill(pi0costheta, weight);
                        //double pi0theta = TMath::RadToDeg() * TMath::ACos(pi0costheta);
                        AnaIO::hRecoPi0ThetaHist->Fill(pi0theta, weight);

                      //}

                      
                      /*if(interactingE_reco > 700 && interactingE_reco < 750) {
                        plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_EVT_Medium, pi0KineticE*1000, evtXStype);
                      }
                      if(interactingE_reco > 650 && interactingE_reco < 700) {
                        plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_EVT_Low, pi0KineticE*1000, evtXStype);
                      }*/

                    }
                  //}

                    
                  
                }
                else{ // Not pass beam cuts or not pass event topology or both
                  if(isFakeData) reco_notPionCEXevt++;
                  // Fill missing response matrix
                  //if(!isFakeData) {
                    //if(interactingE != -999) 
                    uf.response_SliceID_Int.Miss(interactingE, weight);
                    //cout << "Not pass beam cuts interactingE: " << interactingE << endl;               
                  //}
                  //if(/*!isFakeData && */interactingE > 750 && interactingE < 800) {
                  if(/*!isFakeData && */interactingE > 650 && interactingE < 800) {
                    uf.response_SliceID_Pi0KE.Miss(LeadingPiZeroKE*1000, weight);
                    uf.response_SliceID_Pi0CosTheta.Miss(LeadingPiZeroCosTheta, weight);
                    uf.response_SliceID_Pi0Theta.Miss(LeadingPiZeroTheta, weight);
                    AnaIO::pi0theta->Fill(LeadingPiZeroTheta, weight);

                    //cout << "Not pass beam cuts LeadingPiZeroCosTheta: " << LeadingPiZeroCosTheta << endl;
                    //cout << "Not pass beam cuts interactingE: " << interactingE << endl;               
                    //cout << "Not pass beam cuts nPi0: " << AnaIO::true_daughter_nPi0 << endl;
                    //cout << "Not pass beam cuts nPi+: " << AnaIO::true_daughter_nPiPlus << endl;               
                    //cout << "Not pass beam cuts nPi-: " << AnaIO::true_daughter_nPiMinus << endl;  


                    if(isFakeData) {recoMatch_avaDiffCEXevt++;reco_notDiffCEXevt++;}
                  }
                }
              }
            }
            else{ // reco KE beam size = 0
              if(isFakeData) {reco_missedBeam++;}
            
              if(!isFakeData){
                uf.response_SliceID_Ini.Miss((*AnaIO::true_beam_incidentEnergies)[1], weight);
                if(interactingE!=-999) uf.response_SliceID_BeamInt.Miss(interactingE, weight);
                else uf.response_SliceID_BeamInt.Miss(AnaIO::true_beam_incidentEnergies->back(), weight);

                for(unsigned int i = 0; i < AnaIO::true_beam_incidentEnergies->size(); i++){
                  uf.response_SliceID_Inc.Miss((*AnaIO::true_beam_incidentEnergies)[i]);
                }
              }  
              if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" &&  AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0){

                //if(!isFakeData) {
                  //if(interactingE != -999) 
                  uf.response_SliceID_Int.Miss(interactingE, weight);
                  //cout << "Not 0 beam size interactingE: " << interactingE << endl;               
                //}
                
                if(isFakeData) {recoMatch_avaPionCEXevt++; reco_missedevt++;}
                //if(/*!isFakeData && */interactingE > 750 && interactingE < 800) {
                if(/*!isFakeData && */interactingE > 650 && interactingE < 800) {
                  uf.response_SliceID_Pi0KE.Miss(LeadingPiZeroKE*1000, weight);
                  uf.response_SliceID_Pi0CosTheta.Miss(LeadingPiZeroCosTheta, weight);
                  uf.response_SliceID_Pi0Theta.Miss(LeadingPiZeroTheta, weight);
                      AnaIO::pi0theta->Fill(LeadingPiZeroTheta, weight);

                  //cout << "Not 0 beam size LeadingPiZeroCosTheta: " << LeadingPiZeroCosTheta << endl;               
                  //cout << "Not 0 beam size interactingE: " << interactingE << endl;
                  //cout << "Not 0 beam size nPi0: " << AnaIO::true_daughter_nPi0 << endl;
                  //cout << "Not 0 beam size nPi+: " << AnaIO::true_daughter_nPiPlus << endl;               
                  //cout << "Not 0 beam size nPi-: " << AnaIO::true_daughter_nPiMinus << endl;               
                  
                  if(isFakeData) {recoMatch_avaDiffCEXevt++; reco_missedDiffevt++;}
                }
              }
            }
          } // End of if(trackLen != -999)...

          // Unphysical reco track length (also need to be used to fill the missing response matrix)
          else{
            if(isFakeData) {recoMatch_avaPionBeam++; reco_missedBeam++;}
            
            if(!isFakeData){
              uf.response_SliceID_Ini.Miss((*AnaIO::true_beam_incidentEnergies)[1], weight);
              if(interactingE!=-999) uf.response_SliceID_BeamInt.Miss(interactingE, weight);
              else uf.response_SliceID_BeamInt.Miss(AnaIO::true_beam_incidentEnergies->back(), weight);

              for(unsigned int i = 0; i < AnaIO::true_beam_incidentEnergies->size(); i++){
                uf.response_SliceID_Inc.Miss((*AnaIO::true_beam_incidentEnergies)[i]);
              }
            }  
            if((*AnaIO::true_beam_endProcess) == "pi+Inelastic" &&  AnaIO::true_daughter_nPi0 == 1 && AnaIO::true_daughter_nPiPlus == 0 &&  AnaIO::true_daughter_nPiMinus == 0){
              //if(!isFakeData) {
                //if(interactingE != -999) 
                uf.response_SliceID_Int.Miss(interactingE, weight);
                //cout << "Unphysical reco track length interactingE: " << interactingE << endl;               

              //}
              
              if(isFakeData) {recoMatch_avaPionCEXevt++; reco_missedevt++;}
              //if(/*!isFakeData && */interactingE > 750 && interactingE < 800) {
              if(/*!isFakeData && */interactingE > 650 && interactingE < 800) {
                uf.response_SliceID_Pi0KE.Miss(LeadingPiZeroKE*1000, weight);
                uf.response_SliceID_Pi0CosTheta.Miss(LeadingPiZeroCosTheta, weight);
                uf.response_SliceID_Pi0Theta.Miss(LeadingPiZeroTheta, weight);
                      AnaIO::pi0theta->Fill(LeadingPiZeroTheta, weight);

                //cout << "Unphysical reco track length LeadingPiZeroCosTheta: " << LeadingPiZeroCosTheta << endl;               
                //cout << "Unphysical reco track length interactingE: " << interactingE << endl;
                //cout << "Unphysical reco track length nPi0: " << AnaIO::true_daughter_nPi0 << endl;
                //cout << "Unphysical reco track length nPi+: " << AnaIO::true_daughter_nPiPlus << endl;               
                //cout << "Unphysical reco track length nPi-: " << AnaIO::true_daughter_nPiMinus << endl;     

                if(isFakeData) {recoMatch_avaDiffCEXevt++; reco_missedDiffevt++;}
              }
            }
          } // End of Unphysical reco track length ...
        }
      } // End of if(AnaIO::true_beam_PDG == 211)...

    } // End of if(kMC) ...

    //For data need to think more about the structure
    if(!kMC && doXS){

      //cout << "finName: " << finName << endl;
      /*for(int ii=0; ii<lout->GetSize(); ii++){
        const TString tag = lout->At(ii)->GetName();
        if(tag.Contains("Data")) cout << "data tag name: " << tag << endl;
      }*/
      // fake data test!
      //const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
      //const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
      //const int evtXStype = kMC ? anaUtils.GetFillXSEventType() : anaUtils.gkXSBmBkg;

      // Get the reco track length from reco_beam_calo_X,Y,Z
      double trackLen = anaUtils.GetRecoTrackLength();
      vector<double> trackLenAccum = anaUtils.GetRecoTrackLengthAccumVect();
      double pi0KineticE, pi0costheta;
      //cout << "data trackLen: " << trackLen << endl;
      if(trackLen != -999) {
        // Calculate the reco beam energy at the entry point 
        double beam_inst_KE = sqrt(pow(AnaIO::beam_inst_P,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
        //double ff_energy_reco = beam_inst_KE*1000 - 12.74;
        double ff_energy_reco = beam_inst_KE*1000 - 2.65229; //FIXME upper Syst; 

        double initialE_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[0]); 
        // Based on the reco beam accumulated length get the reco beam incident energy
        // Refill the AnaIO::reco_beam_incidentEnergies (analog to AnaIO::true_beam_traj_KE in true version)
        AnaIO::reco_beam_incidentEnergies->clear();
        for(unsigned int idx = 0; idx < trackLenAccum.size(); idx++){
          double energy_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[idx]);
          AnaIO::reco_beam_incidentEnergies->push_back(energy_reco);
        }
        // Reco beam incident energy (vector) and interaction energy calculation (use reco_beam_calo_Z and reco_beam_incidentEnergies) 
        AnaIO::reco_beam_new_incidentEnergies->clear();
        double interactingE_reco = anaUtils.MakeRecoIncidentEnergies(AnaIO::reco_beam_calo_Z, AnaIO::reco_beam_incidentEnergies, AnaIO::reco_beam_new_incidentEnergies);

        //int size_reco = AnaIO::reco_beam_new_incidentEnergies->size();

        if(anaCut.CutBeamAllInOne(false)){ // Pass the beam cuts
          // Fill the reco incident histograms
          //for(int i = 0; i < size_reco; i++){
          //  double energy_reco = (*AnaIO::reco_beam_new_incidentEnergies)[i];
          //  //if(parType == anaUtils.gkBeamPiPlus) 
          //  AnaIO::hRecoIncidentHistData->Fill(energy_reco,0.889673);
          //}
          //plotUtils.FillHist(AnaIO::hRecPiPlusEnergy_OVERLAY_After_EVTXS, interactingE_reco, parType);

          if(interactingE_reco==-999) interactingE_reco = AnaIO::reco_beam_new_incidentEnergies->back();
 
          //double iniwt = anaUtils.CalBeamIniWeight((*AnaIO::reco_beam_new_incidentEnergies)[0]);
          double iniwt = anaUtils.CalBeamIniWeight(initialE_reco);
          double intwt = anaUtils.CalBeamIntWeight(interactingE_reco);

          //cout << "interactingE_reco: " << interactingE_reco << "intwt: " << intwt << endl;
          //cout << "initialE_reco: " << (*AnaIO::reco_beam_new_incidentEnergies)[0] << "iniwt: " << iniwt << endl;

          //AnaIO::hRecoBeamInitialHistData->Fill((*AnaIO::reco_beam_new_incidentEnergies)[0],iniwt);
          AnaIO::hRecoBeamInitialHistData->Fill(initialE_reco,iniwt);
          AnaIO::hRecoBeamInteractingHistData->Fill(interactingE_reco,intwt);
        }

        if(anaCut.CutBeamAllInOne(false) && anaCut.CutTopology(false, pi0KineticE, pi0costheta)) {  // Pass the event topology
          //0.638783
          double intcexwt = anaUtils.CalCEXIntWeight(interactingE_reco);
          plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyBckSub, interactingE_reco, anaUtils.gkXSBmBkg, intcexwt);
          AnaIO::hRecoInteractingHistData->Fill(interactingE_reco,intcexwt);
          double pi0KEweight = 1.0;//anaUtils.CalCEXPi0KEWeight(pi0KineticE*1000);
          if(interactingE_reco > 650 && interactingE_reco < 800) {
            //if(evtXStype == anaUtils.gkXSSignal && parType == anaUtils.gkBeamPiPlus) 
            AnaIO::hRecoPi0KEHistData->Fill(pi0KineticE*1000, pi0KEweight);
            AnaIO::hRecoPi0CosThetaHistData->Fill(pi0costheta, pi0KEweight);
            double pi0theta = TMath::RadToDeg() * TMath::ACos(pi0costheta);
            AnaIO::hRecoPi0ThetaHistData->Fill(pi0theta, pi0KEweight);

          }
          //plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_EVTXS, pi0KineticE*1000, evtXStype);
        }

      }
    }
 

    //====================== Do Beam Cuts (both MC and data) ======================//
    // Do the beam cut and fill histogram
    if(!anaCut.CutBeamAllInOne(kMC, true)) continue;

    //====================== Selected Pion Beam Events Below ====================//
    
    // Get particle and event type
    const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
    //const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
    const int evtXStype = kMC ? anaUtils.GetFillXSEventType() : anaUtils.gkXSBmBkg;

    // Calculate event weight
    double weight = anaUtils.CalWeight(kMC);
    // No need to include bck weight for now (always = 1.0)
    double bckweight = anaUtils.CalBckWeight(kMC);

    // Count beam after beam cut before other cuts
    //BeamCount++;
    BeamCount += 1.0*weight*bckweight;

    // Get final state particles vector in this event and fill histogram
    anaUtils.GetFSParticlesTruth(true);
    // Get final state pi0 decay daughter in this event and fill histogram
    if(anaUtils.GetNParticles()[2] > 0) anaUtils.GetFSPiZeroDecayDaughterTruth(true);

    // Fill beam kinematics
    anaUtils.FillBeamKinematics(kMC);

    // Get the reco track length from reco_beam_calo_X,Y,Z
    double trackLen = anaUtils.GetRecoTrackLength();
    vector<double> trackLenAccum = anaUtils.GetRecoTrackLengthAccumVect();
    
    double interactingE_reco = -999;
    // Only use the event where reco pion beam track length is physical
    if(trackLen != -999) {
      double beam_inst_P_smearing = AnaIO::beam_inst_P;
      // Calculate the beam inst KE
      double beam_inst_KE = sqrt(pow(beam_inst_P_smearing,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
      if(kMC){
        // Smearing the MC beam momentum 
        beam_inst_P_smearing = AnaIO::beam_inst_P + radGaus;//grdm->Gaus(-0.00955164,0.0176993);//radGaus;//grdm->Gaus(-0.00985,0.017756843);
        // Get the beam inst (not smeared) and front-face KE (upstream energy loss study)
        double true_ffKE = -999;
        anaUtils.SetBeamInstKEandFrontFaceKE(beam_inst_KE,true_ffKE,false);
        // Need to reset beam_inst_KE to use the smeared beam momentum for MC
        beam_inst_KE = sqrt(pow(beam_inst_P_smearing,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
        // Calculate the upstream energy loss
        double UpStreamELoss = beam_inst_KE*1000 - true_ffKE;
        // Fill the upstream energy loss after smearing
        if(beam_inst_KE > 0.7 && beam_inst_KE < 0.75) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,0,weight);
        if(beam_inst_KE > 0.75 && beam_inst_KE < 0.8) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,1,weight);
        if(beam_inst_KE > 0.8 && beam_inst_KE < 0.85) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,2,weight);
        if(beam_inst_KE > 0.85 && beam_inst_KE < 0.9) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,3,weight);
        if(beam_inst_KE > 0.9 && beam_inst_KE < 0.95) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,4,weight);
        if(beam_inst_KE > 0.95 && beam_inst_KE < 1.0) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,5,weight);
        if(beam_inst_KE > 1.0 && beam_inst_KE < 1.05) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,6,weight);
        if(beam_inst_KE > 1.05 && beam_inst_KE < 1.1) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,7,weight);
      }
      // Get the energy dependent energy loss (not used for now)
      //double Eloss = anaUtils.GetUpStreamEnergyLoss(kMC, beam_inst_KE);
      // Substract the constant energy loss after smearing
      double ff_energy_reco = beam_inst_KE*1000 - 2.65229; // - 12.74; //FIXME upper Syst
      
      double initialE_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[0]); 

      AnaIO::reco_beam_incidentEnergies->clear();
      for(unsigned int idx = 0; idx < trackLenAccum.size(); idx++){
        double energy_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[idx]);
        //cout << "idx: " << idx << " energy_reco: " << energy_reco << endl;
        AnaIO::reco_beam_incidentEnergies->push_back(energy_reco);
      }
      // Reco beam incident energy (vector) and interaction energy calculation (use reco_beam_calo_Z and reco_beam_incidentEnergies) 
      AnaIO::reco_beam_new_incidentEnergies->clear();
      // Get the interacting energy
      interactingE_reco = anaUtils.MakeRecoIncidentEnergies(AnaIO::reco_beam_calo_Z, AnaIO::reco_beam_incidentEnergies, AnaIO::reco_beam_new_incidentEnergies);
      // If beam ends after APA3 use the last slide as the interacting energy
      if(interactingE_reco==-999) interactingE_reco = AnaIO::reco_beam_new_incidentEnergies->back();
      // Direct intial energy obtained from Jake's method
      double initialE_reco_Jake = (*AnaIO::reco_beam_new_incidentEnergies)[0];

      // Fill some histograms
      plotUtils.FillHist(AnaIO::hRecPiPlusEnergy_OVERLAY_After_EVTXS, interactingE_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergy, interactingE_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInstMomentum, beam_inst_P_smearing*1000, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusFrontFaceEnergy, ff_energy_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInitialEnergy, initialE_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInitialEnergyPosZCut, initialE_reco_Jake, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusTrackLength, trackLen, parType, weight*bckweight);
    } // End of if(trackLen != -999)

    //====================== Do Event Topology Cuts (both MC and data) ======================//

    // Declare pi0 related variables
    double pi0KineticE, pi0costheta;
    if(!anaCut.CutTopology(kMC,pi0KineticE,pi0costheta,true)) continue;
    // Count the number of CEX event after topology cuts
    CEXEventCount++;
    // Get the weight for evt bck reweight
    double XSevtweight = anaUtils.CalXSEvtWeight(kMC, interactingE_reco, evtXStype);
    // Get the signal fraction for data bck sub
    double intcexwt = anaUtils.CalCEXIntWeight(interactingE_reco);
    // Fill the beam interacting energy decomposed by event and particle type
    plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyEvt, interactingE_reco, evtXStype, weight*bckweight*XSevtweight);
    plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyPar, interactingE_reco, parType, weight*bckweight*XSevtweight);

    // Fill signal only MC beam int 
    if(kMC && evtXStype == anaUtils.gkXSSignal) plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyBckSubCheck, interactingE_reco, evtXStype, weight*bckweight);
    // Fill bck substracted data beam int
    if(!kMC) plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyBckSubCheck, interactingE_reco, evtXStype, intcexwt);
    // Fill the pi0 kinetic energy for all pion beams
    plotUtils.FillHist(AnaIO::hRecPiZeroKineticEnergyEvt, pi0KineticE*1000, evtXStype, XSevtweight*weight);
    // Fill the pi0 kinetic energy for selected pion int energy regions
    if(interactingE_reco> 750 && interactingE_reco < 800){
      plotUtils.FillHist(AnaIO::hRecPiZeroSliceKineticEnergyEvt, pi0KineticE*1000, evtXStype, weight);
      if(kMC) AnaIO::hRecPiZeroSliceKineticEnergyEvtMC->Fill(pi0KineticE*1000,weight);
      if(!kMC) AnaIO::hRecPiZeroSliceKineticEnergyEvtData->Fill(pi0KineticE*1000,weight);
    }
    if(interactingE_reco> 650 && interactingE_reco < 700) plotUtils.FillHist(AnaIO::hRecPiZeroSlice2KineticEnergyEvt, pi0KineticE*1000, evtXStype, weight);
    if(interactingE_reco> 800 && interactingE_reco < 850) plotUtils.FillHist(AnaIO::hRecPiZeroSlice3KineticEnergyEvt, pi0KineticE*1000, evtXStype, weight);
    if(interactingE_reco> 650 && interactingE_reco < 800) {
      plotUtils.FillHist(AnaIO::hRecPiZeroRangeKineticEnergyEvt, pi0KineticE*1000, evtXStype, XSevtweight*weight);
      plotUtils.FillHist(AnaIO::hRecPiZeroRangeCosThetaEvt, pi0costheta, evtXStype, XSevtweight*weight);
      double pi0theta = TMath::RadToDeg() * TMath::ACos(pi0costheta);
      plotUtils.FillHist(AnaIO::hRecPiZeroRangeThetaEvt, pi0theta, evtXStype, XSevtweight*weight);

    }
    if(interactingE_reco> 650 && interactingE_reco < 800) plotUtils.FillHist(AnaIO::hRecPiZeroRangeKineticEnergyEvtNoWeight, pi0KineticE*1000, evtXStype);
    if(interactingE_reco> 650 && interactingE_reco < 800) plotUtils.FillHist(AnaIO::hRecPiZeroRangeKineticEnergyEvtOneWeight, pi0KineticE*1000, evtXStype, XSevtweight);
    if(interactingE_reco> 650 && interactingE_reco < 800) {
      plotUtils.FillHist(AnaIO::hRecPiZeroRangeKineticEnergyEvtAnotherWeight, pi0KineticE*1000, evtXStype, weight);
      plotUtils.FillHist(AnaIO::hRecPiZeroRangeCosThetaEvtAnotherWeight, pi0costheta, evtXStype, weight);
      double pi0theta = TMath::RadToDeg() * TMath::ACos(pi0costheta);
      plotUtils.FillHist(AnaIO::hRecPiZeroRangeThetaEvtAnotherWeight, pi0theta, evtXStype, weight);

    }
    // Fill hist related to the effiency calculation
    const double beamEndZ = AnaIO::reco_beam_calo_endZ;
    const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
    plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_afterEvtTop, beamEndZ, channelType);
    
    // Do TKI calculation 
    //anaUtils.TruthMatchingTKI(anaUtils.RecPi0LTVet,anaUtils.RecProtonLTVet,anaUtils.TruthPi0LTVet,anaUtils.TruthProtonLTVet,kMC,anaUtils.GoodTruthMatch);
  
    // Fill output tree
    tout->Fill();
  } // End of while loop
/*
  if(kMC){
    TH1D * IncidentHist = plotUtils.GetIncidentHist(AnaIO::hTruthBeamInitialHist,AnaIO::hTruthBeamInteractingHist);
    cout << "ana IncidentHist integral: " << IncidentHist->Integral(0,100000) << endl;

    const Int_t x0_c = IncidentHist->GetXaxis()->GetFirst();
    const Int_t x1_c = IncidentHist->GetXaxis()->GetLast();

    for(Int_t ix=x0_c; ix<=x1_c; ix++){
      double bin = IncidentHist->GetBinContent(ix);
      double error = IncidentHist->GetBinError(ix);

      AnaIO::hTruthCalcIncidentHist->SetBinContent(ix,bin);
      AnaIO::hTruthCalcIncidentHist->SetBinError(ix,error);

      AnaIO::hTruthBeamCalcIncidentHist->SetBinContent(ix,bin);
      AnaIO::hTruthBeamCalcIncidentHist->SetBinError(ix,error);

    }

    TH1D * NewIncidentHist = plotUtils.GetIncidentHist(AnaIO::hNewTruthBeamInitialHist,AnaIO::hNewTruthBeamInteractingHist);
    cout << "ana NewIncidentHist integral: " << NewIncidentHist->Integral(0,100000) << endl;

    const Int_t x0_new = NewIncidentHist->GetXaxis()->GetFirst();
    const Int_t x1_new = NewIncidentHist->GetXaxis()->GetLast();

    for(Int_t ix=x0_new; ix<=x1_new; ix++){
      double bin = NewIncidentHist->GetBinContent(ix);
      double error = NewIncidentHist->GetBinError(ix);

      AnaIO::hNewTruthBeamCalcIncidentHist->SetBinContent(ix,bin);
      AnaIO::hNewTruthBeamCalcIncidentHist->SetBinError(ix,error);

    }

    TH1D * IncidentHist_50MeVbin = plotUtils.GetIncidentHist(AnaIO::hTruthBeamInitialHist_50MeVbin,AnaIO::hTruthBeamInteractingHist_50MeVbin);
    cout << "ana IncidentHist integral: " << IncidentHist_50MeVbin->Integral(0,100000) << endl;

    const Int_t x0_c_50MeVbin = IncidentHist_50MeVbin->GetXaxis()->GetFirst();
    const Int_t x1_c_50MeVbin = IncidentHist_50MeVbin->GetXaxis()->GetLast();

    for(Int_t ix=x0_c_50MeVbin; ix<=x1_c_50MeVbin; ix++){
      double bin = IncidentHist_50MeVbin->GetBinContent(ix);
      double error = IncidentHist_50MeVbin->GetBinError(ix);

      //AnaIO::hTruthCalcIncidentHist_50MeVbin->SetBinContent(ix,bin);
      //AnaIO::hTruthCalcIncidentHist_50MeVbin->SetBinError(ix,error);

      AnaIO::hTruthBeamCalcIncidentHist_50MeVbin->SetBinContent(ix,bin);
      AnaIO::hTruthBeamCalcIncidentHist_50MeVbin->SetBinError(ix,error);

    }

    cout << "int overall eff: " << uf.eff_num_Int->Integral(0,10000)/uf.eff_den_Int->Integral(0,10000) << endl;
    cout << "beam int overall eff: " << uf.eff_num_BeamInt->Integral(0,10000)/uf.eff_den_BeamInt->Integral(0,10000) << endl;
    cout << "ini overall eff: " << uf.eff_num_Ini->Integral(0,10000)/uf.eff_den_Ini->Integral(0,10000) << endl;
    
    cout << "Pi0KE overall eff: " << uf.eff_num_Pi0KE->Integral(0,10000)/uf.eff_den_Pi0KE->Integral(0,10000) << endl;
    cout << "Pi0CosTheta overall eff: " << uf.eff_num_Pi0CosTheta->Integral(0,10000)/uf.eff_den_Pi0CosTheta->Integral(0,10000) << endl;
    cout << "Pi0Theta overall eff: " << uf.eff_num_Pi0Theta->Integral(0,10000)/uf.eff_den_Pi0Theta->Integral(0,10000) << endl;

  }
*/
  // Print unfolding info
  cout << "--------------- Unfold beam initial ----------------" << endl;
  plotUtils.PrintXSUnfoldingInfo(true_avaPionBeam,recoMatch_avaPionBeam,reco_selectedPionBeam,reco_notPionBeam,reco_missedBeam);
  cout << "--------------- Unfold beam int ----------------" << endl;
  plotUtils.PrintXSUnfoldingInfo(true_avaPionCEXevt,recoMatch_avaPionCEXevt,reco_selectedPionCEXevt,reco_notPionCEXevt,reco_missedevt);
  cout << "--------------- Unfold daughter pi0 ----------------" << endl;
  plotUtils.PrintXSUnfoldingInfo(true_avaDiffCEXevt,recoMatch_avaDiffCEXevt,reco_selectedDiffCEXevt,reco_notDiffCEXevt,reco_missedDiffevt);
  
  uf.SaveHistograms();

  // Print info
  cout << "All entries: " << ientry << endl;
  cout << "BeamCount: " << BeamCount << endl;
  cout << "CEXEventCount: " << CEXEventCount << endl;
  // Kinematic Fitting for Pi0 shower
  //if(kMC) anaUtils.DoKinematicFitting();

  // Add response matrix and effeciency histograms to the list
  if(kMC && doXS) {
    lout->Add(uf.eff_Int);
    lout->Add(uf.eff_Inc);
    lout->Add(uf.eff_Ini);
    lout->Add(uf.eff_BeamInt);
    lout->Add(uf.eff_Pi0KE);
    lout->Add(uf.eff_Pi0CosTheta);
    lout->Add(uf.eff_Pi0Theta);

    lout->Add(uf.response_SliceID_Ini.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_BeamInt.HresponseNoOverflow());

    lout->Add(uf.response_SliceID_Inc.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Int.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Pi0KE.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Pi0CosTheta.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Pi0Theta.HresponseNoOverflow());

    lout->Add(&uf.response_SliceID_Ini);
    lout->Add(&uf.response_SliceID_BeamInt);
    lout->Add(&uf.response_SliceID_Inc);
    lout->Add(&uf.response_SliceID_Int);
    lout->Add(&uf.response_SliceID_Pi0KE);
    lout->Add(&uf.response_SliceID_Pi0CosTheta);
    lout->Add(&uf.response_SliceID_Pi0Theta);
  }
/*
  if(kMC && doXS) {

    TH1D *hini = (TH1D*)AnaIO::hRecoInitialHist->Clone("hini");
    TH1D *hbeamint = (TH1D*)AnaIO::hRecoBeamInteractingHist->Clone("hbeamint");

    TH1D *hinc = (TH1D*)AnaIO::hRecoIncidentHist->Clone("hinc");
    TH1D *hint = (TH1D*)AnaIO::hRecoInteractingHist->Clone("hint");

    TH1D *hpi0KE =(TH1D*)AnaIO::hRecoPi0KEHist->Clone("hpi0KE");
    TH1D *hpi0CosTheta =(TH1D*)AnaIO::hRecoPi0CosThetaHist->Clone("hpi0CosTheta");
    TH1D *hpi0Theta =(TH1D*)AnaIO::hRecoPi0ThetaHist->Clone("hpi0Theta");

    //uf.pur_Inc->SetMaximum(1.0);
    //lout->Add(uf.pur_Inc);
    //lout->Add(uf.pur_Int);
    lout->Add(uf.eff_Int);
    lout->Add(uf.eff_Inc);
    lout->Add(uf.eff_Ini);
    lout->Add(uf.eff_BeamInt);
    lout->Add(uf.eff_Pi0KE);
    lout->Add(uf.eff_Pi0CosTheta);
    lout->Add(uf.eff_Pi0Theta);

    
    //hinc->Multiply(uf.pur_Inc);
    //hint->Multiply(uf.pur_Int);

    //TH1D *hinc = (TH1D*)AnaIO::hTruthMatchedIncidentHist->Clone("hinc");  
    //TH1D *hint = (TH1D*)AnaIO::hTruthMatchedInteractingHist->Clone("hint");

    

    //TH1D *hback = (TH1D*)AnaIO::hNonPionBeamInteractingHist->Clone("hback");
    //TH1D *hback_inc = (TH1D*)AnaIO::hNonPionBeamIncidentHist->Clone("hback");


    //hint->Multiply(uf.pur_Int);
    //hint->Add(hback,-1);
    //hinc->Add(hback_inc,-1);

    //cout << "hint: " << hint->Integral(0,10000) << endl;
    //cout << "hint bin: " << hint->GetBinContent(10) << endl;

    //cout << "hback: " << hback->Integral(0,10000) << endl;
    //cout << "hback bin: " << hback->GetBinContent(10) << endl;

    RooUnfoldBayes unfold_Ini (&uf.response_SliceID_Ini, hini, 4);
    RooUnfoldBayes unfold_BeamInt (&uf.response_SliceID_BeamInt, hbeamint, 4);

    RooUnfoldBayes unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
    RooUnfoldBayes unfold_Int (&uf.response_SliceID_Int, hint, 4);

    RooUnfoldBayes unfold_Pi0KE (&uf.response_SliceID_Pi0KE, hpi0KE, 4);
    RooUnfoldBayes unfold_Pi0CosTheta (&uf.response_SliceID_Pi0CosTheta, hpi0CosTheta, 4);
    RooUnfoldBayes unfold_Pi0Theta (&uf.response_SliceID_Pi0Theta, hpi0Theta, 4);

    //uf.response_SliceID_Inc.SetNameTitleDefault("response_Inc");
    //uf.response_SliceID_Int.SetNameTitleDefault("response_Int");

    lout->Add(uf.response_SliceID_Ini.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_BeamInt.HresponseNoOverflow());

    lout->Add(uf.response_SliceID_Inc.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Int.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Pi0KE.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Pi0CosTheta.HresponseNoOverflow());
    lout->Add(uf.response_SliceID_Pi0Theta.HresponseNoOverflow());

    lout->Add(&uf.response_SliceID_Ini);
    lout->Add(&uf.response_SliceID_BeamInt);
    lout->Add(&uf.response_SliceID_Inc);
    lout->Add(&uf.response_SliceID_Int);
    lout->Add(&uf.response_SliceID_Pi0KE);
    lout->Add(&uf.response_SliceID_Pi0CosTheta);
    lout->Add(&uf.response_SliceID_Pi0Theta);


    TH1D *hini_uf = (TH1D*)unfold_Ini.Hreco();

    cout << "hini_uf: " << hini_uf->Integral(0,10000) << endl;
    const Int_t x0_ci = hini_uf->GetXaxis()->GetFirst();
    const Int_t x1_ci = hini_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0_ci; ix<=x1_ci; ix++){
      double bin = hini_uf->GetBinContent(ix);
      double error = hini_uf->GetBinError(ix);

      AnaIO::hUnFoldedInitialHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedInitialHist->SetBinError(ix,error);

    }

    TH1D *hbeamint_uf = (TH1D*)unfold_BeamInt.Hreco();

    cout << "hbeamint_uf: " << hbeamint_uf->Integral(0,10000) << endl;
    const Int_t x0_cib = hbeamint_uf->GetXaxis()->GetFirst();
    const Int_t x1_cib = hbeamint_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0_cib; ix<=x1_cib; ix++){
      double bin = hbeamint_uf->GetBinContent(ix);
      double error = hbeamint_uf->GetBinError(ix);

      AnaIO::hUnFoldedBeamInteractingHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedBeamInteractingHist->SetBinError(ix,error);

    }

    TH1D *hinc_uf = (TH1D*)unfold_Inc.Hreco();

    cout << "hinc_uf: " << hinc_uf->Integral(0,10000) << endl;
    const Int_t x0_c = hinc_uf->GetXaxis()->GetFirst();
    const Int_t x1_c = hinc_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0_c; ix<=x1_c; ix++){
      double bin = hinc_uf->GetBinContent(ix);
      double error = hinc_uf->GetBinError(ix);

      AnaIO::hUnFoldedIncidentHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedIncidentHist->SetBinError(ix,error);

    }

    //AnaIO::hUnFoldedIncidentHist->Add(hback_inc,-1);


    TH1D *hint_uf = (TH1D*)unfold_Int.Hreco();

    cout << "hint_uf: " << hint_uf->Integral(0,10000) << endl;
    const Int_t x0 = hint_uf->GetXaxis()->GetFirst();
    const Int_t x1 = hint_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0; ix<=x1; ix++){
      double bin = hint_uf->GetBinContent(ix);
      double error = hint_uf->GetBinError(ix);

      AnaIO::hUnFoldedInteractingHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedInteractingHist->SetBinError(ix,error);

    }
    //AnaIO::hUnFoldedInteractingHist->Add(hback,-1);

    TH1D *hpi0KE_uf = (TH1D*)unfold_Pi0KE.Hreco();
    cout << "hpi0KE_uf: " << hpi0KE_uf->Integral(0,10000) << endl;
    const Int_t x0_pi0 = hpi0KE_uf->GetXaxis()->GetFirst();
    const Int_t x1_pi0 = hpi0KE_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0_pi0; ix<=x1_pi0; ix++){
      double bin = hpi0KE_uf->GetBinContent(ix);
      double error = hpi0KE_uf->GetBinError(ix);

      AnaIO::hUnFoldedPi0KEHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedPi0KEHist->SetBinError(ix,error);

    }

    TH1D *hpi0CosTheta_uf = (TH1D*)unfold_Pi0CosTheta.Hreco();
    cout << "hpi0CosTheta_uf: " << hpi0CosTheta_uf->Integral(0,10000) << endl;
    const Int_t x0_pi0costheta = hpi0CosTheta_uf->GetXaxis()->GetFirst();
    const Int_t x1_pi0costheta = hpi0CosTheta_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0_pi0costheta; ix<=x1_pi0costheta; ix++){
      double bin = hpi0CosTheta_uf->GetBinContent(ix);
      double error = hpi0CosTheta_uf->GetBinError(ix);

      AnaIO::hUnFoldedPi0CosThetaHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedPi0CosThetaHist->SetBinError(ix,error);

    }


    TH1D *hpi0Theta_uf = (TH1D*)unfold_Pi0Theta.Hreco();
    cout << "hpi0Theta_uf: " << hpi0Theta_uf->Integral(0,10000) << endl;
    const Int_t x0_pi0theta = hpi0Theta_uf->GetXaxis()->GetFirst();
    const Int_t x1_pi0theta = hpi0Theta_uf->GetXaxis()->GetLast();

    for(Int_t ix=x0_pi0theta; ix<=x1_pi0theta; ix++){
      double bin = hpi0Theta_uf->GetBinContent(ix);
      double error = hpi0Theta_uf->GetBinError(ix);

      AnaIO::hUnFoldedPi0ThetaHist->SetBinContent(ix,bin);
      AnaIO::hUnFoldedPi0ThetaHist->SetBinError(ix,error);

    }



    TH1D * RecoIncidentHist = plotUtils.GetIncidentHist(AnaIO::hUnFoldedInitialHist,AnaIO::hUnFoldedBeamInteractingHist);
    cout << "ana RecoIncidentHist integral: " << RecoIncidentHist->Integral(0,100000) << endl;

    const Int_t x0_co = RecoIncidentHist->GetXaxis()->GetFirst();
    const Int_t x1_co = RecoIncidentHist->GetXaxis()->GetLast();

    for(Int_t ix=x0_co; ix<=x1_co; ix++){
      double bin = RecoIncidentHist->GetBinContent(ix);
      double error = RecoIncidentHist->GetBinError(ix);

      AnaIO::hUnfoldedBeamIncidentHist->SetBinContent(ix,bin);
      AnaIO::hUnfoldedBeamIncidentHist->SetBinError(ix,error);

    }


  } // End of if(kMC && doXS) ...

  if(!kMC && doXS){
    
    TH1D *hiniData = (TH1D*)AnaIO::hRecoBeamInitialHistData->Clone("hiniData");
    TH1D *hintData = (TH1D*)AnaIO::hRecoInteractingHistData->Clone("hintData");
    TH1D *hbeamintData = (TH1D*)AnaIO::hRecoBeamInteractingHistData->Clone("hbeamintData");

    TH1D *hpi0KEData =(TH1D*)AnaIO::hRecoPi0KEHistData->Clone("hpi0KEData");

    
    RooUnfoldBayes unfold_IniData (&uf.response_SliceID_Ini, hiniData, 4);
    RooUnfoldBayes unfold_BeamIntData (&uf.response_SliceID_BeamInt, hbeamintData, 4);
    RooUnfoldBayes unfold_IntData (&uf.response_SliceID_Int, hintData, 4);
    RooUnfoldBayes unfold_Pi0KEData (&uf.response_SliceID_Pi0KE, hpi0KEData, 4);

    
    // Data
    TH1D *hini_ufData = (TH1D*)unfold_IniData.Hreco();

    cout << "hini_ufData: " << hini_ufData->Integral(0,10000) << endl;
    const Int_t x0_cData = hini_ufData->GetXaxis()->GetFirst();
    const Int_t x1_cData = hini_ufData->GetXaxis()->GetLast();

    for(Int_t ix=x0_cData; ix<=x1_cData; ix++){
      double bin = hini_ufData->GetBinContent(ix);
      double error = hini_ufData->GetBinError(ix);

      AnaIO::hUnFoldedBeamInitialHistData->SetBinContent(ix,bin);
      AnaIO::hUnFoldedBeamInitialHistData->SetBinError(ix,error);

    }

    //AnaIO::hUnFoldedIncidentHist->Add(hback_inc,-1);


    TH1D *hbeamint_ufData = (TH1D*)unfold_BeamIntData.Hreco();

    cout << "hbeamint_ufData: " << hbeamint_ufData->Integral(0,10000) << endl;
    const Int_t x0beamData = hbeamint_ufData->GetXaxis()->GetFirst();
    const Int_t x1beamData = hbeamint_ufData->GetXaxis()->GetLast();

    for(Int_t ix=x0beamData; ix<=x1beamData; ix++){
      double bin = hbeamint_ufData->GetBinContent(ix);
      double error = hbeamint_ufData->GetBinError(ix);

      AnaIO::hUnFoldedBeamInteractingHistData->SetBinContent(ix,bin);
      AnaIO::hUnFoldedBeamInteractingHistData->SetBinError(ix,error);

    }

    TH1D *hint_ufData = (TH1D*)unfold_IntData.Hreco();

    cout << "hint_ufData: " << hint_ufData->Integral(0,10000) << endl;
    const Int_t x0Data = hint_ufData->GetXaxis()->GetFirst();
    const Int_t x1Data = hint_ufData->GetXaxis()->GetLast();

    for(Int_t ix=x0Data; ix<=x1Data; ix++){
      double bin = hint_ufData->GetBinContent(ix);
      double error = hint_ufData->GetBinError(ix);

      AnaIO::hUnFoldedInteractingHistData->SetBinContent(ix,bin);
      AnaIO::hUnFoldedInteractingHistData->SetBinError(ix,error);
    }
    

    TH1D *hpi0KE_ufData = (TH1D*)unfold_Pi0KEData.Hreco();
    cout << "hpi0KE_ufData: " << hpi0KE_ufData->Integral(0,10000) << endl;
    const Int_t x0_pi0Data = hpi0KE_ufData->GetXaxis()->GetFirst();
    const Int_t x1_pi0Data = hpi0KE_ufData->GetXaxis()->GetLast();

    for(Int_t ix=x0_pi0Data; ix<=x1_pi0Data; ix++){
      double bin = hpi0KE_ufData->GetBinContent(ix);
      double error = hpi0KE_ufData->GetBinError(ix);

      AnaIO::hUnFoldedPi0KEHistData->SetBinContent(ix,bin);
      AnaIO::hUnFoldedPi0KEHistData->SetBinError(ix,error);

    }
    

    TH1D * RecoIncidentHistData = plotUtils.GetIncidentHist(AnaIO::hUnFoldedBeamInitialHistData,AnaIO::hUnFoldedBeamInteractingHistData);
    cout << "ana RecoIncidentHistData integral: " << RecoIncidentHistData->Integral(0,100000) << endl;

    const Int_t x0_co = RecoIncidentHistData->GetXaxis()->GetFirst();
    const Int_t x1_co = RecoIncidentHistData->GetXaxis()->GetLast();

    for(Int_t ix=x0_co; ix<=x1_co; ix++){
      double bin = RecoIncidentHistData->GetBinContent(ix);
      double error = RecoIncidentHistData->GetBinError(ix);

      AnaIO::hUnfoldedBeamIncidentHistData->SetBinContent(ix,bin);
      AnaIO::hUnfoldedBeamIncidentHistData->SetBinError(ix,error);

    }

  }
*/

  // Print cut flow statistics
  int icut = 0;
  double nsel = -999;

  // Beam cuts
  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam PDG",  icut++), AnaIO::hCutBeamPDGPass, 1, 1, ientry);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Pandora slice",  icut++), AnaIO::hCutPandoraSlicePass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Calo size",  icut++), AnaIO::hCutCaloSizePass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam Quality",  icut++), AnaIO::hCutBeamQualityPass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. APA3 endZ",  icut++), AnaIO::hCutAPA3EndZPass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Michel score",  icut++), AnaIO::hCutMichelScorePass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Chi2/DOF",  icut++), AnaIO::hCutProtonChi2Pass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam Scraper",  icut++), AnaIO::hCutBeamScraperPass, 1, 1, nsel);

  // FS particle cuts
  //nsel = plotUtils.PrintStat(tag+Form(" %d. Nproton",  icut++), AnaIO::hCutnproton, 1, 1, nsel); // TKI proton selection
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nshower",  icut++), AnaIO::hCutnshower, 2, 100000, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npiplus",  icut++), AnaIO::hCutnpiplus, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nmichel",  icut++), AnaIO::hCutnmichel, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npi0",  icut++), AnaIO::hCutnpi0, 1, 1, nsel);
  //nsel = plotUtils.PrintStat(tag+Form(" %d. KF Pass",  icut++), AnaIO::hCutKFPass, 1, 1, nsel);

  cout << "Proton Cuts: " << endl;
  int icut_proton = 0;
  double nsel_proton = -999;
  //nsel_proton = plotUtils.PrintStat(tag+Form(" %d. proton Pandora",  icut_proton++), AnaIO::hCutDaughterPandoraprotonPass, 1, 1, nsel_proton);
  nsel_proton = plotUtils.PrintStat(tag+Form(" %d. proton track score",  icut_proton++), AnaIO::hCutDaughterTrackScorePass, 1, 1, nsel_proton);
  nsel_proton = plotUtils.PrintStat(tag+Form(" %d. proton nhits",  icut_proton++), AnaIO::hCutDaughterTracknHitsPass, 1, 1, nsel_proton);
  nsel_proton = plotUtils.PrintStat(tag+Form(" %d. proton SubPID",  icut_proton++), AnaIO::hCutDaughterProtonSubPIDPass, 1, 1, nsel_proton);
  
  cout << "Shower Cuts: " << endl;
  int icut_shower = 0;
  double nsel_shower = -999;
  //nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower Pandora",  icut_shower++), AnaIO::hCutDaughterPandoraShowerPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower EM score",  icut_shower++), AnaIO::hCutDaughterShowerScorePass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower nhits",  icut_shower++), AnaIO::hCutDaughterShowernHitsPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower startZ",  icut_shower++), AnaIO::hCutDaughterShowerStartZPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower non-empty E",  icut_shower++), AnaIO::hCutDaughterShowerNonEmptyEPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower distance",  icut_shower++), AnaIO::hCutDaughterShowerDistPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower IP",  icut_shower++), AnaIO::hCutDaughterShowerIPPass, 1, 1, nsel_shower);

  cout << "Pi0 Cuts: " << endl;
  int icut_pi0 = 0;
  double nsel_pi0 = -999;
  nsel_pi0 = plotUtils.PrintStat(tag+Form(" %d. Pi0 No cuts",  icut_pi0++), AnaIO::hCutDaughterPi0NOCutsPass, 1, 1, nsel_pi0);
  nsel_pi0 = plotUtils.PrintStat(tag+Form(" %d. Pi0 Mass",  icut_pi0++), AnaIO::hCutDaughterPi0MassPass, 1, 1, nsel_pi0);
  nsel_pi0 = plotUtils.PrintStat(tag+Form(" %d. Pi0 OA",  icut_pi0++), AnaIO::hCutDaughterPi0OAPass, 1, 1, nsel_pi0);


  cout << "\n------------- Beam Purity and Efficiency--------------------" << endl;
  int icut_beam_pe = 0;
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam PDG",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutPDG);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam Pandora",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutPandora);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam CaloSize",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutCaloSize);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam Quality",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutBeamQuality);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam APA3",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutAPA3);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam Michel",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutMichelScore);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam Chi2",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutChi2DOF);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Beam Scraper",  icut_beam_pe++), AnaIO::hBeamEndZ_TrueAvailable, AnaIO::hBeamEndZ_CutBeamScraper);
  
  cout << "\n------------- Proton Purity and Efficiency--------------------" << endl;
  int icut_proton_pe = 0;
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Proton track score",  icut_proton_pe++), AnaIO::hTruthProtonP, AnaIO::hProtonMom_CutTrackScore);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Proton nhits",  icut_proton_pe++), AnaIO::hTruthProtonP, AnaIO::hProtonMom_CutnHits);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Proton SubPID",  icut_proton_pe++), AnaIO::hTruthProtonP, AnaIO::hProtonMom_CutSubPID);

  cout << "\n------------- Shower Purity and Efficiency--------------------" << endl;
  int icut_shower_pe = 0;
  if(kMC) plotUtils.PrintShowerPurityandEff(tag+Form(" %d. Shower EM score",  icut_shower_pe++), AnaIO::hTruthLeadingPi0GammaP,AnaIO::hTruthSubLeadingPi0GammaP, AnaIO::hShowerEnergy_CutEMscore);
  if(kMC) plotUtils.PrintShowerPurityandEff(tag+Form(" %d. Shower nhits",  icut_shower_pe++), AnaIO::hTruthLeadingPi0GammaP,AnaIO::hTruthSubLeadingPi0GammaP, AnaIO::hShowerEnergy_CutnHits);
  if(kMC) plotUtils.PrintShowerPurityandEff(tag+Form(" %d. Shower startZ",  icut_shower_pe++), AnaIO::hTruthLeadingPi0GammaP,AnaIO::hTruthSubLeadingPi0GammaP, AnaIO::hShowerEnergy_CutStartZ);
  if(kMC) plotUtils.PrintShowerPurityandEff(tag+Form(" %d. Shower distance",  icut_shower_pe++), AnaIO::hTruthLeadingPi0GammaP,AnaIO::hTruthSubLeadingPi0GammaP, AnaIO::hShowerEnergy_CutDistance);
  if(kMC) plotUtils.PrintShowerPurityandEff(tag+Form(" %d. Shower IP",  icut_shower_pe++), AnaIO::hTruthLeadingPi0GammaP,AnaIO::hTruthSubLeadingPi0GammaP, AnaIO::hShowerEnergy_CutIP);

  cout << "\n------------- Pi0 Purity and Efficiency--------------------" << endl;
  int icut_pi0_pe = 0;
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. No Cut",  icut_pi0_pe++), AnaIO::hTruthLeadingPiZeroE, AnaIO::hPi0Energy_NoCut);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Pi0 Mass",  icut_pi0_pe++), AnaIO::hTruthLeadingPiZeroE, AnaIO::hPi0Energy_CutMass);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Pi0 OA",  icut_pi0_pe++), AnaIO::hTruthLeadingPiZeroE, AnaIO::hPi0Energy_CutOA);

  cout << "\n------------- Exclusive Channel Purity and Efficiency (1p0n only)--------------------" << endl;
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Exclusive 1p0n",  0), AnaIO::hTruthDalphat1p0n, AnaIO::hRecdalphat);

  cout << "\n------------- Exclusive Channel Purity and Efficiency (Signal)--------------------" << endl;
  if(kMC) plotUtils.PrintExcPurityandEff(tag+Form(" %d. Exclusive signal",  0), AnaIO::hTruthDalphat1p0n, AnaIO::hTruthDalphat1pMn, AnaIO::hTruthDalphatNp0n, AnaIO::hTruthDalphatNpMn, AnaIO::hRecdalphat);

  cout << "\n------------- Charge Exchange Channel Purity and Efficiency --------------------" << endl;
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. Charge Exchange",  0), AnaIO::hBeamEndZ_ChannelsTrueSignal, AnaIO::hBeamEndZ_Channels_afterEvtTop);

  cout << "\n------------- Charge Exchange Channel Purity and Efficiency new Def--------------------" << endl;
  int icut_chex_pe = 0;
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam PDG",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamPDG);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam Pandora",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamPandora);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam Calo",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamCalo);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam Qaulity",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamQuality);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam APA3",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamAPA3);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam Michel",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamMichel);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam Chi2",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamChi2);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Beam Scraper",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_BeamScraper);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx Two Showers",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_TwoShowers);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx No PiPlus",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_NoPiPlus);
  if(kMC) plotUtils.PrintPi0PurityandEff(tag+Form(" %d. ChEx One Pi0",  icut_chex_pe++), AnaIO::hBeamEndZ_ChargeExcChannel, AnaIO::hBeamEndZ_Channels_OnePi0);



  printf("End of %d cuts: %.1f selected\n", icut, nsel);

  printf("End of %d shower cuts: %.1f selected\n", icut_shower, nsel_shower);

  cout << "beam count: " << BeamCount << endl;
  cout << "1p0n: " << AnaIO::hTruthDalphat1p0n->Integral(0,10000)/BeamCount << endl;
  cout << "total: " << (AnaIO::hTruthDalphat1p0n->Integral(0,10000) + AnaIO::hTruthDalphat1pMn->Integral(0,10000) + AnaIO::hTruthDalphatNp0n->Integral(0,10000) +AnaIO::hTruthDalphatNpMn->Integral(0,10000) )/BeamCount << endl;
/*
  // Print signal/background info MC
  const double nsig = AnaIO::hTruthSignal->GetBinContent(2);
  const double nbk = AnaIO::hTruthSignal->GetBinContent(1);
  const double nall = nsig+nbk;
  const double purity = nsig/nall;

  cout << "nAll: " << nall << " nSignal: " << nsig << " nBackground: " << nbk << endl;
  cout << "purity: " << purity << endl;
*/
  const double ntotal = AnaIO::hPi0Total->Integral(0,10000);
  const double nselected = AnaIO::hPi0Selected->Integral(0,10000);
  const double neff = AnaIO::hTruthLeadingPiZeroE->Integral(0,10000);
  
  const double eff = nselected/neff;
  const double purity = nselected/ntotal;

  cout << "nAll: " << ntotal << " nSignal: " << nselected << endl;
  cout << "purity: " << purity << "eff: " << eff << endl;

  return BeamCount;
  
} // End of anaRec

int main(int argc, char * argv[])
{
  // Initialise the input file name
  const TString mcfinName = "input/protoDUNE_mc_reco_flattree_prod4a_ntuple.root";
  const TString datafinName  = "input/protoDUNE_data_reco_flattree_prod4_ntuple.root";

  // Declare output list
  TList * mclout = 0x0;
  TList * datalout = 0x0;

  mclout = new TList;
  datalout = new TList;

  //mclout->Add(g_cex_600MeV);
  
  // Run reco loop for both mc and data
  int nEntryToStop = -1; // Default is looping over all entries
  // Get the input break entries if provided 
  if(argc!=1) nEntryToStop = atoi(argv[1]);
  // bb and uf
  BetheBloch bb(211);
  Unfold uf(20, 0, 1000);
  // MC Analysis
  double mcBeamCount = anaRec(mcfinName,mclout,"mc", nEntryToStop, bb, uf);
  cout << "test intgral mc: " <<  AnaIO::hRecoBeamInitialHistData->Integral(0,10000) << endl;
  // Data Analysis
  double dataBeamCount = anaRec(datafinName,datalout,"data", nEntryToStop, bb, uf);
  cout << "test intgral data: " <<  AnaIO::hRecoBeamInitialHistData->Integral(0,10000) << endl;
  
  // Process histograms
  PlotUtils plotUtils;
  plotUtils.ProcessHist(mclout,true);
  plotUtils.ProcessHist(datalout,false);  

  // Declare output root file
  TFile * fout = new TFile("output/outana.root","recreate");
  // Create mc subdirectory
  TDirectory * ld_mc = gDirectory->mkdir("mc");
  ld_mc->cd();
  // Write list to the root file
  mclout->Write();
  gDirectory->cd("../");
  // Create data subdirectory
  TDirectory * ld_data = gDirectory->mkdir("data");
  ld_data->cd();
  // Write list to the root file
  datalout->Write();
  gDirectory->cd("../");
  // Save the info
  fout->Save();
  fout->Close();
  // Calculate the Data/MC normalisation constant
  double plotScale = dataBeamCount/mcBeamCount;
  cout << "plotScale: " << plotScale << endl;
  // Draw all histograms
  plotUtils.DrawHist(mclout,plotScale,datalout,"output");

}
