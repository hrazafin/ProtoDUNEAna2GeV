#include "include/AnaIO.h"
#include "src/PlotUtils.cxx"
#include "src/AnaUtils.cxx"
#include "src/AnaCut.cxx"
//#include "src/BetheBloch.cxx"
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

  double mc_i = 0;
  double data_i = 0;

  
  int xloop = 0;
  int xloopin = 0;
  int xloopout = 0;
  int xloopoutbadtrack = 0;
  int xcexloop = 0;

  int trueloop = 0;
  int truecexloop = 0;

  int truecex800loop = 0;
  int cex800loopin = 0;
  int cex800loopout = 0;

  int cex800loopfill = 0;
  int cex800loopfillnonzero = 0;

  double beampiplus = 0;
  double totalbeam = 0;

  double cex = 0;
  double totalcex = 0;

  double pi0cex = 0;
  double pi0totalcex = 0;

  double badevt = 0;
  double goodevt = 0;

  bool doit = false;

  // Loop over TTree
  while(tree->GetEntry(ientry)){

    // Break 
    if(nEntryToStop > 0 && ientry>=nEntryToStop){
      cout << "Break the loop after " << nEntryToStop << " entries!" << endl;
      break;
    }
    // update ientry after each loop
    ientry++;
    //fake-data control
    //if(tag.Contains("data")){
    //  bool isFakeData = anaUtils.IsFakeData();
    //  if(isFakeData) continue;
    //}
    double radGaus = 0.0;
    //double radGausffe = 0.0;

    //if(kMC) radGaus = grdm->Gaus(-0.00726,0.022311326);//grdm->Gaus(-0.00985,0.017756843);
    if(kMC && doit){
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
      anaUtils.FillXSTrueHistograms();

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
    // == Binnings for cross sections
    int N_binning_100MeV = 20;
    vector<double> binning_100MeV = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000.};
    // Thin-slice method, incident and interaction histograms calcualtion and filling
    if(kMC && doit){
      // Divide the MC sample into half according to the event number
      bool isFakeData = anaUtils.IsFakeData();
      // Get the MC weight for each event
      double weight = anaUtils.CalWeight(kMC);
      //double bckweight = anaUtils.CalBckWeight(kMC);

      // Test truth matched beam is true pion beam (and yes it is)
      if(AnaIO::reco_beam_true_byHits_PDG == 211){
        AnaIO::hTruthPDGCode_PionMatched->Fill(AnaIO::true_beam_PDG);
      }
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
          vector<TLorentzVector> vecFSParticle = anaUtils.GetFSParticlesTruth(false);
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
            if(isFakeData) xloop++;
            // Calculate the reco beam energy at the entry point
            //radGaus = grdm->Gaus(-0.00985,0.017756843);
            double beam_inst_P_smearing = AnaIO::beam_inst_P + radGaus; //grdm->Gaus(-0.00985,0.017756843);

            double beam_inst_KE = sqrt(pow(beam_inst_P_smearing,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
            // Taking into account the front face energy lost 12.74 MeV is used here
            //double Eloss = anaUtils.GetUpStreamEnergyLoss(kMC,beam_inst_KE);
            //double ff_energy_reco = beam_inst_KE*1000 - 12.74; 2.65229
            double ff_energy_reco = beam_inst_KE*1000 - 2.65229 + 25.0; //FIXME upper Syst

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
                mc_i++;
                // Fake data used to fill incident histogram
                if(isFakeData){
                  // Fill initial histogram
                  //AnaIO::hRecoInitialHist->Fill((*AnaIO::reco_beam_new_incidentEnergies)[0]);
                  AnaIO::hRecoInitialHist->Fill(initialE_reco, weight);
                  //cout << "(*AnaIO::reco_beam_new_incidentEnergies)[0]: " << (*AnaIO::reco_beam_new_incidentEnergies)[0] << " initialE_reco: " << initialE_reco << endl;
                  if(interactingE_reco!=-999) AnaIO::hRecoBeamInteractingHist->Fill(interactingE_reco, weight);
                  else AnaIO::hRecoBeamInteractingHist->Fill(AnaIO::reco_beam_new_incidentEnergies->back(), weight);

                  if(interactingE_reco!=-999) interactingE_reco = AnaIO::reco_beam_new_incidentEnergies->back();
                  // Sungbin's new method
                  anaUtils.FillEsliceHistograms(AnaIO::hNewRecoInitialHist, AnaIO::hNewRecoBeamInteractingHist, AnaIO::hNewRecoIncidentHist, AnaIO::hNewRecoInteractingHist, initialE_reco, interactingE_reco, -1, weight, binning_100MeV, N_binning_100MeV, false);

                  xloopin++;
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
                if(isFakeData) xloopout++;
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
                if(isFakeData) xcexloop++;
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
                    //>!
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
                      cex800loopin++;
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
                        cex800loopfill++;
                        plotUtils.FillHist(AnaIO::hRecPi0Energy_OVERLAY_After_EVT_High, pi0KineticE*1000, evtXStype, weight);
                        AnaIO::hRecoPi0KEHist->Fill(pi0KineticE*1000, weight);
                        AnaIO::hRecoPi0CosThetaHist->Fill(pi0costheta, weight);
                        //double pi0theta = TMath::RadToDeg() * TMath::ACos(pi0costheta);
                        AnaIO::hRecoPi0ThetaHist->Fill(pi0theta, weight);
                        if(pi0KineticE != -999) {
                          //cout << "pi0KineticE*1000: " << pi0KineticE*1000 << endl;
                          cex800loopfillnonzero++;
                        }
                      //}

                      /*else{
                        //if(interactingE > 750 && interactingE < 800) {
                          uf.response_SliceID_Pi0KE.Miss(LeadingPiZeroKE*1000, weight);
                          uf.response_SliceID_Pi0CosTheta.Miss(LeadingPiZeroCosTheta, weight);
                          uf.response_SliceID_Pi0Theta.Miss(LeadingPiZeroTheta, weight);
                          cex800loopout++;
                        //}
                      }*/
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


                    cex800loopout++;
                  }
                }
              }
            }
            else{ // reco KE beam size = 0
              //cout << "reco KE beam size = 0" << endl;
              if(isFakeData) {xloop++; xloopoutbadtrack++;}
            
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
                
                if(isFakeData) xcexloop++;
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
                  
                  cex800loopout++;
                }
              }
            }
          } // End of if(trackLen != -999)...

          // Unphysical reco track length (also need to be used to fill the missing response matrix)
          else{
            if(isFakeData) {xloop++; xloopoutbadtrack++;}
            
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
              
              if(isFakeData) xcexloop++;
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

                cex800loopout++;
              }
            }
          } // End of Unphysical reco track length ...
        }
      } // End of if(AnaIO::true_beam_PDG == 211)...

    } // End of if(kMC) ...

    //For data need to think more about the structure
    if(!kMC && doit){

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
        double ff_energy_reco = beam_inst_KE*1000 - 2.65229 + 25.0; //FIXME upper Syst; 

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
          data_i++;
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
 

    //====================== Do cuts (both MC and data) ======================//
    // Do the beam cut
    if(!anaCut.CutBeamAllInOne(kMC, true)) continue;

    double weight = anaUtils.CalWeight(kMC);
    double bckweight = anaUtils.CalBckWeight(kMC);

    // Count beam after beam cut before other cuts
    BeamCount++;
    //BeamCount += 1.0*weight*bckweight;

    anaUtils.GetFSParticlesTruth();
    if(anaUtils.GetNParticles()[2] > 0) anaUtils.GetFSPiZeroDecayDaughterTruth();

    // Fill beam info
    anaUtils.FillBeamKinematics(kMC);

    const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
    //const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
    const int evtXStype = kMC ? anaUtils.GetFillXSEventType() : anaUtils.gkXSBmBkg;

    // Get the reco track length from reco_beam_calo_X,Y,Z
    double trackLen = anaUtils.GetRecoTrackLength();
    vector<double> trackLenAccum = anaUtils.GetRecoTrackLengthAccumVect();
    

    double interactingE_reco = -999;
    //double pi0KineticE;
    //cout << "data trackLen: " << trackLen << endl;
    if(trackLen != -999) {
      double beam_inst_P_smearing = AnaIO::beam_inst_P;
      //if(radGaus == 0) radGaus = grdm->Gaus(-0.00985,0.017756843);
      if(kMC) beam_inst_P_smearing = AnaIO::beam_inst_P + radGaus;//grdm->Gaus(-0.00955164,0.0176993);//radGaus;//grdm->Gaus(-0.00985,0.017756843);
      // Calculate the reco beam energy at the entry point 
      double beam_inst_KE = sqrt(pow(beam_inst_P_smearing,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();

      if(kMC){
        int start_idx = -1; double true_ffKE = -999;
        for (unsigned int i=0; i<AnaIO::true_beam_traj_Z->size(); i++){
          if ((*AnaIO::true_beam_traj_Z)[i] >= 0){
            start_idx = i-1; // the trajectory point before entering the TPC
            if (start_idx < 0) start_idx = -1;
            break;
          }
        }
        if(start_idx >= 0) true_ffKE = (*AnaIO::true_beam_traj_KE)[start_idx];

        double UpStreamELoss = beam_inst_KE*1000 - true_ffKE;
  
        if(beam_inst_KE > 0.7 && beam_inst_KE < 0.75) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,0,weight);
        if(beam_inst_KE > 0.75 && beam_inst_KE < 0.8) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,1,weight);
        if(beam_inst_KE > 0.8 && beam_inst_KE < 0.85) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,2,weight);
        if(beam_inst_KE > 0.85 && beam_inst_KE < 0.9) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,3,weight);
        if(beam_inst_KE > 0.9 && beam_inst_KE < 0.95) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,4,weight);
        if(beam_inst_KE > 0.95 && beam_inst_KE < 1.0) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,5,weight);
        if(beam_inst_KE > 1.0 && beam_inst_KE < 1.05) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,6,weight);
        if(beam_inst_KE > 1.05 && beam_inst_KE < 1.1) plotUtils.FillHist(AnaIO::hUpStreamELossAfterSmearingAndWeight,UpStreamELoss,7,weight);
      }

      //double Eloss = anaUtils.GetUpStreamEnergyLoss(kMC, beam_inst_KE);
      //double ff_energy_reco = beam_inst_KE*1000 - 12.74;
      double ff_energy_reco = beam_inst_KE*1000 - 2.65229 + 25.0; //FIXME upper Syst
      
      //ff_energy_reco = ff_energy_reco + grdm->Gaus(0.00541435,0.0367387);
      double initialE_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[0]); 

      AnaIO::reco_beam_incidentEnergies->clear();
      for(unsigned int idx = 0; idx < trackLenAccum.size(); idx++){
        double energy_reco = bb.KEAtLength(ff_energy_reco, trackLenAccum[idx]);
        //cout << "idx: " << idx << " energy_reco: " << energy_reco << endl;
        AnaIO::reco_beam_incidentEnergies->push_back(energy_reco);
      }
      // Reco beam incident energy (vector) and interaction energy calculation (use reco_beam_calo_Z and reco_beam_incidentEnergies) 
      AnaIO::reco_beam_new_incidentEnergies->clear();
      interactingE_reco = anaUtils.MakeRecoIncidentEnergies(AnaIO::reco_beam_calo_Z, AnaIO::reco_beam_incidentEnergies, AnaIO::reco_beam_new_incidentEnergies);
      if(interactingE_reco==-999) interactingE_reco = AnaIO::reco_beam_new_incidentEnergies->back();
      
      double initialE_reco_Jake = (*AnaIO::reco_beam_new_incidentEnergies)[0];

      int size_reco = AnaIO::reco_beam_new_incidentEnergies->size();
      int caloz_reco = AnaIO::reco_beam_calo_Z->size();


      //if(anaCut.CutBeamAllInOne(kMC)){ // Pass the beam cuts
      // Pass the beam cuts
      // Fill the reco incident histograms
      for(int i = 0; i < size_reco; i++){ // size_reco -1 is the interactingE_reco
        double energy_reco = (*AnaIO::reco_beam_new_incidentEnergies)[i];
        
        //if(parType == anaUtils.gkBeamPiPlus) 
        //AnaIO::hRecoIncidentHistData->Fill(energy_reco,0.889673);
        plotUtils.FillHist(AnaIO::hRecPiPlusIncidentEnergy, energy_reco, parType, weight*bckweight);
      }

      // Fill the reco incident histograms
      for(int i = 0; i < caloz_reco; i++){ // size_reco -1 is the interactingE_reco
        double caloz = (*AnaIO::reco_beam_calo_Z)[i];
        
        //if(parType == anaUtils.gkBeamPiPlus) 
        //AnaIO::hRecoIncidentHistData->Fill(caloz,0.889673);
        plotUtils.FillHist(AnaIO::hRecPiPlusCaloZ, caloz, parType, weight*bckweight);
      }

      
      
      plotUtils.FillHist(AnaIO::hRecPiPlusEnergy_OVERLAY_After_EVTXS, interactingE_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergy, interactingE_reco, parType, weight*bckweight);

      plotUtils.FillHist(AnaIO::hRecPiPlusInstMomentum, beam_inst_P_smearing*1000, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusFrontFaceEnergy, ff_energy_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInitialEnergy, initialE_reco, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusInitialEnergyPosZCut, initialE_reco_Jake, parType, weight*bckweight);
      plotUtils.FillHist(AnaIO::hRecPiPlusTrackLength, trackLen, parType, weight*bckweight);

      // Calculate the incident entry
      double interval = initialE_reco - interactingE_reco; // 10MeV bin
      double nbin = interval/1.; // Number of bins 
      double inibin = initialE_reco/1.;
      double intbin = interactingE_reco/1.;


      int a = (int)inibin;
      int b = (int)intbin;

      //double x = nbin + 0.5 - (x<0); // x is now 55.499999...
      int y = (int)nbin + 1;
      //double intbin = interactingE_reco/10.;
      int y0 = (int)intbin + 1;
      int y1 = y0 + y;

      if(a != b){
        
        //if(y0 == y) y1 = 0;

        //cout << "nbin: " << nbin << " y: " << y << endl;
        for(int idy = y0; idy <= y1; idy++){
          double incE = idy * 1.0 - 0.5;
          plotUtils.FillHist(AnaIO::hRecPiPlusIncidentEnergyNew, incE, parType, weight*bckweight);
        }
      }

      else {
        for(int idy = y0; idy < y1; idy++){
          double incE = idy * 1.0 - 0.5;
          plotUtils.FillHist(AnaIO::hRecPiPlusIncidentEnergyNew, incE, parType, weight*bckweight);
        }
      }
      //}
    }

  
    // Do event topology cut
    double pi0KineticE, pi0costheta;
    if(!anaCut.CutTopology(kMC,pi0KineticE,pi0costheta,true)) continue;

    CEXEventCount++;
    double XSevtweight = anaUtils.CalXSEvtWeight(kMC, interactingE_reco, evtXStype);

    double intcexwt = anaUtils.CalCEXIntWeight(interactingE_reco);

    plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyEvt, interactingE_reco, evtXStype, weight*bckweight*XSevtweight);
    plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyPar, interactingE_reco, parType, weight*bckweight*XSevtweight);

    
    if(kMC && evtXStype == anaUtils.gkXSSignal) plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyBckSubCheck, interactingE_reco, evtXStype, weight*bckweight);
    if(!kMC) plotUtils.FillHist(AnaIO::hRecPiPlusInteractingEnergyBckSubCheck, interactingE_reco, evtXStype, intcexwt);

    plotUtils.FillHist(AnaIO::hRecPiZeroKineticEnergyEvt, pi0KineticE*1000, evtXStype, XSevtweight*weight);

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
    
    const double beamEndZ = AnaIO::reco_beam_calo_endZ;
    const int channelType = kMC ? anaUtils.GetChannelType((*AnaIO::true_beam_endProcess),AnaIO::true_beam_PDG, AnaIO::true_daughter_nPi0, AnaIO::true_daughter_nPiPlus, AnaIO::true_daughter_nPiMinus) : anaUtils.gkBackground;
    plotUtils.FillHist(AnaIO::hBeamEndZ_Channels_afterEvtTop, beamEndZ, channelType);
    
    // Do TKI calculation 
    //anaUtils.TruthMatchingTKI(anaUtils.RecPi0LTVet,anaUtils.RecProtonLTVet,anaUtils.TruthPi0LTVet,anaUtils.TruthProtonLTVet,kMC,anaUtils.GoodTruthMatch);
  
    // Fill output tree
    tout->Fill();
  } // End of while loop

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



  cout << "xloop: " << xloop << endl; 
  cout << "trueloop: " << trueloop << endl; 
  cout << "xloopin: " << xloopin << "xloopout: " << xloopout << "xloopoutbadtrack: " << xloopoutbadtrack << endl;
  cout << "xlooptotal: " << xloopin + xloopout + xloopoutbadtrack << endl;

  cout << "xcexloop: " << xcexloop << endl; 
  cout << "truecexloop: " << truecexloop << endl;

  cout << "truecex800loop: " << truecex800loop << endl;
  cout << "cex800loopin: " << cex800loopin << endl; 
  cout << "cex800loopout: " << cex800loopout << endl; 
  cout << "cex800looptotal: " << cex800loopin + cex800loopout << endl; 
  cout << "cex800looptotalhalf: " << (cex800loopin + cex800loopout)/2.0 << endl; 

  cout << "cex800loopfill: " << cex800loopfill << endl;
  cout << "cex800loopfillnonzero: " << cex800loopfillnonzero << endl;

  cout << "beam back: " << beampiplus/(totalbeam+beampiplus) << endl;
  cout << "cex back: " << cex/(totalcex+cex) << endl;

  cout << "pi0cex back: " << pi0cex/(pi0totalcex+pi0cex) << endl;

  cout << "good: " << goodevt/(goodevt+badevt) << endl;


  uf.SaveHistograms();

  // Print info
  cout << "All entries: " << ientry << endl;
  cout << "BeamCount: " << BeamCount << endl;
  cout << "CEXEventCount: " << CEXEventCount << endl;
  // Kinematic Fitting for Pi0 shower
  //if(kMC) anaUtils.DoKinematicFitting();



  if(kMC && doit) {

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


  } // End of if(kMC && doit) ...

  if(!kMC && doit){
    
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
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nshower",  icut++), AnaIO::hCutnshower, 2, 100000, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npiplus",  icut++), AnaIO::hCutnpiplus, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nmichel",  icut++), AnaIO::hCutnmichel, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npi0",  icut++), AnaIO::hCutnpi0, 1, 1, nsel);

  /*nsel = plotUtils.PrintStat(tag+Form(" %d. Nproton",  icut++), AnaIO::hCutnproton, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nshower",  icut++), AnaIO::hCutnshower, 2, 100000, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npiplus",  icut++), AnaIO::hCutnpiplus, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nmichel",  icut++), AnaIO::hCutnmichel, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npi0",  icut++), AnaIO::hCutnpi0, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. KF Pass",  icut++), AnaIO::hCutKFPass, 1, 1, nsel);
*/

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

  //if(kMC) BeamCount = mc_i++;
  //else BeamCount = data_i++;

  return BeamCount;
  
} // End of anaRec

int main(int argc, char * argv[])
{
  // Initialise the input file name
  const TString mcfinName = "input/protoDUNE_mc_reco_flattree_prod4a_ntuple.root";
  //const TString mcfinName = "input/protoDUNE_data_reco_flattree_prod4_ntuple_fakeSCEoff.root";
  //const TString mcfinName = "input/protoDUNE_mc_reco_flattree_prod4a_ntuple_SCEoff.root";
  const TString datafinName  = "input/protoDUNE_data_reco_flattree_prod4_ntuple.root";
  //const TString datafinName  = "input/protoDUNE_mc_reco_flattree_prod4a_ntuple.root";

  //const TString datafinName  = "input/protoDUNE_fake_data_reco_flattree_prod4_ntuple.root";

  //const TString datafinName  = "input/protoDUNE_data_reco_flattree_prod4_ntuple_fakeSCEoff.root";

  const TString xcesfinName = "input/LAr_PiPlus_exclusive_cross_section.root";
  TFile *f_CrossSection = TFile::Open(xcesfinName);
  if(!f_CrossSection->IsOpen()){
    cout << "xsce file not open" << endl;
    exit(1);
  }
  //TGraph *g_inel;
  TGraph *g_el;
  TGraph *g_cex;
  TGraph *g_totalinel;


  //f_CrossSection->GetObject("inel_KE",g_inel);
  f_CrossSection->GetObject("cex_KE",g_cex);
  f_CrossSection->GetObject("el_KE",g_el);
  f_CrossSection->GetObject("total_inel_KE",g_totalinel);

  f_CrossSection->Close();


  //const TString DiffxcesfinName = "input/cross_section_out_new1GeV.root";
  // New G4 file with costheta info
  const TString DiffxcesfinName = "input/cross_section_out_withNewVar.root";
  //const TString DiffxcesfinName = "input/cross_section_out_withNewVarLargeBin.root";
  TFile *f_DiffCrossSection = TFile::Open(DiffxcesfinName);
  //TFile *f_DiffCrossSection = new TFile(DiffxcesfinName);
  if(!f_DiffCrossSection->IsOpen()) {
    cout << "Diffxsce file not open" << endl;
    exit(1);
  }
  
  //TH1D *g_cex_400MeV;
  TH1D *g_cex_675MeV;
  TH1D *g_cex_775MeV;
  TH1D *g_cex_875MeV;

  TH1D *g_cexTheta_675MeV;
  TH1D *g_cexTheta_775MeV;
  TH1D *g_cexTheta_875MeV;

  TH1D *g_cexCosTheta_675MeV;
  TH1D *g_cexCosTheta_775MeV;
  TH1D *g_cexCosTheta_875MeV;

  //f_DiffCrossSection->GetObject("inel_cex_1dKEpi0400_MeV",g_cex_400MeV);
  //f_DiffCrossSection->GetObject("inel_cex_1dKEpi0600_MeV",g_cex_600MeV);
  //f_DiffCrossSection->GetObject("inel_cex_1dKEpi0800_MeV",g_cex_800MeV);

  g_cex_675MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dKEpi0675_MeV");
  g_cex_775MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dKEpi0775_MeV");
  g_cex_875MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dKEpi0875_MeV");

  g_cexTheta_675MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dThetapi0675_MeV");
  g_cexTheta_775MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dThetapi0775_MeV");
  g_cexTheta_875MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dThetapi0875_MeV");

  g_cexCosTheta_675MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dcosThetapi0675_MeV");
  g_cexCosTheta_775MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dcosThetapi0775_MeV");
  g_cexCosTheta_875MeV = (TH1D*)f_DiffCrossSection->Get("inel_cex_1dcosThetapi0875_MeV");

  //f_DiffCrossSection->Close();

  //g_cexTheta_675MeV->SetName("");
  //g_cexTheta_775MeV->SetName("");
  //g_cexTheta_875MeV->SetName("");

  //cout << "Name: " << g_cex_675MeV->GetName() << endl;
  const Int_t x0_diff = g_cex_675MeV->GetXaxis()->GetFirst();
  const Int_t x1_diff = g_cex_675MeV->GetXaxis()->GetLast();
  double x_675[x1_diff], y_675[x1_diff];
  for(Int_t ix=x0_diff; ix<=x1_diff; ix++){
    //cout << "Bin 675: "  << ix << " " << g_cex_675MeV->GetBinContent(ix) << endl;
    x_675[ix] = g_cex_675MeV->GetBinCenter(ix);
    y_675[ix] = g_cex_675MeV->GetBinContent(ix);
  }

  auto g_675 = new TGraph(40,x_675,y_675);
/*
  // Theta case
  cout << "Name: " << g_cexTheta_675MeV->GetName() << endl;
  const Int_t x0_difftheta = g_cexTheta_675MeV->GetXaxis()->GetFirst();
  const Int_t x1_difftheta = g_cexTheta_675MeV->GetXaxis()->GetLast();
  double x_675theta[x1_difftheta], y_675theta[x1_difftheta];
  for(Int_t ix=x0_difftheta; ix<=x1_difftheta; ix++){
    //cout << "Bin 675: "  << ix << " " << g_cexTheta_675MeV->GetBinContent(ix) << endl;
    x_675theta[ix] = g_cexTheta_675MeV->GetBinCenter(ix);
    y_675theta[ix] = g_cexTheta_675MeV->GetBinContent(ix);
  }
*/
  //auto g_675theta = new TGraph(40,x_675theta,y_675theta);
  auto g_675theta = new TGraph();
  for(int i = 1; i <= g_cexTheta_675MeV->GetNbinsX(); i++ ) g_675theta->SetPoint(i+1, g_cexTheta_675MeV->GetBinCenter(i), g_cexTheta_675MeV->GetBinContent(i));
/*
  // Cos Theta case
  cout << "Name: " << g_cexCosTheta_675MeV->GetName() << endl;
  const Int_t x0_diffcostheta675 = g_cexCosTheta_675MeV->GetXaxis()->GetFirst();
  const Int_t x1_diffcostheta675 = g_cexCosTheta_675MeV->GetXaxis()->GetLast();
  double x_675costheta[x1_diffcostheta675], y_675costheta[x1_diffcostheta675];
  for(Int_t ix=x0_diffcostheta675; ix<=x1_diffcostheta675; ix++){
    //cout << "Bin 675: "  << ix << " " << g_cexCosTheta_675MeV->GetBinContent(ix) << endl;
    x_675costheta[ix] = g_cexCosTheta_675MeV->GetBinCenter(ix);
    y_675costheta[ix] = g_cexCosTheta_675MeV->GetBinContent(ix);
  }

  auto g_675costheta = new TGraph(40,x_675costheta,y_675costheta);
*/
  auto g_675costheta = new TGraph();
  for(int i = 1; i <= g_cexCosTheta_675MeV->GetNbinsX(); i++ ) g_675costheta->SetPoint(i+1, g_cexCosTheta_675MeV->GetBinCenter(i), g_cexCosTheta_675MeV->GetBinContent(i));
  g_675costheta->RemovePoint(0);
  g_675costheta->RemovePoint(1);
/*
  const int N_points = g_675costheta->GetN();
  double x, y;
  for(int i=0; i < N_points; i++)
  {
    g_675costheta->GetPoint(i, x, y); 
    cout << "f(" << x << ") =" << y <<endl;
  }
*/
  //cout << "Name 775: " << g_cex_775MeV->GetName() << endl;
  const Int_t x0_diff775 = g_cex_775MeV->GetXaxis()->GetFirst();
  const Int_t x1_diff775 = g_cex_775MeV->GetXaxis()->GetLast();
  double x_775[x1_diff775], y_775[x1_diff775];
  for(Int_t ix=x0_diff775; ix<=x1_diff775; ix++){
    //cout << "Bin 775: "  << ix << " " << g_cex_775MeV->GetBinContent(ix) << endl;
    x_775[ix] = g_cex_775MeV->GetBinCenter(ix);
    y_775[ix] = g_cex_775MeV->GetBinContent(ix);
  }

  auto g_775 = new TGraph(40,x_775,y_775);
/*
  // Theta case
  cout << "Name: " << g_cexTheta_775MeV->GetName() << endl;
  const Int_t x0_difftheta775 = g_cexTheta_775MeV->GetXaxis()->GetFirst();
  const Int_t x1_difftheta775 = g_cexTheta_775MeV->GetXaxis()->GetLast();
  double x_775theta[x1_difftheta775], y_775theta[x1_difftheta775];
  for(Int_t ix=x0_difftheta775; ix<=x1_difftheta775; ix++){
    //cout << "Bin 775: "  << ix << " " << g_cexTheta_775MeV->GetBinContent(ix) << endl;
    x_775theta[ix] = g_cexTheta_775MeV->GetBinCenter(ix);
    y_775theta[ix] = g_cexTheta_775MeV->GetBinContent(ix);
  }
*/
  //auto g_775theta = new TGraph(40,x_775theta,y_775theta);
  auto g_775theta = new TGraph();
  for(int i = 1; i <= g_cexTheta_775MeV->GetNbinsX(); i++ ) g_775theta->SetPoint(i+1, g_cexTheta_775MeV->GetBinCenter(i), g_cexTheta_775MeV->GetBinContent(i));
/*
  // Cos Theta case
  cout << "Name: " << g_cexCosTheta_775MeV->GetName() << endl;
  const Int_t x0_diffcostheta775 = g_cexCosTheta_775MeV->GetXaxis()->GetFirst();
  const Int_t x1_diffcostheta775 = g_cexCosTheta_775MeV->GetXaxis()->GetLast();
  double x_775costheta[x1_diffcostheta775], y_775costheta[x1_diffcostheta775];
  for(Int_t ix=x0_diffcostheta775; ix<=x1_diffcostheta775; ix++){
    //cout << "Bin 775: "  << ix << " " << g_cexCosTheta_775MeV->GetBinContent(ix) << endl;
    x_775costheta[ix] = g_cexCosTheta_775MeV->GetBinCenter(ix);
    y_775costheta[ix] = g_cexCosTheta_775MeV->GetBinContent(ix);
  }

  auto g_775costheta = new TGraph(40,x_775costheta,y_775costheta);
*/

  auto g_775costheta = new TGraph();
  for(int i = 1; i <= g_cexCosTheta_775MeV->GetNbinsX(); i++ ) g_775costheta->SetPoint(i+1, g_cexCosTheta_775MeV->GetBinCenter(i), g_cexCosTheta_775MeV->GetBinContent(i));
  g_775costheta->RemovePoint(0);
  g_775costheta->RemovePoint(1);

  //cout << "Name 875: " << g_cex_875MeV->GetName() << endl;
  const Int_t x0_diff875 = g_cex_875MeV->GetXaxis()->GetFirst();
  const Int_t x1_diff875 = g_cex_875MeV->GetXaxis()->GetLast();
  double x_875[x1_diff875], y_875[x1_diff875];
  for(Int_t ix=x0_diff875; ix<=x1_diff875; ix++){
    //cout << "Bin 875: "  << ix << " " << g_cex_875MeV->GetBinContent(ix) << endl;
    x_875[ix] = g_cex_875MeV->GetBinCenter(ix);
    y_875[ix] = g_cex_875MeV->GetBinContent(ix);
  }

  auto g_875 = new TGraph(40,x_875,y_875);
/*
  // Theta case
  cout << "Name: " << g_cexTheta_875MeV->GetName() << endl;
  const Int_t x0_difftheta875 = g_cexTheta_875MeV->GetXaxis()->GetFirst();
  const Int_t x1_difftheta875 = g_cexTheta_875MeV->GetXaxis()->GetLast();
  double x_875theta[x1_difftheta875], y_875theta[x1_difftheta875];
  for(Int_t ix=x0_difftheta875; ix<=x1_difftheta875; ix++){
    //cout << "Bin 875: "  << ix << " " << g_cexTheta_875MeV->GetBinContent(ix) << endl;
    x_875theta[ix] = g_cexTheta_875MeV->GetBinCenter(ix);
    y_875theta[ix] = g_cexTheta_875MeV->GetBinContent(ix);
  }
*/
  //auto g_875theta = new TGraph(40,x_875theta,y_875theta);
  auto g_875theta = new TGraph();
  for(int i = 1; i <= g_cexTheta_875MeV->GetNbinsX(); i++ ) g_875theta->SetPoint(i+1, g_cexTheta_875MeV->GetBinCenter(i), g_cexTheta_875MeV->GetBinContent(i));
/*
  // Cos Theta case
  cout << "Name: " << g_cexCosTheta_875MeV->GetName() << endl;
  const Int_t x0_diffcostheta875 = g_cexCosTheta_875MeV->GetXaxis()->GetFirst();
  const Int_t x1_diffcostheta875 = g_cexCosTheta_875MeV->GetXaxis()->GetLast();
  double x_875costheta[x1_diffcostheta875], y_875costheta[x1_diffcostheta875];
  for(Int_t ix=x0_diffcostheta875; ix<=x1_diffcostheta875; ix++){
    //cout << "\nBin 875: "  << ix << " " << g_cexCosTheta_875MeV->GetBinContent(ix) << endl;
    //cout << "Bin 875 center: "  << ix << " " << g_cexCosTheta_875MeV->GetBinCenter(ix) << endl;
    x_875costheta[ix] = g_cexCosTheta_875MeV->GetBinCenter(ix);
    y_875costheta[ix] = g_cexCosTheta_875MeV->GetBinContent(ix);
  }

  auto g_875costheta = new TGraph(40,x_875costheta,y_875costheta);
*/
  auto g_875costheta = new TGraph();
  for(int i = 1; i <= g_cexCosTheta_875MeV->GetNbinsX(); i++ ) g_875costheta->SetPoint(i+1, g_cexCosTheta_875MeV->GetBinCenter(i), g_cexCosTheta_875MeV->GetBinContent(i));
  g_775costheta->RemovePoint(0);
  g_775costheta->RemovePoint(1);
  
  //f_DiffCrossSection->Close();
  cout << "Name: " << g_cexTheta_875MeV->GetName() << endl;

  delete g_cex_675MeV;
  delete g_cex_775MeV;
  delete g_cex_875MeV;

  delete g_cexTheta_675MeV;
  delete g_cexTheta_775MeV;
  delete g_cexTheta_875MeV;

  delete g_cexCosTheta_675MeV;
  delete g_cexCosTheta_775MeV;
  delete g_cexCosTheta_875MeV;

  f_DiffCrossSection->Close();

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
  plotUtils.DrawHist(mclout,plotScale,datalout,"output",g_totalinel,g_cex,g_675,g_775,g_875,g_675theta,g_775theta,g_875theta,g_675costheta,g_775costheta,g_875costheta);

  //f_CrossSection->Close();
  //f_DiffCrossSection->Close();
}
