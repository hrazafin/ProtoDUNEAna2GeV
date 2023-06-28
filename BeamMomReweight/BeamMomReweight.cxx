#include "include/AnaIO.h"
#include "src/AnaUtils.cxx"
#include "src/AnaCut.cxx"
#include "src/Unfold.cxx"
#include "src/PlotUtils.cxx"


double SliceTruthBeamEnergy(vector<double> *true_beam_traj_Z, vector<double> *true_beam_traj_KE, vector<double> *true_beam_inc_KE);
double SliceRecoBeamEnergy(const double & ffe, vector<double> *reco_beam_inc_KE, vector<double> *reco_beam_inc_KEweight);
TH1D * TotalCEXXSCal(TH1D * hh, TH1D * InteractingHist, std::map< int, BetheBloch* > map_BB);
void DrawTotalXSPlots(TH1D *hinc, TH1D *hint, TGraph* g_cex, std::map< int, BetheBloch* > map_BB);
void SetTitleFormat(TH1 * hh);
void DrawDataMCRatio(TH1D * hratio, bool xsec, const TString tag);
bool IsFakeData();
double GetTrueTrackLength(double & true_ffKE, double & int_energy_true, std::map< int, BetheBloch* > map_BB);
double Get_true_ffKE(double KE_in_TPC, double length_to_ff, std::map< int, BetheBloch* > map_BB);
double GetLoss(double InstE);
double GetUpStreamEnergyLoss(const bool & kMC, const double & InstE);


int anaRec(const TString finName, TList *lout, const TString tag, const int nEntryToStop, double &beamCount_longTrack)
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

  std::map< int, BetheBloch* > map_BB;
  map_BB[13] = new BetheBloch(13);
  map_BB[211] = new BetheBloch(211);
  map_BB[321] = new BetheBloch(321);
  map_BB[2212] = new BetheBloch(2212);

  // Get the TTree from input file
  TTree  * tree = AnaIO::GetInputTree(fin, "pduneana/beamana", tag);
  TTree  * tout = AnaIO::GetOutputTree(lout, tag);

  // Initialise reco histograms
  AnaIO::IniHist(lout, kMC);
  AnaUtils anaUtils;
  AnaCut anaCut;
  PlotUtils plotUtils;
  Unfold uf(24, 0, 1200);
  gStyle->SetOptStat(0);

  // Initialise entry and counters
  double ientry = 0;
  double beamCount = 0;
  beamCount_longTrack = 0;

  // Loop over TTree
  while(tree->GetEntry(ientry)){

    // Break 
    if(nEntryToStop > 0 && ientry>=nEntryToStop){
      cout << "Break the loop after " << nEntryToStop << " entries!" << endl;
      break;
    }
    // update ientry after each loop
    ientry++;

    // ================ Truth Stopping Muon Sample ============= //
    if(kMC){
      if(AnaIO::true_beam_PDG == -13){
        double true_ffKE = -999, true_int_energy = -999;
        double trkLen_true = GetTrueTrackLength(true_ffKE,true_int_energy,map_BB);
        double byRange_len_true = map_BB[13]->RangeFromKE(true_ffKE);
        AnaIO::hTrackLenRatio->Fill(trkLen_true/byRange_len_true);
        if(trkLen_true/byRange_len_true > 0.9 && AnaIO::true_beam_endZ > 300.0) AnaIO::hTrueBeamMom->Fill(AnaIO::true_beam_startP);

      }
    }


    // ================= Long Track Muon ================ //
    if(anaCut.CutBeamLongTrack(kMC)){
      beamCount_longTrack++;
      const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
      double beam_trackLength = anaUtils.GetRecoTrackLength();
      plotUtils.FillHist(AnaIO::hBeamTrackLength_STK, beam_trackLength, parType, 1.0);
      if(beam_trackLength > 150.0) {
        plotUtils.FillHist(AnaIO::hBeamLongTrackLength_STK, beam_trackLength, parType, 1.0);
        if(kMC){
          if(parType == anaUtils.gkBeamMuon) AnaIO::hBeamLongTrackLength_Muon->Fill(beam_trackLength);
          else AnaIO::hBeamLongTrackLength_Other->Fill(beam_trackLength);
        }
      }
    }


    // ===================== Stop Muon ===============//
    const bool kStopMuon = true;
    if(!anaCut.CutBeamAllInOne(kMC,kStopMuon)) continue;
    
    const int parType = kMC ? anaUtils.GetBeamParticleType(AnaIO::reco_beam_true_byHits_PDG) : anaUtils.gkBeamOthers;
    double beam_trackLength = anaUtils.GetRecoTrackLength();

    double beam_inst_KE = sqrt(pow(AnaIO::beam_inst_P,2)+pow(AnaFunctions::PionMass(),2)) - AnaFunctions::PionMass();
    // Taking into account the front face energy lost
    double Eloss = anaUtils.GetUpStreamEnergyLoss(kMC,beam_inst_KE);
    double ff_energy_reco = beam_inst_KE*1000 - Eloss;//12.74;

    double byRange_len = map_BB[13]->RangeFromKE(ff_energy_reco);

    plotUtils.FillHist(AnaIO::hBeamTrackLength, beam_trackLength, parType, 1.0);

    AnaIO::hTrackLenRatio_Reco->Fill(beam_trackLength/byRange_len);
  
    plotUtils.FillHist(AnaIO::hTrackLenRatio_Reco_STK, beam_trackLength/byRange_len, parType, 1.0);

    // Select the stop muon samples
    if(beam_trackLength/byRange_len > 0.9){
      beamCount++;
      plotUtils.FillHist(AnaIO::hStopMuonInstP_STK, AnaIO::beam_inst_P, parType, 1.0);

      double ffE_muon = map_BB[13]->KEFromRangeSpline(beam_trackLength);

      plotUtils.FillHist(AnaIO::hStopMuonKEff_STK, ffE_muon, parType, 1.0);
      
      AnaIO::hTrueBeamMom_matched->Fill(AnaIO::true_beam_startP);

      AnaIO::trueBeamMom = AnaIO::true_beam_startP;
      AnaIO::recoFFKE = ffE_muon;
      AnaIO::instP = AnaIO::beam_inst_P;

      tout->Fill();

    }

    // Fill output tree
    //tout->Fill();
  } // End of while loop


  return beamCount;
  
} // End of anaRec

double GetUpStreamEnergyLoss(const bool & kMC, const double & InstE){
  
  double delta_E = 0.0;
  //if(kMC && InstE > 0){
  if(InstE > 0){

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

double GetLoss(double InstE){

  TF1 *fpCor = new TF1("fpCor","pol2",-400,1400);
  fpCor->SetParameter(0,170.777);
  fpCor->SetParameter(1,-0.595308);
  fpCor->SetParameter(2,0.000456932);
  const double factor = fpCor->Eval(InstE);
  return factor;
}


bool IsFakeData(){
  // Add fake data type for xsec study
  const int event = AnaIO::event;
  if(event%2) return false;
  else return true;
}


TH1D * TotalCEXXSCal(TH1D * hh, TH1D * InteractingHist, std::map< int, BetheBloch* > map_BB){
  TH1D *xsec = (TH1D*)hh->Clone();
  // Scale down it to zero for later calculation
  xsec->Scale(0);

  double avogadro_constant = 6.02214076e23;  // 1 / mol
  double argon_molar_mass = 39.95;           // g / mol
  double liquid_argon_density = 1.39;        // g / cm^3
  double fiducial_thickness = 1.0;//slice_width;//0.479; 

  double sigma_factor = argon_molar_mass / (avogadro_constant * liquid_argon_density * fiducial_thickness); // pre-factor

  // Define start and end bins
  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();

  // Loop over each bin and combine the incident and interacting histograms to get xsec
  for(Int_t ix=x0; ix<=x1; ix++){
    // Read off the entry for each bin of the two histograms
    double incE = hh->GetBinContent(ix);
    double intE = InteractingHist->GetBinContent(ix);
    // Calculate the ratio in the full log form
    double ratio = 0;
    if(incE > intE) ratio = intE/incE;//log(incE/(incE-intE));//*meandEdx_cal[ix-1]*sigma_factor*1e27;

    double dedx = map_BB[211]->meandEdx(hh->GetBinCenter(ix));

    //cout << "ix: " << ix << " hh->GetBinCenter(ix): " << hh->GetBinCenter(ix) << " dedx : " << dedx << "meandEdx_cal[ix-1]: " << endl;//<< meandEdx_cal[ix-1] << endl;

    if(intE > 0) ratio = ratio * dedx*sigma_factor*1e27;

    //if(incE != 0) cout << "ix: "<< ix << "incE*fiducial_thickness: " << incE*fiducial_thickness << "intE: " << intE << " ratio: " << ratio << endl;//"meandEdx_cal[ix]: " << meandEdx_cal[ix-1] << endl;

    // If the incE entry is not zero set the bin content
    if(incE != 0) xsec->SetBinContent(ix,ratio);
    //if( Eslice && ratio > 150) xsec->SetBinContent(ix,9999);
    //if( Eslice && ratio < 50) xsec->SetBinContent(ix,9999);

    // Error propagation
    double einc = hh->GetBinError(ix);
    double eint = InteractingHist->GetBinError(ix);
    double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)));

    // If the ratio is not zero set the error
    if(ratio != 0 ) xsec->SetBinError(ix,error);
    
  }
  // The xsec histogram entry is now set
  xsec->SetMaximum(250);
  xsec->SetMinimum(0);
  return xsec;

}

void DrawTotalXSPlots(TH1D *hinc, TH1D *hint, TGraph* g_cex, std::map< int, BetheBloch* > map_BB){
  
  TLatex tt;
  tt.SetNDC();

  TCanvas * c1 = new TCanvas("c1", "c1", 1200, 800);

  TH1D * hcex = TotalCEXXSCal(hinc,hint,map_BB); 

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 0.98);
  // Set 0.01 will hide axis label
  pad1->SetBottomMargin(0.01);
  //  pad1->SetGridx();
  //  pad1->SetGridy();
  pad1->Draw();
  pad1->cd();

  hcex->SetMarkerStyle(22);
  hcex->SetMarkerSize(1.5);
  hcex->SetMarkerColor(kBlue);
  hcex->SetLineColor(kBlue);
  hcex->SetLineWidth(1);
  //hcex->SetMaximum(0.4);
  hcex->SetMinimum(0.);
  hcex->GetXaxis()->SetRangeUser(0, 1000);
  //hcex->GetXaxis()->SetRangeUser(400, 1000);
  //hcex->GetYaxis()->SetRangeUser(500, 700);
  SetTitleFormat(hcex);
  //hcex->Draw("E1 sames");
  hcex->Draw("E1");

  g_cex->SetLineColor(kRed);
  g_cex->SetLineWidth(2);
  g_cex->Draw("sames C");

  TLegend * legend = new TLegend(0.15, 0.7, 0.38, 0.85);
  legend->AddEntry(g_cex, "Geant4 Input", "l");
  legend->AddEntry(hcex, "Truth Level", "lep");
  //TString lheader("1GeV Pion Data");
  //legend->SetHeader(lheader);
  legend->SetBorderSize(0);
  legend->Draw("same");

  c1->Update();
  c1->cd();

  // ======= pad2 =======//
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.2);
    
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.42);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  // ratio with data
  TH1D *hratio = (TH1D*)hcex->Clone("hcex");

  const Int_t x0_x = hcex->GetXaxis()->GetFirst();
  const Int_t x1_x = hcex->GetXaxis()->GetLast();

  hratio->Scale(0);

  for(Int_t ix=x0_x; ix<=x1_x; ix++){
    double bin_center = hcex->GetBinCenter(ix);
    
    double MC = g_cex->Eval(bin_center);
    
    double data = hcex->GetBinContent(ix);
    double edata = hcex->GetBinError(ix);

    if(MC != 0){
      double ratio = data/MC;
      //cout << "MC ratio: " << ratio << endl;
      hratio->SetBinContent(ix,ratio);
      double error = sqrt(ratio*ratio*(pow(edata/data,2)));
      hratio->SetBinError(ix,error);
    }
  }
  hratio->GetXaxis()->SetRangeUser(0, 1000);
  //hratio->GetXaxis()->SetRangeUser(400, 1000);
  
  hratio->SetTitle(" ");
  DrawDataMCRatio(hratio,true,"");
  c1->cd();

  tt.SetTextSize(0.035);
  tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
  tt.DrawLatex(0.725,0.925,"1GeV/c Pion Beam");
  tt.DrawLatex(0.685,0.775,"Truth Validation");


  c1->Print("output/TruthCEX_XSvalidation.png");

}

void SetTitleFormat(TH1 * hh){
  hh->SetTitle(" ");
  hh->GetYaxis()->CenterTitle();
  hh->GetYaxis()->SetTitleFont(22);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetYaxis()->SetTitleOffset(0.92);
  hh->GetXaxis()->CenterTitle();
  hh->GetXaxis()->SetTitleFont(22);
  hh->GetXaxis()->SetTitleSize(0.05);
  hh->GetXaxis()->SetTitleOffset(0.9);
  hh->GetYaxis()->SetTitle("#sigma_{CEX} (mb)");

}

void DrawDataMCRatio(TH1D * hratio, bool xsec, const TString tag){
  
  hratio->GetYaxis()->SetTitle("Data/MC");
  if(xsec) hratio->GetXaxis()->SetTitle("T_{#pi^{+}} (MeV)");
  if(!xsec && tag.Contains("CosTheta")) hratio->GetXaxis()->SetTitle("cos#theta_{#pi^{0}}");
  if(!xsec && tag.Contains("Theta") && !tag.Contains("Cos")) hratio->GetXaxis()->SetTitle("#theta_{#pi^{0}} (deg.)");
  
  hratio->GetXaxis()->SetLabelSize(0.15);
  hratio->GetXaxis()->SetTitleSize(0.15);
  hratio->GetXaxis()->SetTitleOffset(1.);
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetYaxis()->SetTitleOffset(.3);
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->GetYaxis()->SetRangeUser(0,2.1);
  hratio->GetYaxis()->SetRangeUser(0.5,1.6);
  hratio->GetXaxis()->CenterTitle();
  hratio->GetYaxis()->CenterTitle();

  hratio->GetXaxis()->SetTitleFont(22);
  hratio->GetYaxis()->SetTitleFont(22);
  hratio->GetXaxis()->SetTitleOffset(0.7);
  hratio->GetXaxis()->SetTitleSize(0.25);

  hratio->Draw();
}

void DrawHist(TList *lout)
{
  // Loop over all histograms inside Tlist
  for(int ii=0; ii<lout->GetSize(); ii++){
    // Get the name of the histogram
    const TString tag = lout->At(ii)->GetName();
    // Create new canvas for drawing histogram 
    TCanvas * c1 = new TCanvas(Form("c1_%s",tag.Data()), "", 1200, 800);

    if(tag.Contains("tree")) continue;
    TH1D *hh = (TH1D*)lout->FindObject(tag);
    hh->SetMaximum(hh->GetMaximum()*1.2);
    hh->SetLineWidth(2);
    hh->Draw("hist");

    if(tag.Contains("a00") && !tag.Contains("_uf")){
      hh->SetTitle("Incident Unfold Closure Test");
      if(tag.Contains("Interacting")) hh->SetTitle("Interacting Unfold Closure Test");
      hh->GetYaxis()->SetTitle("Candidates");

      TH1D *hoverlay = (TH1D*)lout->FindObject(tag+"_uf");
      hoverlay->SetLineColor(kRed);
      hoverlay->SetLineWidth(2);
      hoverlay->SetMarkerStyle(8);
      hoverlay->SetMarkerColor(kRed);
      hoverlay->Draw("E1 sames");

      TLegend * legend = new TLegend(0.15, 0.68, 0.48, 0.88);
      legend->AddEntry(hh, "MC Truth", "l");
      legend->AddEntry(hoverlay, "MC Reco. - Unfolded", "lep");
      legend->SetBorderSize(0);
      legend->Draw("sames");
    }
    
    c1->Print("output/"+tag+".png");
  }
}

double SliceTruthBeamEnergy(vector<double> *true_beam_traj_Z, vector<double> *true_beam_traj_KE, vector<double> *true_beam_inc_KE)
{ 
  // Clear the incident KE vetor before filling it
  true_beam_inc_KE->clear();
  // Only include trajectory points starting in the active volume
  double fTrajZStart = -0.49375;
  // ProtoDUNE TPC energy deposit [MeV]
  double fDetalE = 1.0; 
  // APA3 cut
  double endZ_APA = 220.0;

  double tmp_intE = -999, truth_intE = -999;
  int startix = -1;

  for (size_t j = 1; j < true_beam_traj_Z->size(); j++) {
    double z = true_beam_traj_Z->at(j);
    // check if then beam enters the TPC
    if (z > fTrajZStart){
      startix = j; break;
    }
  }
  
  // Only use the beam with initial slice within the TPC
  if(startix >= 0){
    double ffe = true_beam_traj_KE->at(startix-1);
    
    for (size_t i = startix; i < true_beam_traj_KE->size(); i++) {
      // beam z of this slice
      double z = true_beam_traj_Z->at(i-1);
      // Check if this beam endZ less than APA3 and get the interacting energy
      if(z < endZ_APA){
        if((*AnaIO::true_beam_traj_KE)[i] > 0.1){
          tmp_intE = (*AnaIO::true_beam_traj_KE)[i];
          truth_intE = tmp_intE;
        }
      }
      // Beam ends beyond APA3 then get the element at APA3 point
      else{
        // In this case the truth_intE is not set
        tmp_intE = (*AnaIO::true_beam_traj_KE)[i]; break;
      }
    }
    // Set incident energy as the front face energy to start
    double incE = ffe;
    // Slicing the beam KE and fill the incident histogram
    // Even though the beam ends after the APA3 it will still contribute to the incident flux
    while(incE > tmp_intE){
      // Fill the incident histogram
      true_beam_inc_KE->push_back(incE);
      // Decrement the incident energy by fDetalE
      incE -= fDetalE;
    }
  }

  return truth_intE;
  
}

double SliceRecoBeamEnergy(const double & ffe, vector<double> *reco_beam_inc_KE, vector<double> *reco_beam_inc_KEweight)
{ 
  // Clear the incident KE vetor before filling it
  reco_beam_inc_KE->clear();
  reco_beam_inc_KEweight->clear();
  // APA3 cut
  double endZ_APA = 220.0;
  // Declare calorimetry information 
  const vector<double> * trackPitch_SCE = AnaIO::reco_beam_TrkPitch_SCE;
  const vector<double> * dEdx_SCE = AnaIO::reco_beam_calibrated_dEdX_SCE;
  const vector<double> * caloZ = AnaIO::reco_beam_calo_Z;
  
  double incE = ffe;
  double intE = -999;
  if(dEdx_SCE->size() > 0){
    for(unsigned int ii = 0; ii < trackPitch_SCE->size()-1; ii++){ //-1 to not count the last slice
      double trkpit = (*trackPitch_SCE)[ii];
      double dEdx = (*dEdx_SCE)[ii];
      double z = (*caloZ)[ii];
      if(dEdx < 0) continue;
      // Energy slice
      double DeltaE = dEdx*trkpit;
      if(DeltaE > 1.0) DeltaE = 1.0; //cout << "DeltaE: " << DeltaE << endl;

      if(z < endZ_APA){
        reco_beam_inc_KE->push_back(incE);
        reco_beam_inc_KEweight->push_back(1.0);
        incE -= DeltaE;
        intE = incE;
      }
      else{
        intE = -999; break;
      }
    }
  }
  return intE;

}

double GetTrueTrackLength(double & true_ffKE, double & int_energy_true, std::map< int, BetheBloch* > map_BB)
{
  int start_idx = -1;
  vector<double> true_trklen_accum;
  //true_trklen_accum.clear();
  int size = AnaIO::true_beam_traj_Z->size();
  true_trklen_accum.reserve(size); // initialize true_trklen_accum

  for (int i=0; i<size; i++){
    if ((*AnaIO::true_beam_traj_Z)[i] >= 0){
      start_idx = i-1; // the trajectory point before entering the TPC
      if (start_idx < 0) start_idx = -1;
      break;
    }
    true_trklen_accum[i] = 0.; // initialize true_trklen_accum
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
      //true_trklen_accum.push_back(true_trklen);
      true_trklen_accum[i] = true_trklen;

    }

    true_ffKE = -999;

    //true_ffKE = (*AnaIO::true_beam_traj_KE)[start_idx];// + 2.18*(true_trklen_accum)[start_idx+1];

    true_ffKE = Get_true_ffKE((*AnaIO::true_beam_traj_KE)[start_idx+1], (true_trklen_accum)[start_idx+1], map_BB);
    

    int_energy_true = -999;
    /*int traj_max = AnaIO::true_beam_traj_Z->size()-1;
    if ((*AnaIO::true_beam_traj_KE)[traj_max] != 0) {
      int_energy_true = (*AnaIO::true_beam_traj_KE)[traj_max];
    }
    else {
      int temp = traj_max-1;
      while ((*AnaIO::true_beam_traj_KE)[temp] == 0) temp--;
      //int_energy_true = (*AnaIO::true_beam_traj_KE)[temp] - 2.1*((true_trklen_accum)[traj_max]-(true_trklen_accum)[temp]); // 2.1 MeV/cm
      int_energy_true = Get_true_intE((*AnaIO::true_beam_traj_KE)[temp], ((true_trklen_accum)[traj_max]-(true_trklen_accum)[temp]), map_BB);

    }
    //cout << "\nint_energy_true 1: " << int_energy_true << endl;
    */
    // Sunbin's method
    for(unsigned int i = 0; i < AnaIO::true_beam_traj_KE->size(); i++){
      if((*AnaIO::true_beam_traj_KE)[i] > 0.1){
        int_energy_true = (*AnaIO::true_beam_traj_KE)[i];
      }

      //cout << "[true_beam_traj_KE.at(" << i << ") : " << (*AnaIO::true_beam_traj_KE)->at(i) << endl;
    }

    //cout << "int_energy_true 2: " << int_energy_true << endl;

  }

  return true_trklen;
}


double Get_true_ffKE(double KE_in_TPC, double length_to_ff, std::map< int, BetheBloch* > map_BB){
  double this_dEdx = map_BB[13] -> meandEdx(KE_in_TPC);
  return KE_in_TPC + this_dEdx * length_to_ff;
}


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

  // Run reco loop for both mc and data
  int nEntryToStop = -1; // Default is looping over all entries
  // Get the input break entries if provided 
  if(argc!=1) nEntryToStop = atoi(argv[1]);

  // MC Analysis
  double beamCount_longTrack_mc = 0;
  double beamCount_longTrack_data = 0;

  double mcBeamCount = anaRec(mcfinName,mclout,"mc", nEntryToStop, beamCount_longTrack_mc);
  double dataBeamCount = anaRec(datafinName,datalout,"data", nEntryToStop, beamCount_longTrack_data);

  // Draw histograms
  //DrawHist(mclout);

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

  double plotScale = dataBeamCount/mcBeamCount;
  cout << "plotScale: " << plotScale << endl;

  double plotScale_longtrk = beamCount_longTrack_data/beamCount_longTrack_mc;
  cout << "plotScale_longtrk: " << plotScale_longtrk << endl;

  // Draw all histograms
  plotUtils.DrawHist(mclout,plotScale,plotScale_longtrk,datalout,"output");
}
