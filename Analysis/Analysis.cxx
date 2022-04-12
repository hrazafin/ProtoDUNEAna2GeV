#include "include/AnaIO.h"
#include "src/PlotUtils.cxx"
#include "src/AnaUtils.cxx"
#include "src/AnaCut.cxx"

int anaRec(const TString finName, TList *lout, const TString tag, const int nEntryToStop)
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
  
  // Loop over TTree
  while(tree->GetEntry(ientry)){

    // Break 
    if(nEntryToStop > 0 && ientry>=nEntryToStop){
      cout << "Break the loop after " << nEntryToStop << " entries!" << endl;
      break;
    }
    // update ientry after each loop
    ientry++;

    //====================== Do cuts (both MC and data) ======================//
    // Do the beam cut
    if(!anaCut.CutBeamAllInOne(kMC)) continue;
    // Count beam after beam cut before other cuts
    BeamCount++;
    
    //====================== Extract truth information (MC only)======================//
    if(kMC){
      // Get true beam particle type from it's pdg code
      int TruthBeamType = anaUtils.GetParticleType(AnaIO::true_beam_PDG); 
      // Fill the TruthBeamType histogram
      AnaIO::hTruthBeamType->Fill(TruthBeamType); 
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
        // Calcualte the TKI in truth level
        anaUtils.DoTruthTKICalculation();
      }
    }
    //====================== End truth information (MC only)======================//
    //if(kMC) cout << "Purity num: " << anaUtils.selected << endl;
    //if(kMC) cout << "Purity demo: " << anaUtils.total << endl;

    // Fill beam info
    anaUtils.FillBeamKinematics(kMC);
    
    // Do event topology cut
    if(!anaCut.CutTopology(kMC)) continue;

    
    // Do TKI calculation 
    anaUtils.TruthMatchingTKI(anaUtils.RecPi0LTVet,anaUtils.RecProtonLTVet,anaUtils.TruthPi0LTVet,anaUtils.TruthProtonLTVet,kMC);
    
    // Fill output tree
    tout->Fill();
  } // End of while loop

  // Print info
  cout << "All entries: " << ientry << endl;
  cout << "BeamCount: " << BeamCount << endl;
  
  // Kinematic Fitting for Pi0 shower
  //if(kMC) anaUtils.DoKinematicFitting();

 
  // Print cut flow statistics
  int icut = 0;
  double nsel = -999;

  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam PDG",  icut++), AnaIO::hCutBeamPDGPass, 1, 1, ientry);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Pandora slice",  icut++), AnaIO::hCutPandoraSlicePass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Calo size",  icut++), AnaIO::hCutCaloSizePass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam Quality",  icut++), AnaIO::hCutBeamQualityPass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. APA3 endZ",  icut++), AnaIO::hCutAPA3EndZPass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Michel score",  icut++), AnaIO::hCutMichelScorePass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Median dEdx",  icut++), AnaIO::hCutProtonChi2Pass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nshower",  icut++), AnaIO::hCutnshower, 2, 100000, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npi0",  icut++), AnaIO::hCutnpi0, 1, 100000, nsel);
  //nsel = plotUtils.PrintStat(tag+Form(" %d. KF Pass",  icut++), AnaIO::hKFPassRate, 1, 1, nsel);

  cout << "Shower Cuts: " << endl;
  int icut_shower = 0;
  double nsel_shower = -999;
  //nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower Pandora",  icut_shower++), AnaIO::hCutDaughterPandoraShowerPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower EM score",  icut_shower++), AnaIO::hCutDaughterShowerScorePass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower nhits",  icut_shower++), AnaIO::hCutDaughterShowernHitsPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower non-empty E",  icut_shower++), AnaIO::hCutDaughterShowerNonEmptyEPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower distance",  icut_shower++), AnaIO::hCutDaughterShowerDistPass, 1, 1, nsel_shower);
  nsel_shower = plotUtils.PrintStat(tag+Form(" %d. Shower IP",  icut_shower++), AnaIO::hCutDaughterShowerIPPass, 1, 1, nsel_shower);

  printf("End of %d cuts: %.1f selected\n", icut, nsel);

  printf("End of %d shower cuts: %.1f selected\n", icut_shower, nsel_shower);
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
  const double purity = nselected/ntotal;

  cout << "nAll: " << ntotal << " nSignal: " << nselected << endl;
  cout << "purity: " << purity << endl;

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
  // Run reco loop for both mc and data
  int nEntryToStop = -1; // Default is looping over all entries
  // Get the input break entries if provided 
  if(argc!=1) nEntryToStop = atoi(argv[1]);
  // MC Analysis
  double mcBeamCount = anaRec(mcfinName,mclout,"mc", nEntryToStop);
  // Data Analysis
  double dataBeamCount = anaRec(datafinName,datalout,"data", nEntryToStop);
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
