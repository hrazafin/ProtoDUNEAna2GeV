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
        anaUtils.DoTruthTKICalculation();
      }
    }

    //====================== Do cuts (both MC and data) ======================//
    // Do the beam cut
    if(!anaCut.CutBeamAllInOne(kMC)) continue;
    // Count beam after beam cut before other cuts
    BeamCount++; 
    // Fill beam info
    anaUtils.FillBeamKinematics(kMC);
    
    // Do event topology cut
    if(!anaCut.CutTopology(kMC)) continue;
    

    // Fill output tree
    tout->Fill();
  } // End of while loop

  cout << "All entries: " << ientry << endl;
  cout << "BeamCount: " << BeamCount << endl;
  
  // Print cut flow statistics
  int icut = 0;
  double nsel = -999;

  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam ID",  icut++), AnaIO::hCutBeamIDPass, 1, 1, ientry);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Pandora beam type",  icut++), AnaIO::hCutBeamType, 13, 13, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Beam Pos",  icut++), AnaIO::hCutBeamPosPass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. APA3",  icut++), AnaIO::hCutBeamEndZPass, 1, 1, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nproton", icut++), AnaIO::hCutnproton, 1, 100000, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nshower",  icut++), AnaIO::hCutnshower, 2, 100000, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Npiplus",  icut++), AnaIO::hCutnpiplus, 0, 0, nsel);
  nsel = plotUtils.PrintStat(tag+Form(" %d. Nmichel",  icut++), AnaIO::hCutnmichel, 0, 0, nsel);
  
  printf("End of %d cuts: %.1f selected\n", icut, nsel);

  const double nsig = AnaIO::hTruthSignal->GetBinContent(2);
  const double nbk = AnaIO::hTruthSignal->GetBinContent(1);
  const double nall = nsig+nbk;
  const double purity = nsig/nall;

  cout << "nAll: " << nall << " nSignal: " << nsig << " nBackground: " << nbk << endl;
  cout << "purity: " << purity << endl;
  return BeamCount;
  
} // End of anaRec

int main(int argc, char * argv[])
{
  // Initialise the input file name
  const TString mcfinName = "input/protoDUNE_mc_reco_flattree_prod4a.root";
  const TString datafinName  = "input/protoDUNE_data_reco_flattree_prod4a.root";
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
