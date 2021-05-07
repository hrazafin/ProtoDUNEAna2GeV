#include "include/AnaIO.h"
#include "src/PlotUtils.cxx"
#include "src/AnaUtils.cxx"
#include "src/AnaCut.cxx"

int anaRec(const TString finName, TList *lout, const TString tag, const int nEntryToStop)
{
  //_____________________________________________________ basic settings _____________________________________________________ 
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
  TTree  * tree = AnaIO::GetInputTree(fin, "pionana/beamana", tag);
  // Initialise reco histograms
  AnaIO::IniHist(lout, kMC);

  AnaUtils anaUtils;
  AnaCut anaCut;
  // Initialise entry and beam counter
  int ientry = 0;
  int BeamCount = 0;
  // Loop over TTree
  while(tree->GetEntry(ientry)){
    AnaIO::hEvent->Fill(AnaIO::event); 
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
    //anaCut.CountPFP(kMC,true);
  } // End of while loop
  return BeamCount;
} // End of anaRec

int main(int argc, char * argv[])
{
  // Initialise the input file name
  const TString mcfinName = "input/protoDUNE_mc_reco_flattree.root";
  const TString datafinName  = "input/protoDUNE_data_reco_flattree.root";
  // Declare output list
  TList * mclout = 0x0;
  TList * datalout = 0x0;
  mclout = new TList;
  datalout = new TList;
  // Run reco loop for both mc and data
  int nEntryToStop = -1; // Default is looping over all entries
  // Get the input break entries
  if(argc!=1) nEntryToStop = atoi(argv[1]);
  double mcBeamCount = anaRec(mcfinName,mclout,"mc", nEntryToStop);
  double dataBeamCount = anaRec(datafinName,datalout,"data", nEntryToStop);
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
  plotUtils.DrawHist(mclout,plotScale,datalout,"output");

}
