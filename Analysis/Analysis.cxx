#include "include/AnaIO.h"

void anaRec(TString finName, TList *lout)
{
  //_____________________________________________________ basic settings _____________________________________________________ 
  cout << "Input file:" << endl;
  gSystem->Exec(Form("readlink -f %s", finName.Data()));
  cout << endl;
  // Open the input file
  TFile * fin = new TFile(finName);
  if(!fin->IsOpen()){
    cout << "fin not open!" << endl;
    exit(1);
  }
  else cout << "fin is open!" << endl;

  // Get the TTree from input file
  TTree  * tree = AnaIO::GetInputTree(fin, "pionana/beamana");
  // Initialise reco histograms
  AnaIO::IniRecHist(lout);

  // Initialise entry
  int ientry = 0;
  // Loop over TTree
  while(tree->GetEntry(ientry)){
    
    AnaIO::hevent->Fill(AnaIO::event); 
    // Break 
    if(ientry>=100){
      cout << "Breaking after " << ientry << endl;
      break;
    }
    // update ientry after each loop
    ientry++;
  } // End of while loop

} // End of anaRec

int main()
{
  // Initialise the input file name
  const TString mcfinName = "input/protoDUNE_mc_reco_flattree.root";
  
  // Declare mc list
  TList * mclout = 0x0;
  mclout = new TList;
  // Run reco loop
  anaRec(mcfinName,mclout);
  // Declare output root file
  TFile * fout = new TFile("output/outana.root","recreate");
  // Create mc subdirectory
  TDirectory * ld = gDirectory->mkdir("mc");
  ld->cd();
  // Write list to the root file
  mclout->Write();
  gDirectory->cd("../");
  // Save the info
  fout->Save();
  fout->Close();
}
