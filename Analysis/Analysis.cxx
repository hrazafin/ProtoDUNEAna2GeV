#include "include/AnaIO.h"
#include "src/PlotUtils.cxx"

void anaRec(const TString finName, TList *lout, const TString tag)
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
  AnaIO::IniRecHist(lout,tag);

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
  const TString datafinName  = "input/protoDUNE_data_reco_flattree.root";
  // Declare output list
  TList * mclout = 0x0;
  TList * datalout = 0x0;
  mclout = new TList;
  datalout = new TList;
  // Run reco loop for both mc and data
  anaRec(mcfinName,mclout,"mc");
  anaRec(datafinName,datalout,"data");
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
  // Draw all histograms in mclout and datalout
  PlotUtils PlotUtils;
  PlotUtils.DrawHist(mclout,"output");
  PlotUtils.DrawHist(datalout,"output");
}
