#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "PlotUtils.h"

namespace AnaIO
{
  // Declare variables
  int event;
  // Declare histograms
  TH1D * hevent = 0x0;
  // Get input tree
  TTree * GetInputTree(TFile * fin, const TString tname)
  {
    TTree * tree=(TTree*) fin->Get(tname);
    if(!tree){
      cout<<"no tree!"<<endl;
      gDirectory->ls();
      exit(1);
    }
    else cout << "Successuflly get tree with name: " << tname << endl;

    tree->SetBranchAddress("event", &event); 
    return tree;
  } // End of GetInputTree
  
  // Initialise reco histograms
  void IniRecHist(TList * lout, const TString tag)
  {
    hevent = new TH1D("hevent_"+tag, "", 10000, 0, 100000); lout->Add(hevent);
  }// End of IniRecHist

} // End of namespace
#endif
