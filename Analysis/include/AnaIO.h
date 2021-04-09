#ifndef _ANAIO_H_
#define _ANAIO_H_

#include "stdio.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"
#include <iostream>      //std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
using namespace std;

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
  void IniRecHist(TList * lout)
  {
    hevent = new TH1D("hevent", "", 1000, 0, 1000); lout->Add(hevent);
  }// End of IniRecHist

} // End of namespace
#endif
