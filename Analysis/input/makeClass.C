void makeClass(){
  TFile* f = new TFile("protoDUNE_mc_reco_flattree_prod4a_ntuple_6GeV.root");
  TTree *MyTree;
  f->GetObject("pduneana/beamana",MyTree);
  MyTree->MakeClass("MC_Prod4a_ntuple_6GeV");
}
