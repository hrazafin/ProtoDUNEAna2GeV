void makeClass(){
  TFile* f = new TFile("protoDUNE_mc_reco_flattree_prod4a_whole.root");
  TTree *MyTree;
  f->GetObject("pduneana/beamana",MyTree);
  MyTree->MakeClass("MC_Prod4a");
}
