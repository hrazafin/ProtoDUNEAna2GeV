#include "../include/Unfold.h"

Unfold::Unfold(int nb, double xlo, double xhi)
  : response_Int(nb, xlo, xhi,"response_Int")
  , response_BeamInt(nb, xlo, xhi,"response_BeamInt")
  , response_Inc(nb, xlo, xhi,"response_Inc")
  , response_Ini(nb, xlo, xhi,"response_Ini")
  , response_Pi0KE(nb, xlo, xhi,"response_Pi0KE")
  , response_Pi0CosTheta(10, -1., 1.,"response_Pi0CosTheta") 
  , response_Pi0Theta(10, 0., 180.,"response_Pi0Theta")

{

  response_Int.UseOverflow(false);
  response_BeamInt.UseOverflow(false);
  response_Inc.UseOverflow(false);
  response_Ini.UseOverflow(false);
  response_Pi0KE.UseOverflow(false);
  response_Pi0CosTheta.UseOverflow(false);

  eff_num_Int = new TH1D("eff_num_Int", "eff_num_Int", nb, xlo, xhi);
  eff_den_Int = new TH1D("eff_den_Int", "eff_den_Int", nb, xlo, xhi);
  eff_num_Inc = new TH1D("eff_num_Inc", "eff_num_Inc", nb, xlo, xhi);
  eff_den_Inc = new TH1D("eff_den_Inc", "eff_den_Inc", nb, xlo, xhi);
  eff_num_Ini = new TH1D("eff_num_Ini", "eff_num_Ini", nb, xlo, xhi);
  eff_den_Ini = new TH1D("eff_den_Ini", "eff_den_Ini", nb, xlo, xhi);

  eff_num_BeamInt = new TH1D("eff_num_BeamInt", "eff_num_BeamInt", nb, xlo, xhi);
  eff_den_BeamInt = new TH1D("eff_den_BeamInt", "eff_den_BeamInt", nb, xlo, xhi);

  eff_num_Pi0KE = new TH1D("eff_num_Pi0KE", "eff_num_Pi0KE", nb, xlo, xhi);
  eff_den_Pi0KE = new TH1D("eff_den_Pi0KE", "eff_den_Pi0KE", nb, xlo, xhi);

  eff_num_Pi0CosTheta = new TH1D("eff_num_Pi0CosTheta", "eff_num_Pi0CosTheta", 10, -1., 1.);
  eff_den_Pi0CosTheta = new TH1D("eff_den_Pi0CosTheta", "eff_den_Pi0CosTheta", 10, -1., 1.);

  eff_num_Pi0Theta = new TH1D("eff_num_Pi0Theta", "eff_num_Pi0Theta", 9, 0., 180.);
  eff_den_Pi0Theta = new TH1D("eff_den_Pi0Theta", "eff_den_Pi0Theta", 9, 0., 180.);


  pur_num_Int = new TH1D("pur_num_Int", "pur_num_Int", nb, xlo, xhi);
  pur_num_Inc = new TH1D("pur_num_Inc", "pur_num_Inc", nb, xlo, xhi);
  pur_den_Int = new TH1D("pur_den_Int", "pur_den_Int", nb, xlo, xhi);
  pur_den_Inc = new TH1D("pur_den_Inc", "pur_den_Inc", nb, xlo, xhi);


  eff_num_Int->Sumw2();
  eff_den_Int->Sumw2();
  eff_num_Inc->Sumw2();
  eff_den_Inc->Sumw2();
  pur_num_Int->Sumw2();
  pur_num_Inc->Sumw2();
  pur_den_Int->Sumw2();
  pur_den_Inc->Sumw2();

  eff_num_BeamInt->Sumw2();
  eff_den_BeamInt->Sumw2();
  eff_num_Ini->Sumw2();
  eff_den_Ini->Sumw2();
  eff_num_Pi0KE->Sumw2();
  eff_den_Pi0KE->Sumw2();
  eff_num_Pi0CosTheta->Sumw2();
  eff_den_Pi0CosTheta->Sumw2();
  eff_num_Pi0Theta->Sumw2();
  eff_den_Pi0Theta->Sumw2();
}  

void Unfold::SaveHistograms(){

  //eff_num_Int->Write("eff_num_Int");
  //eff_den_Int->Write("eff_den_Int");
  //eff_num_Inc->Write("eff_num_Inc");
  //eff_den_Inc->Write("eff_den_Inc");
  //pur_num_Int->Write("pur_num_Int");
  //pur_num_Inc->Write("pur_num_Inc");
  //pur_den->Write("pur_den");

  eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
  eff_Int->Divide(eff_den_Int);
  //eff_Int->Write("eff_Int");
  eff_Int->SetName("i910eff_Int");

  eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
  eff_Inc->Divide(eff_den_Inc);
  //eff_Inc->Write("eff_Inc");
  eff_Inc->SetName("i911eff_Inc");

  eff_Ini = (TH1D*)eff_num_Ini->Clone("eff_Ini");
  eff_Ini->Divide(eff_den_Ini);
  //eff_Inc->Write("eff_Inc");
  eff_Ini->SetName("i912eff_Ini");

  eff_BeamInt = (TH1D*)eff_num_BeamInt->Clone("eff_BeamInt");
  eff_BeamInt->Divide(eff_den_BeamInt);
  //eff_Int->Write("eff_Int");
  eff_BeamInt->SetName("i913eff_BeamInt");

  eff_Pi0KE = (TH1D*)eff_num_Pi0KE->Clone("eff_Pi0KE");
  eff_Pi0KE->Divide(eff_den_Pi0KE);
  //eff_Int->Write("eff_Int");
  eff_Pi0KE->SetName("i914eff_Pi0KE");


  eff_Pi0CosTheta = (TH1D*)eff_num_Pi0CosTheta->Clone("eff_Pi0CosTheta");
  eff_Pi0CosTheta->Divide(eff_den_Pi0CosTheta);
  //eff_Int->Write("eff_Int");
  eff_Pi0CosTheta->SetName("i915eff_Pi0CosTheta");

  eff_Pi0Theta = (TH1D*)eff_num_Pi0Theta->Clone("eff_Pi0Theta");
  eff_Pi0Theta->Divide(eff_den_Pi0Theta);
  //eff_Int->Write("eff_Int");
  eff_Pi0Theta->SetName("i916eff_Pi0Theta");


  pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
  pur_Int->Divide(pur_den_Int);
  //pur_Int->Write("pur_Int");

  pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
  pur_Inc->Divide(pur_den_Inc);
  //pur_Inc->Write("pur_Inc");

  TH2D *hint = (TH2D*)response_Int.Hresponse();
  hint->SetTitle("Interactions;Reco Pion Energy (MeV);True Pion Energy (MeV)");
  hint->SetName("i901response_Int");
  TH2D *hbeamint = (TH2D*)response_BeamInt.Hresponse();
  hbeamint->SetTitle("Interactions;Reco Pion Energy (MeV);True Pion Energy (MeV)");
  hbeamint->SetName("i902response_BeamInt");
  //hbeamint->RebinX(50);
  //hbeamint->RebinY(50);

  //hint->Write("hresponse_Int");
  TH2D *hinc = (TH2D*)response_Inc.Hresponse();
  hinc->SetTitle("Incidents; Reco Pion Energy (MeV); True Pion Energy (MeV)");
  hinc->SetName("i903response_Inc");

  //hinc->Write("hresponse_Inc");
  TH2D *hini = (TH2D*)response_Ini.Hresponse();
  hini->SetTitle("Initial; Reco Pion Energy (MeV); True Pion Energy (MeV)");
  hinc->SetName("i904response_Ini");
  //hini->Write("hresponse_Ini");

  //response_Int.Write("response_Int");
  //response_Inc.Write("response_Inc");
  //response_Ini.Write("response_Ini");

  TH2D *hpi0KE = (TH2D*)response_Pi0KE.Hresponse();
  hpi0KE->SetTitle("Pi0 KE;Reco Pion KE (MeV);True Pion KE (MeV)");
  hpi0KE->SetName("i905response_pi0KE");
  //hint->Write("hresponse_Int");
  TH2D *hpi0CosTheta = (TH2D*)response_Pi0CosTheta.Hresponse();
  hpi0CosTheta->SetTitle("Pi0 CosTheta;Reco Pion CosTheta;True Pion CosTheta");
  hpi0CosTheta->SetName("i906response_pi0CosTheta");

  TH2D *hpi0Theta = (TH2D*)response_Pi0Theta.Hresponse();
  hpi0Theta->SetTitle("Pi0 Theta;Reco Pion Theta;True Pion Theta");
  hpi0Theta->SetName("i907response_pi0Theta");
}
