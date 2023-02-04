#include "SliceParams.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TemplateFitter.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "json/json.h"
#include <fstream>
#include <iostream>
#include <string>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TH2D.h"

static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name <<" <option(s)>\n"
            << "Options:\n"
            << "\t-h,--help\t\tShow this help message.\n"
            << "\t-c config.json\t\tSpecify configuration file."
            << std::endl;
}

int main(int argc, char** argv){

  //bool fitfakedata = false;

  bool found_config = false;

  string config_file;

  //for (int i = 1; i < argc; ++i) {
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     config_file = argv[++iArg];
     found_config = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      show_usage(argv[0]);
      return 1;
    }
  }

  if (!found_config){
    show_usage(argv[0]);
    return 1;
  }

  Json::Value root;
  ifstream file(config_file);
  file >> root;
  cout<<root<<endl;

  TFile *fdata = TFile::Open(root["datafile"].asString().c_str());
  TFile *fmc   = TFile::Open(root["mcfile"].asString().c_str());
  TFile *fbkg  = TFile::Open(root["bkgfile"].asString().c_str());
  TFile *fout  = TFile::Open(root["outfile"].asString().c_str(), "recreate");

  //TVectorD sf(2); sf[0] = 1.; sf[1] = 0.;
  TVectorD *sf_mu = (TVectorD*)fbkg->Get("sf_mu");
  TVectorD *sf_p = (TVectorD*)fbkg->Get("sf_p");
  TVectorD *sf_spi = (TVectorD*)fbkg->Get("sf_spi");

  cout<<"Muon scaling factor: "<<(*sf_mu)[0]<<"+-"<<(*sf_mu)[1]<<endl;
  cout<<"Proton scaling factor: "<<(*sf_p)[0]<<"+-"<<(*sf_p)[1]<<endl;
  cout<<"Pion scaling factor: "<<(*sf_spi)[0]<<"+-"<<(*sf_spi)[1]<<endl;

  TH1D *hsliceID[pi::nIntTypes+1];
  TH1D *hdata = new TH1D("hdata","Data;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1); // h_recosliceid_allevts_cuts (hreco_sliceID_6_0)
  TH1D *hproton = new TH1D("hproton","Proton background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hmu = new TH1D("hmu","Muon background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hspi = new TH1D("hspi","Secondary pion background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hpiel = new TH1D("hpiel","Pion elastic;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hother = new TH1D("hother","Other backgrounds;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  hdata->Sumw2();
  hproton->Sumw2();
  hmu->Sumw2();
  hspi->Sumw2();
  hpiel->Sumw2();
  hother->Sumw2();

  // first loop to get total number of events
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    if (i==0){
      hsliceID[i] = (TH1D*)fdata->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
    }
    else {
      hsliceID[i] = (TH1D*)fmc->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
    }
  }
  TH1D *hmc = new TH1D("hmc","MC;Slice ID;Events",hsliceID[0]->GetNbinsX(),-1,hsliceID[0]->GetNbinsX()-1);
  hmc->Sumw2();
  for (int j = 0; j < pi::nthinslices+2; ++j){
    int bin = hsliceID[0]->FindBin(j-0.5);
    hdata->SetBinContent(j+1, hsliceID[0]->GetBinContent(bin));
    hdata->SetBinError(j+1, hsliceID[0]->GetBinError(bin));
    double nmc = 0;
    for (int i = 1; i < pi::nIntTypes+1; ++i){
      nmc += hsliceID[i]->GetBinContent(bin);
    }
    hmc->SetBinContent(j+1, nmc);
    hmc->SetBinError(j+1, sqrt(nmc));
  }

  // second loop to fill histograms
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    if (i==0){
      hsliceID[i] = (TH1D*)fdata->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
    }
    else {
      hsliceID[i] = (TH1D*)fmc->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
      //hsliceID[i]->Scale(totaldata/totalmc);
      hsliceID[i]->Multiply(hsliceID[0]);
      hsliceID[i]->Divide(hmc);
    }
    for (int j = 0; j < pi::nthinslices+2; ++j){
      int bin = hsliceID[i]->FindBin(j-0.5);
      if (i!=0){
        if (i == pi::kMuon || i == pi::kMIDmu){
          double binc = hmu->GetBinContent(bin);
          double bine = hmu->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin)*(*sf_mu)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin)*(*sf_mu)[0],2)
                      + pow(hsliceID[i]->GetBinContent(bin)*(*sf_mu)[1],2));
          hmu->SetBinContent(j+1, binc);
          hmu->SetBinError(j+1, bine);
        }
        else if (i == pi::kMIDp){
          double binc = hproton->GetBinContent(bin);
          double bine = hproton->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin)*(*sf_p)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin)*(*sf_p)[0],2)
                      + pow(hsliceID[i]->GetBinContent(bin)*(*sf_p)[1],2));
          hproton->SetBinContent(j+1, binc);
          hproton->SetBinError(j+1, bine);
        }
        else if (i == pi::kMIDpi){
          double binc = hspi->GetBinContent(bin);
          double bine = hspi->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin)*(*sf_spi)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin)*(*sf_spi)[0],2)
                      + pow(hsliceID[i]->GetBinContent(bin)*(*sf_spi)[1],2));
          hspi->SetBinContent(j+1, binc);
          hspi->SetBinError(j+1, bine);
        }
        else if (i == pi::kPiElas){
          double binc = hpiel->GetBinContent(bin);
          double bine = hpiel->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin),2));
          hpiel->SetBinContent(j+1, binc);
          hpiel->SetBinError(j+1, bine);
        }
        else if (i != pi::kPiInel){
          double binc = hother->GetBinContent(bin);
          double bine = hother->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin),2));
          hother->SetBinContent(j+1, binc);
          hother->SetBinError(j+1, bine);
        }
      }
    }
  }

  TH1D *hsiginc = (TH1D*)hdata->Clone("hsiginc");
  hsiginc->SetTitle("Pion incident signal;Slice ID;Events");
  hsiginc->Add(hmu,-1);
  hsiginc->Add(hproton,-1);
  hsiginc->Add(hspi,-1);
  hsiginc->Add(hother,-1);
  TH1D *hsignal = (TH1D*)hsiginc->Clone("hsignal");
  hsignal->SetTitle("Pion interaction signal;Slice ID;Events");
  hsignal->Add(hpiel,-1);
  
  // unfolding
  RooUnfoldResponse *response_SliceID_Inc = (RooUnfoldResponse*)fmc->Get("response_SliceID_Inc");
  RooUnfoldBayes unfold_Inc (response_SliceID_Inc, hsiginc, 4);
  RooUnfoldResponse *response_SliceID_Int = (RooUnfoldResponse*)fmc->Get("response_SliceID_Int");
  RooUnfoldBayes unfold_Int (response_SliceID_Int, hsignal, 4);
  
  TH1D *hsiginc_uf;
  TH1D *hsignal_uf;
  //hsignal_uf->Sumw2();
  hsiginc_uf = (TH1D*)unfold_Inc.Hreco();
  hsiginc_uf->SetNameTitle("hsiginc_uf", "Unfolded incident signal;Slice ID;Events");
  hsignal_uf = (TH1D*)unfold_Int.Hreco();
  hsignal_uf->SetNameTitle("hsignal_uf", "Unfolded interaction signal;Slice ID;Events");
  
  // get Ninc and Nint
  double Ninc[pi::nthinslices] = {0};
  double Nint[pi::nthinslices] = {0};
  double err_inc[pi::nthinslices] = {0};
  double err_int[pi::nthinslices] = {0};
  double SliceID[pi::nthinslices] = {0};

  for (int i = 0; i<pi::nthinslices; ++i){
    SliceID[i] = i+1;
    Nint[i] = hsignal_uf->GetBinContent(i+2);
    err_int[i] = hsignal_uf->GetBinError(i+2);
    for (int j = i; j<=pi::nthinslices; ++j){
      Ninc[i] += hsiginc_uf->GetBinContent(j+2);
      err_inc[i] += pow(hsiginc_uf->GetBinError(j+2),2);
    }
    err_inc[i] = sqrt(err_inc[i]);
  }
  TGraphErrors *gr_inc = new TGraphErrors(pi::nthinslices, SliceID, Ninc, 0, err_inc);
  gr_inc->SetNameTitle("gr_inc", "Incident number;Slice ID;Events");
  gr_inc->Write();
  TGraphErrors *gr_int = new TGraphErrors(pi::nthinslices, SliceID, Nint, 0, err_int);
  gr_int->SetNameTitle("gr_int", "Interaction number;Slice ID;Events");
  gr_int->Write();
  TGraphErrors *gr_trueincE = (TGraphErrors*)fmc->Get("gr_trueincE");
  gr_trueincE->SetNameTitle("gr_trueincE", "True incident energy;Slice ID;Energy (MeV)");
  gr_trueincE->Write();
  TGraphErrors *gr_recoincE = (TGraphErrors*)fdata->Get("gr_recoincE");
  gr_recoincE->SetNameTitle("gr_recoincE", "Reco incident energy;Slice ID;Energy (MeV)");
  gr_recoincE->Write();
  
  // Calculate cross-section
  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // g/cm^3
  double xs[pi::nthinslices] = {0};
  double err_xs[pi::nthinslices] = {0};
  double incE[pi::nthinslices] = {0};
  double err_incE[pi::nthinslices] = {0};
  for (int i = 0; i<pi::nthinslices; ++i){
    xs[i] = MAr/(Density*NA*pi::thinslicewidth)*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
    err_xs[i] = MAr/(Density*NA*pi::thinslicewidth)*1e27*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2));
    incE[i] = gr_trueincE->GetPointY(i);
    err_incE[i] = gr_trueincE->GetErrorY(i);
  }
  TGraphErrors *gr_recoxs = new TGraphErrors(pi::nthinslices, incE, xs, err_incE, err_xs);
  gr_recoxs->SetNameTitle("gr_recoxs", "Reco cross-section;Energy (MeV); Cross-section (mb)");
  gr_recoxs->Write();
  /*TFile f2("../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
  TCanvas *c1 = new TCanvas("c1", "c1");
  gr_recoxs->SetTitle("Pion Inelastic Cross Section");
  gr_recoxs->GetXaxis()->SetTitle("Pion Kinetic Energy (MeV)");
  gr_recoxs->GetXaxis()->SetRangeUser(360, 900);
  gr_recoxs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
  gr_recoxs->GetYaxis()->SetRangeUser(400, 900);
  gr_recoxs->SetLineWidth(2);
  gr_recoxs->Draw("ape");
  total_inel_KE->SetLineColor(2);
  total_inel_KE->Draw("c");
  c1->Print("recoxs.png");*/
  
  // test sample validation
  TH1D *hval_siginc_reco = (TH1D*)fmc->Get("h_recosliceid_pion_cuts");
  hval_siginc_reco->Write("hval_siginc_reco");
  TH1D *hval_trueinc = (TH1D*)fmc->Get("h_truesliceid_pion_all");
  hval_trueinc->Write("hval_trueinc");
  TH1D *hval_signal_reco = (TH1D*)fmc->Get("h_recosliceid_pioninelastic_cuts");
  hval_signal_reco->Write("hval_signal_reco");
  TH1D *hval_trueint = (TH1D*)fmc->Get("h_truesliceid_pioninelastic_all");
  hval_trueint->Write("hval_trueint");
  double Ninc_t[pi::nthinslices] = {0};
  double Nint_t[pi::nthinslices] = {0};
  double err_inc_t[pi::nthinslices] = {0};
  double err_int_t[pi::nthinslices] = {0};
  for (int i = 0; i<pi::nthinslices; ++i){
    Nint_t[i] = hval_trueint->GetBinContent(i+2);
    err_int_t[i] = hval_trueint->GetBinError(i+2);
    for (int j = i; j<=pi::nthinslices; ++j){
      Ninc_t[i] += hval_trueinc->GetBinContent(j+2);
      err_inc_t[i] += pow(hval_trueinc->GetBinError(j+2),2);
    }
    err_inc_t[i] = sqrt(err_inc_t[i]);
  }
  double xs_t[pi::nthinslices] = {0};
  double err_xs_t[pi::nthinslices] = {0};
  for (int i = 0; i<pi::nthinslices; ++i){
    xs_t[i] = MAr/(Density*NA*pi::thinslicewidth)*log(Ninc_t[i]/(Ninc_t[i]-Nint_t[i]))*1e27;
    err_xs_t[i] = MAr/(Density*NA*pi::thinslicewidth)*1e27*sqrt(pow(Nint_t[i]*err_inc_t[i]/Ninc_t[i]/(Ninc_t[i]-Nint_t[i]),2)+pow(err_int_t[i]/(Ninc_t[i]-Nint_t[i]),2));
  }
  TGraphErrors *gr_truexs = new TGraphErrors(pi::nthinslices, incE, xs_t, 0, err_xs_t);
  gr_truexs->SetNameTitle("gr_truexs", "Reco cross-section;Energy (MeV); Cross-section (mb)");
  gr_truexs->Write();
  TGraphErrors *gr_truexs_allMC = (TGraphErrors*)fmc->Get("gr_truexs");
  gr_truexs_allMC->Write("gr_truexs_allMC");
  
  fout->Write();
  fout->Close();
  return 0;
}
