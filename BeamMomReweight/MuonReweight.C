#include "include/TemplateFitter.h"
#include "include/TemplateFitter.cxx"

const int typenum = 8;
const int boundbin = 22; // frontier of first TPC bin
const int uppbound = 50; // frontier of all interested TPC bin
const double xbins[boundbin+2] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,10.*uppbound};

double GetChi2(TH1D * hdata, TH1D * hmc){
    int nx = hdata->GetNbinsX();
    double chisq = 0;
    for(int ix=1; ix<=nx; ix++){
      double x = hdata->GetBinContent(ix);
      double y = hmc->GetBinContent(ix);
        
      double ex = hdata->GetBinError(ix);
      double ey = hmc->GetBinError(ix);

      if (!ex) ex = 1;
      if (!ey) ey = 1;
      chisq += pow(x-y,2)/(pow(ex,2)+pow(ey,2));
    }
    return chisq;
}


double calchi2(double muweight, TH1D *data_dist, TH1D *mc_dist_sep[typenum], TH1D *hdata, TH1D *hmc, TH1D *hmc0, bool save_plot=false){

  const double typeweight[typenum] = {1, muweight, 1,1,1,1,1,1};
  //double mc_inte = mc_dist->Integral(0, uppbound);
  double data_inte = data_dist->Integral(0, uppbound);

  double mc_inte_sep;
  double mc_inte_sep0;
  double mc_inte_sep_sq;
  for (int j=0; j<typenum; ++j){
    mc_inte_sep += mc_dist_sep[j]->Integral(0, uppbound)*typeweight[j];
    mc_inte_sep0 += mc_dist_sep[j]->Integral(0, uppbound);
    mc_inte_sep_sq += mc_dist_sep[j]->Integral(0, uppbound)*pow(typeweight[j],2);
  }
  //cout<<mc_inte<<"\t"<<mc_inte_sep<<"\t"<<data_inte<<endl;
  double mcnorm = data_inte/mc_inte_sep;
  double mcnorm0 = data_inte/mc_inte_sep0;
  double sfactor = mc_inte_sep_sq/mc_inte_sep; //sfactor = 1.0;
  double fom = 0;
  double chi = 0;
  double mcbin;
  double mcbin0;
  int initi = 16;
  for (int i=initi; i<=boundbin; ++i){
    //mcbin = mc_dist->GetBinContent(i+1) * mcnorm;
    mcbin = 0;
    mcbin0 = 0;
    for (int j=0; j<typenum; ++j){
      mcbin += mc_dist_sep[j]->GetBinContent(i+1)*typeweight[j];
      mcbin0 += mc_dist_sep[j]->GetBinContent(i+1);
    }
    mcbin *= mcnorm;
    mcbin0 *= mcnorm0;
    if(data_dist->GetBinContent(i+1) != 0){
        chi = (data_dist->GetBinContent(i+1) - mcbin)/sqrt(sfactor*data_dist->GetBinContent(i+1) + mcbin);
        fom += chi*chi;
    }
    

    //cout<<i<<endl;
    //cout<<mcbin<<endl;
    //cout<<data_dist->GetBinError(i)<<endl;
    if (save_plot) {
      hdata->SetBinContent(i, data_dist->GetBinContent(i+1)/10);
      hdata->SetBinError(i, data_dist->GetBinError(i+1)/10);
      hmc->SetBinContent(i, mcbin/10);
      hmc0->SetBinContent(i, mcbin0/10);
    }
  }
  double over_data = data_dist->Integral(23, uppbound);
  //mcbin = mc_dist->Integral(23, uppbound) * mcnorm;
  mcbin = 0;
  mcbin0 = 0;
  for (int j=0; j<typenum; ++j){
    mcbin += mc_dist_sep[j]->Integral(23, uppbound)*typeweight[j];
    mcbin0 += mc_dist_sep[j]->Integral(23, uppbound);
  }
  mcbin *= mcnorm;
  mcbin0 *= mcnorm0;
  if(over_data != 0){
    chi = (over_data - mcbin)/sqrt(over_data + mcbin);
    fom += chi*chi;
    fom /= (boundbin-initi+2);
  }
  if (save_plot) {
    cout<<"\nWeight = "<<muweight<<"\t FOM = "<<fom<<endl;

     cout << "error line: " << fom+1./(boundbin+1) << endl;
    
    hdata->SetBinContent(boundbin+1, over_data/280);
    hdata->SetBinError(boundbin+1, sqrt(over_data)/280);
    hmc->SetBinContent(boundbin+1, mcbin/280);
    hmc0->SetBinContent(boundbin+1, mcbin0/280);
    
    //TCanvas* c1 = new TCanvas("c1","c1");
    TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
    hmc->SetLineColor(kRed);
    hmc0->SetLineColor(kBlue);
    hdata->GetYaxis()->SetRangeUser(0, 1.3*hmc0->GetBinContent(hmc0->GetMaximumBin()));

    hdata->Draw();
    hdata->GetXaxis()->SetTitle("Track length (cm)");
    hdata->GetYaxis()->SetTitle("Bin Normalised");
    hdata->SetTitle(" ");
    hdata->GetYaxis()->CenterTitle();
    hdata->GetYaxis()->SetTitleFont(22);
    hdata->GetYaxis()->SetTitleSize(0.05);
    hdata->GetYaxis()->SetTitleOffset(0.9);
    hdata->GetXaxis()->CenterTitle();
    hdata->GetXaxis()->SetTitleFont(22);
    hdata->GetXaxis()->SetTitleSize(0.05);
    hdata->GetXaxis()->SetTitleOffset(0.9);
    hdata->SetLineColor(kBlack);
    hdata->SetFillStyle(0);
    hdata->SetMarkerStyle(8);
    hdata->SetMarkerSize(1);
    hdata->SetMarkerColor(kBlack);
    hdata->Draw("e1");

    hmc0->SetLineColor(kGreen+3);
    hmc0->SetLineStyle(7);
    hmc0->SetLineWidth(2);

    hmc->SetLineColor(kRed);
    hmc->SetLineWidth(2);
    hmc->SetFillStyle(0);

    hmc->Draw("hist same");
    hmc0->Draw("hist same");
    hdata->SetTitle(" ");

    auto lgfit = new TLegend(0.58,0.60,0.88,0.88);
    //lgfit->SetHeader(Form("#chi^{2}: %f",GetChi2(h0,hmc_fit)));
    lgfit->SetHeader(Form("#chi^{2}_{def}: %.2f, #chi^{2}_{fit}: %.2f", GetChi2(hdata,hmc0), GetChi2(hdata,hmc)));
    lgfit->AddEntry(hdata,"data","lp");
    lgfit->AddEntry(hmc0,"MC Default","l");
    lgfit->AddEntry(hmc,"MC Reweighted","l");
    lgfit->SetBorderSize(0);
    lgfit->Draw("sames");

    TLatex tt;
    tt.SetNDC();

    tt.DrawLatex(0.135,0.825,Form("#color[4]{#bf{#it{#alpha_{#mu} = %.2f#pm0.02}}}",muweight));

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    tt.DrawLatex(0.705,0.925,"Long Track Sample");

    c1->Print("output/muon_reweight.pdf");
  }
  return fom;
}


void MuonReweight(){

    double plotScale = 0.698653;
    gStyle->SetOptStat(0);

    TString finName = "output/outana.root";

    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
        cout << "file file not open" << endl;
        exit(1);
    }

    // Get Histograms
    TH2D * h2d_LongtrkLen_mc = (TH2D*) file->Get("mc/a005hBeamTrackLength_STK;1");
    TH2D * h2d_LongtrkLen_data = (TH2D*) file->Get("data/a005hBeamTrackLength_STK;1");

    const double LongTrkBin[] = {-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 500};
    TH1D * h1d_LongtrkLen_mc_muon = new TH1D("h1d_LongtrkLen_mc_muon",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    TH1D * h1d_LongtrkLen_mc_other = new TH1D("h1d_LongtrkLen_mc_other",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    TH1D * h1d_LongtrkLen_data = new TH1D("h1d_LongtrkLen_data",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);

    int nx = h2d_LongtrkLen_mc->GetNbinsX();
    int ny = h2d_LongtrkLen_mc->GetNbinsY();
    cout << "nx: " << nx << endl;
    for(int ix=0; ix<=nx+1; ix++){

        double binwidth = 10;
        if(ix == 24) binwidth = 280;
        

        // Data
        double ientry_data = h2d_LongtrkLen_data->GetBinContent(ix, ny);
        double ientry_dataerr = h2d_LongtrkLen_data->GetBinError(ix, ny);

        // MC component
        double ientry_mc_muon = h2d_LongtrkLen_mc->GetBinContent(ix, 1);
        double ientry_mc_other = h2d_LongtrkLen_mc->GetBinContent(ix, 0) + h2d_LongtrkLen_mc->GetBinContent(ix, 2) + h2d_LongtrkLen_mc->GetBinContent(ix, 3) + h2d_LongtrkLen_mc->GetBinContent(ix, 4) + h2d_LongtrkLen_mc->GetBinContent(ix, 5);
    
        // MC component error
        double ientry_mc_muonerr = h2d_LongtrkLen_mc->GetBinError(ix, 1);
        double ientry_mc_othererr = sqrt(pow(h2d_LongtrkLen_mc->GetBinError(ix, 0),2) + pow(h2d_LongtrkLen_mc->GetBinError(ix, 2),2) + pow(h2d_LongtrkLen_mc->GetBinError(ix, 3),2) + pow(h2d_LongtrkLen_mc->GetBinError(ix, 4),2) + pow(h2d_LongtrkLen_mc->GetBinError(ix, 5),2));

        h1d_LongtrkLen_data->SetBinContent(ix, ientry_data/binwidth);
        h1d_LongtrkLen_data->SetBinError(ix, ientry_dataerr/binwidth);

        h1d_LongtrkLen_mc_muon->SetBinContent(ix, ientry_mc_muon/binwidth);
        h1d_LongtrkLen_mc_muon->SetBinError(ix, ientry_mc_muonerr/binwidth);

        h1d_LongtrkLen_mc_other->SetBinContent(ix, ientry_mc_other/binwidth);
        h1d_LongtrkLen_mc_other->SetBinError(ix, ientry_mc_othererr/binwidth);

    }

    TemplateFitter fitter;   
    h1d_LongtrkLen_mc_muon->Scale(plotScale); h1d_LongtrkLen_mc_other->Scale(plotScale);

    fitter.SetHistograms(h1d_LongtrkLen_data, h1d_LongtrkLen_mc_other, h1d_LongtrkLen_mc_muon);
    fitter.SetFitRange(1, nx+1);
    fitter.Fit();
    double par = fitter.GetPar();
    double parerr = fitter.GetParError();
    cout << "===== Muon Reweight ==== " << endl;
    cout << par << " " << parerr << endl;


    TCanvas * c11 = new TCanvas("c11", "", 1200, 800);
    //h1d_LongtrkLen_mc_muon->Scale(plotScale);
    //h1d_LongtrkLen_mc_other->Scale(plotScale);

    //h1d_LongtrkLen_mc_muon->Add(h1d_LongtrkLen_mc_other);
    //h1d_LongtrkLen_mc_muon->Draw("hists");

    //h1d_LongtrkLen_mc_other->Draw("hists sames");


    h1d_LongtrkLen_data->SetMarkerStyle(25);
    h1d_LongtrkLen_data->SetMarkerSize(1);
    h1d_LongtrkLen_data->SetMarkerColor(kBlack);
    h1d_LongtrkLen_data->SetLineColor(kBlack);
    h1d_LongtrkLen_data->SetTitle(" ");
    h1d_LongtrkLen_data->GetYaxis()->CenterTitle();
    h1d_LongtrkLen_data->GetYaxis()->SetTitleFont(22);
    h1d_LongtrkLen_data->GetYaxis()->SetTitleSize(0.05);
    h1d_LongtrkLen_data->GetYaxis()->SetTitleOffset(0.9);
    h1d_LongtrkLen_data->GetXaxis()->CenterTitle();
    h1d_LongtrkLen_data->GetXaxis()->SetTitleFont(22);
    h1d_LongtrkLen_data->GetXaxis()->SetTitleSize(0.05);
    h1d_LongtrkLen_data->GetXaxis()->SetTitleOffset(0.9);
    h1d_LongtrkLen_data->SetMaximum(h1d_LongtrkLen_data->GetMaximum()*1.5);

    h1d_LongtrkLen_data->Draw("e1");

    TH1D * hmc_fit = (TH1D*)h1d_LongtrkLen_mc_other->Clone();

    hmc_fit->Add(h1d_LongtrkLen_mc_muon,1.71);
    
    hmc_fit->SetLineColor(kRed);
    hmc_fit->SetLineWidth(3);
    
    hmc_fit->Draw("hists sames");

    TH1D * hmc_def = (TH1D*)h1d_LongtrkLen_mc_other->Clone();
    hmc_def->Add(h1d_LongtrkLen_mc_muon,1.0);

    hmc_def->SetLineColor(kGreen+3);
    hmc_def->SetLineStyle(7);
    hmc_def->SetLineWidth(3);
    
    hmc_def->Draw("hists sames");

    //c11->Print("output/hLongtrkLen.pdf");


    // ================ Muon reweight =============== //

    TH1D *data_dist = (TH1D*) h2d_LongtrkLen_data->ProjectionX(Form("hdata_trklen_%d", ny),ny,ny);  
    
    //TH2D * h2d_LongtrkLen_mc = (TH2D*) file->Get("mc/a005hBeamTrackLength_STK;1");
    cout << "data_dist->Integral(): " << data_dist->Integral() << endl;

    TH1D *mc_dist_sep[typenum];
    for (int j=0; j<typenum; ++j){
        mc_dist_sep[j] = h2d_LongtrkLen_mc->ProjectionX(Form("hreco_trklen_%d", j+1),j+1,j+1);
        cout << "mc_dist_sep[j]->Integral(): " << mc_dist_sep[j]->Integral() << endl;
    }

    /*const double LongTrkBin[] = {-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 500};
    TH1D * hmc0 = new TH1D("hmc0",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    TH1D * hmc = new TH1D("hmc_muon",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    TH1D * hdata = new TH1D("hdata",";Track Length (cm); Candidates", sizeof(LongTrkBin)/sizeof(double)-1, LongTrkBin);
    */
    TH1D *hdata = new TH1D("hdata", "hdata", boundbin+1, xbins);
    TH1D *hmc = new TH1D("hmc", "hmc", boundbin+1, xbins);
    TH1D *hmc0 = new TH1D("hmc0", "hmc0", boundbin+1, xbins);

    vector<double> fom_list;
    double selw = 1.;
    double minfom = 999.;
    TGraph* cur = new TGraph();
    int pi = 0;
    double minw = 1.3;
    double maxw = 1.9;
    double stepw = 0.01;
    for (double w = minw; w <= maxw; w+=stepw) {
        double fom = calchi2(w, data_dist, mc_dist_sep, hdata, hmc, hmc0);
        fom_list.push_back(fom);
        //cout<<fom<<", "; // python list format for plotting
        cur->SetPoint(pi, w, fom);
        ++pi;
        cout<<"weight = "<<w<<"\t fom = "<<fom<<endl;
        if (fom<minfom) {
            minfom = fom;
            selw = w;
        }
    }
    calchi2(selw, data_dist, mc_dist_sep, hdata, hmc, hmc0, true);
  
    TCanvas* c2 = new TCanvas("c2","c2");
    cur->GetXaxis()->SetTitle("Muon weight");
    cur->GetYaxis()->SetTitle("Chi2/Ndf");
    cur->Draw();
    TLine *minline = new TLine(minw, minfom, maxw, minfom);
    minline->SetLineColor(kRed);
    minline->Draw("same");
    TLine *errline = new TLine(minw, minfom+1./(boundbin+1), maxw, minfom+1./(boundbin+1));
    errline->SetLineColor(kRed);
    errline->SetLineStyle(2);
    errline->Draw("same");
    c2->Print("output/muon_reweight_curve.pdf");


}
