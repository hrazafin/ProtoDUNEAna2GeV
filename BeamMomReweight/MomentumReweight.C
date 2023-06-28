
double CalWeight(const bool kMC, const double & trueBeamMom){
    
  double weight = 1.;
  //return weight;
  
  //double mom_mu0 = 1.0033;
  //double mom_sigma0 = 0.0609;
  double mom_mu0 = 1.00362;
  double mom_sigma0 = 0.0595377;
  //double mom_mu = 1.01656;//1.01642;//1.01818;
  //double mom_sigma = 0.0711996;//0.07192;

  double mom_mu = 1.01604;//1.01642;//1.01818;
  double mom_sigma = 0.0708354;//0.07192;

  double wlimit = 3.5;

  // 10 - 13.6069
  // 9 - 13.5196
  // 8 - 13.436
  // 7 - 13.4201
  // 6 - 13.5564
  // 5 - 13.9625
  // 4 - 14.7017
  // 3 - 16.7721
  if(kMC){
    // momentum reweight (outlier weights are set to 1e-5)
    double deno = exp(-pow((trueBeamMom-mom_mu0)/mom_sigma0,2)/2);
    //if (deno < wlimit) deno = wlimit;
    double numo = exp(-pow((trueBeamMom-mom_mu)/mom_sigma,2)/2);
    //if (numo < wlimit) numo = wlimit;
    weight *= numo;
    weight /= deno;
    if (weight>wlimit) weight=wlimit;
    if (weight<1./wlimit) weight=1./wlimit;
  }

  return weight;
}

double GetChi2(TH1D * hdata, TH1D * hmc){
    int nx = hdata->GetNbinsX();
    double chisq = 0;
    for(int ix=1; ix<=nx; ix++){
      double x = hdata->GetBinContent(ix);
      double y = hmc->GetBinContent(ix);
        
      double ex = hdata->GetBinError(ix);
      double ey = hmc->GetBinError(ix);

      if (!ex) ex = 1;
      chisq += pow(x-y,2)/(pow(ex,2)+pow(ey,2));
    }
    return chisq;
}

// Define a function that takes as input the parameters that you want to vary in the MC distribution
void MCfunction(int &npars, double *grad, double &value, double *par, int flag){
    
    // Generate a new MC histogram with the updated parameters
    TH1D *mcHist = new TH1D("mcHist", "New MC Histogram", 20, 600, 1300);
    double mom_mu0 = 1.00362;
    double mom_sigma0 = 0.0595377;
    double wlimit = 3.5; 

    //TString finName = "input_file/outana_treeBin1.root";
    TString finName = "output/outana.root";

    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
        cout << "file file not open" << endl;
        exit(1);
    }

    TH1D * dataHist = (TH1D*) file->Get("data/a004hBeamStopMuonKEff_STK_sum;1");

    TTree * tree=(TTree*) file->Get("mc/tree");
    double ientry = 0;

    double trueBeamMom, recoFFKE;
    tree->SetBranchAddress("trueBeamMom", &trueBeamMom);
    tree->SetBranchAddress("recoFFKE", &recoFFKE);
    cout << "par[0]: " <<par[0] << "par[1]: " << par[1] << endl;
    // Loop over TTree
    while(tree->GetEntry(ientry)){
        ientry++;
        
        double deno = exp(-pow((trueBeamMom-mom_mu0)/mom_sigma0,2)/2);
        double numo = exp(-pow((trueBeamMom-par[0])/par[1],2)/2);

        double weight = numo/deno;// calculate the weight for this event using the ratio of Gaussian functions
        if (weight>wlimit) weight=wlimit;
        if (weight<1./wlimit) weight=1./wlimit;
        
        mcHist->Fill(recoFFKE, weight);

    }

    double mcIntegral = mcHist->Integral();
    double dataIntegral = dataHist->Integral();

    mcHist->Scale(1/mcIntegral);
    dataHist->Scale(1/dataIntegral);
/*
    double chi2 = 0;
    int ndf = 0;

    for (int i = 1; i <= mcHist->GetNbinsX(); i++) {
        double dataContent = dataHist->GetBinContent(i);
        double dataError = dataHist->GetBinError(i);
        double mcContent = mcHist->GetBinContent(i);
        double mcError = mcHist->GetBinError(i);
        
        //cout << "dataError: " << dataError << endl;
        //cout << "mcError: " << mcError << endl;

        if (dataContent > 0 && mcContent > 0) {
            double temp = (dataContent - mcContent);// / (dataError*dataError);// + mcError*mcError);
            chi2 += temp*temp;
            ndf++;
        }
    }

    cout << "chi2 : " << chi2 << endl;
    //cout << "ndf : " << ndf << endl;

    //cout << "chi2 / ndf : " << chi2 / ndf << endl;
*/


    // Calculate the chi-square between the data and new MC histograms
    double chi2 = dataHist->Chi2Test(mcHist, "WWCHI2");
    
    cout << "chi2 : " << chi2 << endl;

    //delete mcHist;

    value = chi2;

}


void MomentumReweight(){

    double plotScale = 1.63792;
    gStyle->SetOptStat(0);

    //TString finName = "input_file/outana_Total.root";
    TString finName = "output/outana.root";

    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
        cout << "file file not open" << endl;
        exit(1);
    }

    // Get Histograms
    TH1D * h1d_trueMom = (TH1D*) file->Get("mc/b003hTrueBeamMom_matched;1");
    
    //TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
    h1d_trueMom->GetXaxis()->SetRangeUser(800,1300);

    TF1 *f1 = new TF1("f1","gaus",800,1500);
    h1d_trueMom->Fit("f1");
    h1d_trueMom->SetLineColor(kBlue);
    h1d_trueMom->SetLineWidth(2);
    h1d_trueMom->Draw("E");
    f1->Draw("sames");

    cout << "true pars: " << "1. Mean: " << f1->GetParameter(1) << endl;
    cout << "true pars: " << "1. Sigma: " << f1->GetParameter(2) << endl;

    //c1->Print("output/hStopMuontrueMom.pdf");

    TH2D * h2d_ffE = (TH2D*) file->Get("mc/a004hBeamStopMuonKEff_STK;1");
    TH1D * h1d_ffE_reco = h2d_ffE->ProjectionX("hStopMuonrecffE_projX");
    TH1D * h1d_ffE_data = (TH1D*) file->Get("data/a004hBeamStopMuonKEff_STK_sum;1");

    TH2D * h2d_instP = (TH2D*) file->Get("mc/a003hBeamStopMuonInstP_STK;1");
    TH1D * h1d_instP_reco = h2d_instP->ProjectionX("hStopMuonrecinstP_projX");
    TH1D * h1d_instP_data = (TH1D*) file->Get("data/a003hBeamStopMuonInstP_STK_sum;1");

    cout << "done" << endl;



    TTree * tree=(TTree*) file->Get("mc/tree");

    double ientry = 0;
    TH1D * hmuonffKE_recoTree = new TH1D("z001hmuonffKE_recoTree",";Reco FF KE (MeV); Candidates", 20, 600, 1300);


    //TH1D * hmuonffKE_dataTree = new TH1D("z002hmuonffKE_dataTree",";Data FF KE (MeV); Candidates", 30, 700, 1300);
    TH1D * hmuonstartP_trueTree = new TH1D("z003hmuonstartP_trueTree",";True Start P (GeV); Candidates", 20, 0.7, 1.3);
    TH1D * hmuonstartP_trueTree_wt = new TH1D("z004hmuonstartP_trueTree_wt",";True Start P (GeV); Candidates", 20, 0.7, 1.3);

    TH1D * hmuon_instP_recoTree = new TH1D("z002hmuon_instP_recoTree",";Reco Inst P (MeV); Candidates", 20, 0.6, 1.3);


    double trueBeamMom;
    tree->SetBranchAddress("trueBeamMom", &trueBeamMom);

    double recoFFKE;
    tree->SetBranchAddress("recoFFKE", &recoFFKE);

    double instP;
    tree->SetBranchAddress("instP", &instP);
    
    // Loop over TTree
    while(tree->GetEntry(ientry)){
        ientry++;
        
        double weight = CalWeight(true,trueBeamMom);

        hmuonstartP_trueTree->Fill(trueBeamMom);
        hmuonstartP_trueTree_wt->Fill(trueBeamMom,weight);

        hmuonffKE_recoTree->Fill(recoFFKE,weight);
        hmuon_instP_recoTree->Fill(instP,weight);
        
    }

    TCanvas * c11 = new TCanvas("c11", "", 1200, 800);
    //hmuonstartP_trueTree->GetXaxis()->SetRangeUser(800,1300);

    TF1 *f11 = new TF1("f11","gaus",800,1500);
    hmuonstartP_trueTree->Fit("f11");
    hmuonstartP_trueTree->SetLineColor(kBlue);
    hmuonstartP_trueTree->SetLineWidth(2);
    hmuonstartP_trueTree->Draw("E");

    hmuonstartP_trueTree_wt->SetLineColor(kRed);
    hmuonstartP_trueTree_wt->SetLineWidth(2);
    hmuonstartP_trueTree_wt->Draw("hists sames");

    f11->Draw("sames");

    cout << "1 true pars: " << "1. Mean: " << f1->GetParameter(1) << endl;
    cout << "1 true pars: " << "1. Sigma: " << f1->GetParameter(2) << endl;

    c11->Print("output/hStopMuontrueMom.pdf");

    //TCanvas * c22 = new TCanvas("c22", "", 1200, 800);
    TF1 *f12_reco = new TF1("f12_reco","gaus",800,1500);

    hmuonffKE_recoTree->SetLineColor(kBlue);
    hmuonffKE_recoTree->SetLineWidth(2);
    hmuonffKE_recoTree->SetMaximum(hmuonffKE_recoTree->GetMaximum()*1.5);
    hmuonffKE_recoTree->SetFillStyle(0);
    hmuonffKE_recoTree->Scale(1/hmuonffKE_recoTree->Integral());
    hmuonffKE_recoTree->Fit("f12_reco","0"); // Only fit not draw

    hmuonffKE_recoTree->Draw("hists");

    //c22->Print("output/hmuonffKE_recoTree.pdf");

    TCanvas * c2 = new TCanvas("c2", "", 1200, 800);
    TF1 *f1_reco = new TF1("f1_reco","gaus",800,1500);

    h1d_ffE_reco->SetTitle(" ");
    h1d_ffE_reco->GetYaxis()->CenterTitle();
    h1d_ffE_reco->GetYaxis()->SetTitleFont(22);
    h1d_ffE_reco->GetYaxis()->SetTitleSize(0.05);
    h1d_ffE_reco->GetYaxis()->SetTitleOffset(0.9);
    h1d_ffE_reco->GetXaxis()->CenterTitle();
    h1d_ffE_reco->GetXaxis()->SetTitleFont(22);
    h1d_ffE_reco->GetXaxis()->SetTitleSize(0.05);
    h1d_ffE_reco->GetXaxis()->SetTitleOffset(0.9);

    h1d_ffE_reco->GetYaxis()->SetTitle("Area Normalised");

    h1d_ffE_reco->Draw();

    h1d_ffE_reco->SetLineColor(kGreen+3);
    h1d_ffE_reco->SetLineStyle(7);
    h1d_ffE_reco->SetLineWidth(2);

    //h1d_ffE_reco->SetLineColor(kBlue);
    //h1d_ffE_reco->SetLineWidth(2);
    h1d_ffE_reco->SetMaximum(h1d_ffE_reco->GetMaximum()*1.5);
    h1d_ffE_reco->SetFillStyle(0);
    //h1d_ffE_reco->Scale(plotScale);
    h1d_ffE_reco->Scale(1/h1d_ffE_reco->Integral());
    h1d_ffE_reco->Fit("f1_reco","0"); // Only fit not draw

    h1d_ffE_reco->Draw("hists");

    cout << "reco pars: " << "1. Mean: " << f1_reco->GetParameter(1) << endl;
    cout << "reco pars: " << "1. Sigma: " << f1_reco->GetParameter(2) << endl;

    h1d_ffE_data->SetLineColor(kBlack);
    h1d_ffE_data->SetFillStyle(0);
    h1d_ffE_data->SetMarkerStyle(8);
    h1d_ffE_data->SetMarkerSize(1);
    h1d_ffE_data->SetMarkerColor(kBlack);
    //h1d_ffE_data->SetLineWidth(2);
    h1d_ffE_data->Scale(1/h1d_ffE_data->Integral());
    TF1 *f1_data = new TF1("f1_data","gaus",800,1500);
    h1d_ffE_data->Fit("f1_data","0");

    h1d_ffE_data->Draw("E1 sames");

    cout << "data pars: " << "1. Mean: " << f1_data->GetParameter(1) << endl;
    cout << "data pars: " << "1. Sigma: " << f1_data->GetParameter(2) << endl;
    
    hmuonffKE_recoTree->SetLineColor(kRed);
    hmuonffKE_recoTree->SetLineWidth(2);
    hmuonffKE_recoTree->SetMaximum(hmuonffKE_recoTree->GetMaximum()*1.5);
    hmuonffKE_recoTree->SetFillStyle(0);
    hmuonffKE_recoTree->Fit("f12_reco","0"); // Only fit not draw

    hmuonffKE_recoTree->Draw("hists sames");
/*
    auto lg = new TLegend(0.68,0.68,0.88,0.88);
    lg->AddEntry(h1d_ffE_reco,"MC","lp");
    lg->AddEntry(hmuonffKE_recoTree,"MC weighted","lp");
    lg->AddEntry(h1d_ffE_data,"data","lep");
    lg->SetBorderSize(0);
    lg->Draw("sames");
*/
    auto lgfit = new TLegend(0.58,0.60,0.88,0.88);
    //lgfit->SetHeader(Form("#chi^{2}: %f",GetChi2(h0,hmc_fit)));
    lgfit->SetHeader(Form("#chi^{2}_{def}: %.2f, #chi^{2}_{fit}: %.2f", GetChi2(h1d_ffE_data,h1d_ffE_reco), GetChi2(h1d_ffE_data,hmuonffKE_recoTree)));
    lgfit->AddEntry(h1d_ffE_data,"data","lp");
    lgfit->AddEntry(h1d_ffE_reco,"MC Default","l");
    lgfit->AddEntry(hmuonffKE_recoTree,"MC Reweighted","l");
    lgfit->SetBorderSize(0);
    lgfit->Draw("sames");

    TLatex tt;
    tt.SetNDC();

    tt.DrawLatex(0.135,0.825,"#color[4]{#bf{#it{#mu_{fit}} = 1.016}}");
    tt.DrawLatex(0.135,0.755,"#color[4]{#bf{#it{#sigma_{fit}} = 0.071}}");


    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    tt.DrawLatex(0.685,0.925,"Stopping Beam Muon");

    c2->Print("output/hStopMuonrecffE.pdf");


    TCanvas * c3 = new TCanvas("c3", "", 1200, 800);
    //double newplotScale = h1d_instP_data->Integral()/hmuon_instP_recoTree->Integral();
    cout << "old plotScale: " << h1d_instP_data->Integral()/h1d_instP_reco->Integral() << endl;
    cout << "new plotScale: " << h1d_instP_data->Integral()/hmuon_instP_recoTree->Integral() << endl;


 

    h1d_instP_reco->Draw();

    h1d_instP_reco->SetLineColor(kBlue);
    h1d_instP_reco->SetLineWidth(2);
    h1d_instP_reco->SetMaximum(h1d_instP_reco->GetMaximum()*1.5);
    h1d_instP_reco->SetFillStyle(0);
    //h1d_instP_reco->Scale(plotScale);
    h1d_instP_reco->Scale(1/h1d_instP_reco->Integral());

    h1d_instP_reco->Draw("hists");


    h1d_instP_data->SetLineColor(kBlack);
    h1d_instP_data->SetFillStyle(0);
    h1d_instP_data->SetMarkerStyle(8);
    h1d_instP_data->SetMarkerSize(1);
    h1d_instP_data->SetLineWidth(2);
    h1d_instP_data->Scale(1/h1d_instP_data->Integral());

    h1d_instP_data->Draw("E sames");

    hmuon_instP_recoTree->SetLineColor(kRed);
    hmuon_instP_recoTree->SetLineWidth(2);
    //hmuon_instP_recoTree->SetMaximum(hmuon_instP_recoTree->GetMaximum()*1.5);
    //hmuon_instP_recoTree->Scale(newplotScale);
    hmuon_instP_recoTree->Scale(1/hmuon_instP_recoTree->Integral());

    hmuon_instP_recoTree->SetFillStyle(0);
    hmuon_instP_recoTree->Draw("hists sames");

    auto lg3 = new TLegend(0.68,0.68,0.88,0.88);
    lg3->AddEntry(h1d_instP_reco,"MC","lp");
    lg3->AddEntry(hmuon_instP_recoTree,"MC weighted","lp");
    lg3->AddEntry(h1d_instP_data,"data","lep");
    lg3->SetBorderSize(0);
    lg3->Draw("sames");

    //c3->Print("output/hStopMuonrecInstP.pdf");


    // Minimize the chi-square by varying the parameters of the MC distribution using a minimization algorithm
    TMinuit *fitter = new TMinuit(2);
    cout << "MCfunction" << endl;
    fitter->SetFCN(MCfunction);
    fitter->DefineParameter(0, "mean", 1.00000, 0.01, 0.9, 1.1);
    fitter->DefineParameter(1, "sigma", 0.07000, 0.01, 0.06, 0.08);

    int flag = fitter->Command("MIGRAD");
    cout << "flag: " << flag << endl;
    double mean_fit, errmean_fit;
    fitter->GetParameter(0,mean_fit,errmean_fit);
    cout << "mean_fit: " << mean_fit << endl;

    double sigma_fit, errsigma_fit;
    fitter->GetParameter(1,sigma_fit,errsigma_fit);
    cout << "sigma_fit: " << sigma_fit << endl;

    Double_t chi2 = h1d_ffE_data->Chi2Test(hmuonffKE_recoTree, "WWCHI2");
    cout << "Final Chi-Square: " << chi2 << endl;

    double chi2_def = GetChi2(h1d_ffE_data,h1d_ffE_reco);
    cout << "chi2_def: " << chi2_def << endl;
    double chi2_fit = GetChi2(h1d_ffE_data,hmuonffKE_recoTree);
    cout << "chi2_fit: " << chi2_fit << endl;

    file->Close();

}

