TGraphErrors * ConvertTH1DtoGraph(TH1D * hist, const TString tag){
    int n = hist->GetNbinsX(); 

    double x[n], y[n], ex[n], ey[n]; // define arrays to store the x and y values, as well as the error bars

    for (int i=1; i<=n; i++) {
        x[i-1] = hist->GetBinCenter(i); // get the bin center for x
        y[i-1] = hist->GetBinContent(i); // get the bin content for y
        if(tag.Contains("data")) ex[i-1] = hist->GetBinWidth(i)*0.0; // set the error bar for x to half the bin width
        else ex[i-1] = hist->GetBinWidth(i)/2.0; // set the error bar for x to half the bin width
        //if(tag.Contains("data")) ey[i-1] = sqrt(y[i-1]); // set the error bar for y to the square root of the bin content hist->GetBinError(i);//
        //else ey[i-1] = 
        ey[i-1] = hist->GetBinError(i);//
        //ey[i-1] = sqrt(y[i-1]); 
    }

    TGraphErrors *graph = new TGraphErrors(n,x,y,ex,ey);
    

    return graph;
}
void RemovePointsWithZero(TGraphErrors *graph, const TString tag){
    Int_t npoints = graph->GetN();

    for (Int_t i=0; i<npoints; i++) {
        Double_t x, y;
        graph->GetPoint(i, x, y);
        if (y == 0) {
            graph->RemovePoint(i);
            npoints--;
            i--;
        }
    }
}

void DrawSystError(TGraphErrors *graph, const TString tag){

    double sys_err[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 35.574,25.271,16.661,17.593,15.142,13.621,14.135,14.484,21.067,18.678};
    double sys_err_Pi0KE[] = {0.0483,0.0243,0.0322,0.0222,0.0150,0.0123,0.0107,0.0022};
    double sys_err_Pi0Theta[] = {0.0834,0.1510,0.1499,0.1339,0.0806,0.0874,0.0600,0.0235,0.0606,0};
    double sys_err_Pi0CosTheta[] = {3.387,1.846,2.605,3.773,6.091,3.878,6.375,9.912,12.07,21.95};

    for (int i = 0; i < graph->GetN(); i++) {
        double stat_err = graph->GetErrorY(i);
        double sys_err_i = sys_err[i];
        if(tag.Contains("Pi0KE")) sys_err_i = sys_err_Pi0KE[i];
        if(tag.Contains("Pi0Theta")) sys_err_i = sys_err_Pi0Theta[i];
        if(tag.Contains("Pi0CosTheta")) sys_err_i = sys_err_Pi0CosTheta[i];

        double total_err = std::sqrt(stat_err * stat_err + sys_err_i * sys_err_i);
        if(tag.Contains("Pi0CosTheta")) {
            //cout << "i: " << i << " stat_err: " << stat_err << endl;
            //cout << "i: " << i << " sys_err_i: " << sys_err_i << endl;
            cout << "i: " << i << " total_err: " << total_err << endl;
        }
        graph->SetPointError(i, 0, total_err);
        
    }
    if(tag.Contains("Pi0Theta")) cout << graph->RemovePoint(graph->GetN()-1);
    //cout <<"last: " << graph->GetPointY(graph->GetN()-1) << endl; //cout << graph->RemovePoint(graph->GetN());


    graph->Draw("P sames");

}

void DrawHist(TGraphErrors* graph_dataTotal, TGraphErrors* graph_recoTotal, TH1D * h1d_recoTotal, TGraph* g_cex, const TString tag){

    TLatex tt;
    tt.SetNDC();
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(3);
    auto tmpGraph = (TGraphErrors*)graph_dataTotal->Clone();
    tmpGraph->SetTitle("");
    tmpGraph->SetMarkerSize(0);
    if(tag.Contains("Total")) tmpGraph->SetMaximum(300);
    else if(tag.Contains("Pi0Theta")) tmpGraph->SetMaximum(2);
    else if(tag.Contains("Pi0CosTheta")) tmpGraph->SetMaximum(250);
    else tmpGraph->SetMaximum(0.4);
    if(tag.Contains("Total")) tmpGraph->GetXaxis()->SetRangeUser(450.0,950.0);
    else if(tag.Contains("Pi0Theta")) tmpGraph->GetXaxis()->SetRangeUser(0,180);
    else if(tag.Contains("Pi0CosTheta")) tmpGraph->GetXaxis()->SetRangeUser(-1,1);
    else tmpGraph->GetXaxis()->SetRangeUser(0,800);

    tmpGraph->GetXaxis()->SetTitle("T_{#pi^{+}} (MeV)");
    tmpGraph->GetYaxis()->SetTitle("#sigma_{CEX} (mb)");
    if(tag.Contains("Pi0KE")) {
        tmpGraph->GetXaxis()->SetTitle("T_{#pi^{0}} (MeV)");
        tmpGraph->GetYaxis()->SetTitle("d#sigma/dT_{#pi^{0}} (mb/MeV)");
    }
    else if(tag.Contains("Pi0Theta")) {
        tmpGraph->GetXaxis()->SetTitle("#theta_{#pi^{0}} (deg.)");
        tmpGraph->GetYaxis()->SetTitle("d#sigma/d#theta_{#pi^{0}} (mb/deg.)");
    }
    else if(tag.Contains("Pi0CosTheta")) {
        tmpGraph->GetXaxis()->SetTitle("cos#theta_{#pi^{0}}");
        tmpGraph->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi^{0}} (mb)");
    }

    tmpGraph->GetYaxis()->CenterTitle();
    tmpGraph->GetYaxis()->SetTitleFont(22);
    tmpGraph->GetYaxis()->SetTitleSize(0.05);
    tmpGraph->GetYaxis()->SetTitleOffset(0.92);
    tmpGraph->GetXaxis()->CenterTitle();
    tmpGraph->GetXaxis()->SetTitleFont(22);
    tmpGraph->GetXaxis()->SetTitleSize(0.05);
    tmpGraph->GetXaxis()->SetTitleOffset(0.9);
    tmpGraph->Draw("AP sames");


    RemovePointsWithZero(graph_dataTotal,tag);
    graph_dataTotal->SetMarkerColor(kBlack);
    graph_dataTotal->SetMarkerStyle(8);
    tmpGraph->SetLineWidth(2);
    graph_dataTotal->SetLineWidth(2);
    //graph_dataTotal->Draw("P same");
    //if(tag.Contains("Total")) DrawSystError(tmpGraph,"Total");

    RemovePointsWithZero(graph_recoTotal,tag);
    graph_recoTotal->SetLineWidth(2);
    graph_recoTotal->SetLineColor(kGray+2);
    graph_recoTotal->SetMarkerSize(0);
    graph_recoTotal->SetFillColorAlpha(kGray, 0.4);
    graph_recoTotal->Draw("e2 sames");
    /*auto graph_recoTotal_overlay = (TGraphErrors*)graph_recoTotal->Clone("graph_recoTotal");
    graph_recoTotal_overlay->SetMarkerSize(0);
    graph_recoTotal_overlay->SetLineWidth(3);
    graph_recoTotal_overlay->SetLineColor(kGray+2);
    for (int i = 0; i < graph_recoTotal_overlay->GetN(); i++) {
        graph_recoTotal_overlay->SetPointError(i, 25, 0);
    }   
    graph_recoTotal_overlay->Draw("Z sames");
    */
    h1d_recoTotal->SetFillStyle(0);
    h1d_recoTotal->SetLineColor(kGray+2);
    h1d_recoTotal->SetLineWidth(3);
    h1d_recoTotal->Draw("HIST sames");
    
    g_cex->SetLineColor(kRed);
    g_cex->SetLineWidth(2);
    g_cex->Draw("sames C");

    graph_dataTotal->Draw("P same");
    DrawSystError(tmpGraph,tag);


    TLegend * legend = new TLegend(0.15, 0.58, 0.48, 0.88);
    if(tag.Contains("Pi0KE") || tag.Contains("Pi0Theta")) legend = new TLegend(0.55, 0.58, 0.88, 0.88);

    legend->AddEntry(g_cex, "Geant4 v4.10.6", "l");
    TLegendEntry *entry_reco = legend->AddEntry(graph_recoTotal, "Simulation (stats.)", "lep");
    TLegendEntry *entry_data = legend->AddEntry(graph_dataTotal, "Data", "lep");

    entry_data->SetOption("pl");
    entry_data->SetMarkerSize(1.2);
    entry_data->SetMarkerColor(graph_dataTotal->GetMarkerColor());
    entry_data->SetLineColor(graph_dataTotal->GetLineColor());
    entry_data->SetLineWidth(graph_dataTotal->GetLineWidth());
    entry_data->SetFillColor(graph_dataTotal->GetFillColor());
    entry_data->SetFillStyle(graph_dataTotal->GetFillStyle());

    entry_reco->SetOption("lf");
    entry_reco->SetFillColorAlpha(kGray, 0.4);

    //TString lheader("1GeV Pion Data");
    //legend->SetHeader(lheader);
    legend->SetBorderSize(0);
    legend->Draw("same");

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    
    if(tag.Contains("Total")) tt.DrawLatex(0.725,0.925,"1GeV/c Pion Data");
    else tt.DrawLatex(0.600,0.925,"Beam T_{#pi^{+}} = 650 - 800 MeV Data");
    if(tag.Contains("Total") || tag.Contains("Pi0CosTheta")){
        tt.DrawLatex(0.635,0.825,"#bf{#it{Preliminary}}");
        tt.DrawLatex(0.635,0.775,"#bf{#it{(Stats. + Syst. Error)}}");
        tt.DrawLatex(0.635,0.725,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");

    }
    else{
        tt.DrawLatex(0.185,0.825,"#bf{#it{Preliminary}}");
        tt.DrawLatex(0.185,0.775,"#bf{#it{(Stats. + Syst. Error)}}");
        tt.DrawLatex(0.185,0.725,"#color[4]{#pi^{+} #bf{+ Ar} #rightarrow #pi^{0} #bf{+ (nucleons)}}");

    }

}

void macro_XSplots(){

    TString finName = "output/outXShists.root";

    TFile *file = TFile::Open(finName);

    if(!file->IsOpen()){
        cout << "file file not open" << endl;
        exit(1);
    }

    // ======= GEANT4 Input XS ======== //
    const TString xcesfinName = "../Analysis/input/LAr_PiPlus_exclusive_cross_section.root";
    TFile *f_CrossSection = TFile::Open(xcesfinName);
    if(!f_CrossSection->IsOpen()){
        cout << "xsce file not open" << endl;
        exit(1);
    }

    TGraph *g_cex;

    f_CrossSection->GetObject("cex_KE",g_cex);

    const TString DiffxcesfinName = "../Analysis/input/cross_section_out_withNewVar.root";
    TFile *f_DiffCrossSection = TFile::Open(DiffxcesfinName);
    if(!f_DiffCrossSection->IsOpen()) {
        cout << "Diffxsce file not open" << endl;
        exit(1);
    }

   const TString finNameXS = "input/cross_section_out.root";
    TFile *fileXS = TFile::Open(finNameXS);

    TH1D * g_cex_650to800MeV = (TH1D*)fileXS->Get("inel_cex_1dKEpi0650to800_MeV");

    const Int_t x0_diff650to800 = g_cex_650to800MeV->GetXaxis()->GetFirst();
    const Int_t x1_diff650to800 = g_cex_650to800MeV->GetXaxis()->GetLast();
    double x_650to800[x1_diff650to800], y_650to800[x1_diff650to800];
    for(Int_t ix=x0_diff650to800; ix<=x1_diff650to800; ix++){
        //cout << "Bin 650to800: "  << ix << " " << g_cex_650to800MeV->GetBinContent(ix) << endl;
        x_650to800[ix] = g_cex_650to800MeV->GetBinCenter(ix);
        y_650to800[ix] = g_cex_650to800MeV->GetBinContent(ix);
    }

    TGraph *g_650to800 = new TGraph(40,x_650to800,y_650to800);



    // Get Histograms
    TH1D * h1d_dataTotal = (TH1D*) file->Get("data/data_TotalCEXxsec;1");
    TH1D * h1d_recoTotal = (TH1D*) file->Get("data/reco_TotalCEXxsec;1");
    //TH1D * h1d_truthTotal = (TH1D*) file->Get("data/truth_TotalCEXxsec;1");

    auto graph_dataTotal = ConvertTH1DtoGraph(h1d_dataTotal,"data");
    auto graph_recoTotal = ConvertTH1DtoGraph(h1d_recoTotal,"reco");

    TCanvas * cTotal = new TCanvas("cTotal", "", 1200, 800);
    DrawHist(graph_dataTotal,graph_recoTotal,h1d_recoTotal,g_cex,"Total");
    cTotal->Print("output/1GeV_Pion_TotalCEX_XS.pdf");


    // Get Histograms
    TH1D * h1d_dataPi0KE = (TH1D*) file->Get("data/data_DiffCEXxsec_Pi0KE;1");
    TH1D * h1d_recoPi0KE = (TH1D*) file->Get("data/reco_DiffCEXxsec_Pi0KE;1");
    //TH1D * h1d_truthPi0KE = (TH1D*) file->Get("data/truth_Pi0KECEXxsec;1");

    auto graph_dataPi0KE = ConvertTH1DtoGraph(h1d_dataPi0KE,"data");
    auto graph_recoPi0KE = ConvertTH1DtoGraph(h1d_recoPi0KE,"reco");

    TCanvas * cPi0KE = new TCanvas("cPi0KE", "", 1200, 800);
    DrawHist(graph_dataPi0KE,graph_recoPi0KE,h1d_recoPi0KE,g_650to800,"Pi0KE");
    //graph_dataPi0KE->Draw("AP");
    cPi0KE->Print("output/Pi0KE_DifferentialCEX_XS.pdf");


    TH1D * g_cex_650to800MeV_Theta = (TH1D*)fileXS->Get("inel_cex_1dThetapi0650to800_MeV");
    const Int_t x0_diff650to800_Theta = g_cex_650to800MeV_Theta->GetXaxis()->GetFirst();
    const Int_t x1_diff650to800_Theta = g_cex_650to800MeV_Theta->GetXaxis()->GetLast();
    double x_650to800_Theta[x1_diff650to800_Theta], y_650to800_Theta[x1_diff650to800_Theta];
    for(Int_t ix=x0_diff650to800_Theta; ix<=x1_diff650to800_Theta; ix++){
        //cout << "Bin 650to800: "  << ix << " " << g_cex_650to800MeV_Theta->GetBinContent(ix) << endl;
        x_650to800_Theta[ix] = g_cex_650to800MeV_Theta->GetBinCenter(ix);
        y_650to800_Theta[ix] = g_cex_650to800MeV_Theta->GetBinContent(ix);
    }

    TGraph *g_650to800_Theta = new TGraph(40,x_650to800_Theta,y_650to800_Theta);
  
    // Get Histograms
    TH1D * h1d_dataPi0Theta = (TH1D*) file->Get("data/data_DiffCEXxsec_Pi0Theta;1");
    TH1D * h1d_recoPi0Theta = (TH1D*) file->Get("data/reco_DiffCEXxsec_Pi0Theta;1");
    //TH1D * h1d_truthPi0Theta = (TH1D*) file->Get("data/truth_Pi0ThetaCEXxsec;1");

    auto graph_dataPi0Theta = ConvertTH1DtoGraph(h1d_dataPi0Theta,"data");
    auto graph_recoPi0Theta = ConvertTH1DtoGraph(h1d_recoPi0Theta,"reco");

    TCanvas * cPi0Theta = new TCanvas("cPi0Theta", "", 1200, 800);
    DrawHist(graph_dataPi0Theta,graph_recoPi0Theta,h1d_recoPi0Theta,g_650to800_Theta,"Pi0Theta");
    //graph_dataPi0Theta->Draw("AP");
    cPi0Theta->Print("output/Pi0Theta_DifferentialCEX_XS.pdf");

    TH1D * g_cex_650to800MeV_CosTheta = (TH1D*)fileXS->Get("inel_cex_1dcosThetapi0650to800_MeV");
    const Int_t x0_diff650to800_CosTheta = g_cex_650to800MeV_CosTheta->GetXaxis()->GetFirst();
    const Int_t x1_diff650to800_CosTheta = g_cex_650to800MeV_CosTheta->GetXaxis()->GetLast();
    double x_650to800_CosTheta[x1_diff650to800_CosTheta], y_650to800_CosTheta[x1_diff650to800_CosTheta];
    for(Int_t ix=x0_diff650to800_CosTheta; ix<=x1_diff650to800_CosTheta; ix++){
        //cout << "Bin 650to800: "  << ix << " " << g_cex_650to800MeV_CosTheta->GetBinContent(ix) << endl;
        x_650to800_CosTheta[ix] = g_cex_650to800MeV_CosTheta->GetBinCenter(ix);
        y_650to800_CosTheta[ix] = g_cex_650to800MeV_CosTheta->GetBinContent(ix);
    }

    TGraph *g_650to800_CosTheta = new TGraph(40,x_650to800_CosTheta,y_650to800_CosTheta);
    g_650to800_CosTheta->RemovePoint(0);
    g_650to800_CosTheta->RemovePoint(1);

    // Get Histograms
    TH1D * h1d_dataPi0CosTheta = (TH1D*) file->Get("data/data_DiffCEXxsec_Pi0CosTheta;1");
    TH1D * h1d_recoPi0CosTheta = (TH1D*) file->Get("data/reco_DiffCEXxsec_Pi0CosTheta;1");
    //TH1D * h1d_truthPi0CosTheta = (TH1D*) file->Get("data/truth_Pi0CosThetaCEXxsec;1");

    auto graph_dataPi0CosTheta = ConvertTH1DtoGraph(h1d_dataPi0CosTheta,"data");
    auto graph_recoPi0CosTheta = ConvertTH1DtoGraph(h1d_recoPi0CosTheta,"reco");

    TCanvas * cPi0CosTheta = new TCanvas("cPi0CosTheta", "", 1200, 800);
    DrawHist(graph_dataPi0CosTheta,graph_recoPi0CosTheta,h1d_recoPi0CosTheta,g_650to800_CosTheta,"Pi0CosTheta");
    cPi0CosTheta->Print("output/Pi0CosTheta_DifferentialCEX_XS.pdf");



    cout << "done" << endl;


}

