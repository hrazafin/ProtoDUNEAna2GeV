#include "../include/PlotUtils.h"

#ifndef EPSILON
#define EPSILON 1E-12
#endif

double PlotUtils::PrintStat(const TString tag, TH1 *hh, const double val0, const double val1, const double oldsel)
{
  const double newall = hh->Integral(0,100000);
  if(oldsel!=-999){
    if( fabs(newall-oldsel)>EPSILON ){
      printf("style::PrintStat newall != oldsel %f %f\n", newall, oldsel); exit(1);
    }
  }

  double nsel = -999;
  TH2 * h2d = dynamic_cast<TH2*>(hh);
  if(h2d){
    TAxis *ax = h2d->GetXaxis();
    const int xbin0 = ax->FindBin(val0);
    const int xbin1 = ax->FindBin(val1);
    nsel = h2d->Integral(xbin0, xbin1, 0, 100000);
  }
  else{
    const int xbin0 = hh->FindBin(val0);
    const int xbin1 = hh->FindBin(val1);
    nsel = hh->Integral(xbin0, xbin1);
  }

  printf("%-50s: all %5.1f selected %5.1f fraction %.1f%% histogram %s \n", tag.Data(), newall, nsel, nsel/newall*100, hh->GetName());
  return nsel;
}


void PlotUtils::FillHist(TH1 * hh,  double xx, const double yy)
{
  // Get histogram name
  const TString tag = hh->GetName();
  // Get X bin size
  const int nbx = hh->GetNbinsX();
  // Get min and max value
  const double xmin = hh->GetXaxis()->GetBinLowEdge(1);
  const double xmax = hh->GetXaxis()->GetBinUpEdge(nbx);  
  
  if(xx<xmin){
    xx = hh->GetXaxis()->GetBinCenter(1);
  }
  if(xx>=xmax){
    xx = hh->GetXaxis()->GetBinCenter(nbx);
  }

  TH2 * h2 = dynamic_cast<TH2*>(hh);
  if(h2){
    hh->Fill(xx, yy);
  }
  else{
    hh->Fill(xx);
  }

  double tmpoverflow = -999;
  if(h2){
    tmpoverflow = h2->Integral(nbx+1, nbx+1, 0, 10000);
  }
  else{
    tmpoverflow = hh->GetBinContent(nbx+1);
  }
  if(tmpoverflow>1E-12){
    printf("PlotUtils::FillInRange still overflow! %s %f %d %f %f\n", tag.Data(), tmpoverflow, nbx, xx, yy); exit(1);
  }
}

void PlotUtils::ProcessHist(TList *lout, const bool kMC)
{
  // Get the size of TList
  const int nhist = lout->GetSize();
  // Loop over each histogram
  for(int ii=0; ii<nhist; ii++){
    TH2D * htmp = dynamic_cast<TH2D *>( lout->At(ii) );
    if(htmp){
      const TString tag = htmp->GetName();
      // Check if the histogram name contians stack tags
      if(tag.Contains("STK") || tag.Contains("OVERLAY") || tag.Contains("COMPOSE")){
        // Convert this tmp 2D histogram to a stack histogram
        THStack * stk = ConvertToStack(htmp,kMC,typeMaps);
        // Check if the stk exist 
        if(stk){
          // Add to lout for MC
          if(kMC) {
            lout->Add(stk);
            // Normalize the stk histogram
            THStack * snor = NormalizeStack(stk);
            lout->Add(snor);
          } 
          // Only need the sum of this stack histogram in data 
          else{
            // Get the sum of this hitogram
            TH1D * hsum = GetStackedSum(stk);
            lout->Add(hsum);
          }
        }
      } // End of stack tag

      // Check if the histogram name contians RES (resolution)
      else if(tag.Contains("RES") && kMC){
        // Column normalise each bin to easily see the maximum
        TH1D * hmean = 0x0;
        TH1D * hcdf = 0x0;
        TH2D * hnor = NormalHist(htmp, 5, true, hmean, hcdf);
        hnor->SetTitle(tag);
        lout->Add(hnor);
        // Get the 1D histogram resolution plot using projectY
        TH1D * hprojY = htmp->ProjectionY(tag+"_projY");
	      hprojY->SetStats(1);
        lout->Add(hprojY);
        TH1D * hprojX = htmp->ProjectionX(tag+"_projX");
        hprojX->SetStats(1);
        lout->Add(hprojX);
      }
      // Check if the histogram name contians REG (regular)
      else if(tag.Contains("REG") && kMC){
        // Column normalise each bin to easily see the maximum
        TH1D * hmean = 0x0;
        TH1D * hcdf = 0x0;
        TH2D * hnor = NormalHist(htmp, 5, true, hmean, hcdf);
        hnor->SetTitle(tag);
        lout->Add(hnor);
        lout->Add(hmean);
        lout->Add(hcdf);
        TH1D * hprojX = htmp->ProjectionX(tag+"_projX");
        hprojX->SetStats(1);
        lout->Add(hprojX);
        // Do the corrections
        if (tag.Contains("Correction") && kMC){
          xSlicedEnergyCorrection(htmp);
        }
      }
      else if(tag.Contains("sigma")){   
        xSlicedSigma(htmp,tag);
      }
      // Do nothing (You can add more else if to process more tags)
      else {}
    } // End of if(htmp)

    // Stack histogram
    THStack * hstk = dynamic_cast<THStack *>(lout->At(ii));
    if(hstk){
      const TString tag = hstk->GetName();
      // TKI related variables
      if(tag.Contains("stkTruth") && kMC){
        // Get the name of raw shower histogram
        TRegexp re("stk");
        TString tmp = tag;
        tmp(re) = "h";
        // Get the histogram name
        TString name1p0n = tmp+"1p0n";
        TString nameNp0n = tmp+"Np0n";
        TString name1pMn = tmp+"1pMn";
        TString nameNpMn = tmp+"NpMn";

        // Find histograms in the list
        TH1D *hh1p0n = (TH1D*)lout->FindObject(name1p0n);
        TH1D *hhNp0n = (TH1D*)lout->FindObject(nameNp0n); 
        TH1D *hh1pMn = (TH1D*)lout->FindObject(name1pMn); 
        TH1D *hhNpMn = (TH1D*)lout->FindObject(nameNpMn); 
        // Get the total number of entries
        double ntotall = -999;
        double n1 = hh1p0n->Integral(0,hh1p0n->GetNbinsX()+1); 
        double n2 = hhNp0n->Integral(0,hhNp0n->GetNbinsX()+1); 
        double n3 = hh1pMn->Integral(0,hh1pMn->GetNbinsX()+1); 
        double n4 = hhNpMn->Integral(0,hhNpMn->GetNbinsX()+1);
        ntotall = n1 + n2 + n3 + n4;
        // Add to stk
        hh1p0n->SetFillColor(GetColor(1014));
        hstk->Add(hh1p0n);
        hhNp0n->SetFillColor(GetColor(1011));
        hstk->Add(hhNp0n);
        hh1pMn->SetFillColor(GetColor(1007));
        hstk->Add(hh1pMn);
        hhNpMn->SetFillColor(GetColor(kOrange));
        hstk->Add(hhNpMn);
        hstk->SetTitle(Form("Total %.0f events", ntotall));
      }
      else{}
    } // End of if(hstk)

    // 1D histogram
    TH1 * hh = dynamic_cast<TH1*> (lout->At(ii));
    if(hh){
      const TString tag = hh->GetName();
      // Pi0 rec. efficiency
      if(tag.Contains("hMatchedTruthPi0Momentum") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e008hTruthLeadingPiZeroE";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Leading photon rec. efficiency
      if(tag.Contains("hMatchedTruthldShowerEnergy") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e012hTruthLeadingPi0GammaP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Leading photon rec. efficiency
      if(tag.Contains("hMatchedTruthslShowerEnergy") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e013hTruthSubLeadingPi0GammaP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // OA eff.
      // Pi0 rec. efficiency
      if(tag.Contains("hMatchedTruthPi0OA") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e021hTruthPi0OA";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Leading photon rec. efficiency
      if(tag.Contains("hMatchedTruthldShowerOA") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e022hTruthLeadingPi0GammaOA";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Leading photon rec. efficiency
      if(tag.Contains("hMatchedTruthslShowerOA") && kMC && !tag.Contains("_ratio")){
        const TString TM_tmp = "e023hTruthSubLeadingPi0GammaOA";
        TH1D *htrue = (TH1D*)lout->FindObject(TM_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
    }
  } // End of for loop
}

void PlotUtils::DrawHist(TList *lout, const double plotscale, TList * overlayList, const TString outdir)
{
  TLatex tt;
  tt.SetNDC();
  // Setup gStyle
  gStyleSetup();
  // Loop over all histograms inside Tlist
  for(int ii=0; ii<lout->GetSize(); ii++){
    // Create new canvas for drawing histogram
    TCanvas * c1 = new TCanvas(Form("c1_%i", ii), "", 1200, 800);
    // Setup Pad
    PadSetup(c1);
    // Get the name of the histogram
    const TString tag = lout->At(ii)->GetName();
    // Skip TTree object in the list
    if(tag == "tree") continue;
    // Create histograms and stack histogram if they exist
    // 1D histogram
    TH1 * hh = dynamic_cast<TH1*> (lout->At(ii));
    // Stack histogram
    THStack * hstk = dynamic_cast<THStack *>(lout->At(ii));
    // 2D histogram
    TH2 * h2d = 0x0;
    // Data overlay histogram
    TH1D *holay = (TH1D*)overlayList->FindObject(tag);

    //=========================== 1D and 2D histograms ===========================//
    if(hh){
      // Cast to 2D histogram
      h2d = dynamic_cast<TH2 *>(hh);
      // Check if h2d exist
      if(h2d){
        SetTitleFormat(h2d);
        // Draw 2D histogram from its name (MC only)
        if(!tag.Contains("proj") && (tag.Contains("RES") || tag.Contains("REG"))) {
          //h2d->Draw("colz");
          if(tag.Contains("DIAG")){
            // Draw the diagonal line for rec. VS true histogram
            double max = h2d->GetXaxis()->GetBinUpEdge(h2d->GetXaxis()->GetLast());
            double min = h2d->GetXaxis()->GetBinLowEdge(h2d->GetXaxis()->GetFirst());
            TLine *line = new TLine(min,min,max,max);
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->Draw();
            h2d->Draw("colz");
          }
          if((tag.Contains("REG") || tag.Contains("RES")) && tag.Contains("_nor")){
            //TPad *h2_pad = new TPad("h2_pad", "h2_pad",0.0,0.0,1.0,1.0);
            //h2_pad->Draw();

            //TPad *proj1_pad = new TPad("proj1_pad", "proj1_pad",0.0,0.0,1.0,1.0);
            //proj1_pad->Draw();
            //h2_pad->cd();
            //h2d->Draw("COL");
            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.3, 1, 0.98);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.02);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();
            h2d->GetXaxis()->SetLabelOffset(999);
            h2d->GetXaxis()->SetLabelSize(0);
            h2d->GetYaxis()->SetTitleSize(0.07);
            h2d->GetYaxis()->SetTitleOffset(0.6);
            h2d->Draw("Colz");
            c1->Update();
            c1->cd();

            TRegexp re("_nor");
            TString name = tag;
            name(re) = "_projX";
            cout << "name: " << name << endl;
            TH1D *hcdf = (TH1D*)lout->FindObject(name);

            TPad *pad2 = new TPad(Form("pad2_%d",ii), Form("pad2_%d",ii), 0, 0, 1, 0.3);
             
            pad2->SetTopMargin(0.03);
            pad2->SetBottomMargin(0.42);
            pad2->Draw();
            pad2->cd();

            //Float_t rightmax = 1.1*hcdf->GetMaximum();
            //Float_t scale = abs(h2d->GetYaxis()->GetXmax()) + abs(h2d->GetYaxis()->GetXmin());
            //cout << "scale: " << scale << endl;
   
            //hcdf->SetLineColor(kBlue);
            //hcdf->SetFillStyle(3004);
            hcdf->GetXaxis()->SetTitle(h2d->GetXaxis()->GetTitle());
            hcdf->GetYaxis()->SetTitle("Event fraction");
            hcdf->GetXaxis()->SetLabelSize(0.1);
            hcdf->GetXaxis()->SetTitleSize(0.15);
            hcdf->GetYaxis()->SetLabelSize(0.1);
            hcdf->GetYaxis()->SetTitleSize(0.05);
            hcdf->GetYaxis()->SetTitleSize(0.1);
            hcdf->GetYaxis()->SetTitleOffset(0.4);
            hcdf->GetYaxis()->SetNdivisions(505);
            hcdf->GetXaxis()->CenterTitle();
            hcdf->GetYaxis()->CenterTitle();

            hcdf->GetXaxis()->SetTitleFont(22);
            hcdf->GetYaxis()->SetTitleFont(22);
            hcdf->GetXaxis()->SetTitleOffset(0.9);
            hcdf->Scale(1/hcdf->Integral(0,10000));
            //const Int_t x0 = hh->GetXaxis()->GetFirst();
            //const Int_t x1 = hh->GetXaxis()->GetLast();

            //for(Int_t ix=x0; ix<=x1; ix++){
              //double binCont = hcdf->GetBinContent(ix) - abs(h2d->GetYaxis()->GetXmin());
              //hcdf->SetBinContent(ix, binCont);
            //}
            
            //proj1_pad->cd();
            //proj1_pad->SetFillStyle(0);
            //hcdf->SetFillColorAlpha(kGreen,0.3);
            //hcdf->Draw("bar Y+");
            hcdf->SetFillColorAlpha(kBlue,0.3);
            //hcdf->Draw("bar Y+");
            hcdf->Draw("hist");
            c1->cd();
          }
        }
        if(!tag.Contains("proj") && tag.Contains("Bin")){
          h2d->SetMarkerSize(2);
          gStyle->SetPaintTextFormat("4.3f");
          h2d->Draw("colz text");
        }

        if(tag.Contains("sigma")){
          h2d->SetMarkerSize(2);
          //gStyle->SetPaintTextFormat("4.3f");
          h2d->Draw("colz text");
        }
      }
      // There is no h2d, we only have 1D histogram
      else {
        SetTitleFormat(hh);
        // Check if we have data overlay histogram
        if(holay){
          // Fracitonal plots
          if(tag.Contains("FCN")){
            auto lg = new TLegend(0.7,0.7,0.85,0.88);
            TF1 *fitfuncMC = 0x0;
            TF1 *fitfuncData = 0x0;
            if(tag.Contains("Start")){
              fitfuncMC = new TF1("fMC","gaus",-100,100);
              fitfuncMC->SetLineColor(kRed);
              fitfuncData = new TF1("fData","gaus",-100,100);
              fitfuncData->SetLineColor(kBlack);
            }
            else{
              fitfuncMC = new TF1("fMC",CauchyDens,-100,100,3);
              //TF1 *fitfuncMC = new TF1("fMC", "TMath::Voigt(x - [0], [1], [2], 4)", -1000,1000);
              Double_t parMC[3];
              parMC[0] = hh->GetMean(); 
              parMC[1] = hh->GetRMS();
              parMC[2] = hh->GetMaximum();
              fitfuncMC->SetParameters(parMC);
              fitfuncMC->SetParNames("Mean","FWHM","Constant");

              fitfuncMC->SetLineColor(kRed);

              fitfuncData = new TF1("fData",CauchyDens,-100,100,3);
              //TF1 *fitfuncData = new TF1("fData", "TMath::Voigt(x - [0], [1], [2], 4)", -1000,1000);
              Double_t parData[3];
              parData[0] = holay->GetMean(); 
              parData[1] = holay->GetRMS();
              parData[2] = holay->GetMaximum();
              fitfuncData->SetParameters(parData);
              fitfuncData->SetParNames("Mean","FWHM","Constant");

              fitfuncData->SetLineColor(kBlack);
              gPad->Update();
            }

            hh->Scale(1/hh->Integral(0,1000));
            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kRed);
            hh->SetLineColor(kRed);
            hh->SetLineWidth(1);

            hh->Fit("fMC");
            hh->Draw("E");
/*            
            gPad->Update();
            TPaveStats *st = (TPaveStats*)hh->GetListOfFunctions()->FindObject("stats");
            //st->SetOptFit(1);
            gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

            st->SetY1NDC(0.48); 
            st->SetY2NDC(0.68);
            st->SetX1NDC(0.7); 
            st->SetX2NDC(0.93);
            st->SetTextColor(kRed);
*/
            holay->Scale(1/holay->Integral(0,1000));
            
            DrawOverlay(holay);
            holay->Fit("fData");
            
            //holay->Draw("sames E");
/*
            gPad->Update();
            TPaveStats *st1 = (TPaveStats*)holay->GetListOfFunctions()->FindObject("stats");
            //st1->SetOptFit(1);
            gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

            st1->SetY1NDC(0.28); 
            st1->SetY2NDC(0.48);
            st1->SetX1NDC(0.7); 
            st1->SetX2NDC(0.93);
            st1->SetTextColor(kBlack);
*/
            TLegendEntry *le = lg->AddEntry(hh,"MC","lp");
            le->SetTextColor(kRed);
            lg->AddEntry(holay,"DATA","lp");
            lg->Draw("same");

            c1->Update();
          }
          else{
            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 1.);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.02);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();

            // Scale data to MC
            hh->Scale(plotscale);
            hh->GetYaxis()->CenterTitle();
            hh->GetYaxis()->SetTitleFont(22);
            hh->GetYaxis()->SetTitleSize(0.06);
            hh->GetYaxis()->SetTitleOffset(0.8);

            //hh->Draw("hist text");
            hh->Draw("hist");
            DrawOverlay(holay);
            c1->Update();
            c1->cd();
            TPad *pad2 = new TPad(Form("pad2_%d",ii), Form("pad2_%d",ii), 0, 0, 1, 0.2);
             
            pad2->SetTopMargin(0.03);
            pad2->SetBottomMargin(0.42);
            pad2->SetGridx();
            pad2->SetGridy();
            pad2->Draw();
            pad2->cd();

            TH1D *hratio = (TH1D*)holay->Clone(Form("hratio_%d",ii));
            hratio->SetTitle(" ");
            hratio->Divide(hh);
            DrawDataMCRatio(hratio);
            c1->cd();
          }

        }
        // No data overlay histogram
        else{

          // Draw overlay for raw and corrected shower resolution
          if(tag.Contains("Pi0MomNorm")){
            // Get the name of raw shower histogram
            TRegexp re("Norm");
            TString name = tag;
            name(re) = "EOA";
            TH1D *holayRaw = (TH1D*)lout->FindObject(name);    
            TString fit = tag;
            fit(re) = "Asym";
            TH1D *holayFit = (TH1D*)lout->FindObject(fit);
            if(holayRaw){
              if(holayRaw->GetMaximum() > hh->GetMaximum()) hh->SetMaximum(holayRaw->GetMaximum()*1.2);
              else hh->SetMaximum(hh->GetMaximum()*1.2);
              hh->SetFillStyle(4050);
              hh->SetLineColor(kRed);
              hh->SetFillColorAlpha(kRed, 0.35);
              hh->Draw("hist");
              c1->Update();
	            holayRaw->SetFillStyle(4050);
              holayRaw->SetLineColor(kGreen);
              holayRaw->SetFillColorAlpha(kGreen, 0.35);
              holayRaw->Draw("SAMES hist"); 
              c1->Update();
              if(holayFit){
                holayFit->SetLineColor(kMagenta);
                holayFit->SetLineWidth(2);
                holayFit->Draw("SAMES hist");
                c1->Update();
                cout << "E1E2 mean: " << hh->GetMean() << " RMS: " << hh->GetRMS() << endl;
                cout << "E1OA mean: " << holayRaw->GetMean() << " RMS: " << holayRaw->GetRMS() << endl;
                cout << "Asym mean: " << holayFit->GetMean() << " RMS: " << holayFit->GetRMS() << endl;
              }
              auto lg = new TLegend(0.6,0.5,0.85,0.88);

              lg->AddEntry(hh,"E1 + E2","f");
              TLegendEntry* l1 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",hh->GetMean(),hh->GetRMS()), "");
              l1->SetTextSize(0.03);
              l1->SetTextColor(kRed);
              lg->AddEntry(holayRaw,"E1 + #theta","f");
              TLegendEntry* l2 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",holayRaw->GetMean(),holayRaw->GetRMS()), "");
              l2->SetTextSize(0.03);
              l2->SetTextColor(kGreen);
              lg->AddEntry(holayFit,"Asymmetry","f");
              TLegendEntry* l3 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",holayFit->GetMean(),holayFit->GetRMS()), "");
              l3->SetTextSize(0.03);
              l3->SetTextColor(kMagenta);
              lg->Draw("same");

            }
            else hh->Draw("hist");
          }
          
          else if(tag.Contains("FitRes") && tag.Contains("projY")){
            TF1 *fCauchy = new TF1("fCauchy",CauchyDens,-1000,1000,3);
            Double_t par[3];
            par[0] = hh->GetMean(); 
            par[1] = hh->GetRMS();
            par[2] = hh->GetMaximum();
            fCauchy->SetParameters(par);
            fCauchy->SetParNames("Mean","FWHM","Constant");
            fCauchy->SetLineColor(kRed);
            hh->Fit("fCauchy");
            hh->SetMaximum(hh->GetMaximum()*1.2);
            hh->Draw("hist");
            fCauchy->Draw("same");
          }
          
          // Compare fitted results
          else if(tag.Contains("Compare")){
            // Get the name of raw shower histogram
            TRegexp re("Compare");
            TString name = tag;
            name(re) = "ComparePost";
            TH1D *holay = (TH1D*)lout->FindObject(name);
            if(holay){
            if(holay->GetMaximum() > hh->GetMaximum()) hh->SetMaximum(holay->GetMaximum()*1.2);
            else hh->SetMaximum(hh->GetMaximum()*1.2);
            //if(tag.Contains("E1")) hh->SetMaximum(350);
            //if(tag.Contains("E2")) hh->SetMaximum(200);
            //if(tag.Contains("OA")) hh->SetMaximum(200);
            auto lg = new TLegend(0.6,0.5,0.85,0.88);
            //hh->SetStats(0);
            hh->SetFillStyle(4050);
            hh->SetFillColor(24);
            hh->SetLineColor(24);
            hh->Draw("hist");
            //hh->SetName("Pre Fit");
            gPad->Update();
/*        
            TPaveStats *st = (TPaveStats*)hh->GetListOfFunctions()->FindObject("stats");
            st->SetY1NDC(0.7); 
            st->SetY2NDC(0.9);
            st->SetX1NDC(0.7); 
            st->SetX2NDC(0.93);
*/    
            holay->SetFillStyle(3001);
            holay->SetFillColor(46);
            holay->SetLineColor(46);
            //holay->SetStats(0);
            holay->Draw("SAMES hist");
            //holay->SetName("Post Fit");
            gPad->Update();
/*            
            TPaveStats *st_holay = (TPaveStats*)holay->GetListOfFunctions()->FindObject("stats");
            
            st_holay->SetY1NDC(0.47); 
            st_holay->SetY2NDC(0.67);
            st_holay->SetX1NDC(0.7); 
            st_holay->SetX2NDC(0.93);
*/

            lg->AddEntry(hh,"Before Fitting","f");
            TLegendEntry* l1 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",hh->GetMean(),hh->GetRMS()), "");
            l1->SetTextSize(0.03);
            lg->AddEntry(holay,"After Fitting","f");
            TLegendEntry* l2 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",holay->GetMean(),holay->GetRMS()), "");
            l2->SetTextSize(0.03);
            l2->SetTextColor(46);
            lg->Draw("same");
            }
          }
          
          else if(tag.Contains("CorMean")){
            TF1 *fCorMean = 0x0;
            if(tag.Contains("ProtonMomentum")){
              fCorMean = new TF1("fCorMean",CorrectionFCN,0,2,4);
              fCorMean->SetParameters(-178558,8.639592,0.08331107,-0.00658262);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }
            
            else if(tag.Contains("ShowerE")){
              fCorMean = new TF1("fCorMean",CorrectionFCN,0,1,4);
              fCorMean->SetParameters(-0.58,5.01,0.09,-0.18);
              fCorMean->SetLineColor(kRed);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }
            else if(tag.Contains("ShowerTheta")){
              TF1 *fCorMean_0 = new TF1("fCorMean_0","gaus",0,100);
              TF1 *fCorMean_1 = new TF1("fCorMean_1","gaus",100,180);
              
              fCorMean = new TF1("fCorMean","gaus(0)+gaus(3)",0,180);

              Double_t par[7];

              hh->Fit(fCorMean_0,"R");
              hh->Fit(fCorMean_1,"R+");

              fCorMean_0->GetParameters(&par[0]);
              fCorMean_1->GetParameters(&par[3]);

              fCorMean->SetParameters(par);

              hh->Fit(fCorMean,"R+");
              hh->Draw("hist");
              fCorMean->Draw("same");              
            }
            else if(tag.Contains("ShowerPhi")){
              fCorMean = new TF1("fCorMean","pol7",-280,280);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");              
            }

            else if(tag.Contains("ProtonTheta")){
              TF1 *fCorMean_0 = new TF1("fCorMean_0","gaus",0,100);
              TF1 *fCorMean_1 = new TF1("fCorMean_1","gaus",100,180);
              
              fCorMean = new TF1("fCorMean","gaus(0)+gaus(3)",0,180);

              Double_t par[7];

              hh->Fit(fCorMean_0,"R");
              hh->Fit(fCorMean_1,"R+");

              fCorMean_0->GetParameters(&par[0]);
              fCorMean_1->GetParameters(&par[3]);

              fCorMean->SetParameters(par);

              hh->Fit(fCorMean,"R+");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }
            else if(tag.Contains("ProtonPhi")){
              fCorMean = new TF1("fCorMean","pol5",-280,280);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }
            else if(tag.Contains("ShowerLDE")){
              fCorMean = new TF1("fCorMean",CorrectionFCN,0,1,4);
              fCorMean->SetParameters(-0.58,5.01,0.09,-0.18);
              fCorMean->SetLineColor(kRed);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }

            else if(tag.Contains("ShowerSLE")){
              fCorMean = new TF1("fCorMean",CorrectionFCN,0,1,4);
              fCorMean->SetParameters(-0.58,5.01,0.09,-0.18);
              fCorMean->SetLineColor(kRed);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }

            else {
              cout << "not found!!" << endl;
            }

          }

          else if(tag.Contains("_ratio")){

            hh->GetXaxis()->SetTitle("Energy (GeV)");
            hh->GetYaxis()->SetTitle("Efficiency");
 
            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->SetMinimum(0.0);
            hh->SetMaximum(1);
            //hh->Draw("E");
            hh->DrawCopy("HIST c");
            hh->SetFillColor(kBlack);
            hh->SetFillStyle(3001);
            hh->SetFillColorAlpha(kBlack, 0.35);
            hh->Draw("e2same");

            if(tag.Contains("m010hMatchedTruthPi0Momentum")){
              TString ld_name = "m011hMatchedTruthldShowerEnergy_ratio";
              TString sl_name = "m012hMatchedTruthslShowerEnergy_ratio";
              TH1D *holay_ld = (TH1D*)lout->FindObject(ld_name);
              TH1D *holay_sl = (TH1D*)lout->FindObject(sl_name);
              holay_ld->SetMarkerStyle(8);
              holay_ld->SetMarkerSize(1);
              holay_ld->SetMarkerColor(kRed);
              holay_ld->SetLineColor(kRed);
              holay_ld->SetLineWidth(1);
              holay_ld->SetMinimum(0.0);
              holay_ld->SetMaximum(1);
              //holay_ld->Draw("SAME E");
              holay_ld->DrawCopy("SAME HIST c");
              holay_ld->SetFillColor(kRed);
              holay_ld->SetFillStyle(3001);
              holay_ld->SetFillColorAlpha(kRed, 0.35);
              holay_ld->Draw("e2same");

              holay_sl->SetMarkerStyle(8);
              holay_sl->SetMarkerSize(1);
              holay_sl->SetMarkerColor(kBlue);
              holay_sl->SetLineColor(kBlue);
              holay_sl->SetLineWidth(1);
              holay_sl->SetMinimum(0.0);
              holay_sl->SetMaximum(1);
              //holay_sl->Draw("SAME E");

              holay_sl->DrawCopy("SAME HIST c");
              holay_sl->SetFillColor(kBlue);
              holay_sl->SetFillStyle(3001);
              holay_sl->SetFillColorAlpha(kBlue, 0.35);
              holay_sl->Draw("e2same");

              auto lg = new TLegend(0.54,0.54,0.85,0.88);
              lg->AddEntry(hh,"#pi^{0}-particle","l");
              lg->AddEntry(holay_ld,"Leanding Photon","l");
              lg->AddEntry(holay_sl,"SubLeanding Photon","l");
              lg->Draw("same");
            }
          }
          else if (tag.Contains("_mean") && tag.Contains("PurityVS")){
            // Get the name of raw shower histogram
            TRegexp re("Purity");
            TString name = tag;
            name(re) = "Completeness";
            TH1D *holay = (TH1D*)lout->FindObject(name);
            hh->SetMinimum(0);
            hh->SetMaximum(1.3);
            if(tag.Contains("IP")) hh->SetMaximum(1.6);
            hh->SetLineColor(kRed);
            hh->SetLineWidth(2);
            hh->DrawCopy("HIST c");
            hh->SetFillColor(kRed);
            hh->SetFillStyle(3001);
            hh->SetFillColorAlpha(kRed, 0.35);
            hh->Draw("e2same");
            
            holay->SetLineColor(kBlue);
            holay->SetLineWidth(2);
            holay->DrawCopy("SAME HIST c"); 
            holay->SetFillColor(kBlue);
            holay->SetFillStyle(3001);
            holay->SetFillColorAlpha(kBlue, 0.35);
            holay->Draw("e2same");
            
            auto lg = new TLegend(0.65,0.2,0.85,0.38);
            if(tag.Contains("IP")) lg = new TLegend(0.65,0.7,0.85,0.88);
            lg->AddEntry(hh,"Purity","l");
            lg->AddEntry(holay,"Completeness","l");
            lg->Draw("same");

          }

          else if(tag.Contains("sigma11") && tag.Contains("Mean")){
            // Get the name of raw shower histogram
            TRegexp re("sigma11");
            TString name = tag;
            name(re) = "sigma22";
            TH1D *holay = (TH1D*)lout->FindObject(name);
            hh->SetMinimum(0);
            hh->SetMaximum(0.5);
            hh->SetLineColor(kRed);
            hh->SetLineWidth(2);
            hh->DrawCopy("HIST c");
            hh->SetFillColor(kRed);
            hh->SetFillStyle(3001);
            hh->SetFillColorAlpha(kRed, 0.35);
            hh->Draw("e2same");
            
            holay->SetLineColor(kBlue);
            holay->SetLineWidth(2);
            holay->DrawCopy("SAME HIST c"); 
            holay->SetFillColor(kBlue);
            holay->SetFillStyle(3001);
            holay->SetFillColorAlpha(kBlue, 0.35);
            holay->Draw("e2same");
            
            auto lg = new TLegend(0.65,0.7,0.85,0.88);
            lg->AddEntry(hh,"Leading #gamma","l");
            lg->AddEntry(holay,"Subleading #gamma","l");
            lg->Draw("same");
          }
          else {
            hh->Draw("hist");
          }
        }
      }
    } // End of if(hh)

    ////=========================== stack histogram ===========================//
    else if (hstk) {
      // Get the data overlay histogram for stack
      holay = (TH1D*)overlayList->FindObject(tag+"_sum");
      hstk->Draw("hist");
      hstk->SetTitle("");
      // Draw legend if needed
      if(tag.Contains("stkTruth")){
        const TString tag = hstk->GetName();
        hstk->GetYaxis()->SetTitle("Candidates");
        hstk->GetXaxis()->SetTitle("#delta#alpha_{T} (deg)");
        if(tag.Contains("stkTruthPn")){
          hstk->GetXaxis()->SetTitle("p_{n} (GeV/c)");
        }
        
        // Get the name of raw shower histogram
        TRegexp re("stk");
        TString tmp = tag;
        tmp(re) = "h";
        // Get the histogram name
        TString name1p0n = tmp+"1p0n";
        TString nameNp0n = tmp+"Np0n";
        TString name1pMn = tmp+"1pMn";
        TString nameNpMn = tmp+"NpMn";
        // Find histograms in the list
        TH1D *hh1p0n = (TH1D*)lout->FindObject(name1p0n);
        TH1D *hhNp0n = (TH1D*)lout->FindObject(nameNp0n); 
        TH1D *hh1pMn = (TH1D*)lout->FindObject(name1pMn); 
        TH1D *hhNpMn = (TH1D*)lout->FindObject(nameNpMn); 
        SetTitleFormat(hstk,false);
        double lxoff = 0.45;
        if(tmp.Contains("Dalphat") || tmp.Contains("MomIniPi") || tmp.Contains("ThetaIniPi")) lxoff = 0.03;
        TLegend * lg = new TLegend(lxoff+0.15, 0.58, lxoff+0.38, 0.88);
        TString lheader("0.45<p_{p}<1, p^{s.l.}_{p}<0.45");
        lg->SetHeader(lheader);
        lg->AddEntry(hh1p0n, "1p0n", "f");
        lg->AddEntry(hhNp0n, "Np0n", "f");
        lg->AddEntry(hh1pMn, "1pMn", "f");
        lg->AddEntry(hhNpMn, "NpMn", "f");
        lg->Draw("same");
      }

      if(holay){
        const TString name = hstk->GetName();
        /*
        const TString name = hstk->GetName();
        cout << "Test name: " << name << endl;
        for(auto ie : typeMaps[name]){
          cout << "Test!! : typeMaps: " << ie << endl;
        }
        */
        if(tag.Contains("COMPOSE") && !tag.Contains("LOG")){
          TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 1.);
          pad1->SetBottomMargin(0.02);
          //  pad1->SetGridx();
          //  pad1->SetGridy();
          pad1->Draw();
          pad1->cd();
        
          // Scale data to MC
          ScaleStack(hstk, plotscale);

          hstk->GetYaxis()->SetTitle(holay->GetYaxis()->GetTitle());
          SetTitleFormat(hstk);

          TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
          hstk->SetMaximum(hsum->GetMaximum()*1.2);
          if(tag.Contains("Pi0Energy_COMPOSE")) hstk->SetMaximum(hsum->GetMaximum()*1.7);
          hstk->Draw("nostack HIST");
          DrawOverlay(holay);
          if(tag.Contains("Mass")){
            double ymax = hsum->GetMaximum()*1.2;
            TLine *line = new TLine(0.134977,0,0.134977,ymax);
            line->SetLineColor(kGray+3);
            line->SetLineStyle(2);
            line->Draw("sames");
          }

          hsum->SetLineColor(9);
          hsum->SetLineWidth(3);
          hsum->SetFillStyle(0);
          hsum->Draw("same HIST");
          auto* legend = new TLegend(0.55, 0.55, 0.89, 0.88);
          // Cheat Legend
          TH1D * h1 = 0x0;
          TH1D * h2 = 0x0;
          TH1D * h3 = 0x0;
          h1 = new TH1D(Form("h1%s", tag.Data()),  "", 20, -0.5, 19.5); 
          h2 = new TH1D(Form("h2%s", tag.Data()),  "", 20, -0.5, 19.5); 
          h3 = new TH1D(Form("h3%s", tag.Data()),  "", 20, -0.5, 19.5); 
          h1->SetFillStyle(3004);h2->SetFillStyle(3004);h3->SetFillStyle(3004);
          h1->SetFillColor(1505);h2->SetFillColor(1509);h1->SetLineColor(1505);h2->SetLineColor(1509);
          h3->SetFillColor(1502);h3->SetLineColor(1502);

          double tot = 0;
          for(auto ie : typeMaps[name]){
            tot += ie;
          }
          legend->AddEntry(h1, Form("Two Gammas same Pi0 (%.1f%s)",(typeMaps[name][0])/tot*100,"%"), "f");
          legend->AddEntry(h2, Form("Two Gammas diff Pi0 (%.1f%s)",(typeMaps[name][1])/tot*100,"%"), "f");
          legend->AddEntry(h3, Form("One Gamma (%.2f%s)",(typeMaps[name][2])/tot*100,"%"), "f");
          legend->AddEntry(hsum, "MC", "f");
          legend->Draw("same");

          c1->Update();
          c1->cd();
          TPad *pad2 = new TPad(Form("pad2_%d",ii), Form("pad2_%d",ii), 0, 0, 1, 0.2);
           
          pad2->SetTopMargin(0.03);
          pad2->SetBottomMargin(0.42);
          pad2->SetGridx();
          pad2->SetGridy();
          pad2->Draw();
          pad2->cd();

          TH1D *hratio = (TH1D*)holay->Clone(Form("hratio_%d",ii));
          hratio->SetTitle(" ");
          hratio->Divide(hsum);
          DrawDataMCRatio(hratio);
          c1->cd();

        }

        else if(tag.Contains("COMPOSE") && tag.Contains("LOG")){
          
          TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 1.);
          pad1->SetBottomMargin(0.02);
          //  pad1->SetGridx();
          //  pad1->SetGridy();
          pad1->Draw();
          pad1->cd();
          pad1->SetLogy();
          // Scale data to MC
          ScaleStack(hstk, plotscale);

          hstk->GetYaxis()->SetTitle(holay->GetYaxis()->GetTitle());
          SetTitleFormat(hstk);

          TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
          hstk->SetMaximum(hsum->GetMaximum()*10); 
          hstk->Draw("nostack HIST c");
         
          DrawOverlay(holay);

          hsum->SetLineColor(9);
          hsum->SetLineWidth(3);
          hsum->SetFillStyle(0);
          hsum->Draw("same HIST c");

          // Cheat Legend
          vector<TString> parType = FillLegendType("parType",hstk->GetName());

          vector<TString> htype = FillLegendStyle(1,"parType");//need to matcy parType

          const int mrks[]={1,1,1,1, 1,1,1,1,1,1, 1,6};
          int *cols=GetColorArray(parType.size());
          cols[parType.size()-1]= kBlack;
          TLegend * lg = 0x0;
          lg = DrawLegend(parType, htype, tag, cols, mrks, 4);
          lg->Draw("same");

          c1->Update();
          c1->cd();
          TPad *pad2 = new TPad(Form("pad2_%d",ii), Form("pad2_%d",ii), 0, 0, 1, 0.2);
           
          pad2->SetTopMargin(0.03);
          pad2->SetBottomMargin(0.42);
          pad2->SetGridx();
          pad2->SetGridy();
          pad2->Draw();
          pad2->cd();

          TH1D *hratio = (TH1D*)holay->Clone(Form("hratio_%d",ii));
          hratio->SetTitle(" ");
          hratio->Divide(hsum);
          DrawDataMCRatio(hratio);
          c1->cd();
        }

        else{
          TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 1.);
          pad1->SetBottomMargin(0.02);
          //  pad1->SetGridx();
          //  pad1->SetGridy();
          pad1->Draw();
          pad1->cd();

          // Scale data to MC
          ScaleStack(hstk, plotscale);
          
          hstk->GetYaxis()->SetTitle(holay->GetYaxis()->GetTitle());
          SetTitleFormat(hstk);

          TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
          hstk->SetMaximum(holay->GetMaximum()*1.2);
          hstk->Draw("hist");
          DrawOverlay(holay);
          //c1->Update();
          
          const int overlayColor = kBlack;

          // Event type legend
          if(hstk->GetNhists() <= 3){
            vector<TString> evtType = FillLegendType("evtType",hstk->GetName());
            vector<TString> htype = FillLegendStyle(2,"evtType");

            int *cols=GetColorArray(evtType.size());
            cols[3]=overlayColor;
            const int mrks[]={1,1,1,6};
            TLegend * lg = 0x0;
            lg = DrawLegend(evtType, htype, tag, cols, mrks);
            lg->Draw("same");
          }
          // Beam particle type
          else if(hstk->GetNhists() == 6){
            vector<TString> beamType = FillLegendType("beamType",hstk->GetName());
            vector<TString> htype = FillLegendStyle(2,"beamType");

            int *cols=GetColorArray(beamType.size());
            const int mrks[]={1,1,1,1,1,1,6};
            cols[beamType.size()-1]=overlayColor;

            TLegend * lg = 0x0;
            lg = DrawLegend(beamType, htype, tag, cols, mrks, 2);
            lg->Draw("same");
          }
          // Particle type legend
          else{
            vector<TString> parType = FillLegendType("parType",hstk->GetName());

            vector<TString> htype = FillLegendStyle(2,"parType");//need to matcy parType         

            const int mrks[]={1,1,1,1, 1,1,1,1,1,1, 1,6};
            int *cols=GetColorArray(parType.size());
            cols[parType.size()-1]=overlayColor;
            TLegend * lg = 0x0;
            lg = DrawLegend(parType, htype, tag, cols, mrks, 4);
            lg->Draw("same");

          }
          c1->Update();
          c1->cd();
          TPad *pad2 = new TPad(Form("pad2_%d",ii), Form("pad2_%d",ii), 0, 0, 1, 0.2);
           
          pad2->SetTopMargin(0.03);
          pad2->SetBottomMargin(0.42);
          pad2->SetGridx();
          pad2->SetGridy();
          pad2->Draw();
          pad2->cd();

          TH1D *hratio = (TH1D*)holay->Clone(Form("hratio_%d",ii));
          hratio->SetTitle(" ");
          hratio->Divide(hsum);
          DrawDataMCRatio(hratio);
          c1->cd();
        }
  
      }

      // Spectial case
      else if(tag.Contains("OVERLAY")){
        TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
        hstk->SetMaximum(hsum->GetMaximum()*1.2);

        hstk->GetXaxis()->SetTitle("Energy (GeV)");
        hstk->GetYaxis()->SetTitle("Candidates");
        SetTitleFormat(hstk,false);

        hstk->Draw("nostack");
        hsum->SetLineColor(kBlack);
        hsum->SetLineWidth(3);
        hsum->SetFillStyle(0);
        hsum->Draw("same");
        const TString overlayName = "e008hTruthLeadingPiZeroE";
        TH1D *holay_Pi0 = (TH1D*)lout->FindObject(overlayName);
        holay_Pi0->SetLineColor(kRed);
        holay_Pi0->SetLineWidth(3);
        holay_Pi0->Draw("same");
        auto* legend = new TLegend(0.55, 0.55, 0.85, 0.85);
        // Cheat Legend
        TH1D * E1 = 0x0;
        TH1D * E2 = 0x0;
        E1 = new TH1D(Form("E1%s", tag.Data()), "", 20, -0.5, 19.5); 
        E2 = new TH1D(Form("E2%s", tag.Data()),  "", 20, -0.5, 19.5); 
        E1->SetFillStyle(3004);E2->SetFillStyle(3004);
        E1->SetFillColor(1505);E2->SetFillColor(1509);E1->SetLineColor(1505);E2->SetLineColor(1509);

        legend->AddEntry(E1, "Leading photon", "f");
        legend->AddEntry(E2, "SubLeading photon", "f");
        legend->AddEntry(hsum, "Photon Spectrum", "f");
        legend->AddEntry(holay_Pi0, "#pi^{0} Spectrum", "f");
        legend->Draw("same");
        
      }  
    }
    else cout << "PlotUtils::DrawHist not found correct histogram!" << " name: " << tag << endl;

    tt.SetTextSize(0.035);
    tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP");
    c1->Print(outdir+"/"+tag+".png");
    
  } // End of for loop
}

THStack * PlotUtils::ConvertToStack(const TH2D * hh, const bool kMC, std::map<TString, vector<double>> &typeMaps)
{
  TString tag = hh->GetName();
  cout << "tag: " << tag << endl;
  const TString tit = hh->GetTitle();
  const char* Xtitle = hh->GetXaxis()->GetTitle();
  const char* Ytitle = hh->GetYaxis()->GetTitle();

  TString typ = "DATA";
  if(kMC) typ = "MC";
  const int nx = hh->GetNbinsX();
  const double xmin = hh->GetXaxis()->GetBinLowEdge(1);
  const double xmax = hh->GetXaxis()->GetBinUpEdge(nx);

  const int ny = hh->GetNbinsY();

  const double oldintegral = hh->Integral(0, 10000, 0, 1000000);

  const int y0 = 1;
  const int y1 = ny;

  double newintegral = 0;
  //THStack * stk = new THStack(tag+"_stack", tag);
  THStack * stk = new THStack(tag,tag);
  // Get color for each y value
  const int *col=GetColorArray(ny);

  // Save info for each component
  vector<double> TruthTypeVect;
  //need to take into account of overflow in y
  for(int iy = y0; iy<=y1; iy++){
    const double toty = hh->Integral(0, 100000, iy, iy);
    TruthTypeVect.push_back(toty);
    //cout << "iy: " << iy << " toty: " << toty << endl;
    if(toty<1E-12){
      continue;
    }

    TH1D * htmp = new TH1D(Form("%s%sy%d", typ.Data(), tag.Data(), iy), Form("My Stack;%s;%s", Xtitle, Ytitle), nx, xmin, xmax);
    for(int ix=0; ix<=nx+1; ix++){
      const double ientry = hh->GetBinContent(ix, iy);
      newintegral += ientry;

      htmp->SetBinContent(ix, ientry);
    }

    const int icol = GetColor(col[iy-y0]);//need constant map between y and color
    htmp->SetFillColor(icol);
    htmp->SetLineColor(kBlack);
    htmp->SetLineWidth(1);
    htmp->SetMarkerSize(2);
    if(tag.Contains("OVERLAY")){
      htmp->SetLineColor(icol);
      htmp->SetFillStyle(3004);
    }
    if(tag.Contains("COMPOSE")){
      htmp->SetLineColor(icol);
      htmp->SetFillStyle(3004);
      if(tag.Contains("LOG")){
        htmp->SetFillStyle(0);
        htmp->SetLineWidth(2);
      }
    }
    printf("PlotUtils::ConvertToStack %s adding y %f with color %d\n", tag.Data(), hh->GetYaxis()->GetBinCenter(iy), icol);
    stk->Add(htmp);
  }
  //cout << "newintegral: " << newintegral << endl;
  for(double in : TruthTypeVect){
    cout << "Type percent: " << in/newintegral*100 << "%" << endl;
  }
  if(oldintegral!=newintegral){
    printf("PlotUtils::ConvertToStack integral not matched! %s old %f new %f\n", tag.Data(), oldintegral, newintegral); exit(1);
  }
  if(kMC) typeMaps[tag] = TruthTypeVect;
  
  return stk;
}

TH1D * PlotUtils::GetStackedSum(THStack *stk)
{
  const TList * ll = stk->GetHists();
  const TString tag = stk->GetName();

  TH1D * hout = 0x0;
  if(stk->GetHists()){
    hout = (TH1D*)ll->At(0)->Clone(tag);
    hout->SetName(tag+"_sum");
    hout->SetTitle(tag);
    hout->SetDirectory(0);
    for(Int_t ii=1; ii<ll->GetEntries(); ii++){
      hout->Add((TH1D*)ll->At(ii));
    }

    hout->SetEntries(hout->Integral(0,10000));
 }
  return hout;
}
void PlotUtils::ScaleStack(THStack *stk, const double scale)
{
  const TList * ll = stk->GetHists();
  for(Int_t ii=0; ii<ll->GetEntries(); ii++){
    TH1D *hh = (TH1D*)ll->At(ii);
    hh->Scale(scale);
  }
}

TH2D * PlotUtils::NormalHist(const TH2D *hraw, const Double_t thres, const Bool_t kmax, TH1D * &hmean, TH1D * &hcdf)
{
  hcdf = GetCDF(hraw, Form("%s_cdf",hraw->GetName()));

  TH2D *hh=(TH2D*)hraw->Clone(Form("%s_nor",hraw->GetName()));
  hh->Scale(0);

  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  Double_t hmax = -1e10;
  Double_t hmin = 1e10;
  Double_t nent = 0;

  TH1D * projX = hraw->ProjectionX(Form("tmpnormalhist%s", hh->GetName()));
  hmean = (TH1D*)projX->Clone(Form("%s_mean",hraw->GetName()));
  hmean->Scale(0);
  for(Int_t ix=x0; ix<=x1; ix++){

    TH1D * sliceh = hraw->ProjectionY(Form("tmpnormalhist%sx%d", hh->GetName(), ix), ix, ix, "oe");
    const Double_t tot = sliceh->GetEntries();
    const Double_t mean = sliceh->GetMean();
    const Double_t RMS = sliceh->GetRMS();
    //cout << "name: " << hh->GetName() << " " << ix << " mean: " << mean << endl;
    const TString tag = hh->GetName();
    //if(tag.Contains("VSnHits")){
    hmean->SetBinContent(ix,mean);
    hmean->SetBinError(ix,RMS);
    //}
    TH1D * pdfh=0x0;


    if(tot>1E-12){
      nent += tot;

      Double_t imax = -999;

      imax = sliceh->GetBinContent(sliceh->GetMaximumBin());

      for(Int_t iy=y0; iy<=y1; iy++){
        const Double_t cont = kmax ? sliceh->GetBinContent(iy)/imax : pdfh->GetBinContent(iy);
        const Double_t ierr = kmax ? sliceh->GetBinError(iy)/imax   : pdfh->GetBinError(iy); 
        if(tot>thres && cont>0){
          hh->SetBinContent(ix, iy, cont);
          hh->SetBinError(ix,iy, ierr);
          if(cont>hmax) hmax = cont;
          if(cont<hmin) hmin = cont;
        }
      }
    }
    delete pdfh;
    delete sliceh;
  }

  hh->SetEntries(nent);
  hh->SetMinimum(0.99*hmin);
  hh->SetMaximum(1.1*hmax);

  TString xtit(hraw->GetXaxis()->GetTitle());
  if(xtit.Contains("(")){
    xtit=xtit(0, xtit.First('('));
  }

  TString ytit(hraw->GetYaxis()->GetTitle());
  if(ytit.Contains("(")){
    ytit=ytit(0, ytit.First('('));
  }

  hh->SetTitle(Form("f(%s|%s) %s", ytit.Data(), xtit.Data(), hraw->GetTitle()));

  return hh; 
}

THStack * PlotUtils::NormalizeStack(THStack * hstk)
{ 
  const TString tag = hstk->GetName();
  THStack * hout = new THStack(tag+"_normalized", tag);
  
  const TH1D * hsum = GetStackedSum(hstk);
  if(hsum){
    const double sumintegral = hsum->Integral(0,1000000);
    if(sumintegral<EPSILON){
      printf("PlotUtils::NormalizeStack sum integral 0 %e\n", sumintegral); exit(1);
    }
    const TList * ll = hstk->GetHists();
    for(int ii=0; ii<ll->GetEntries(); ii++){
      const TH1D * hold=(TH1D*) ll->At(ii);
      TH1D * htmp=(TH1D*)hold->Clone(Form("%scopy", hold->GetName()));
      
      htmp->Sumw2();
      htmp->Divide(hsum);
      hout->Add(htmp);
    }

    delete hsum;
  }
  return hout;
}

int PlotUtils::GetColor(const int col)
{
  if(col>=1000){
    return 1500+col-1000;
  }
  else{
    return col;
  }
}

int * PlotUtils::GetColorArray(const int minsize)
{
  const int col[]={
                   1005, 1009, 1002, kOrange,
                   1014, 1007, 1003, 1015,
                   1008, 1004, 1006, 1010,
                   1012, 1013, 1011, kGreen+3,
                   1008, 1009, 1002, 1011, 1014, 1007, 1003, 1015, 1005, 1008, 1009, 1002, 1011, 1014, 1007, 1003, 1015, 1005, 1008, 1009, 1002, 1011, 1014, 1007, 1003, 1015, 1005, 1008, 1009, 1002, 1011, 1014, 1007, 1003, 1015, 1005};
  
  const int nc = sizeof(col)/sizeof(int);
  if(nc<minsize){
    printf("PlotUtils::GetColorArray too small size %d %d\n", nc, minsize); exit(1);
  }
  int *outcl = new int[nc];
  for(int ii=0; ii<nc; ii++){
    outcl[ii] = col[ii];
  }

  return outcl;
}

void PlotUtils::PadSetup(TPad *currentPad, const Double_t currentLeft, const Double_t currentTop, const Double_t currentRight, const Double_t currentBottom)
{
  currentPad->SetTicks(1,1);
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);

  currentPad->SetFillColor(0);
}

void PlotUtils::IniColorCB()
{ 
  static bool kset = false;
  if(kset){
    printf("PlotUtils::IniColorCB arleady set\n");
    return;
  }
  else{
    printf("PlotUtils::IniColorCB creating new color\n");
  }
  const int fgkColorBase = 1500;  
  //http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
  Int_t id=fgkColorBase+1;
  new TColor(id++, 0./255., 0./255., 0./255., "CB1_Black",1.0);
  new TColor(id++, 0./255., 73./255., 73./255., "CB2_Forest",1.0);
  new TColor(id++, 0./255., 146./255., 146./255., "CB3_Teal",1.0);
  new TColor(id++, 255./255., 109./255., 182./255., "CB4_HotPink",1.0);
  new TColor(id++, 255./255., 182./255., 119./255., "CB5_BabyPink",1.0);
  new TColor(id++, 73./255., 0./255., 146./255., "CB6_Purple",1.0);
  new TColor(id++, 0./255., 109./255., 219./255., "CB7_RoyalBlue",1.0);
  new TColor(id++, 182./255., 109./255., 255./255., "CB8_Lilac",1.0);
  new TColor(id++, 109./255., 182./255., 255./255., "CB9_BlueGrey",1.0);
  new TColor(id++, 182./255., 219./255., 255./255., "CB10_SpaceWolves",1.0);
  new TColor(id++, 146./255., 0./255., 0./255., "CB11_Maroon",1.0);
  new TColor(id++, 146./255., 73./255., 0./255., "CB12_Tan",1.0);
  new TColor(id++, 219./255., 209./255., 0./255., "CB13_Orange",1.0);
  new TColor(id++, 36./255., 255./255., 36./255., "CB14_DayGleen",1.0);
  new TColor(id++, 255./255., 255./255., 109./255., "CB15_SunFlower",1.0);
  
  kset = true;
}

void PlotUtils::SetColor()
{
  gStyle->SetHistFillColor(0);
  //gStyle->SetFillColor(0);//it conflicts with color palette
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFillColor(0);

  gStyle->SetPalette(56);//only 56 available

  return;
  //---

  const Int_t nRGBs = 5;
  const Int_t nCont = 100;
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  const Int_t cgc = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);

  const Int_t nc = nCont;
  Int_t colors[nc];
  gStyle->SetNumberContours(nc);
  for(Int_t ii=0; ii<nc; ii++){
    colors[ii]=cgc+ii;
  }
  gStyle->SetPalette(nc, colors);
}


void PlotUtils::gStyleSetup()
{
  IniColorCB();
  //TGaxis::SetMaxDigits(2);
  gStyle->SetOptFit(1);
  // Set overall stats box
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetStatY(0.87);
  gStyle->SetStatH(0.12);
  gStyle->SetStatX(0.83);
  gStyle->SetStatW(0.12);
  //gStyle->SetStatColor(0);
  gStyle->SetStatStyle(1001);
  gStyle->SetStatBorderSize(2);
  gStyle->SetTitleX(0.55);
  
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(-1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(-1);
  gStyle->SetLegendBorderSize(-1);
  gStyle->SetPalette(1,0);  
  // Set color blind color
  SetColor();
  
  gROOT->ForceStyle();
}

void PlotUtils::getProfileFit(TH2D * h2d)
{
  TH1D *hprof = h2d->ProfileX();
  auto cctmp = new TCanvas("cctmp","cctmp",800,600);
  hprof->SetLineColor(kRed);
  hprof->SetLineWidth(3);
  hprof->Fit("pol1");
  hprof->SetStats(1);
  hprof->Draw();  
  cctmp->Print("output/hprof_withP.png");
}

TH1D * PlotUtils::GetCDF(const TH2D *hraw, const TString hname){

  TH1D * hcdf = hraw->ProjectionX(hname);
  TH1D * hold = (TH1D*) hcdf->Clone("hold");

  hcdf->Scale(0);

  const Int_t n2 = hold->GetNbinsX()+1;
  const Double_t ntot = hold->Integral(0, n2);

  for(Int_t ii=0; ii<= n2; ii++){
    const Double_t num = hold->Integral(0, ii);
    if(num<EPSILON){
      continue;
    }

    const Double_t nee = TMath::Sqrt(num);

    const Double_t fra = num/(ntot+EPSILON);
    const Double_t fee = nee/(ntot+EPSILON);

    hcdf->SetBinContent(ii, fra);
    hcdf->SetBinError(ii, fee);
  }

  delete hold;

  return hcdf;
}

void PlotUtils::DrawOverlay(TH1D *holay)
{
  holay->SetMarkerStyle(8);
  holay->SetMarkerSize(1);
  holay->SetMarkerColor(kBlack);
  holay->SetLineColor(kBlack);
  holay->SetLineWidth(1);
  holay->Draw("sames E");
}

void PlotUtils::getSliceXDrawY(TH2D * h2d)
{ 
  auto cc = new TCanvas("cc","cc",1600,1200);
  cc->Divide(4,5,0,0);
  
  for(int ii=1; ii<=h2d->GetNbinsX(); ii++){
    cc->cd(ii);
    TF1 *f1 = new TF1(Form("f1%d",ii),"gaus",-.5,.5);
    //TF1 *f1 = new TF1(Form("f1%d",ii),"[0]*TMath::Landau(x,[1],[2])",0,1);
    TH1D * hpj= h2d->ProjectionY(Form("f1%d",ii), ii, ii);
    if (hpj->GetEntries() != 0) {
      f1->SetParameters(hpj->GetMaximum(), hpj->GetMean(), hpj->GetRMS() );
      hpj->Fit(Form("f1%d",ii));
    }
    //hpj->SetStats(0);
    hpj->Draw();
  }  
    cc->Print("output/hSliceXDrawY.png");
}


void PlotUtils::xSlicedSigma(TH2D * h2d, TString tag)
{ 
  const TString name = h2d->GetName();
  for(int ii=1; ii<=h2d->GetNbinsX(); ii++){
    TH1D * hpj= h2d->ProjectionY(Form("ftmp%d",ii), ii, ii);
    if (hpj->GetEntries() != 0) {
      //double bin_center = hpj->GetBinCenter(ii);
      //cout << "bin_center: " << bin_center << " mean: " << mean << endl;
      if(tag.Contains("11")){
        AnaIO::hsigma11Mean->SetBinContent(ii,hpj->GetMean());
        AnaIO::hsigma11Mean->SetBinError(ii,hpj->GetMeanError());
      }
      if(tag.Contains("22")){
        AnaIO::hsigma22Mean->SetBinContent(ii,hpj->GetMean());
        AnaIO::hsigma22Mean->SetBinError(ii,hpj->GetMeanError());
      }
      
    }
    
  }
}

void PlotUtils::xSlicedEnergyCorrection(TH2D * h2d)
{ 
  auto cc = new TCanvas("cc","cc",1600,1200);
  cc->Divide(3,4,0,0);
  const TString name = h2d->GetName();
  for(int ii=1; ii<=h2d->GetNbinsX(); ii++){
    cc->cd(ii);
    TF1 *fitfunc = new TF1(Form("fitfunc%d",ii),CauchyDens,-0.6,0.6,3);
    //TF1 *fitfunc = new TF1("f", "TMath::Voigt(x - [0], [1], [2], 4)", -1, 1);
    TH1D * hpj= h2d->ProjectionY(Form("ftmp%d",ii), ii, ii);
    if (hpj->GetEntries() != 0) {
      // Sets initial values and parameter names
      Double_t par[3];
      par[0] = hpj->GetMean(); 
      par[1] = hpj->GetRMS();
      par[2] = hpj->GetMaximum();
      fitfunc->SetParameters(par);
      fitfunc->SetParNames("Mean","FWHM","Constant");

      hpj->SetMarkerStyle(kFullCircle);
      hpj->SetMarkerSize(0.5);
      hpj->SetMarkerColor(kBlack);
      hpj->SetLineColor(kBlack);
      hpj->SetLineWidth(1);
      
      hpj->Fit(Form("fitfunc%d",ii));
      fitfunc->Draw("same");
   
      if(name.Contains("protonPT")){
        AnaIO::hMeanPMomT->SetMinimum(-0.1);
        AnaIO::hMeanPMomT->SetMaximum(0.1);
        AnaIO::hMeanPMomT->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanPMomT->SetBinError(ii,fitfunc->GetParError(0));
      }
      if(name.Contains("protonPAll")){
        AnaIO::hMeanPMom->SetMinimum(-0.1);
        AnaIO::hMeanPMom->SetMaximum(0.1);
        AnaIO::hMeanPMom->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanPMom->SetBinError(ii,fitfunc->GetParError(0));
        cout << "proton all bin : " << AnaIO::hMeanPMom->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "proton all: " << fitfunc->GetParameter(0) << endl;
      }
      if(name.Contains("showerE")){
        AnaIO::hMeanShowerE->SetMinimum(-1);
        AnaIO::hMeanShowerE->SetMaximum(1);
        AnaIO::hMeanShowerE->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanShowerE->SetBinError(ii,fitfunc->GetParError(0));
        cout << "shower E: " << fitfunc->GetParameter(0) << endl;
      }
      if(name.Contains("showerTheta")){
        AnaIO::hMeanShowerTheta->SetMinimum(-30);
        AnaIO::hMeanShowerTheta->SetMaximum(30);
        if(ii == 13) AnaIO::hMeanShowerTheta->SetBinContent(ii,0);
        else if(ii == 14) AnaIO::hMeanShowerTheta->SetBinContent(ii,0);
        else {
          AnaIO::hMeanShowerTheta->SetBinContent(ii,fitfunc->GetParameter(0));
          AnaIO::hMeanShowerTheta->SetBinError(ii,fitfunc->GetParError(0));
        }
        cout << "shower Theta bin : " << AnaIO::hMeanShowerTheta->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "shower Theta: " << fitfunc->GetParameter(0) << endl;
      }
      if(name.Contains("showerPhi")){
        AnaIO::hMeanShowerPhi->SetMinimum(-30);
        AnaIO::hMeanShowerPhi->SetMaximum(30);
        AnaIO::hMeanShowerPhi->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanShowerPhi->SetBinError(ii,fitfunc->GetParError(0));
        cout << "shower Phi bin : "  << AnaIO::hMeanShowerPhi->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "shower Phi: " << fitfunc->GetParameter(0) << endl;
      }
      if(name.Contains("ProtonTheta")){
        AnaIO::hMeanProtonTheta->SetMinimum(-30);
        AnaIO::hMeanProtonTheta->SetMaximum(30);
        if(ii == 13) AnaIO::hMeanProtonTheta->SetBinContent(ii,0);
        else if(ii == 14) AnaIO::hMeanProtonTheta->SetBinContent(ii,0);
        else {
          AnaIO::hMeanProtonTheta->SetBinContent(ii,fitfunc->GetParameter(0));
          AnaIO::hMeanProtonTheta->SetBinError(ii,fitfunc->GetParError(0));
        }
        cout << "Proton Theta bin : " << AnaIO::hMeanProtonTheta->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "Proton Theta: " << fitfunc->GetParameter(0) << endl;
      }
      if(name.Contains("ProtonPhi")){
        AnaIO::hMeanProtonPhi->SetMinimum(-30);
        AnaIO::hMeanProtonPhi->SetMaximum(30);
        AnaIO::hMeanProtonPhi->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanProtonPhi->SetBinError(ii,fitfunc->GetParError(0));
        cout << "Proton Phi bin : "  << AnaIO::hMeanProtonPhi->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "Proton Phi: " << fitfunc->GetParameter(0) << endl;
      }

      if(name.Contains("showerLDE")){
        AnaIO::hMeanLDShowerE->SetMinimum(-1);
        AnaIO::hMeanLDShowerE->SetMaximum(1);
        AnaIO::hMeanLDShowerE->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanLDShowerE->SetBinError(ii,fitfunc->GetParError(0));
        cout << "shower LD E: " << fitfunc->GetParameter(0) << endl;
      }

      if(name.Contains("showerSLE")){
        AnaIO::hMeanSLShowerE->SetMinimum(-1);
        AnaIO::hMeanSLShowerE->SetMaximum(1);
        AnaIO::hMeanSLShowerE->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanSLShowerE->SetBinError(ii,fitfunc->GetParError(0));
        cout << "shower SL E: " << fitfunc->GetParameter(0) << endl;
      }
      

    }
    //hpj->SetStats(0);
    hpj->Draw("E");
  }

  cc->Print("output/"+name+"_sliced.png");
}

TLegend *PlotUtils::DrawLegend(const vector<TString> &entries, const vector<TString>& htype, const TString tag, const int *tmpcol, const int * tmpmkr, const int ncol)
{
  const int *defcol=GetColorArray();
  const int * cols=0x0;
  if(tmpcol){
    cols=tmpcol;
  }
  else{
    cols=defcol;
    printf("PlotUtils::DrawLegend using default color\n");
  }

  const int defmkr[]={20, 24, 21, 25, 22, 26, 23, 32, 34, 28, 29, 30, 20, 24, 21, 25, 22, 26, 23, 32, 34, 28, 29, 30};
  const int * mkrs=0x0;
  if(tmpmkr){
    mkrs=tmpmkr;
  }
  else{
    mkrs=defmkr;
    printf("PlotUtils::DrawLegend using default maker\n");
  }

  const int nent = entries.size();

  //SetGlobalStyle();
  TLegend * lg = 0x0;
  // Larger legend box
  //lg = new TLegend(0.55, 0.58, 0.88, 0.88);
  lg = new TLegend(0.65, 0.68, 0.88, 0.88);
  if(tag.Contains("COMPOSE")){
    lg = new TLegend(0.4, 0.55, 0.88, 0.88);
    if(tag.Contains("LOG")) lg = new TLegend(0.4, 0.6, 0.88, 0.9);
  }

  for(int ii=0; ii<nent; ii++){
    TH1D * hh=new TH1D(Form("h%d%s",ii,tag.Data()),"",1,0,1);
    const int col = GetColor(cols[ii]);
    hh->SetFillColor(col);
    if(tag.Contains("COMPOSE")){
      hh->SetFillColor(kWhite);
      hh->SetLineWidth(6);
    }
    hh->SetLineColor(col);
    hh->SetMarkerStyle(mkrs[ii]);
    hh->SetMarkerSize(3);
    hh->SetMarkerColor(col);
    lg->AddEntry(hh, entries[ii], htype[ii]);
  }

  lg->SetNColumns(ncol);

  return lg;
}

TH1D * PlotUtils::GetRecEfficiency(TH1 * hh, TH1D * htrue, const TString tag){

  TH1D *hEff = (TH1D*)hh->Clone();
  hEff->Scale(0);
  const int nx = hh->GetNbinsX();
  cout << tag.Data() << " rec Eff: " << hh->Integral(0,10000)/htrue->Integral(0,10000) << endl;
  cout << "num: " << hh->Integral(0,10000) << endl;
  cout << "demo: " << htrue->Integral(0,10000) << endl;
  for(int ix=0; ix<=nx+1; ix++){
    const double ientry = hh->GetBinContent(ix);
    const double error = hh->GetBinError(ix);
    const double ientry_demo = htrue->GetBinContent(ix);
    const double error_demo = htrue->GetBinError(ix);
    if(ientry_demo > 0 && ientry > 0){
      const double ratio = ientry/ientry_demo;
      const double error_bin = sqrt(ratio*ratio*(pow(error/ientry,2)+pow(error_demo/ientry_demo,2)));
      hEff->SetBinContent(ix, ratio);
      hEff->SetBinError(ix, error_bin);
      //cout << "ratio: " << ratio << endl;
      //cout << "error_bin: " << error_bin << endl;
    }
  }
  hEff->SetName(tag+"_ratio");
  return hEff;
}

void PlotUtils::SetTitleFormat(TH1 * hh){
  hh->SetTitle(" ");
  hh->GetYaxis()->CenterTitle();
  hh->GetYaxis()->SetTitleFont(22);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetYaxis()->SetTitleOffset(0.9);
  hh->GetXaxis()->CenterTitle();
  hh->GetXaxis()->SetTitleFont(22);
  hh->GetXaxis()->SetTitleSize(0.05);
  hh->GetXaxis()->SetTitleOffset(0.9);
  //hh->GetXaxis()->SetLabelOffset(999);
  //hh->GetXaxis()->SetLabelSize(0);
}

void PlotUtils::SetTitleFormat(TH2 * h2d){
  h2d->SetTitle(" ");
  h2d->GetYaxis()->CenterTitle();
  h2d->GetYaxis()->SetTitleFont(22);
  h2d->GetYaxis()->SetTitleSize(0.05);
  h2d->GetYaxis()->SetTitleOffset(0.9);
  h2d->GetXaxis()->CenterTitle();
  h2d->GetXaxis()->SetTitleFont(22);
  h2d->GetXaxis()->SetTitleSize(0.05);
  h2d->GetXaxis()->SetTitleOffset(0.9);
}

void PlotUtils::SetTitleFormat(THStack * stk, bool offSet){
  stk->SetTitle(" ");
  stk->GetYaxis()->CenterTitle();
  stk->GetYaxis()->SetTitleFont(22);
  stk->GetYaxis()->SetTitleSize(0.05);
  stk->GetYaxis()->SetTitleOffset(0.9);
  stk->GetXaxis()->CenterTitle();
  stk->GetXaxis()->SetTitleFont(22);
  stk->GetXaxis()->SetTitleSize(0.05);
  stk->GetXaxis()->SetTitleOffset(0.9);
  if(offSet) {
    stk->GetXaxis()->SetLabelOffset(999);
    stk->GetXaxis()->SetLabelSize(0);
  }
}

void PlotUtils::DrawDataMCRatio(TH1D * hratio){
  
  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetXaxis()->SetLabelSize(0.15);
  hratio->GetXaxis()->SetTitleSize(0.15);
  hratio->GetXaxis()->SetTitleOffset(1.);
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetYaxis()->SetTitleOffset(.3);
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->GetYaxis()->SetRangeUser(0,2);
  hratio->GetXaxis()->CenterTitle();
  hratio->GetYaxis()->CenterTitle();

  hratio->GetXaxis()->SetTitleFont(22);
  hratio->GetYaxis()->SetTitleFont(22);
  hratio->GetXaxis()->SetTitleOffset(0.7);
  hratio->GetXaxis()->SetTitleSize(0.25);

  hratio->Draw();
}

vector<TString> PlotUtils::FillLegendType(TString tag, TString name){
  vector<TString> tmpType;
  double tot = 0;
  for(auto ie : typeMaps[name]){
    tot += ie;
  }
  // Particle type
  if(tag.Contains("parType")){
    tmpType.push_back(Form("p (%.1f%s)",(typeMaps[name][0])/tot*100,"%"));
    tmpType.push_back(Form("#pi^{+} (%.1f%s)",(typeMaps[name][1])/tot*100,"%"));
    tmpType.push_back(Form("#pi^{#minus} (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back(Form("#gamma (%.1f%s)",(typeMaps[name][3])/tot*100,"%"));

    tmpType.push_back(Form("2ry p (%.1f%s)",(typeMaps[name][4])/tot*100,"%"));
    tmpType.push_back(Form("2ry #pi^{+} (%.1f%s)",(typeMaps[name][5])/tot*100,"%"));
    tmpType.push_back(Form("2ry #pi^{#minus} (%.1f%s)",(typeMaps[name][6])/tot*100,"%"));
    tmpType.push_back(Form("2ry #gamma (%.1f%s)",(typeMaps[name][7])/tot*100,"%"));
    tmpType.push_back(Form("2ry e^{#pm} (%.1f%s)",(typeMaps[name][8])/tot*100,"%"));
    tmpType.push_back(Form("2ry #mu^{#pm} (%.1f%s)",(typeMaps[name][9])/tot*100,"%"));

    tmpType.push_back(Form("others (%.1f%s)",(typeMaps[name][10])/tot*100,"%"));
    tmpType.push_back("data");
  }

  if(tag.Contains("evtType")){
    tmpType.push_back(Form("signal (%.1f%s)",(typeMaps[name][0])/tot*100,"%"));
    tmpType.push_back(Form("background (%.1f%s)",(typeMaps[name][1])/tot*100,"%"));
    tmpType.push_back(Form("non-#pi^{+} beam (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back("data");
  }

  if(tag.Contains("beamType")){
    tmpType.push_back(Form("p beam (%.1f%s)",(typeMaps[name][0])/tot*100,"%"));
    tmpType.push_back(Form("#pi^{+} (%.1f%s)",(typeMaps[name][1])/tot*100,"%"));
    tmpType.push_back(Form("#mu^{#pm} (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back(Form("e/#gamma (%.1f%s)",(typeMaps[name][3])/tot*100,"%"));
    tmpType.push_back(Form("#pi^{-} (%.1f%s)",(typeMaps[name][4])/tot*100,"%"));
    tmpType.push_back(Form("others (%.1f%s)",(typeMaps[name][5])/tot*100,"%"));
    tmpType.push_back("data");
  }

  return tmpType;
}

vector<TString> PlotUtils::FillLegendStyle(int opt, TString tag){
  vector<TString> tmpType;
  // Particle type
  if(tag.Contains("parType") && opt == 1){
    tmpType.push_back("l");
    tmpType.push_back("l");
    tmpType.push_back("l");
    tmpType.push_back("l");

    tmpType.push_back("l");
    tmpType.push_back("l");
    tmpType.push_back("l");
    tmpType.push_back("l");
    tmpType.push_back("l");
    tmpType.push_back("l");

    tmpType.push_back("l");
    tmpType.push_back("pl");
  }

  if(tag.Contains("parType") && opt == 2){
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");

    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");

    tmpType.push_back("f");
    tmpType.push_back("pl");
  }
  if(tag.Contains("evtType") && opt == 2){
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("ple");
  }

  if(tag.Contains("beamType") && opt == 2){
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("ple");
  }

  return tmpType;
}


