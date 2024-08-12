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


void PlotUtils::FillHist(TH1 * hh,  double xx, const double yy, const double & weight)
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
    hh->Fill(xx, weight);
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

void PlotUtils::FillHist(TH2 * hh,  double xx, const double yy, const double & weight)
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
    hh->Fill(xx, yy, weight);
  }
  else{
    hh->Fill(xx, weight);
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
        //TH2D * hnor = NormalHist(htmp, 5, true, hmean, hcdf);
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
        //TH2D * hnor = NormalHist(htmp, 5, true, hmean, hcdf);
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
      // CVM things not used
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
        hh1p0n->SetFillColor(GetColor(1014,true));
        hh1p0n->SetLineColor(GetColor(1014,true));

        hstk->Add(hh1p0n);
        hhNp0n->SetFillColor(GetColor(1011,true));
        hhNp0n->SetLineColor(GetColor(1011,true));

        hstk->Add(hhNp0n);
        hh1pMn->SetFillColor(GetColor(1007,true));
        hh1pMn->SetLineColor(GetColor(1007,true));

        hstk->Add(hh1pMn);
        hhNpMn->SetFillColorAlpha(kOrange,0.35);
        hhNpMn->SetLineColor(GetColor(kOrange,true));

        hstk->Add(hhNpMn);
        hstk->SetTitle(Form("Total %.0f events", ntotall));
        cout << "Total event: " << ntotall << endl;

      }
      else{}
    } // End of if(hstk)

    // 1D histogram
    TH1 * hh = dynamic_cast<TH1*> (lout->At(ii));
    if(hh){
      const TString tag = hh->GetName();

      // Beam rec. efficiency
      if(tag.Contains("hMatchedTruthBeamMomentum") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e000hTruthBeamMomentum";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Proton rec. efficiency
      if(tag.Contains("hMatchedTruthProtonMomentum") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e006hTruthProtonP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }

      // Leading Proton rec. efficiency
      if(tag.Contains("hMatchedTruthLeadingProtonMomentum") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e006hTruthLeadingProtonP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }

      // PiPlus rec. efficiency
      if(tag.Contains("hMatchedTruthPiPlusMomentum") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e028hTruthPiPlusP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }

      // Leading PiPlus rec. efficiency
      if(tag.Contains("hMatchedTruthLeadingPiPlusMomentum") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e037hTruthLeadingPiPlusP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }

      // Pi0 rec. efficiency
      if(tag.Contains("hMatchedTruthPi0Energy") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e008hTruthLeadingPiZeroE";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Leading photon rec. efficiency
      if(tag.Contains("hMatchedTruthLeadingShowerEnergy") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e012hTruthLeadingPi0GammaPBin1";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Leading pi0 photon rec. efficiency
      if(tag.Contains("hMatchedTruthldShowerEnergy") && kMC && !tag.Contains("_ratio")){
        const TString truth_tmp = "e012hTruthLeadingPi0GammaP";
        TH1D *htrue = (TH1D*)lout->FindObject(truth_tmp);
        TH1D *hEff = GetRecEfficiency(hh,htrue,tag);
        lout->Add(hEff);
      }
      // Subleading pi0 photon rec. efficiency
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

void PlotUtils::DrawHist(TList *lout, const double plotscale, const double plotscale_wt, TList * overlayList, const TString outdir)
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
          h2d->Draw("colz");
          /*if(tag.Contains("DIAG")){
            // Draw the diagonal line for rec. VS true histogram
            double max = h2d->GetXaxis()->GetBinUpEdge(h2d->GetXaxis()->GetLast());
            double min = h2d->GetXaxis()->GetBinLowEdge(h2d->GetXaxis()->GetFirst());
            TLine *line = new TLine(min,min,max,max);
            line->SetLineColor(kRed);
            line->SetLineWidth(2);
            line->Draw("sames");
            //h2d->Draw("colz");
          }*/
          if((tag.Contains("REG") || tag.Contains("RES")) && tag.Contains("_nor")){

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

            gPad->Update();

            TPaletteAxis* palette = (TPaletteAxis*)h2d->GetListOfFunctions()->FindObject("palette");
            gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn
            if(palette){
              palette->SetX1NDC(0.910-0.005);
              palette->SetY1NDC(0.02);
              palette->SetX2NDC(0.960-0.005);
              palette->SetY2NDC(0.90);
              palette->Draw();
            }
            c1->Update();

            if(tag.Contains("DIAG")){
              // Draw the diagonal line for rec. VS true histogram
              double max = h2d->GetXaxis()->GetBinUpEdge(h2d->GetXaxis()->GetLast());
              double min = h2d->GetXaxis()->GetBinLowEdge(h2d->GetXaxis()->GetFirst());
              TLine *line = new TLine(min,min,max,max);
              line->SetLineColor(1505);
              line->SetLineWidth(2);
              line->Draw("sames");
            }
            c1->Update();
            c1->cd();

            TRegexp re("_nor");
            TString name = tag;
            name(re) = "_projX";
            //cout << "name: " << name << endl;
            TH1D *htmp = (TH1D*)lout->FindObject(name);
            TH1D *hcdf = (TH1D*)htmp->Clone(Form("hcdf_%d",ii));

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
          else if(tag.Contains("InstXY") && !tag.Contains("nor") && !tag.Contains("mis")){

            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0.01, 0, 0.5, 0.98);
            pad1->Draw();
            pad1->cd();
            h2d->GetYaxis()->SetTitleOffset(1.01);
            h2d->Draw("colz");
            auto lgMC = new TLegend(0.7,0.7,0.85,0.88);
            lgMC->SetHeader("MC","C");
            lgMC->Draw("same");
            c1->Update();
            c1->cd();

            TPad *pad2 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0.51, 0, 1, 0.98);
            pad2->Draw();
            pad2->cd();
            SetTitleFormat(holay);
            holay->GetYaxis()->SetTitleOffset(1.01);
            holay->Draw("colz");
            auto lgDATA = new TLegend(0.6,0.7,0.85,0.88);
            lgDATA->SetHeader("DATA","C");
            lgDATA->Draw("same");
            c1->cd();
          }
          else{
            h2d->Draw("colz");
          }
        }
        else if(!tag.Contains("proj") && tag.Contains("Bin")){
          h2d->SetMarkerSize(2);
          gStyle->SetPaintTextFormat("4.3f");
          h2d->Draw("colz text");
        }

        else if(tag.Contains("sigma")){
          h2d->SetMarkerSize(2);
          //gStyle->SetPaintTextFormat("4.3f");
          h2d->Draw("colz text");
        }

        else if(tag.Contains("MAP")){
          //h2d->Divide(holay);
          //h2d->SetMarkerSize(2);
          TH2D *hh=(TH2D*)h2d->Clone(Form("%s_clone",h2d->GetName()));
          h2d->Scale(0);
          hh->Scale(plotscale);

          const Int_t x0 = h2d->GetXaxis()->GetFirst();
          const Int_t x1 = h2d->GetXaxis()->GetLast();
          const Int_t y0 = h2d->GetYaxis()->GetFirst();
          const Int_t y1 = h2d->GetYaxis()->GetLast();


          for(Int_t ix=x0; ix<=x1; ix++){
            for(Int_t iy=y0; iy<=y1; iy++){
              double MCbin = hh->GetBinContent(ix,iy);
              double DATAbin = holay->GetBinContent(ix,iy);
              //cout << "ix: " << ix << " iy: " << iy << " MC: " << MCbin << " DATA: " << DATAbin << " ratio: " << DATAbin/(MCbin+1e-10) << endl;
              h2d->SetBinContent(ix,iy,DATAbin/(MCbin+1e-10));
            }
          }

          h2d->GetZaxis()->SetRangeUser(0.5, 1.5);
          h2d->Draw("colz");
        }

        else{
          h2d->Draw("colz");
          if(tag.Contains("response")){
            cout << "tag res: " << tag << endl;
            double max = h2d->GetXaxis()->GetBinUpEdge(h2d->GetXaxis()->GetLast());
            double min = h2d->GetXaxis()->GetBinLowEdge(h2d->GetXaxis()->GetFirst());
            TLine *line = new TLine(min,min,max,max);
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->Draw();
          }
        }

      }
      // There is no h2d, we only have 1D histogram
      else {
        SetTitleFormat(hh);
        // Check if we have data overlay histogram
        if(holay){
          // Fracitonal plots
          if(tag.Contains("FCN")){
            auto lg = new TLegend(0.19,0.7,0.34,0.88);
            if(tag.Contains("BeamThetaZ_FCN") || tag.Contains("BeamDeltaE")) lg = new TLegend(0.63,0.7,0.78,0.88);
            TF1 *fitfuncMC = 0x0;
            TF1 *fitfuncData = 0x0;
            if(tag.Contains("Start") || tag.Contains("Inst") ){
              fitfuncMC = new TF1("fMC","gaus",-100,100);
              fitfuncMC->SetLineColor(kRed);
              fitfuncData = new TF1("fData","gaus",-100,100);
              fitfuncData->SetLineColor(kBlack);
            }
            else if(tag.Contains("BeamDeltaE")){
              fitfuncMC = new TF1("fMC","[0]*TMath::Landau(x,[1],[2])",0,3.5);
              Double_t parMC[3];
              parMC[0] = hh->GetMaximum();
              parMC[1] = hh->GetMean();
              parMC[2] = hh->GetRMS();
              fitfuncMC->SetParameters(parMC);
              fitfuncMC->SetParNames("Constant","MPV","Sigma");

              fitfuncData = new TF1("fData","[0]*TMath::Landau(x,[1],[2])",0,3.5);
              Double_t parData[3];
              parData[0] = holay->GetMaximum();
              parData[1] = holay->GetMean();
              parData[2] = holay->GetRMS();
              fitfuncData->SetParameters(parData);
              fitfuncData->SetParNames("Constant","MPV","Sigma");

              fitfuncData->SetLineColor(kBlack);
              gPad->Update();
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

            gPad->Update();
            TPaveStats *st = (TPaveStats*)hh->GetListOfFunctions()->FindObject("stats");
            //st->SetOptFit(1);
            gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

            st->SetY1NDC(0.48);
            st->SetY2NDC(0.68);
            st->SetX1NDC(0.7);
            st->SetX2NDC(0.93);
            st->SetTextColor(kRed);

            holay->Scale(1/holay->Integral(0,1000));

            DrawOverlay(holay);
            holay->Fit("fData");

            //holay->Draw("sames E");

            gPad->Update();
            TPaveStats *st1 = (TPaveStats*)holay->GetListOfFunctions()->FindObject("stats");
            //st1->SetOptFit(1);
            gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

            st1->SetY1NDC(0.28);
            st1->SetY2NDC(0.48);
            st1->SetX1NDC(0.7);
            st1->SetX2NDC(0.93);
            st1->SetTextColor(kBlack);

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
                //cout << "E1E2 mean: " << hh->GetMean() << " RMS: " << hh->GetRMS() << endl;
                //cout << "E1OA mean: " << holayRaw->GetMean() << " RMS: " << holayRaw->GetRMS() << endl;
                //cout << "Asym mean: " << holayFit->GetMean() << " RMS: " << holayFit->GetRMS() << endl;
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

          else if((tag.Contains("FitRes") || tag.Contains("FitRawRes")) && tag.Contains("projY")){
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

          else if(tag.Contains("FitGaus") && tag.Contains("projY")){
            TF1 *fGaus = new TF1(Form("fGaus%s",tag.Data()),"gaus",-100,100);
            if (hh->GetEntries() != 0) {
              fGaus->SetParameters(hh->GetMaximum(), hh->GetMean(), hh->GetRMS() );
              hh->Fit(Form("fGaus%s",tag.Data()));
            }
            //fGaus->SetParameters(par);
            //fGaus->SetParNames("Mean","Sigma","Constant");
            fGaus->SetLineColor(kRed);
            hh->SetMaximum(hh->GetMaximum()*1.2);
            hh->Draw("hist");
            fGaus->Draw("same");
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
            //hh->SetFillColorAlpha(24, 0.35);
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
            holay->SetFillStyle(4050);
            holay->SetFillColor(46);
            holay->SetLineColor(46);
            holay->SetFillColorAlpha(46, 0.35);
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

            else if(tag.Contains("BeamTheta")){
              fCorMean = new TF1("fCorMean","pol3",-280,280);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }
            else if(tag.Contains("BeamPhi")){
              fCorMean = new TF1("fCorMean","pol3",-280,280);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }

            if(tag.Contains("BeamMomentum")){
              fCorMean = new TF1("fCorMean",CorrectionFCN,0,2,4);
              fCorMean->SetParameters(-178558,8.639592,0.08331107,-0.00658262);
              hh->Fit("fCorMean");
              hh->Draw("hist");
              fCorMean->Draw("same");
            }


            else {
              cout << "not found!!" << endl;
            }

          }

          else if(tag.Contains("_ratio")){

            /*hh->GetXaxis()->SetTitle("Energy (GeV)");
            hh->GetYaxis()->SetTitle("Efficiency");
            if(tag.Contains("Proton")){
              hh->GetXaxis()->SetTitle("Truth Proton Momentum (GeV/c)");
              hh->GetYaxis()->SetTitle("Efficiency");
            }
            if(tag.Contains("PiPlus")){
              hh->GetXaxis()->SetTitle("Truth Proton Momentum (GeV/c)");
              hh->GetYaxis()->SetTitle("Efficiency");
            }*/
            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->SetMinimum(0.0);
            hh->SetMaximum(1.3);
            //hh->Draw("E");
            hh->DrawCopy("HIST c");
            hh->SetFillColor(kBlack);
            hh->SetFillStyle(1001);
            hh->SetFillColorAlpha(kBlack, 0.35);
            // Proton case
            if(tag.Contains("BeamMom")){
              auto lg = new TLegend(0.54,0.70,0.85,0.88);
              lg->AddEntry(hh,"Beam #pi^{+}","pl");
              lg->Draw("same");
            }
            // Proton case
            if(tag.Contains("Proton")){
              hh->SetMarkerColor(kGreen+3);
              hh->SetLineColor(kGreen+3);
              hh->SetFillColor(kGreen+3);
              hh->SetFillStyle(1001);
              hh->SetFillColorAlpha(kGreen+3, 0.35);
              auto lg = new TLegend(0.54,0.70,0.85,0.88);
              lg->AddEntry(hh,"Daughter proton","pl");
              lg->Draw("same");
            }
            // Leading shower case
            if(tag.Contains("LeadingShower")){
              hh->SetMarkerColor(kRed);
              hh->SetLineColor(kRed);
              hh->SetFillColor(kRed);
              hh->SetFillStyle(1001);
              hh->SetFillColorAlpha(kRed, 0.35);
              auto lg = new TLegend(0.54,0.70,0.85,0.88);
              lg->AddEntry(hh,"Daughter #gamma (from #pi^{0} decay)","pl");
              lg->Draw("same");
            }
            // Leading shower case
            if(tag.Contains("PiPlus")){
              hh->SetMarkerColor(kBlue);
              hh->SetLineColor(kBlue);
              hh->SetFillColor(kBlue);
              hh->SetFillStyle(1001);
              hh->SetFillColorAlpha(kBlue, 0.35);
              auto lg = new TLegend(0.54,0.70,0.85,0.88);
              lg->AddEntry(hh,"Daughter #pi^{+}","pl");
              lg->Draw("same");
            }
            hh->Draw("e2same");

            if(tag.Contains("m010hMatchedTruthPi0Energy")){
              hh->SetMaximum(1.0);
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
              holay_ld->SetFillStyle(1001);
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
              holay_sl->SetFillStyle(1001);
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
            hh->SetFillStyle(1001);
            hh->SetFillColorAlpha(kRed, 0.35);
            hh->Draw("e2same");

            holay->SetLineColor(kBlue);
            holay->SetLineWidth(2);
            holay->DrawCopy("SAME HIST c");
            holay->SetFillColor(kBlue);
            holay->SetFillStyle(1001);
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
            hh->SetFillStyle(1001);
            hh->SetFillColorAlpha(kRed, 0.35);
            hh->Draw("e2same");

            holay->SetLineColor(kBlue);
            holay->SetLineWidth(2);
            holay->DrawCopy("SAME HIST c");
            holay->SetFillColor(kBlue);
            holay->SetFillStyle(1001);
            holay->SetFillColorAlpha(kBlue, 0.35);
            holay->Draw("e2same");

            auto lg = new TLegend(0.65,0.7,0.85,0.88);
            lg->AddEntry(hh,"Leading #gamma","l");
            lg->AddEntry(holay,"Subleading #gamma","l");
            lg->Draw("same");
          }

          else if(tag.Contains("UpStreamELoss")){
            TF1 *f1 = new TF1(Form("f1%s",tag.Data()),"gaus",-100,100);
            if (hh->GetEntries() != 0) {
              f1->SetParameters(hh->GetMaximum(), hh->GetMean(), hh->GetRMS() );
              hh->Fit(Form("f1%s",tag.Data()));
            }
            double sigma = f1->GetParameter("Sigma");
            double mean = f1->GetParameter("Mean");
            // 3 sigma value
            double scraperValue = mean+sigma*3.0;
            cout << "tag: " << tag << " scraperValue: " << scraperValue << endl;
            cout << "sigma: " << sigma << " mean: " << mean << endl;

            // Draw the line
            TLine *line = new TLine(scraperValue,0,scraperValue,f1->GetMaximum());
            line->SetLineColor(kBlue);
            line->SetLineStyle(kDashed);
            line->SetLineWidth(2);
            line->Draw();
            //hpj->SetStats(0);
            hh->Draw("sames");
          }
/*
          // Total CEX xsec plots
          // i026hUnFoldedIncidentHist_xsec && i001hTruthIncidentHist_CEXxsec && i001hTruthIncidentHist_xsec (useless)
          else if (tag.Contains("xsec") && !tag.Contains("Total")  && !tag.Contains("Diff")){
            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 0.98);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.01);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();

            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->Draw("sames E1");

            //auto* legend = new TLegend(0.6, 0.6, 0.85, 0.88);
            TLegend * legend = new TLegend(0.15, 0.58, 0.48, 0.88);

            //hh->Draw("hist e");
            if(tag.Contains("CEX") || tag.Contains("UnFold")){
              g_cex->SetLineColor(kRed);
              g_cex->Draw("sames C");

              if(tag.Contains("i001hTruthIncidentHist_CEXxsec")){
                TString truth_Eslice = "i001hTruthIncidentHist_CalCEXxsec";
                TH1D *holay_truth_Eslice = (TH1D*)lout->FindObject(truth_Eslice);

                holay_truth_Eslice->SetMarkerStyle(25);
                holay_truth_Eslice->SetMarkerSize(1);
                holay_truth_Eslice->SetMarkerColor(kBlue);
                holay_truth_Eslice->SetLineColor(kBlue);
                holay_truth_Eslice->SetLineWidth(1);
                holay_truth_Eslice->Draw("sames E1");

                TString truth_Eslice_50MeVbin = "i001hTruthIncidentHist_CalCEXxsec_50MeVbin";
                TH1D *holay_truth_Eslice_50MeVbin = (TH1D*)lout->FindObject(truth_Eslice_50MeVbin);

                holay_truth_Eslice_50MeVbin->SetMarkerStyle(23);
                holay_truth_Eslice_50MeVbin->SetMarkerSize(1);
                holay_truth_Eslice_50MeVbin->SetMarkerColor(kOrange+3);
                holay_truth_Eslice_50MeVbin->SetLineColor(kOrange+3);
                holay_truth_Eslice_50MeVbin->SetLineWidth(1);
                holay_truth_Eslice_50MeVbin->Draw("sames E1");

                TString truth_Eslice_New = "i001hTruthIncidentHist_CalNewCEXxsec";
                TH1D *holay_truth_Eslice_New = (TH1D*)lout->FindObject(truth_Eslice_New);

                holay_truth_Eslice_New->SetMarkerStyle(22);
                holay_truth_Eslice_New->SetMarkerSize(1);
                holay_truth_Eslice_New->SetMarkerColor(kGreen+3);
                holay_truth_Eslice_New->SetLineColor(kGreen+3);
                holay_truth_Eslice_New->SetLineWidth(1);
                holay_truth_Eslice_New->Draw("sames E1");

                legend->AddEntry(g_cex, "Geant4 Prediction", "l");
                legend->AddEntry(hh, "MC Truth (Thin-Slice)", "lep");
                legend->AddEntry(holay_truth_Eslice, "MC Truth (E-Slice)", "lep");
                legend->AddEntry(holay_truth_Eslice_50MeVbin, "MC Truth (E-Slice 50MeV Bin)", "lep");
                legend->AddEntry(holay_truth_Eslice_New, "MC Truth (New Method)", "lep");

                legend->Draw("same");

              }
              else if(tag.Contains("i029")){
                legend->AddEntry(g_cex, "Geant4 Prediction", "l");
                legend->AddEntry(hh, "Data", "lep");
                legend->Draw("same");
              }
              else{
                legend->AddEntry(g_cex, "Geant4 Prediction", "l");
                legend->AddEntry(hh, "MC Reco", "lep");
                legend->Draw("same");
              }
            }
            else{
              g_inel->SetLineColor(kRed);
              g_inel->Draw("sames C");

              //g_inel->SetLineColor(kRed);
              //g_inel->Draw("sames C");

              legend->AddEntry(g_inel, "Geant4 Prediction", "l");
              legend->AddEntry(hh, "MC Reco", "lep");
              legend->Draw("same");
            }

            if(tag.Contains("i026")){
              TString truth = "i001hTruthIncidentHist_CEXxsec";
              TH1D *holay_truth = (TH1D*)lout->FindObject(truth);
              holay_truth->SetMarkerStyle(8);
              holay_truth->SetMarkerSize(1);
              holay_truth->SetMarkerColor(kBlue);
              holay_truth->SetMarkerStyle(kFullSquare);
              holay_truth->SetLineColor(kBlue);
              holay_truth->SetLineWidth(1);
              holay_truth->Draw("sames E1");
              legend->AddEntry(holay_truth, "MC Truth", "lep");
              TString lheader("Fake Data (1/2 MC sample)");
              legend->SetHeader(lheader);
              legend->Draw("same");
            }

            if(tag.Contains("i027")){
              TString truth = "i001hTruthIncidentHist_CalCEXxsec";
              TH1D *holay_truth = (TH1D*)lout->FindObject(truth);
              holay_truth->SetMarkerStyle(8);
              holay_truth->SetMarkerSize(1);
              holay_truth->SetMarkerColor(kBlue);
              holay_truth->SetMarkerStyle(kFullSquare);
              holay_truth->SetLineColor(kBlue);
              holay_truth->SetLineWidth(1);
              holay_truth->Draw("sames E1");
              legend->AddEntry(holay_truth, "MC Truth", "lep");
              TString lheader("Fake Data (1/2 MC sample)");
              legend->SetHeader(lheader);
              legend->Draw("same");
            }

            if(tag.Contains("i029")){
              TString truth = "i001hTruthIncidentHist_CalCEXxsec";
              TH1D *holay_truth = (TH1D*)lout->FindObject(truth);
              holay_truth->SetMarkerStyle(8);
              holay_truth->SetMarkerSize(1);
              holay_truth->SetMarkerColor(kBlue);
              holay_truth->SetMarkerStyle(kFullSquare);
              holay_truth->SetLineColor(kBlue);
              holay_truth->SetLineWidth(1);
              holay_truth->Draw("sames E1");
              legend->AddEntry(holay_truth, "MC Truth", "lep");

              TString reco = "i027hUnFoldedBeamIncidentHist_xsec";
              TH1D *holay_reco = (TH1D*)lout->FindObject(reco);
              holay_reco->SetMarkerStyle(8);
              holay_reco->SetMarkerSize(1);
              holay_reco->SetMarkerColor(kGreen+3);
              holay_reco->SetMarkerStyle(23);
              holay_reco->SetLineColor(kGreen+3);
              holay_reco->SetLineWidth(1);
              holay_reco->Draw("sames E1");
              legend->AddEntry(holay_reco, "MC Reco", "lep");

              TString lheader("Real Data (Full sample)");
              legend->SetHeader(lheader);
              legend->Draw("same");
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

            TH1D *hratio = (TH1D*)hh->Clone(Form("hratio_%d",ii));

            const Int_t x0_x = hh->GetXaxis()->GetFirst();
            const Int_t x1_x = hh->GetXaxis()->GetLast();

            hratio->Scale(0);

            for(Int_t ix=x0_x; ix<=x1_x; ix++){
              double bin_center = hh->GetBinCenter(ix);
              double MC = 0;
              if(tag.Contains("CEX") || tag.Contains("UnFold")) MC = g_cex->Eval(bin_center);
              else  MC = g_inel->Eval(bin_center);

              double data = hh->GetBinContent(ix);
              double edata = hh->GetBinError(ix);
              if(MC != 0){
                double ratio = data/MC;
                hratio->SetBinContent(ix,ratio);
                double error = sqrt(ratio*ratio*(pow(edata/data,2)));
                hratio->SetBinError(ix,error);
              }
              //cout << "bin: " << hratio->GetBinContent(ix) << endl;
            }
            hratio->SetTitle(" ");
            //hratio->Divide(hsum);
            DrawDataMCRatio(hratio,true);
            c1->cd();

          }

          // Three KE truth diff. xsec plots
          else if (tag.Contains("DiffCEXxsec") && !tag.Contains("Theta") && !tag.Contains("Cos")){

            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 0.98);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.01);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();
            // Thin-Slice
            //TString CEXname = "i001hTruthIncidentHist_CEXxsec";
            // E-Slice
            TString CEXname = "i001hTruthIncidentHist_CalCEXxsec";
            TH1D *CEXInteractingHist = (TH1D*)lout->FindObject(CEXname);

            //cout << " no norm: " << hh->Integral("width") << endl;
            hh->Scale(1/50.0);
            SetTitleFormat(hh);
            //cout << " mid norm: " << hh->Integral("width") << endl;
            //hh->Scale(g_cex->Eval(600));
            //cout << "evl: " << g_cex->Eval(600) << endl;
            if(tag.Contains("700")) hh->Scale(CEXInteractingHist->GetBinContent(14));
            else if(tag.Contains("800")) hh->Scale(CEXInteractingHist->GetBinContent(16));
            else if(tag.Contains("900")) hh->Scale(CEXInteractingHist->GetBinContent(18));

            //cout << "true scale(16): " << CEXInteractingHist->GetBinContent(16) << endl;

            //cout << "evl: " << g_cex->Eval(600) << endl;
            //cout << "evl hand: " << CEXInteractingHist->GetBinContent(12) << endl;
            //cout << "norm: " << hh->Integral("width") << endl;
            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->SetMaximum(0.4);
            hh->SetMinimum(0);

            hh->Draw("E");

            auto* legend = new TLegend(0.6, 0.6, 0.85, 0.88);

            if(tag.Contains("700")){
              g_675->SetLineColor(kRed);
              g_675->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 675 MeV", "lep");
              legend->AddEntry(g_675, "Geant4 Prediction", "l");
            }
            else if(tag.Contains("800")){
              g_775->SetLineColor(kRed);
              g_775->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 775 MeV", "lep");
              legend->AddEntry(g_775, "Geant4 Prediction", "l");
            }
            else if(tag.Contains("900")){
              g_875->SetLineColor(kRed);
              g_875->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 875 MeV", "lep");
              legend->AddEntry(g_875, "Geant4 Prediction", "l");
            }

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

            TH1D *hratio = (TH1D*)hh->Clone(Form("hratio_%d",ii));

            const Int_t x0_x = hh->GetXaxis()->GetFirst();
            const Int_t x1_x = hh->GetXaxis()->GetLast();

            hratio->Scale(0);

            for(Int_t ix=x0_x; ix<=x1_x; ix++){
              double bin_center = hh->GetBinCenter(ix);
              double MC = 0;
              if(tag.Contains("700")) MC = g_675->Eval(bin_center);
              else if(tag.Contains("800")) MC = g_775->Eval(bin_center);
              else if(tag.Contains("900")) MC = g_875->Eval(bin_center);
              double data = hh->GetBinContent(ix);
              double edata = hh->GetBinError(ix);
              if(MC != 0){
                double ratio = data/MC;
                hratio->SetBinContent(ix,ratio);
                double error = sqrt(ratio*ratio*(pow(edata/data,2)));
                hratio->SetBinError(ix,error);

              }
              //cout << "bin: " << hratio->GetBinContent(ix) << endl;
            }
            hratio->SetTitle(" ");
            //hratio->Divide(hsum);
            DrawDataMCRatio(hratio,true);

            c1->cd();
          }

          // Three Theta truth diff. xsec plots
          else if (tag.Contains("DiffCEXxsec") && tag.Contains("Theta") && !tag.Contains("Cos")){

            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 0.98);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.01);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();

            // Thin-Slice
            //TString CEXname = "i001hTruthIncidentHist_CEXxsec";
            // E-Slice
            TString CEXname = "i001hTruthIncidentHist_CalCEXxsec";
            TH1D *CEXInteractingHist = (TH1D*)lout->FindObject(CEXname);

            //cout << " no norm: " << hh->Integral("width") << endl;
            hh->Scale(1/20.0);
            //cout << " after norm: " << hh->Integral("width") << endl;

            SetTitleFormat(hh);

            if(tag.Contains("700")) hh->Scale(CEXInteractingHist->GetBinContent(14));
            else if(tag.Contains("800")) hh->Scale(CEXInteractingHist->GetBinContent(16));
            else if(tag.Contains("900")) hh->Scale(CEXInteractingHist->GetBinContent(18));


            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            //hh->SetMaximum(0.4);
            hh->SetMinimum(0);

            hh->Draw("E");

            auto* legend = new TLegend(0.6, 0.6, 0.85, 0.88);

            if(tag.Contains("700")){
              g_675theta->SetLineColor(kRed);
              g_675theta->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 675 MeV", "lep");
              legend->AddEntry(g_675theta, "Geant4 Prediction", "l");
            }
            else if(tag.Contains("800")){
              g_775theta->SetLineColor(kRed);
              g_775theta->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 775 MeV", "lep");
              legend->AddEntry(g_775theta, "Geant4 Prediction", "l");
            }
            else if(tag.Contains("900")){
              g_875theta->SetLineColor(kRed);
              g_875theta->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 875 MeV", "lep");
              legend->AddEntry(g_875theta, "Geant4 Prediction", "l");
            }

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

            TH1D *hratio = (TH1D*)hh->Clone(Form("hratio_%d",ii));

            const Int_t x0_x = hh->GetXaxis()->GetFirst();
            const Int_t x1_x = hh->GetXaxis()->GetLast();

            hratio->Scale(0);

            for(Int_t ix=x0_x; ix<=x1_x; ix++){
              double bin_center = hh->GetBinCenter(ix);
              double MC = 0;
              if(tag.Contains("700")) MC = g_675theta->Eval(bin_center);
              else if(tag.Contains("800")) MC = g_775theta->Eval(bin_center);
              else if(tag.Contains("900")) MC = g_875theta->Eval(bin_center);
              double data = hh->GetBinContent(ix);
              double edata = hh->GetBinError(ix);
              if(MC != 0){
                double ratio = data/MC;
                hratio->SetBinContent(ix,ratio);
                double error = sqrt(ratio*ratio*(pow(edata/data,2)));
                hratio->SetBinError(ix,error);

              }
              //cout << "bin: " << hratio->GetBinContent(ix) << endl;
            }
            hratio->SetTitle(" ");
            //hratio->Divide(hsum);
            DrawDataMCRatio(hratio,true);

            c1->cd();
          }


          // Three Cos Theta truth diff. xsec plots
          else if (tag.Contains("DiffCEXxsec") && tag.Contains("Theta") && tag.Contains("Cos")){

            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 0.98);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.01);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();

            // Thin-Slice
            //TString CEXname = "i001hTruthIncidentHist_CEXxsec";
            // E-Slice
            TString CEXname = "i001hTruthIncidentHist_CalCEXxsec";
            TH1D *CEXInteractingHist = (TH1D*)lout->FindObject(CEXname);

            //cout << " no norm: " << hh->Integral("width") << endl;
            hh->Scale(1/0.2);
            //cout << " after norm: " << hh->Integral("width") << endl;

            SetTitleFormat(hh);

            if(tag.Contains("700")) hh->Scale(CEXInteractingHist->GetBinContent(14));
            else if(tag.Contains("800")) hh->Scale(CEXInteractingHist->GetBinContent(16));
            else if(tag.Contains("900")) hh->Scale(CEXInteractingHist->GetBinContent(18));


            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->SetMaximum(240);
            hh->SetMinimum(0);

            hh->Draw("E");

            auto* legend = new TLegend(0.15, 0.6, 0.4, 0.88);

            if(tag.Contains("700")){
              g_675costheta->SetLineColor(kRed);
              g_675costheta->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 675 MeV", "lep");
              legend->AddEntry(g_675costheta, "Geant4 Prediction", "l");
            }
            else if(tag.Contains("800")){
              g_775costheta->SetLineColor(kRed);
              g_775costheta->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 775 MeV", "lep");
              legend->AddEntry(g_775costheta, "Geant4 Prediction", "l");
            }
            else if(tag.Contains("900")){
              g_875costheta->SetLineColor(kRed);
              g_875costheta->Draw("sames C");
              legend->AddEntry(hh, "Pion KE 875 MeV", "lep");
              legend->AddEntry(g_875costheta, "Geant4 Prediction", "l");
            }

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

            TH1D *hratio = (TH1D*)hh->Clone(Form("hratio_%d",ii));

            const Int_t x0_x = hh->GetXaxis()->GetFirst();
            const Int_t x1_x = hh->GetXaxis()->GetLast();

            hratio->Scale(0);

            for(Int_t ix=x0_x; ix<=x1_x; ix++){
              double bin_center = hh->GetBinCenter(ix);
              double MC = 0;
              if(tag.Contains("700")) MC = g_675costheta->Eval(bin_center);
              else if(tag.Contains("800")) MC = g_775costheta->Eval(bin_center);
              else if(tag.Contains("900")) MC = g_875costheta->Eval(bin_center);
              double data = hh->GetBinContent(ix);
              double edata = hh->GetBinError(ix);
              if(MC != 0){
                double ratio = data/MC;
                hratio->SetBinContent(ix,ratio);
                double error = sqrt(ratio*ratio*(pow(edata/data,2)));
                hratio->SetBinError(ix,error);

              }
              //cout << "bin: " << hratio->GetBinContent(ix) << endl;
            }
            hratio->SetTitle(" ");
            //hratio->Divide(hsum);
            DrawDataMCRatio(hratio,true);

            c1->cd();
          }

          // i031hUnFoldedPi0KEHist_DiffCEXRecoxsec800
          // Rec. plus Truth. final diff xsec plot
          else if (tag.Contains("_DiffCEXRecoxsec800") && !tag.Contains("Data")){

            TPad *pad1 = new TPad(Form("pad1_%d",ii), Form("pad1_%d",ii), 0, 0.2, 1, 0.98);
            // Set 0.01 will hide axis label
            pad1->SetBottomMargin(0.01);
            //  pad1->SetGridx();
            //  pad1->SetGridy();
            pad1->Draw();
            pad1->cd();

            //hh->Scale(1/50.0);
            SetTitleFormat(hh);
            // Thin-slice
            //TString CEXname = "i026hUnFoldedIncidentHist_xsec";
            // E-slice
            TString CEXname = "i027hUnFoldedBeamIncidentHist_xsec";
            TH1D *CEXInteractingHist = (TH1D*)lout->FindObject(CEXname);
            hh->Scale(CEXInteractingHist->GetBinContent(16));
            //cout << "scale(16): " << CEXInteractingHist->GetBinContent(16) << endl;

            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            //hh->SetMaximum(0.4);
            hh->SetMinimum(0.);

            //hh->Draw("E");

            if(!tag.Contains("Theta")){
              hh->Scale(1/50.0);
              hh->SetMaximum(0.4);
              hh->SetMinimum(0.);
              hh->Draw("E");
              TString TruthCEXname = "i001hTruthIncidentHist_DiffCEXxsec800";
              TH1D *holay_truth = (TH1D*)lout->FindObject(TruthCEXname);


              holay_truth->SetMarkerStyle(8);
              holay_truth->SetMarkerSize(1);
              holay_truth->SetMarkerColor(kBlue);
              holay_truth->SetMarkerStyle(kFullSquare);
              holay_truth->SetLineColor(kBlue);
              holay_truth->SetLineWidth(1);
              holay_truth->Draw("sames E1");

              auto* legend = new TLegend(0.6, 0.5, 0.85, 0.88);
              g_775->SetLineColor(kRed);
              g_775->Draw("sames C");
              legend->AddEntry(holay_truth, "MC Truth T_{beam #pi^{+}} = 775 MeV", "lep");
              legend->AddEntry(hh, "MC Reco T_{beam #pi^{+}} = 775 MeV", "lep");
              legend->AddEntry(g_775, "Geant4 Prediction", "l");

              TString lheader("Fake Data (Full MC sample)");
              legend->SetHeader(lheader);
              //legend->SetHeader("#splitline{Fake Data (Full MC sample)}{T_{beam #pi^{+}} = 775 MeV}");
              legend->Draw("same");
            }


            if(tag.Contains("CosTheta")){
              hh->Scale(1/0.2);
              //hh->SetMaximum(0.4);
              hh->SetMinimum(0.);
              hh->Draw("E");
              TString TruthCEXname = "i001hTruthIncidentHist_DiffCEXxsecCosTheta800";
              TH1D *holay_truth = (TH1D*)lout->FindObject(TruthCEXname);

              holay_truth->SetMarkerStyle(8);
              holay_truth->SetMarkerSize(1);
              holay_truth->SetMarkerColor(kBlue);
              holay_truth->SetMarkerStyle(kFullSquare);
              holay_truth->SetLineColor(kBlue);
              holay_truth->SetLineWidth(1);
              holay_truth->Draw("sames E1");

              auto* legend = new TLegend(0.15, 0.5, 0.45, 0.88);
              g_775costheta->SetLineColor(kRed);
              g_775costheta->Draw("sames C");
              legend->AddEntry(holay_truth, "MC Truth T_{beam #pi^{+}} = 775 MeV", "lep");
              legend->AddEntry(hh, "MC Reco T_{beam #pi^{+}} = 775 MeV", "lep");
              legend->AddEntry(g_775costheta, "Geant4 Prediction", "l");

              TString lheader("Fake Data (Full MC sample)");
              legend->SetHeader(lheader);
              //legend->SetHeader("#splitline{Fake Data (Full MC sample)}{T_{beam #pi^{+}} = 775 MeV}");
              legend->Draw("same");
            }

            if(tag.Contains("Theta") && !tag.Contains("Cos")){
              hh->Scale(1/20.0);

              //hh->SetMaximum(0.4);
              hh->SetMinimum(0.);
              hh->Draw("E");
              TString TruthCEXname = "i001hTruthIncidentHist_DiffCEXxsecTheta800";
              TH1D *holay_truth = (TH1D*)lout->FindObject(TruthCEXname);

              holay_truth->SetMarkerStyle(8);
              holay_truth->SetMarkerSize(1);
              holay_truth->SetMarkerColor(kBlue);
              holay_truth->SetMarkerStyle(kFullSquare);
              holay_truth->SetLineColor(kBlue);
              holay_truth->SetLineWidth(1);
              holay_truth->Draw("sames E1");

              auto* legend = new TLegend(0.6, 0.5, 0.85, 0.88);
              g_775theta->SetLineColor(kRed);
              g_775theta->Draw("sames C");
              legend->AddEntry(holay_truth, "MC Truth T_{beam #pi^{+}} = 775 MeV", "lep");
              legend->AddEntry(hh, "MC Reco T_{beam #pi^{+}} = 775 MeV", "lep");
              legend->AddEntry(g_775theta, "Geant4 Prediction", "l");

              TString lheader("Fake Data (Full MC sample)");
              legend->SetHeader(lheader);
              //legend->SetHeader("#splitline{Fake Data (Full MC sample)}{T_{beam #pi^{+}} = 775 MeV}");
              legend->Draw("same");
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

            TH1D *hratio = (TH1D*)hh->Clone(Form("hratio_%d",ii));

            const Int_t x0_x = hh->GetXaxis()->GetFirst();
            const Int_t x1_x = hh->GetXaxis()->GetLast();

            hratio->Scale(0);

            for(Int_t ix=x0_x; ix<=x1_x; ix++){
              double bin_center = hh->GetBinCenter(ix);
              double MC = 0;
              //if(tag.Contains("700")) MC = g_675->Eval(bin_center);
              //else if(tag.Contains("800")) MC = g_775->Eval(bin_center);
              //else if(tag.Contains("900")) MC = g_875->Eval(bin_center);
              if(!tag.Contains("Theta")) MC = g_775->Eval(bin_center);
              if(tag.Contains("CosTheta")) MC = g_775costheta->Eval(bin_center);
              if(tag.Contains("Theta") && !tag.Contains("Cos")) MC = g_775theta->Eval(bin_center);


              double data = hh->GetBinContent(ix);
              double edata = hh->GetBinError(ix);
              if(MC != 0){
                double ratio = data/MC;
                hratio->SetBinContent(ix,ratio);
                double error = sqrt(ratio*ratio*(pow(edata/data,2)));
                hratio->SetBinError(ix,error);

              }
              //cout << "bin: " << hratio->GetBinContent(ix) << endl;
            }
            hratio->SetTitle(" ");
            //hratio->Divide(hsum);
            DrawDataMCRatio(hratio);

            c1->cd();

          }

          else if(tag.Contains("i011hTruthDiffCEXInteractingHist_100MeV")){

            TString Name200MeV = "i010hTruthDiffCEXInteractingHist_200MeV";
            TString Name300MeV = "i009hTruthDiffCEXInteractingHist_300MeV";
            TString Name500MeV = "i008hTruthDiffCEXInteractingHist_500MeV";
            TString Name600MeV = "i005hTruthDiffCEXInteractingHist_700MeV";
            TString Name800MeV = "i006hTruthDiffCEXInteractingHist_800MeV";
            TString Name900MeV = "i007hTruthDiffCEXInteractingHist_900MeV";
            TString NameSum = "i012hTruthDiffCEXInteractingHist";

            TH1D *DiffCEXInteractingHist_900MeV = (TH1D*)lout->FindObject(Name900MeV);
            TH1D *DiffCEXInteractingHist_800MeV = (TH1D*)lout->FindObject(Name800MeV);
            TH1D *DiffCEXInteractingHist_600MeV = (TH1D*)lout->FindObject(Name600MeV);
            TH1D *DiffCEXInteractingHist_500MeV = (TH1D*)lout->FindObject(Name500MeV);
            TH1D *DiffCEXInteractingHist_300MeV = (TH1D*)lout->FindObject(Name300MeV);
            TH1D *DiffCEXInteractingHist_200MeV = (TH1D*)lout->FindObject(Name200MeV);

            TH1D *DiffCEXInteractingHist_Sum = (TH1D*)lout->FindObject(NameSum);


            const int * cols = GetColorArray();

            hh->SetMaximum(500);

            DiffCEXInteractingHist_900MeV->SetLineColor(GetColor(cols[0]));
            hh->Draw("sames C");

            DiffCEXInteractingHist_900MeV->SetLineColor(GetColor(cols[1]));
            DiffCEXInteractingHist_900MeV->Draw("sames C");

            DiffCEXInteractingHist_800MeV->SetLineColor(GetColor(cols[2]));
            DiffCEXInteractingHist_800MeV->Draw("sames C");

            DiffCEXInteractingHist_600MeV->SetLineColor(GetColor(cols[3]));
            DiffCEXInteractingHist_600MeV->Draw("sames C");

            DiffCEXInteractingHist_500MeV->SetLineColor(GetColor(cols[4]));
            DiffCEXInteractingHist_500MeV->Draw("sames C");

            DiffCEXInteractingHist_300MeV->SetLineColor(GetColor(cols[5]));
            DiffCEXInteractingHist_300MeV->Draw("sames C");

            DiffCEXInteractingHist_200MeV->SetLineColor(GetColor(cols[6]));
            DiffCEXInteractingHist_200MeV->Draw("sames C");

            DiffCEXInteractingHist_Sum->SetLineColor(kBlack);
            DiffCEXInteractingHist_Sum->Draw("sames C");


            auto* legend = new TLegend(0.72, 0.55, 0.85, 0.88);

            legend->AddEntry(hh, "100MeV", "f");
            legend->AddEntry(DiffCEXInteractingHist_200MeV, "200MeV", "f");
            legend->AddEntry(DiffCEXInteractingHist_300MeV, "300MeV", "f");
            legend->AddEntry(DiffCEXInteractingHist_500MeV, "500MeV", "f");
            legend->AddEntry(DiffCEXInteractingHist_600MeV, "700MeV", "f");
            legend->AddEntry(DiffCEXInteractingHist_800MeV, "800MeV", "f");
            legend->AddEntry(DiffCEXInteractingHist_900MeV, "900MeV", "f");

            legend->AddEntry(DiffCEXInteractingHist_Sum, "Total", "f");
            legend->Draw("same");

          }

          else if(tag.Contains("i017hRecoIncidentHist")){

            TString truthMatchedInc = "i019hTruthMatchedIncidentHist";
            TH1D *holay = (TH1D*)lout->FindObject(truthMatchedInc);

            TString recoBckSubInc = "i021hRecoBckSubIncidentHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInc);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInc);

            TString BckInc = "i023hNonPionBeamIncidentHist";
            TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->Draw("e1");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kBlue);
            holaybckSub->SetLineColor(kBlue);
            holaybckSub->SetLineWidth(1);
            holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holay->SetLineColor(kRed);
            holay->Draw("hist sames");

          }

          else if(tag.Contains("i018hRecoInteractingHist")){

            TString truthMatchedInt = "i020hTruthMatchedInteractingHist";
            TH1D *holay = (TH1D*)lout->FindObject(truthMatchedInt);

            TString recoBckSubInt = "i022hRecoBckSubInteractingHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInt);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInt);

            TString BckInc = "i024hNonPionBeamInteractingHist";
            TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            hh->Draw("e1");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kBlue);
            holaybckSub->SetLineColor(kBlue);
            holaybckSub->SetLineWidth(1);
            holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holay->SetLineColor(kRed);
            holay->Draw("hist sames");

          }

          else if(tag.Contains("i000hTruthBeamIncidentHist")){

            TString CalcInc = "i001hTruthBeamCalcIncidentHist";
            TH1D *hcalint = (TH1D*)lout->FindObject(CalcInc);

            TString CalcIncOld = "i000hTruthBeamIncidentHistOldM";
            TH1D *hcalintOld = (TH1D*)lout->FindObject(CalcIncOld);



            hh->SetLineColor(kRed);
            //hh->Rebin(50);
            hh->Draw("hist");

            hcalint->SetLineColor(kBlue);
            //hcalint->Rebin(50);
            hcalint->Draw("e1 sames");

            hcalintOld->SetLineColor(kGreen);
            //hcalintOld->Rebin(50);
            hcalintOld->Draw("hist sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "(IniE - IntE)/Eslice", "l");
            lg->AddEntry(hcalint, "Ini-Int Hist Substraction", "l");
            lg->AddEntry(hcalintOld, "Jake's Method", "l");


            lg->Draw("sames");

          }

          // This is the initial histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i000hTruthBeamInitialHist") && !tag.Contains("50MeVbin")){

            TString recoBckSubInc = "i016hRecoInitialHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInc);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            holaybckSub->Rebin(50);
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i024hUnFoldedInitialHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);
            holayunfold->Rebin(50);

            hh->Rebin(50);
            hh->SetMaximum(hh->GetMaximum()*1.5);

            hh->SetLineColor(kRed);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            //holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth", "l");
            lg->AddEntry(holaybckSub, "Bck Sub", "lp");
            lg->AddEntry(holayunfold, "Unfolded", "lp");

            lg->Draw("sames");

          }

          // This is the beam interacting histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i002hTruthBeamInteractingHist") && !tag.Contains("50MeVbin")){

            TString recoBckSubInc = "i016hRecoBeamInteractingHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInc);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            holaybckSub->Rebin(50);
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i024hUnFoldedBeamInteractingHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);
            holayunfold->Rebin(50);

            hh->Rebin(50);
            hh->SetMaximum(hh->GetMaximum()*1.5);

            hh->SetLineColor(kRed);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            //holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth", "l");
            lg->AddEntry(holaybckSub, "Bck Sub", "lp");
            lg->AddEntry(holayunfold, "Unfolded", "lp");

            lg->Draw("sames");

          }

          // This is the beam incident histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i001hTruthBeamCalcIncidentHist")){


            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i027hUnFoldedBeamIncidentHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);
            //cout << "i027hUnfoldedBeamIncidentHist total: " << holayunfold->Integral(0,10000) << endl;

            holayunfold->Rebin(50);

            //hh->Rebin(50);
            hh->SetMaximum(hh->GetMaximum()*1.5);

            hh->SetLineColor(kRed);
            hh->Draw("hist");


            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth", "l");
            //lg->AddEntry(holaybckSub, "Bck Sub", "lp");
            lg->AddEntry(holayunfold, "Unfolded", "lp");

            lg->Draw("sames");

          }

          // Data
          else if(tag.Contains("i029hUnFoldedBeamIncidentHistData")){

            TString truth = "i001hTruthBeamCalcIncidentHist";
            TH1D *htmp = (TH1D*)lout->FindObject(truth);
            TH1D *htruth = (TH1D*)htmp->Clone("hint");
            //htruth->Rebin(50);
            //TH1D *htruth = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i027hUnFoldedBeamIncidentHist";
            TH1D *htmp1 = (TH1D*)lout->FindObject(UnfoldInc);
            TH1D *holayunfold = (TH1D*)htmp1->Clone("hint1");

            //cout << "i027hUnfoldedBeamIncidentHist total: " << holayunfold->Integral(0,10000) << endl;
            double scale = hh->Integral(0,10000)/holayunfold->Integral(0,10000);

            hh->Rebin(50);
            holayunfold->Scale(scale);
            holayunfold->Rebin(50);
            holayunfold->SetMaximum(hh->GetMaximum()*1.5);
            holayunfold->SetMinimum(0.0);

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kGreen+3);
            holayunfold->SetLineColor(kGreen+3);
            holayunfold->SetLineWidth(1);
            holayunfold->Draw("e1");

            htruth->SetLineColor(kRed);
            htruth->Scale(scale);
            htruth->Draw("hist sames");



            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            //hh->Add(holaybck,-1);
            hh->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Real Data (Full sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(holayunfold, "MC Reco", "l");
            lg->AddEntry(htruth, "MC Truth", "l");
            lg->AddEntry(hh, "Data", "lp");

            lg->Draw("sames");

          }

          // Data
          else if(tag.Contains("i027hUnFoldedBeamInitialHistData")){

            TString truth = "i000hTruthBeamInitialHist";
            TH1D *htmp = (TH1D*)lout->FindObject(truth);
            TH1D *htruth = (TH1D*)htmp->Clone("hint");
            //htruth->Rebin(50);
            //TH1D *htruth = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i024hUnFoldedInitialHist";
            TH1D *htmp1 = (TH1D*)lout->FindObject(UnfoldInc);
            TH1D *holayunfold = (TH1D*)htmp1->Clone("hint1");
            //cout << "i027hUnfoldedBeamIncidentHist total: " << holayunfold->Integral(0,10000) << endl;
            double scale = hh->Integral(0,10000)/holayunfold->Integral(0,10000);

            hh->Rebin(50);
            holayunfold->Scale(scale);
            holayunfold->Rebin(50);
            holayunfold->SetMaximum(hh->GetMaximum()*1.5);
            holayunfold->SetMinimum(0.0);

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kGreen+3);
            holayunfold->SetLineColor(kGreen+3);
            holayunfold->SetLineWidth(1);
            holayunfold->Draw("e1");

            htruth->SetLineColor(kRed);
            htruth->Scale(scale);
            htruth->Draw("hist sames");

            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            //hh->Add(holaybck,-1);
            hh->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Real Data (Full sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(holayunfold, "MC Reco", "lp");
            lg->AddEntry(htruth, "MC Truth", "l");
            lg->AddEntry(hh, "Data", "lp");

            lg->Draw("sames");

          }

          // Data
          else if(tag.Contains("i028hUnFoldedBeamInteractingHistData")){

            TString truth = "i002hTruthBeamInteractingHist";
            TH1D *htmp = (TH1D*)lout->FindObject(truth);
            TH1D *htruth = (TH1D*)htmp->Clone("hint");
            //htruth->Rebin(50);
            //TH1D *htruth = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i024hUnFoldedBeamInteractingHist";
            TH1D *htmp1 = (TH1D*)lout->FindObject(UnfoldInc);
            TH1D *holayunfold = (TH1D*)htmp1->Clone("hint1");
            //cout << "i027hUnfoldedBeamIncidentHist total: " << holayunfold->Integral(0,10000) << endl;
            double scale = hh->Integral(0,10000)/holayunfold->Integral(0,10000);

            hh->Rebin(50);
            holayunfold->Scale(scale);
            holayunfold->Rebin(50);
            holayunfold->SetMaximum(hh->GetMaximum()*1.5);
            holayunfold->SetMinimum(0.0);

            //holayunfold->SetLineColor(kGreen+3);
            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kGreen+3);
            holayunfold->SetLineColor(kGreen+3);
            holayunfold->SetLineWidth(1);
            holayunfold->Draw("e1");

            htruth->SetLineColor(kRed);
            htruth->Scale(scale);
            htruth->Draw("hist sames");

            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            //hh->Add(holaybck,-1);
            hh->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Real Data (Full sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(holayunfold, "MC Reco", "lp");
            lg->AddEntry(htruth, "MC Truth", "l");
            lg->AddEntry(hh, "Data", "lp");

            lg->Draw("sames");

          }

          // Data
          else if(tag.Contains("i025hUnFoldedInteractingHistData")){

            TString truth = "i004hTruthCEXInteractingHist";
            TH1D *htmp = (TH1D*)lout->FindObject(truth);
            TH1D *htruth = (TH1D*)htmp->Clone("hint");
            //htruth->Rebin(50);
            //TH1D *htruth = (TH1D*)lout->FindObject(recoBckSubInc);
            TString UnfoldInc = "i025hUnFoldedInteractingHist";
            TH1D *htmp1 = (TH1D*)lout->FindObject(UnfoldInc);
            TH1D *holayunfold = (TH1D*)htmp1->Clone("hint1");
            //cout << "i027hUnfoldedBeamIncidentHist total: " << holayunfold->Integral(0,10000) << endl;
            double scale = hh->Integral(0,10000)/holayunfold->Integral(0,10000);

            holayunfold->SetLineColor(kGreen+3);
            holayunfold->Scale(scale);
            holayunfold->SetMaximum(holayunfold->GetMaximum()*1.5);
            holayunfold->SetMinimum(0.0);

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kGreen+3);
            holayunfold->SetLineColor(kGreen+3);
            holayunfold->SetLineWidth(1);
            holayunfold->Draw("e1");

            htruth->SetLineColor(kRed);
            htruth->Scale(scale);
            htruth->Draw("hist sames");




            hh->SetMarkerStyle(8);
            hh->SetMarkerSize(1);
            hh->SetMarkerColor(kBlack);
            hh->SetLineColor(kBlack);
            hh->SetLineWidth(1);
            //hh->Add(holaybck,-1);
            hh->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Real Data (Full sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(holayunfold, "MC Reco", "lp");
            lg->AddEntry(htruth, "MC Truth", "l");
            lg->AddEntry(hh, "Data", "lp");

            lg->Draw("sames");

          }


          // This is the incident histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i001hTruthIncidentHist")){

            TString CalcInc = "i001hTruthCalcIncidentHist";
            TH1D *hcalint = (TH1D*)lout->FindObject(CalcInc);

            TString recoBckSubInc = "i021hRecoBckSubIncidentHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInc);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInc);

            TString BckInc = "i023hNonPionBeamIncidentHist";
            TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            TString UnfoldInc = "i026hUnFoldedIncidentHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);

            hh->SetLineColor(kRed);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");


            hcalint->SetMarkerStyle(8);
            hcalint->SetMarkerSize(1);
            hcalint->SetMarkerColor(kBlue);
            hcalint->SetLineColor(kBlue);
            hcalint->SetLineWidth(1);
            //hcalint->Add(holaybck,-1);
            hcalint->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth", "l");
            lg->AddEntry(holaybckSub, "Bck Sub", "lp");
            lg->AddEntry(holayunfold, "Unfolded", "lp");
            lg->AddEntry(hcalint, "Calc", "lp");

            lg->Draw("sames");

          }
          // This is the interacting histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i002hTruthInteractingHist")){

            TString recoBckSubInt = "i022hRecoBckSubInteractingHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInt);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInt);

            TString BckInc = "i024hNonPionBeamInteractingHist";
            TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            TString UnfoldInc = "i025hUnFoldedInteractingHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);

            hh->SetLineColor(kRed);
            hh->SetMaximum(hh->GetMaximum()*1.5);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth", "l");
            lg->AddEntry(holaybckSub, "Bck Sub", "lp");
            lg->AddEntry(holayunfold, "Unfolded", "lp");
            lg->Draw("sames");
          }
          // This is the pi0 KE histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i006hTruthDiffCEXInteractingHist_800MeV")){

            TString recoBckSubInt = "i030hRecoPi0KEHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInt);

            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInt);

            //TString BckInc = "i030hRecoPi0KEHist";
            //TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            TString UnfoldInc = "i031hUnFoldedPi0KEHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);

            TString UnfoldIncData = "i031hUnFoldedPi0KEHistData";
            TH1D *holayunfoldData = (TH1D*)lout->FindObject(UnfoldIncData);

            hh->SetLineColor(kRed);
            hh->SetMaximum(hh->GetMaximum()*1.5);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            //holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlue);
            holayunfold->SetLineColor(kBlue);
            holayunfold->SetLineWidth(1);
            //holayunfold->Scale(0.5);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            holayunfoldData->SetMarkerStyle(8);
            holayunfoldData->SetMarkerSize(1);
            holayunfoldData->SetMarkerColor(kBlack);
            holayunfoldData->SetLineColor(kBlack);
            holayunfoldData->SetLineWidth(1);
            //holayunfoldData->Scale(0.5);
            //holayunfoldData->Add(holaybck,-1);
            holayunfoldData->Draw("e1 sames");

            TLegend * lg = new TLegend(0.55, 0.55, 0.85, 0.88);
            TString lheader("Fake Data (Full MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth T_{beam #pi^{+}} = 775 MeV", "l");
            lg->AddEntry(holaybckSub, "Bck Sub T_{beam #pi^{+}} = 775 MeV", "lp");
            lg->AddEntry(holayunfold, "Unfolded T_{beam #pi^{+}} = 775 MeV", "lp");
            lg->AddEntry(holayunfoldData, "Data", "lp");

            lg->Draw("sames");
          }

          // This is the pi0 Costheta histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i666hTruthDiffCEXInteractingHistCosTheta_800MeV")){

            TString recoBckSubInt = "i032hRecoPi0CosThetaHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInt);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInt);

            //TString BckInc = "i030hRecoPi0KEHist";
            //TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            TString UnfoldInc = "i033hUnFoldedPi0CosThetaHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);

            hh->SetLineColor(kRed);
            hh->SetMaximum(hh->GetMaximum()*1.1);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            //holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Scale(0.5);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.55, 0.45, 0.88);
            TString lheader("Fake Data (Full MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth T_{beam #pi^{+}} = 775 MeV", "l");
            lg->AddEntry(holaybckSub, "Bck Sub T_{beam #pi^{+}} = 775 MeV", "lp");
            lg->AddEntry(holayunfold, "Unfolded T_{beam #pi^{+}} = 775 MeV", "lp");
            lg->Draw("sames");
          }


          // This is the pi0 theta histogram plot (bck sub && unfolding && truth)
          else if(tag.Contains("i066hTruthDiffCEXInteractingHistTheta_800MeV")){

            TString recoBckSubInt = "i034hRecoPi0ThetaHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInt);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInt);

            //TString BckInc = "i030hRecoPi0KEHist";
            //TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);

            TString UnfoldInc = "i035hUnFoldedPi0ThetaHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);

            hh->SetLineColor(kRed);
            hh->SetMaximum(hh->GetMaximum()*1.1);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            //holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Scale(0.5);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.55, 0.55, 0.85, 0.88);
            TString lheader("Fake Data (Full MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth T_{beam #pi^{+}} = 775 MeV", "l");
            lg->AddEntry(holaybckSub, "Bck Sub T_{beam #pi^{+}} = 775 MeV", "lp");
            lg->AddEntry(holayunfold, "Unfolded T_{beam #pi^{+}} = 775 MeV", "lp");
            lg->Draw("sames");
          }

          // Same as i002 total Interacting Histogram
          else if(tag.Contains("i004hTruthCEXInteractingHist")){

            TString recoBckSubInt = "i022hRecoBckSubInteractingHist";
            TH1D *htmp = (TH1D*)lout->FindObject(recoBckSubInt);
            TH1D *holaybckSub = (TH1D*)htmp->Clone("hint");
            //TH1D *holaybckSub = (TH1D*)lout->FindObject(recoBckSubInt);
            //holaybckSub->Rebin(50);
            //holaybckSub->Scale(0.5);


            //TString BckInc = "i024hNonPionBeamInteractingHist";
            //TH1D *holaybck = (TH1D*)lout->FindObject(BckInc);


            TString UnfoldInc = "i025hUnFoldedInteractingHist";
            TH1D *holayunfold = (TH1D*)lout->FindObject(UnfoldInc);
            //holayunfold->Rebin(50);
            //holayunfold->Scale(0.5); //== FIXME whole sample



            //hh->Rebin(50);
            hh->SetLineColor(kRed);
            hh->SetMaximum(hh->GetMaximum()*1.5);
            hh->Draw("hist");

            holaybckSub->SetMarkerStyle(8);
            holaybckSub->SetMarkerSize(1);
            holaybckSub->SetMarkerColor(kGreen+3);
            holaybckSub->SetLineColor(kGreen+3);
            holaybckSub->SetLineWidth(1);
            //holaybckSub->Add(holaybck,-1);
            holaybckSub->Draw("e1 sames");

            holayunfold->SetMarkerStyle(8);
            holayunfold->SetMarkerSize(1);
            holayunfold->SetMarkerColor(kBlack);
            holayunfold->SetLineColor(kBlack);
            holayunfold->SetLineWidth(1);
            //holayunfold->Add(holaybck,-1);
            holayunfold->Draw("e1 sames");

            TLegend * lg = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("Fake Data (1/2 MC sample)");
            lg->SetHeader(lheader);
            lg->AddEntry(hh, "MC Truth", "l");
            lg->AddEntry(holaybckSub, "Bck Sub", "lp");
            lg->AddEntry(holayunfold, "Unfolded", "lp");
            lg->Draw("sames");

          }

          else if (tag.Contains("hTruthTestFFEnergyM1")){
            TString M2 = "i000hTruthTestFFEnergyM2";
            TH1D *hm2 = (TH1D*)lout->FindObject(M2);
            hh->SetMaximum(hm2->GetMaximum()*1.3);

            hh->SetLineColor(kRed);
            hh->Draw("hist");

            hm2->SetLineColor(kBlue);
            hm2->Draw("hist sames");

            TLegend * legend = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("True Front-Face Energy Compare");

            legend->SetHeader(lheader);
            legend->AddEntry(hh, "Method1 - EnergySlice", "l");
            legend->AddEntry(hm2, "Method2 - SpaceSlice", "l");
            legend->Draw("same");

          }

          else if (tag.Contains("hTruthTestIntEnergyM1")){
            TString M2 = "i100hTruthTestIntEnergyM2";
            TH1D *hm2 = (TH1D*)lout->FindObject(M2);
            hh->SetMaximum(hm2->GetMaximum()*1.3);

            hh->SetLineColor(kRed);
            hh->Draw("hist");

            hm2->SetLineColor(kBlue);
            hm2->Draw("hist sames");

            TLegend * legend = new TLegend(0.15, 0.58, 0.48, 0.88);
            TString lheader("True Interaction Energy Compare");

            legend->SetHeader(lheader);
            legend->AddEntry(hh, "Method1 - EnergySlice", "l");
            legend->AddEntry(hm2, "Method2 - SpaceSlice", "l");
            legend->Draw("same");

          }
*/
          // Cos theta angle plots
          else if(tag.Contains("OVERLAY") && tag.Contains("CosTheta") && tag.Contains("Low")){

            TString Mid = "i105hTruthPi0DaughtersCosTheta_Middle_OVERLAY";
            TString High = "i104hTruthPi0DaughtersCosTheta_High_OVERLAY";

            TH1D *holaymid = (TH1D*)lout->FindObject(Mid);
            TH1D *holayhigh = (TH1D*)lout->FindObject(High);


            SetTitleFormat(hh);
            hh->SetMaximum(holayhigh->GetMaximum()*1.5);

            hh->SetLineWidth(3);
            holaymid->SetLineWidth(3);
            holayhigh->SetLineWidth(3);

            //hh->SetFillColor(1505);
            hh->SetLineColor(1505);
            hh->Draw("hists");

            //holaymid->SetFillColor(1509);
            holaymid->SetLineColor(1509);
            holaymid->Draw("sames hists");

            //holayhigh->SetFillColor(1502);
            holayhigh->SetLineColor(1502);
            holayhigh->Draw("sames hists");


            //auto* legend = new TLegend(0.55, 0.55, 0.85, 0.88);
            TLegend * legend = new TLegend(0.15, 0.58, 0.48, 0.88);


            TString lheader("T_{beam #pi^{+}} = 650 - 800 MeV");
            //if(tag.Contains("Middle"))  lheader = ("T_{beam #pi^{+}} = 875 MeV, 200 < T_{#pi^{0}} < 400 MeV");
            //if(tag.Contains("Low")) lheader = ("T_{beam #pi^{+}} = 875 MeV, T_{#pi^{0}} < 150 MeV");

            legend->SetHeader(lheader);
            legend->AddEntry(hh, "T_{#pi^{0}} < 150 MeV", "l");
            legend->AddEntry(holaymid, "200 < T_{#pi^{0}} < 400 MeV", "l");
            legend->AddEntry(holayhigh, "T_{#pi^{0}} > 650 MeV", "l");
            legend->Draw("same");

          }

          else if(tag.Contains("_PreFit")){
            TH1D *hh = (TH1D*)lout->FindObject(tag);
            SetTitleFormat(hh);
            TRegexp re("_PreFit");
            TString name = tag;
            name(re) = "_PostFit";
            //cout << "name: " << name << endl;
            TH1D *hpostFit = (TH1D*)lout->FindObject(name);
            if(hpostFit){
              //if(hh->GetMaximum() > hpostFit->GetMaximum()) hh->SetMaximum(hh->GetMaximum()*1.15);
              //else hh->SetMaximum(hpostFit->GetMaximum()*1.15);
              if(tag.Contains("Energy") || tag.Contains("LD")) hh->SetMaximum(hpostFit->GetMaximum()*1.15);
              if(tag.Contains("SL")) hh->SetMaximum(hpostFit->GetMaximum()*1.15);
              if(tag.Contains("OAShower")) hh->SetMaximum(hpostFit->GetMaximum()*1.15);
              //hh->SetStats(1);
              hh->SetFillStyle(4050);
              //hh->SetFillColor(24);
              hh->SetFillColorAlpha(24,0.3);
              hh->SetLineColor(24);
              hh->SetLineWidth(2);

              hh->Draw("hist");
              // Stats box
              /*hh->SetName("Pre Fit");
              gPad->Update();
              TPaveStats *st = (TPaveStats*)hh->GetListOfFunctions()->FindObject("stats");
              gPad->Modified(); gPad->Update();
              st->SetY1NDC(0.7);
              st->SetY2NDC(0.9);
              st->SetX1NDC(0.7);
              st->SetX2NDC(0.93);
              */

              hpostFit->SetFillStyle(4050);
              hpostFit->SetFillColorAlpha(46,0.3);
              hpostFit->SetLineColor(46);
              hpostFit->SetLineWidth(2);

              // Must use SAMES
              hpostFit->Draw("SAMES hist");
              /*hpostFit->SetName("Post Fit");
              gPad->Update();
              TPaveStats *st_hpostFit = (TPaveStats*)hpostFit->GetListOfFunctions()->FindObject("stats");
              gPad->Modified(); gPad->Update();
              st_hpostFit->SetY1NDC(0.47);
              st_hpostFit->SetY2NDC(0.67);
              st_hpostFit->SetX1NDC(0.7);
              st_hpostFit->SetX2NDC(0.93);
              */

              TF1 *fCauchy_PreFit = new TF1("fCauchy_PreFit",CauchyDens,-0.4,0.4,3);
              TF1 *fGaus_PreFit = new TF1("fGaus_PreFit","gaus",-0.4,0.4);

              Double_t par[3];
              par[0] = hh->GetMean();
              par[1] = hh->GetRMS();
              par[2] = hh->GetMaximum();
              fCauchy_PreFit->SetParameters(par);
              fCauchy_PreFit->SetParNames("Mean","FWHM","Constant");
              fCauchy_PreFit->SetLineColor(kRed);
              hh->Fit("fCauchy_PreFit","R0");
              hh->Fit("fGaus_PreFit","R0");

              //fCauchy_PreFit->SetLineColor(kGray);
              //fCauchy_PreFit->Draw("sames C");

              //fGaus_PreFit->SetLineColor(kGray);
              //fGaus_PreFit->SetLineStyle(kDashed);
              //fGaus_PreFit->Draw("sames C");



              TF1 *fCauchy_PostFit = new TF1("fCauchy_PostFit",CauchyDens,-0.4,0.4,3);
              TF1 *fGaus_PostFit = new TF1("fGaus_PostFit","gaus",-0.4,0.4);
              Double_t par1[3];
              par1[0] = hpostFit->GetMean();
              par1[1] = hpostFit->GetRMS();
              par1[2] = hpostFit->GetMaximum();
              fCauchy_PostFit->SetParameters(par1);
              fCauchy_PostFit->SetParNames("Mean","FWHM","Constant");
              fCauchy_PostFit->SetLineColor(kRed);
              hpostFit->Fit("fCauchy_PostFit","R0");
              hpostFit->Fit("fGaus_PostFit","R0");
              //fCauchy_PostFit->Draw("sames C");

              //fGaus_PostFit->SetLineStyle(kDashed);
              //fGaus_PostFit->Draw("sames C");

              cout << "tag: " << tag << endl;
              cout << "PreFit: " << " Mean: " << fCauchy_PreFit->GetParameter(0) << " FWHM: " << fCauchy_PreFit->GetParameter(1) << endl;
              cout << "PostFit: " << " Mean: " << fCauchy_PostFit->GetParameter(0) << " FWHM: " << fCauchy_PostFit->GetParameter(1) << endl;

              cout << "Gaus PreFit: " << " Mean: " << fGaus_PreFit->GetParameter(1) << " FWHM: " << fGaus_PreFit->GetParameter(2) << endl;
              cout << "Gaus PostFit: " << " Mean: " << fGaus_PostFit->GetParameter(1) << " FWHM: " << fGaus_PostFit->GetParameter(2) << endl;


              auto lg = new TLegend(0.62,0.5,0.85,0.88);
              lg->AddEntry(hh,"Before Fitting","f");
              TLegendEntry* l1 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",fCauchy_PreFit->GetParameter(0),abs(fCauchy_PreFit->GetParameter(1))), "");
              l1->SetTextSize(0.03);
              lg->AddEntry(hpostFit,"After Fitting","f");
              TLegendEntry* l2 = lg->AddEntry((TObject*)0, Form("( #mu: %.2f #sigma: %.2f )",fCauchy_PostFit->GetParameter(0),abs(fCauchy_PostFit->GetParameter(1))), "");
              l2->SetTextSize(0.03);
              l2->SetTextColor(46);
              lg->Draw("same");
              TLine *line1 = new TLine(0,0,0,hh->GetMaximum());
              line1->SetLineColor(kBlack);
              line1->SetLineStyle(kDashed);
              line1->SetLineWidth(3);
              line1->Draw("same");
            }
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
      // Draw TKI legend if needed
      if(tag.Contains("stkTruth")){
        const TString tag = hstk->GetName();
        hstk->GetYaxis()->SetTitle("Event / degree");
        hstk->GetXaxis()->SetTitle("p#pi^{0} - #delta#alpha_{T} (deg)");
        double binWidth = 10;
        if(tag.Contains("stkTruthPn")){
          hstk->GetYaxis()->SetTitle("Event / GeV");
          hstk->GetXaxis()->SetTitle("p#pi^{0} - p_{n} (GeV/c)");
          binWidth = 0.05;
        }
        hstk->SetMaximum(hstk->GetMaximum()*1.3);
        ScaleStack(hstk,1/binWidth);
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
        TLegend * lg = new TLegend(lxoff+0.15, 0.58, lxoff+0.38, 0.85);
        //TString lheader("0.45<p_{p}<1, p^{s.l.}_{p}<0.45");
        //lg->SetHeader(lheader);
        lg->SetHeader("#splitline{0.45 < p_{p} < 1 GeV/c}{p^{s.l.}_{p} < 0.45 GeV/c}");

        lg->AddEntry(hh1p0n, "1p0n", "f");
        lg->AddEntry(hhNp0n, "Np0n", "f");
        lg->AddEntry(hh1pMn, "1pMn", "f");
        lg->AddEntry(hhNpMn, "NpMn", "f");
        lg->SetNColumns(2);

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

          pad1->Draw();
          pad1->cd();

          // Scale the beam cuts related hists data to MC no weight
          if(tag.Contains("CutPDG") || tag.Contains("BeamPDG")) ScaleStack(hstk, 0.702223); // CutPDG
          else if(tag.Contains("CutPandora") || tag.Contains("BeamPandora")) ScaleStack(hstk, 0.709511); // CutPandora
          else if(tag.Contains("CutCaloSize") || tag.Contains("BeamCalo")) ScaleStack(hstk, 0.703488); // CutCaloSize
          else if(tag.Contains("CutBeamQuality") || tag.Contains("BeamQuality")) ScaleStack(hstk, 0.698653); // CutBeamQuality
          else if(tag.Contains("CutAPA3") || tag.Contains("BeamAPA3")) ScaleStack(hstk, 0.658395); // CutAPA3
          else if(tag.Contains("CutMichelScore") || tag.Contains("BeamMichel")) ScaleStack(hstk, 0.65585);  // CutMichelScore
          else if(tag.Contains("CutChi2DOF") || tag.Contains("BeamChi2")) ScaleStack(hstk, 0.649126); // CutChi2DOF
          else if(tag.Contains("CutBeamScraper") || tag.Contains("BeamScraper")) ScaleStack(hstk, 0.57804); // CutBeamScraper(final beam scale)

          else if(tag.Contains("hRecPiPlus") && !tag.Contains("b00")) ScaleStack(hstk, plotscale_wt);
          else if(tag.Contains("i087hRecPiZeroKineticEnergyEvt_COMPOSE_XsEvt")) ScaleStack(hstk, plotscale_wt);

          // Scale data to MC
          else ScaleStack(hstk, plotscale);

          hstk->GetYaxis()->SetTitle(holay->GetYaxis()->GetTitle());
          SetTitleFormat(hstk);

          TList* histList = hstk->GetHists(); // Get the list of histograms in the stack
          TIter next(histList);
          TH1* hist;
          while ((hist = dynamic_cast<TH1*>(next()))) {
            hist->SetLineColor(hist->GetFillColor());
            hist->SetLineWidth(3);
          }

          TH1D * hsum_tmp = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
          TH1D *hsum = (TH1D*)hsum_tmp->Clone(Form("hsum_tmp%s",tag.Data()));

          if(hsum->GetMaximum() > holay->GetMaximum()) hstk->SetMaximum(hsum->GetMaximum()*1.2);
          else hstk->SetMaximum(holay->GetMaximum()*1.2);
          if(tag.Contains("Pi0Energy_COMPOSE")) hstk->SetMaximum(hsum->GetMaximum()*1.7);

          if(!tag.Contains("_XsEvt")) hstk->Draw("nostack HIST");
          else hstk->Draw("HIST");
          //hstk->Draw("nostack HIST");

          DrawOverlay(holay);
          if(tag.Contains("Mass")){
            double ymax = hsum->GetMaximum()*1.25;
            TLine *line = new TLine(0.134977,0,0.134977,ymax);
            line->SetLineColor(kGray+3);
            line->SetLineStyle(2);
            line->Draw("sames");
            // MC
            TF1 *f1_MC = new TF1("f1_MC","gaus",0.1,0.2);
            f1_MC->SetParameters(hsum_tmp->GetMaximum(), hsum_tmp->GetMean(), hsum_tmp->GetRMS() );
            // Data
            TF1 *f1_Data = new TF1("f1_Data","gaus",0.1,0.2);
            f1_Data->SetParameters(holay->GetMaximum(), holay->GetMean(), holay->GetRMS() );
            cout << "Tag: " << tag << endl;
            cout << "Fit MC mass: " << endl;
            hsum_tmp->Fit(f1_MC,"R0");
            cout << "Fit Data mass: " << endl;
            holay->Fit(f1_Data,"R0");


          }

          hsum->SetLineColor(9);
          hsum->SetLineWidth(3);
          hsum->SetFillStyle(0);
          hsum->Draw("same HIST");

          auto* legend = new TLegend(0.55, 0.55, 0.89, 0.88);

          const TList * ll = hstk->GetHists();
          // Cheat Legend
          TH1D * h1 = (TH1D*)ll->At(0);
          TH1D * h2 = (TH1D*)ll->At(1);
          TH1D * h3 = (TH1D*)ll->At(2);

          double tot = 0;
          for(auto ie : typeMaps[name]){
            tot += ie;
          }
          if(tag.Contains("RecPi0") && !tag.Contains("_EVT") && !tag.Contains("t0")){
            legend->AddEntry(h1, Form("#pi^{0} Signal (%.1f%s)",(typeMaps[name][0])/tot*100,"%"), "f");
            legend->AddEntry(h2, Form("Two #gamma Background (%.1f%s)",(typeMaps[name][1])/tot*100,"%"), "f");
            legend->AddEntry(h3, Form("One #gamma Background (%.2f%s)",(typeMaps[name][2])/tot*100,"%"), "f");
            //legend->AddEntry(holay, "Fake Data (1/2 MC sample)", "pl");
            legend->AddEntry(holay, "Data", "pl");
            legend->AddEntry(hsum, "MC", "f");
            legend->Draw("same");
          }
          else if(tag.Contains("_EVTXS")) {
            legend = new TLegend(0.15, 0.55, 0.49, 0.88);

            TH1D * h1 = (TH1D*)ll->At(0);
            TH1D * h2 = (TH1D*)ll->At(1);
            TH1D * h3 = (TH1D*)ll->At(2);
            TH1D * h4 = (TH1D*)ll->At(3);
            TH1D * h5 = (TH1D*)ll->At(4);
            TH1D * h6 = (TH1D*)ll->At(5);
            TH1D * h7 = (TH1D*)ll->At(6);
            TH1D * h8 = (TH1D*)ll->At(7);


            /*legend->AddEntry(h1, Form("Signal (Charge Exchange) (%.1f%s)",(typeMaps[name][0])/tot*100,"%"), "f");
            legend->AddEntry(h2, Form("Absorption Backgrounds (%.1f%s)",(typeMaps[name][1])/tot*100,"%"), "f");
            legend->AddEntry(h3, Form("Inelastic Backgrounds (%.2f%s)",(typeMaps[name][2])/tot*100,"%"), "f");
            legend->AddEntry(h4, Form("Sigle #pi^{0} Backgrounds (#pi^{#pm} > 0) (%.2f%s)",(typeMaps[name][3])/tot*100,"%"), "f");
            legend->AddEntry(h5, Form("Multi #pi^{0} Backgrounds (%.2f%s)",(typeMaps[name][4])/tot*100,"%"), "f");
            */
            legend->AddEntry(h1, Form("Beam #pi^{+} (%.1f%s)",(typeMaps[name][0])/tot*100,"%"), "f");
            legend->AddEntry(h2, Form("Beam #mu^{+} (%.1f%s)",(typeMaps[name][1])/tot*100,"%"), "f");
            legend->AddEntry(h3, Form("misID: p (%.1f%s)",(typeMaps[name][2])/tot*100,"%"), "f");
            legend->AddEntry(h4, Form("misID: #pi^{#pm} (%.1f%s)",(typeMaps[name][3])/tot*100,"%"), "f");
            legend->AddEntry(h5, Form("misID: #mu^{+} (%.1f%s)",(typeMaps[name][4])/tot*100,"%"), "f");
            legend->AddEntry(h6, Form("misID: e/#gamma (%.1f%s)",(typeMaps[name][5])/tot*100,"%"), "f");
            legend->AddEntry(h7, Form("Comics (%.1f%s)",(typeMaps[name][6])/tot*100,"%"), "f");
            legend->AddEntry(h8, Form("others (%.1f%s)",(typeMaps[name][7])/tot*100,"%"), "f");

            legend->AddEntry(hsum, "MC", "f");
            legend->Draw("same");
          }
          else{

            if(tag.Contains("i07")) legend = new TLegend(0.15, 0.55, 0.49, 0.88);
            if(tag.Contains("i08")) legend = new TLegend(0.55, 0.55, 0.89, 0.88);
            if(tag.Contains("i09")) legend = new TLegend(0.15, 0.55, 0.49, 0.88);
            if(tag.Contains("i100")) legend = new TLegend(0.15, 0.55, 0.49, 0.88);
            if(tag.Contains("i40")) legend = new TLegend(0.15, 0.55, 0.49, 0.88);

            TH1D * h1 = (TH1D*)ll->At(0);
            TH1D * h2 = (TH1D*)ll->At(1);
            TH1D * h3 = (TH1D*)ll->At(2);
            TH1D * h4 = (TH1D*)ll->At(3);
            TH1D * h5 = (TH1D*)ll->At(4);
            TH1D * h6 = (TH1D*)ll->At(5);
            //TH1D * h7 = (TH1D*)ll->At(6);

            /*legend->AddEntry(h1, Form("Signal (Charge Exchange) (%.1f%s)",(typeMaps[name][0])/tot*100,"%"), "f");
            legend->AddEntry(h2, Form("Absorption Backgrounds (%.1f%s)",(typeMaps[name][1])/tot*100,"%"), "f");
            legend->AddEntry(h3, Form("Inelastic Backgrounds (%.2f%s)",(typeMaps[name][2])/tot*100,"%"), "f");
            legend->AddEntry(h4, Form("Sigle #pi^{0} Backgrounds (#pi^{#pm} > 0) (%.2f%s)",(typeMaps[name][3])/tot*100,"%"), "f");
            legend->AddEntry(h5, Form("Multi #pi^{0} Backgrounds (%.2f%s)",(typeMaps[name][4])/tot*100,"%"), "f");
            */
            legend->AddEntry(h1, Form("Signal (Charge Exchange) (%.1f%s)",(typeMaps[name][0])/tot*100,"%"), "f");
            legend->AddEntry(h2, Form("Absorption Background (%.1f%s)",(typeMaps[name][1])/tot*100,"%"), "f");
            legend->AddEntry(h3, Form("Pion Production Background (#pi^{0} = 0) (%.2f%s)",(typeMaps[name][2])/tot*100,"%"), "f");
            legend->AddEntry(h4, Form("Pion Production Background (#pi^{0} = 1) (%.2f%s)",(typeMaps[name][3])/tot*100,"%"), "f");
            legend->AddEntry(h5, Form("Pion Production Background (#pi^{0} > 1) (%.2f%s)",(typeMaps[name][4])/tot*100,"%"), "f");
            //legend->AddEntry(h6, Form("Other Background (%.2f%s)",(typeMaps[name][5])/tot*100,"%"), "f");
            legend->AddEntry(h6, Form("Beam Background (%.2f%s)",(typeMaps[name][5])/tot*100,"%"), "f");

            //legend->AddEntry(holay, "Fake Data (1/2 MC Sample)", "pl");
            legend->AddEntry(holay, "Data", "pl");
            //legend->AddEntry(h6, Form("Other Backgrounds (%.2f%s)",(typeMaps[name][5])/tot*100,"%"), "f");
            //legend->AddEntry(h7, Form("Beam Backgrounds (%.2f%s)",(typeMaps[name][6])/tot*100,"%"), "f");

            legend->AddEntry(hsum, "MC", "f");
            legend->Draw("same");
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

          const int mrks[]={1,1,1,1, 1,1,1,1,1,1, 1,8};
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
          if(tag.Contains("BeamQualityTheta") || tag.Contains("BeamQualityZ") || tag.Contains("BeamQualityXY")) pad1->SetLogy();
          pad1->SetBottomMargin(0.02);
          //  pad1->SetGridx();
          //  pad1->SetGridy();
          pad1->Draw();
          pad1->cd();

          // Get the correct sacle factor for each beam cut (turn on when needed)
          /*TH1D * hMC = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
          if(tag.Contains("EndZ_Cut")){
            cout << "MC entry: " << hMC->Integral() << endl;
            cout << "Data entry: " << holay->Integral() << endl;
            cout << "name: " << tag << "scale factor: " << holay->Integral()/hMC->Integral() << endl;
          }*/

          // Scale the beam cuts related hists data to MC no weight
          if(tag.Contains("CutPDG") || tag.Contains("BeamPDG")) ScaleStack(hstk, 0.702223); // CutPDG
          else if(tag.Contains("CutPandora") || tag.Contains("BeamPandora")) ScaleStack(hstk, 0.709511); // CutPandora
          else if(tag.Contains("CutCaloSize") || tag.Contains("BeamCalo")) ScaleStack(hstk, 0.703488); // CutCaloSize
          else if(tag.Contains("CutBeamQuality") || tag.Contains("BeamQuality")) ScaleStack(hstk, 0.698653); // CutBeamQuality
          else if(tag.Contains("CutAPA3") || tag.Contains("BeamAPA3")) ScaleStack(hstk, 0.658395); // CutAPA3
          else if(tag.Contains("CutMichelScore") || tag.Contains("BeamMichel")) ScaleStack(hstk, 0.65585);  // CutMichelScore
          else if(tag.Contains("CutChi2DOF") || tag.Contains("BeamChi2")) ScaleStack(hstk, 0.649126); // CutChi2DOF
          else if(tag.Contains("CutBeamScraper") || tag.Contains("BeamScraper")) ScaleStack(hstk, 0.57804); // CutBeamScraper(final beam scale)

          else if(tag.Contains("hRecPiPlus") && !tag.Contains("b00")) ScaleStack(hstk, plotscale_wt);
          else if(tag.Contains("i087hRecPiZeroKineticEnergyEvt_COMPOSE_XsEvt")) ScaleStack(hstk, plotscale_wt);

          // Scale data to MC
          else ScaleStack(hstk, plotscale);

          hstk->GetYaxis()->SetTitle(holay->GetYaxis()->GetTitle());
          SetTitleFormat(hstk);

          TList* histList = hstk->GetHists(); // Get the list of histograms in the stack
          TIter next(histList);
          TH1* hist;
          while ((hist = dynamic_cast<TH1*>(next()))) {
            hist->SetFillColorAlpha(hist->GetFillColor(),1.0);
          }

          TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
          hstk->SetMaximum(holay->GetMaximum()*1.2);
          if(tag.Contains("BeamQualityZ")) hstk->SetMaximum(holay->GetMaximum()*400);
          if(tag.Contains("Recdalphat_STK")) hstk->SetMaximum(holay->GetMaximum()*2);
          hstk->Draw("hist");
          DrawOverlay(holay);
          //c1->Update();

          const int overlayColor = kBlack;

          // Event type legend
          if(hstk->GetNhists() <= 3 && !tag.Contains("Channel")){
            vector<TString> evtType = FillLegendType("evtType",hstk->GetName());
            vector<TString> htype = FillLegendStyle(2,"evtType");

            int *cols=GetColorArray(evtType.size());
            cols[3]=overlayColor;
            const int mrks[]={1,1,1,8};
            TLegend * lg = 0x0;
            lg = DrawLegend(evtType, htype, tag, cols, mrks);
            lg->Draw("same");
          }

          // Beam particle type
          else if((hstk->GetNhists() == 8 && !tag.Contains("Channel")) || tag.Contains("i05") || tag.Contains("i06")){

            vector<TString> beamType = FillLegendType("beamType",hstk->GetName());
            vector<TString> htype = FillLegendStyle(2,"beamType");

            int *cols=GetColorArray(beamType.size());
            const int mrks[]={1,1,1,1,1,1,1,1,8};
            cols[beamType.size()-1]=overlayColor;

            TLegend * lg = 0x0;
            lg = DrawLegend(beamType, htype, tag, cols, mrks, 2);
            lg->Draw("same");
          }

          // Channel type
          else if(tag.Contains("Channel")){
            vector<TString> beamType = FillLegendType("channelType",hstk->GetName());
            vector<TString> htype = FillLegendStyle(2,"channelType");

            int *cols=GetColorArray(beamType.size());
            const int mrks[]={1,1,1,1,1,1,1,8};
            cols[beamType.size()-1]=overlayColor;

            TLegend * lg = 0x0;
            lg = DrawLegend(beamType, htype, tag, cols, mrks, 2);
            lg->Draw("same");
          }

          // TKI particle type
          else if(hstk->GetNhists() == 4){
            vector<TString> beamType = FillLegendType("TKIType",hstk->GetName());
            vector<TString> htype = FillLegendStyle(2,"TKIType");

            int *cols=GetColorArray(beamType.size());
            const int mrks[]={1,1,1,1,1,8};
            cols[beamType.size()-1]=overlayColor;

            TLegend * lg = 0x0;
            lg = DrawLegend(beamType, htype, tag, cols, mrks, 2);
            lg->Draw("same");
          }
          // Particle type legend
          else{
            vector<TString> parType = FillLegendType("parType",hstk->GetName());

            vector<TString> htype = FillLegendStyle(2,"parType");//need to matcy parType

            const int mrks[]={1,1,1,1, 1,1,1,1,1,1, 1,8};
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
      else if(tag.Contains("OVERLAY") && !tag.Contains("CEX")){

        TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
        hstk->SetMaximum(hsum->GetMaximum()*1.2);

        hstk->GetXaxis()->SetTitle("Truth Energy (GeV)");
        hstk->GetYaxis()->SetTitle("Candidates");
        SetTitleFormat(hstk,false);

        TList* histList = hstk->GetHists(); // Get the list of histograms in the stack
        TIter next(histList);
        TH1* hist;
        while ((hist = dynamic_cast<TH1*>(next()))) {
          hist->SetLineColor(hist->GetFillColor());
          hist->SetLineWidth(3);
        }

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
        const TList * ll = hstk->GetHists();

        TH1D * E1 = (TH1D*)ll->At(0);
        TH1D * E2 = (TH1D*)ll->At(1);

        legend->AddEntry(E1, "Leading photon", "f");
        legend->AddEntry(E2, "SubLeading photon", "f");
        legend->AddEntry(hsum, "Photon Spectrum", "f");
        legend->AddEntry(holay_Pi0, "#pi^{0} Spectrum", "f");
        legend->Draw("same");

      }

      else if(tag.Contains("OVERLAY") && tag.Contains("CEX")){

        TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
        hstk->SetMaximum(hsum->GetMaximum()*1.2);

        hstk->GetXaxis()->SetTitle("Number of Particles");
        hstk->GetYaxis()->SetTitle("Candidates");
        SetTitleFormat(hstk,false);

        TList* histList = hstk->GetHists(); // Get the list of histograms in the stack
        TIter next(histList);
        TH1* hist;
        while ((hist = dynamic_cast<TH1*>(next()))) {
          hist->SetFillStyle(0); // Set the fill style for each histogram
        }

        hstk->Draw("nostack L");
        //hstk->Draw("nostack HIST");

        auto* legend = new TLegend(0.5, 0.55, 0.85, 0.88);
        // Cheat Legend
        const TList * ll = hstk->GetHists();

        TH1D * h1 = (TH1D*)ll->At(0);
        TH1D * h2 = (TH1D*)ll->At(1);
        TH1D * h3 = (TH1D*)ll->At(2);

        h1->SetLineWidth(3);
        h2->SetLineWidth(3);
        h3->SetLineWidth(3);

        TString lheader("T_{#pi^{0}} > 650 MeV");
        //legend->SetHeader("#splitline{Fake Data (Full MC sample)}{T_{beam #pi^{+}} = 775 MeV}");
        if(tag.Contains("Middle"))  lheader = ("200 < T_{#pi^{0}} < 400 MeV");
        if(tag.Contains("Low")) lheader = ("T_{#pi^{0}} < 150 MeV");

        legend->SetHeader(lheader);
        legend->AddEntry(h1, "Protons", "l");
        legend->AddEntry(h2, "Neutrons", "l");
        legend->AddEntry(h3, "Nucleus", "l");
        legend->Draw("same");

      }
    }
    else cout << "PlotUtils::DrawHist not found correct histogram!" << " name: " << tag << endl;

    tt.SetTextSize(0.035);
    if(tag.Contains("i001") && !tag.Contains("_")) tt.DrawLatex(0.165,0.925,"DUNE:ProtoDUNE-SP Simulation");
    else tt.DrawLatex(0.125,0.925,"DUNE:ProtoDUNE-SP Simulation");
    if(tag.Contains("i10")) tt.DrawLatex(0.565,0.925,"True Beam T_{#pi^{+}} = 650 - 800 MeV");
    else tt.DrawLatex(0.665,0.925,"2GeV/c Pion Beam");

    c1->Print(outdir+"/"+tag+".pdf");

  } // End of for loop
}

THStack * PlotUtils::ConvertToStack(const TH2D * hh, const bool kMC, std::map<TString, vector<double>> &typeMaps)
{
  TString tag = hh->GetName();
  //cout << "tag: " << tag << endl;
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

    int icol = -999;
    if((tag.Contains("COMPOSE") || tag.Contains("OVERLAY")) && !tag.Contains("LOG")) icol = GetColor(col[iy-y0],true);
    else icol = GetColor(col[iy-y0],false);//need constant map between y and color

    htmp->SetFillColor(icol);
    if((tag.Contains("COMPOSE") || tag.Contains("OVERLAY")) && !tag.Contains("LOG") && icol < 1000) htmp->SetFillColorAlpha(icol, 0.35);

    htmp->SetLineColor(kBlack);
    htmp->SetLineWidth(1);
    htmp->SetMarkerSize(2);
    if(tag.Contains("OVERLAY") && !tag.Contains("CEX")){
      htmp->SetLineColor(icol);
      htmp->SetFillStyle(1001);
    }
    if(tag.Contains("OVERLAY") && tag.Contains("CEX")){
      htmp->SetLineColor(icol);
      htmp->SetFillStyle(1001);
      htmp->SetLineWidth(2);
      //htmp->SetFillColorAlpha(icol,0.3);

    }
    if(tag.Contains("COMPOSE")){
      htmp->SetLineColor(icol);
      htmp->SetFillStyle(1001);
      if(tag.Contains("LOG")){
        htmp->SetFillStyle(0);
        htmp->SetLineWidth(2);
      }
    }

    //printf("PlotUtils::ConvertToStack %s adding y %f with color %d\n", tag.Data(), hh->GetYaxis()->GetBinCenter(iy), icol);
    stk->Add(htmp);
  }
  //cout << "newintegral: " << newintegral << endl;
  //for(double in : TruthTypeVect){
    //cout << "Type percent: " << in/newintegral*100 << "%" << endl;
  //}
  int oldintegralintPart = (int) oldintegral;
  int newintegralintPart = (int) newintegral;

  if(oldintegralintPart!=newintegralintPart){
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

      if(tag.Contains("i054hRecPiPlusInitialEnergy_STK") && ii == 0){
        cout << "i054hRecPiPlusInitialEnergy_STK" << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio: " << htmp->GetBinContent(ix) << endl;
          getWeight_hRecPiPlusInitialEnergy->Fill(htmp->GetBinContent(ix));
        }
      }

      if(tag.Contains("i054hRecPiPlusInitialEnergyPosZCut_STK") && ii == 0){
        cout << "i054hRecPiPlusInitialEnergyPosZCut_STK" << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio: " << htmp->GetBinContent(ix) << endl;
        }
      }

      if(tag.Contains("i050hRecPiPlusInteractingEnergy_STK") && ii == 0){
        cout << "i050hRecPiPlusInteractingEnergy_STK" << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio ShowME : " << htmp->GetBinContent(ix) << endl;
          getWeight_hRecPiPlusInteractingEnergy->Fill(htmp->GetBinContent(ix));
        }
      }

      if(tag.Contains("i076hRecPiPlusInteractingEnergyEvt_COMPOSE") && ii == 0){
        cout << "i076hRecPiPlusInteractingEnergyEvt_COMPOSE" << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio : " << htmp->GetBinContent(ix) << endl;
          getWeight_hRecPiPlusInteractingEnergyEvt->Fill(htmp->GetBinContent(ix));
        }
      }

      if(tag.Contains("i300hRecPiZeroRangeKineticEnergyEvtNoWeight_COMPOSE_XsEvt") && ii == 0){
        cout << "getWeight_hRecPiZeroSliceKineticEnergyEvt" << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio ShowME : " << htmp->GetBinContent(ix) << endl;
          getWeight_hRecPiZeroSliceKineticEnergyEvt->Fill(htmp->GetBinContent(ix));
        }
      }
      // ============ Pi0 KE =========== //
      if(tag.Contains("i30") && ii == 0){
        cout << tag << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio Herilala: " << htmp->GetBinContent(ix) << endl;
        }
      }
      // ============ Pi0 Costheta =========== //
      if(tag.Contains("i40") && ii == 0){
        cout << tag << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio: " << htmp->GetBinContent(ix) << endl;
        }
      }
      // ============ Pi0 Theta =========== //
      if(tag.Contains("i50") && ii == 0){
        cout << tag << endl;
        const Int_t x0 = htmp->GetXaxis()->GetFirst();
        const Int_t x1 = htmp->GetXaxis()->GetLast();
        for(Int_t ix=x0; ix<=x1; ix++){
          cout << "ix: " << ix << "ratio: " << htmp->GetBinContent(ix) << endl;
        }
      }
    }

    delete hsum;
  }
  return hout;
}

int PlotUtils::GetColor(const int col, const bool alpha)
{
  if(col>=1000 && !alpha){
    return 1500+col-1000;
  }
  else if(col>=1000 && alpha){
    return 1700+col-1000;
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

  const int fgkColorBase_solid = 1700;
  //http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
  Int_t ids=fgkColorBase_solid+1;
  new TColor(ids++, 0./255., 0./255., 0./255., "CB1_Black",0.35);
  new TColor(ids++, 0./255., 73./255., 73./255., "CB2_Forest",0.35);
  new TColor(ids++, 0./255., 146./255., 146./255., "CB3_Teal",0.35);
  new TColor(ids++, 255./255., 109./255., 182./255., "CB4_HotPink",0.35);
  new TColor(ids++, 255./255., 182./255., 119./255., "CB5_BabyPink",0.35);
  new TColor(ids++, 73./255., 0./255., 146./255., "CB6_Purple",0.35);
  new TColor(ids++, 0./255., 109./255., 219./255., "CB7_RoyalBlue",0.35);
  new TColor(ids++, 182./255., 109./255., 255./255., "CB8_Lilac",0.35);
  new TColor(ids++, 109./255., 182./255., 255./255., "CB9_BlueGrey",0.35);
  new TColor(ids++, 182./255., 219./255., 255./255., "CB10_SpaceWolves",0.35);
  new TColor(ids++, 146./255., 0./255., 0./255., "CB11_Maroon",0.35);
  new TColor(ids++, 146./255., 73./255., 0./255., "CB12_Tan",0.35);
  new TColor(ids++, 219./255., 209./255., 0./255., "CB13_Orange",0.35);
  new TColor(ids++, 36./255., 255./255., 36./255., "CB14_DayGleen",0.35);
  new TColor(ids++, 255./255., 255./255., 109./255., "CB15_SunFlower",0.35);


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
  cctmp->Print("output/hprof_withP.pdf");
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
    cc->Print("output/hSliceXDrawY.pdf");
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

      // Beam
      if(name.Contains("BeamTheta")){
        AnaIO::hMeanBeamTheta->SetMinimum(-60);
        AnaIO::hMeanBeamTheta->SetMaximum(60);
        AnaIO::hMeanBeamTheta->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanBeamTheta->SetBinError(ii,fitfunc->GetParError(0));

        cout << "Beam Theta bin : " << AnaIO::hMeanBeamTheta->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "Beam Theta: " << fitfunc->GetParameter(0) << endl;
      }
      if(name.Contains("BeamPhi")){
        AnaIO::hMeanBeamPhi->SetMinimum(-60);
        AnaIO::hMeanBeamPhi->SetMaximum(60);
        AnaIO::hMeanBeamPhi->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanBeamPhi->SetBinError(ii,fitfunc->GetParError(0));
        cout << "Beam Phi bin : "  << AnaIO::hMeanBeamPhi->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "Beam Phi: " << fitfunc->GetParameter(0) << endl;
      }

      if(name.Contains("BeamPAll")){
        AnaIO::hMeanBeamPMom->SetMinimum(-0.2);
        AnaIO::hMeanBeamPMom->SetMaximum(0.2);
        AnaIO::hMeanBeamPMom->SetBinContent(ii,fitfunc->GetParameter(0));
        AnaIO::hMeanBeamPMom->SetBinError(ii,fitfunc->GetParError(0));
        cout << "Beam p all bin : " << AnaIO::hMeanPMom->GetXaxis()->GetBinCenter(ii) << endl;
        cout << "Beam p all: " << fitfunc->GetParameter(0) << endl;
      }


    }
    //hpj->SetStats(0);
    hpj->Draw("E");
  }

  cc->Print("output/"+name+"_sliced.pdf");
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
  //lg = new TLegend(0.6, 0.6, 0.88, 0.88);
  //lg = new TLegend(0.4, 0.4, 0.88, 0.88);
  //lg = new TLegend(0.65, 0.68, 0.88, 0.88);
  lg = new TLegend(0.55, 0.55, 0.89, 0.88);

  if(tag.Contains("COMPOSE")){
    lg = new TLegend(0.4, 0.55, 0.88, 0.88);
    if(tag.Contains("LOG")) lg = new TLegend(0.4, 0.61, 0.88, 0.88);
  }
  if(tag.Contains("Channel") || tag.Contains("BeamEndZ")){
    lg = new TLegend(0.5, 0.6, 0.88, 0.88);
  }
  if(tag.Contains("BeamQualityTheta")){
    lg = new TLegend(0.2, 0.6, 0.58, 0.88);
  }
  if(tag.Contains("i05")){
    // small size
    //lg = new TLegend(0.13, 0.6, 0.4, 0.88);
    lg = new TLegend(0.13, 0.55, 0.47, 0.88);
  }
  for(int ii=0; ii<nent; ii++){
    TH1D * hh=new TH1D(Form("h%d%s",ii,tag.Data()),"",1,0,1);
    int col = -999;

    if((tag.Contains("COMPOSE") || tag.Contains("OVERLAY")) && !tag.Contains("LOG")) col = GetColor(cols[ii],true);
    else col = GetColor(cols[ii],false);

    hh->SetFillColor(col);
    if((tag.Contains("COMPOSE") || tag.Contains("OVERLAY")) && !tag.Contains("LOG") && col < 1000) hh->SetFillColorAlpha(col, 0.35);
    if(tag.Contains("COMPOSE") && tag.Contains("LOG")){
      hh->SetLineWidth(6);
    }
    hh->SetLineColor(col);
    hh->SetMarkerStyle(mkrs[ii]);
    hh->SetMarkerSize(1);
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

void PlotUtils::DrawDataMCRatio(TH1D * hratio, bool xsec){

  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetXaxis()->SetLabelSize(0.15);
  hratio->GetXaxis()->SetTitleSize(0.15);
  hratio->GetXaxis()->SetTitleOffset(1.);
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetYaxis()->SetTitleOffset(.3);
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->GetYaxis()->SetRangeUser(0,2.1);
  if(xsec) hratio->GetYaxis()->SetRangeUser(0.5,1.6);
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
    tmpType.push_back(Form("Beam #pi^{+} (%.1f%s)",(typeMaps[name][0])/tot*100,"%"));
    tmpType.push_back(Form("Beam #mu^{+} (%.1f%s)",(typeMaps[name][1])/tot*100,"%"));
    tmpType.push_back(Form("misID: p (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back(Form("misID: #pi^{#pm} (%.1f%s)",(typeMaps[name][3])/tot*100,"%"));
    tmpType.push_back(Form("misID: #mu^{+} (%.1f%s)",(typeMaps[name][4])/tot*100,"%"));
    tmpType.push_back(Form("misID: e/#gamma (%.1f%s)",(typeMaps[name][5])/tot*100,"%"));
    tmpType.push_back(Form("Comisc (%.1f%s)",(typeMaps[name][6])/tot*100,"%"));
    tmpType.push_back(Form("others (%.1f%s)",(typeMaps[name][7])/tot*100,"%"));
    tmpType.push_back("data");
  }

  if(tag.Contains("TKIType")){
    tmpType.push_back(Form("1p0n (%.1f%s)",(typeMaps[name][0])/tot*100,"%"));
    tmpType.push_back(Form("1pMn (%.1f%s)",(typeMaps[name][1])/tot*100,"%"));
    //tmpType.push_back(Form("Np0n (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back(Form("NpMn (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back(Form("background (%.1f%s)",(typeMaps[name][3])/tot*100,"%"));
    tmpType.push_back(Form("non-#pi^{+} beam (%.1f%s)",(typeMaps[name][4])/tot*100,"%"));
    tmpType.push_back("data");
  }
  if(tag.Contains("channelType")){
    tmpType.push_back(Form("Charge Exchange (%.1f%s)",(typeMaps[name][0])/tot*100,"%"));
    tmpType.push_back(Form("Inelastic (%.1f%s)",(typeMaps[name][1])/tot*100,"%"));
    tmpType.push_back(Form("Absorption (%.1f%s)",(typeMaps[name][2])/tot*100,"%"));
    tmpType.push_back(Form("PionDecays (%.1f%s)",(typeMaps[name][3])/tot*100,"%"));
    tmpType.push_back(Form("OtherChannels (%.1f%s)",(typeMaps[name][4])/tot*100,"%"));
    tmpType.push_back(Form("BeamMuons (%.1f%s)",(typeMaps[name][5])/tot*100,"%"));
    tmpType.push_back(Form("Background (%.1f%s)",(typeMaps[name][6])/tot*100,"%"));
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

  if(tag.Contains("TKIType") && opt == 2){
    tmpType.push_back("f");
    tmpType.push_back("f");
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
    tmpType.push_back("f");
    tmpType.push_back("f");
    tmpType.push_back("ple");
  }
  if(tag.Contains("channelType") && opt == 2){
    tmpType.push_back("f");
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

void PlotUtils::PrintShowerPurityandEff(const TString tag, TH1D * h_ldGamma, TH1D * h_slGamma, TH2D * h2d)
{
  TString name = h2d->GetName();
  //cout << "h2d name: " << name << endl;
  std::map<TString,vector<double>> typeMaps;
  ConvertToStack(h2d,true,typeMaps);
  double tot = 0;
  for(auto ie : typeMaps[name]){
    tot += ie;
  }
  double deno = h_ldGamma->Integral(0,10000);// + h_slGamma->Integral(0,10000);
  if(name.Contains("Pi0")) deno = h_ldGamma->Integral(0,10000);

  //cout << "stk name: " << stk->GetName() << endl;
  //cout << "tot: " << tot << "deno: " << deno << endl;
  //cout << "num: " << (typeMaps[name][3]) << endl;
  //cout << "purity: " << (typeMaps[name][3])/tot << " eff: " << (typeMaps[name][3])/deno << endl;

  const double purity = (typeMaps[name][3])/tot;
  const double eff = (typeMaps[name][3])/deno;

  printf("%-50s: purity %5.3f efficiency %5.3f pur*eff %.3f%% \n", tag.Data(), purity, eff, purity*eff);

}

void PlotUtils::PrintPi0PurityandEff(const TString tag, TH1D * h_ldGamma, TH2D * h2d)
{
  TString name = h2d->GetName();
  //cout << "h2d name: " << name << endl;
  std::map<TString,vector<double>> typeMaps;
  ConvertToStack(h2d,true,typeMaps);
  double tot = 0;
  for(auto ie : typeMaps[name]){
    tot += ie;
  }
  double deno = h_ldGamma->Integral(0,10000);

  //cout << "stk name: " << stk->GetName() << endl;
  //cout << "tot: " << tot << "deno: " << deno << endl;
  //cout << "num: " << (typeMaps[name][0]) << endl;
  //cout << "purity: " << (typeMaps[name][0])/tot << " eff: " << (typeMaps[name][0])/deno << endl;

 double purity = (typeMaps[name][0])/tot;
 double eff = (typeMaps[name][0])/deno;

  if(name.Contains("hPiPlusMom")){
   purity = (typeMaps[name][1])/tot;
   eff = (typeMaps[name][1])/deno;
  }

  if(name.Contains("hShowerEnergy")){
   purity = (typeMaps[name][3])/tot;
   eff = (typeMaps[name][3])/deno;
  }

  printf("%-50s: purity %5.3f efficiency %5.3f pur*eff %.3f%% \n", tag.Data(), purity, eff, purity*eff);

}


void PlotUtils::PrintExcPurityandEff(const TString tag, TH1D * h_1, TH1D * h_2, TH1D * h_3, TH1D * h_4, TH2D * h2d)
{
  TString name = h2d->GetName();
  //cout << "h2d name: " << name << endl;
  std::map<TString,vector<double>> typeMaps;
  ConvertToStack(h2d,true,typeMaps);
  double tot = 0;
  for(auto ie : typeMaps[name]){
    tot += ie;
  }
  double deno = h_1->Integral(0,10000) + h_2->Integral(0,10000) + h_3->Integral(0,10000) + h_4->Integral(0,10000);

  //cout << "stk name: " << stk->GetName() << endl;
  //cout << "tot: " << tot << "deno: " << deno << endl;
  //cout << "num: " << (typeMaps[name][0]) << endl;
  //cout << "purity: " << (typeMaps[name][0])/tot << " eff: " << (typeMaps[name][0])/deno << endl;

  const double purity = (typeMaps[name][0]+typeMaps[name][1]+typeMaps[name][2])/tot;
  const double eff = (typeMaps[name][0]+typeMaps[name][1]+typeMaps[name][2])/deno;

  printf("%-50s: purity %5.2f efficiency %5.2f pur*eff %.3f%% \n", tag.Data(), purity, eff, purity*eff);

}

void PlotUtils::PrintXSUnfoldingInfo(const double& i_true, const double& i_reco, const double& i_reco_sel, const double& i_reco_nosel, const double& i_reco_bad){

  cout << "truth entry: " << i_true << endl;
  cout << "reco entry: " << i_reco << endl;
  cout << "------------------------" << endl;
  cout << "reco selected: " << i_reco_sel << " reco not-selected: " << i_reco_nosel << " reco bad: " << i_reco_bad << endl;
  cout << "reco sum entry: " << i_reco_sel + i_reco_nosel + i_reco_bad << endl;
  cout << endl;

}

void PlotUtils::TotalCEXXSCal(TH1 * hh, TH1D * InteractingHist, TH1D * xsec, const bool & Eslice, const bool & widerBin, const bool & newMethod)
{

  // Scale down it to zero for later calculation
  xsec->Scale(0);
  double avogadro_constant = 6.02214076e23;  // 1 / mol
  double argon_molar_mass = 39.95;           // g / mol
  double liquid_argon_density = 1.39;        // g / cm^3
  double fiducial_thickness = 0.479;//slice_width;//0.479;         // cm wire spacing
  if(Eslice) fiducial_thickness = slice_width;
  if(Eslice && newMethod) fiducial_thickness = slice_width*50.0;
  if(Eslice && widerBin) fiducial_thickness = slice_width*50.0;
  //double fiducial_thickness = 20;         // cm wire spacing
  double sigma_factor = argon_molar_mass / (avogadro_constant * liquid_argon_density * fiducial_thickness); // pre-factor
  double meandEdx[] = {5.06921, 2.78604, 2.35879, 2.20657, 2.14137, 2.11321, 2.10322, 2.10293, 2.10805, 2.11628, 2.12627, 2.13724, 2.14869, 2.16033, 2.17195,
                       2.18344, 2.19472, 2.20573, 2.21646, 2.2269};
  // Define start and end bins
  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();

  // Testing slice thickness effect
  //cout << "reco slice * Ninc: " << fiducial_thickness*hh->Integral("width") << endl;
  //cout << "reco Nint: " << InteractingHist->Integral("width") << endl;

  // Loop over each bin and combine the incident and interacting histograms to get xsec
  for(Int_t ix=x0; ix<=x1; ix++){

    // Read off the entry for each bin of the two histograms
    double incE = hh->GetBinContent(ix);
    double intE = InteractingHist->GetBinContent(ix);
    // Calculate the ratio in the full log form
    double ratio = 0;
    if(incE > intE) ratio = log(incE/(incE-intE));//*meandEdx[ix-1]*sigma_factor*1e27;
    if(Eslice) ratio = ratio * meandEdx[ix-1]*sigma_factor*1e27;

    if(incE != 0) cout << "ix: "<< ix << "incE: " << incE << "intE: " << intE << " ratio: " << ratio << "meandEdx[ix]: " << meandEdx[ix-1] << endl;

    // Simple ratio form
    //double ratio = intE/incE;

    // If the incE entry is not zero set the bin content
    if(incE != 0) xsec->SetBinContent(ix,ratio);

    // Error propagation method 1
    //double error = sqrt(intE+pow(intE,2)/incE)/incE;
    // Error propagation method 2
    double einc = hh->GetBinError(ix);
    double eint = InteractingHist->GetBinError(ix);
    double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)));
    // If the ratio is not zero set the error
    if(ratio != 0 ) xsec->SetBinError(ix,error);

  }
  // The xsec histogram entry is now set
  if(!Eslice) xsec->Scale(sigma_factor*1e27);
  xsec->SetMaximum(300);
  TString tag = InteractingHist->GetName();
  if(!tag.Contains("CEX") && !tag.Contains("UnFold")) xsec->SetMaximum(1600);
  xsec->SetMinimum(0);
   /*for(Int_t ix=x0; ix<=x1; ix++){
    cout << "ix: "<< ix << "bin content: " << xsec->GetBinContent(ix) << endl;
   }*/
}


void PlotUtils::DiffCEXXSCal(TH1D * DiffCEXInteractingHist, TH1D * DiffCEXxsec, const double diffInt,  const double diffInterror)
{

  DiffCEXxsec->Scale(0);

  const Int_t x0_diff = DiffCEXInteractingHist->GetXaxis()->GetFirst();
  const Int_t x1_diff = DiffCEXInteractingHist->GetXaxis()->GetLast();

  for(Int_t ix=x0_diff; ix<=x1_diff; ix++){
    double incE = diffInt;
    double intE = DiffCEXInteractingHist->GetBinContent(ix);
    if(intE != 0 ){
      //double ratio = log(incE/(incE-intE));
      double ratio = intE/incE;
      DiffCEXxsec->SetBinContent(ix,ratio);
      //double error = sqrt(intE+pow(intE,2)/incE)/incE;
      double einc = diffInterror;
      double eint = DiffCEXInteractingHist->GetBinError(ix);
      double error = sqrt(ratio*ratio*(pow(einc/incE,2)+pow(eint/intE,2)));
      DiffCEXxsec->SetBinError(ix,error);
      //cout << "error: " << error << endl;
      //cout << "ratio: " << ratio << endl;
    }
  }
}

TH1D * PlotUtils::GetIncidentHist(TH1D * InitialHist, TH1D * InteractingHist)
{
  TH1D *IncidentHist = (TH1D*)InitialHist->Clone();
  IncidentHist->Scale(0);

  //double IniSum = 0, IntSum = 0;
  const Int_t x0 = InteractingHist->GetXaxis()->GetFirst();
  const Int_t x1 = InteractingHist->GetXaxis()->GetLast();
/*  for(Int_t ix=x1; ix>=x0; ix--){
    IniSum += InitialHist->GetBinContent(ix+1);
    //double error = InitialHist->GetBinError(ix);
    IntSum += InteractingHist->GetBinContent(ix+1);
    IncidentHist->SetBinContent(ix, IniSum-IntSum);
  }
*/
  cout << "Name: " << InitialHist->GetName() << endl;
  cout << "InitialHist: " << InitialHist->Integral(0,10000) << endl;
  cout << "InteractingHist: " << InteractingHist->Integral(0,10000) << endl;

  for(Int_t ix=x1; ix>=x0; ix--){
    double iniN = 0;
    double intN = 0;
    for(Int_t jx=x1; jx>=ix; jx--){
      iniN += InitialHist->GetBinContent(jx);
    }
    for(Int_t jx=x1; jx>=ix+1; jx--){
      intN += InteractingHist->GetBinContent(jx);
    }
    /*for(Int_t jx=ix; jx>=x0; jx--){
      intN += InteractingHist->GetBinContent(jx);
    }
    for(Int_t jx=ix-1; jx>=x0; jx--){
      iniN += InitialHist->GetBinContent(jx);
    }*/
    //cout << "iniN: " << iniN << endl;
    //cout << "intN: " << intN << endl;
    //cout << "iniN - intN: " << iniN - intN << endl;

    //cout << "intN: " << intN << endl;
    //cout << "iniN: " << iniN << endl;
    //cout << "intN - iniN: " << intN - iniN << endl;

    //double entry = intN - iniN;
    double entry = iniN - intN;
    //cout << "ix: " << ix << "entry: " << entry << endl;
    if (entry < 0) entry = 0;
    IncidentHist->SetBinContent(ix,entry);
    IncidentHist->SetBinError(ix,InitialHist->GetBinError(ix));

    //hh->SetBinError(ix,error);
  }

  return IncidentHist;
}
