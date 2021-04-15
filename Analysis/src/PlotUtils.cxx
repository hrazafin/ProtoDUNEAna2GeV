#include "../include/PlotUtils.h"

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
      // check STK tag
      if(tag.Contains("STK")){
        THStack * stk = ConvertToStack(htmp);
        if(stk){
          if(kMC) lout->Add(stk);
          else{
            TH1D * hsum = GetStackedSum(stk);
            lout->Add(hsum);
          }
          //THStack * snor = NormalizeStack(stk);
          //lout->Add(snor);
        }
      }
    }
  } // End of for loop
}
void PlotUtils::DrawHist(TList *lout, const double plotscale, TList * overlayList, const TString outdir)
{
  // Create new canvas for drawing histogram
  TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
  // Setup Pad
  PadSetup(c1);
  // Setup gStyle
  gStyleSetup();
  // Loop over all histograms inside Tlist
  for(int ii=0; ii<lout->GetSize(); ii++){
    // Get the name of the histogram
    const TString tag = lout->At(ii)->GetName();
    // Need to convert name from mc to data for finding corresponding data histogram when overlay
    const TRegexp re("_mc");
    const TString name = tag;
    name(re) = "_data";
    // Create histograms and stack
    TH1 * hh = dynamic_cast<TH1*> (lout->At(ii));
    THStack * hstk = dynamic_cast<THStack *>(lout->At(ii));
    TH2 * h2d = 0x0;
    //TH1 * holay = dynamic_cast<TH1*> (overlayList->FindObject(name));
    TH1D *holay = (TH1D*)overlayList->FindObject(name);
    if(hh){
      h2d = dynamic_cast<TH2 *>(hh);
      if(tag.Contains("RES")) {
        h2d->Draw("colz");
      }
      else if(holay){
      if(plotscale!=1) hh->Scale(plotscale);
      hh->Draw("hist");
      holay->SetMarkerStyle(8);
      holay->SetMarkerSize(1);
      holay->SetMarkerColor(kBlack);
      holay->SetLineColor(kBlack);
      holay->SetLineWidth(1);
      holay->Draw("same E");
      c1->Update();
      }
    }
    else if (hstk) {
      hstk->Draw("hist");
      holay = (TH1D*)overlayList->FindObject(name+"_sum");
      if(holay){
      holay->Draw();//to generate statbox
      c1->Update(); 
      if(plotscale!=1) ScaleStack(hstk, plotscale); 
      hstk->Draw("hist");
      holay->SetMarkerStyle(8);
      holay->SetMarkerSize(1);
      holay->SetMarkerColor(kBlack);
      holay->SetLineColor(kBlack);
      holay->SetLineWidth(1);
      holay->Draw("same E");
      c1->Update();
      }
    }
    else cout << "PlotUtils::DrawHist not found correct histogram!" << endl;
    c1->Print(outdir+"/"+tag+".png");
  } // End of for loop
}

THStack * PlotUtils::ConvertToStack(const TH2D * hh)
{
  const TString tag = hh->GetName();
  const TString tit = hh->GetTitle();

  const int nx = hh->GetNbinsX();
  const double xmin = hh->GetXaxis()->GetBinLowEdge(1);
  const double xmax = hh->GetXaxis()->GetBinUpEdge(nx);

  const int ny = hh->GetNbinsY();

  const double oldintegral = hh->Integral(0, 10000, 0, 1000000);

  const int y0 = 1;
  const int y1 = ny;

  double newintegral = 0;
  //THStack * stk = new THStack(tag+"_stack", tag);
  THStack * stk = new THStack(tag, tag);
  // Get color for each y value
  const int *col=GetColorArray(ny);

  //need to take into account of overflow in y
  for(int iy = y0; iy<=y1; iy++){
    const double toty = hh->Integral(0, 100000, iy, iy);
    if(toty<1E-12){
      continue;
    }

    TH1D * htmp = new TH1D(Form("%sy%d", tag.Data(), iy), tit.Data(), nx, xmin, xmax);
    for(int ix=0; ix<=nx+1; ix++){
      const double ientry = hh->GetBinContent(ix, iy);
      newintegral += ientry;

      htmp->SetBinContent(ix, ientry);
    }

    const int icol = GetColor(col[iy-y0]);//need constant map between y and color
    htmp->SetFillColor(icol);
    htmp->SetLineColor(kBlack);
    htmp->SetMarkerSize(2);
    printf("style::ConvertToStack %s adding y %f with color %d\n", tag.Data(), hh->GetYaxis()->GetBinCenter(iy), icol);
    stk->Add(htmp);
  }
  if(oldintegral!=newintegral){
    printf("style::ConvertToStack integral not matched! %s old %f new %f\n", tag.Data(), oldintegral, newintegral); exit(1);
  }

  return stk;
}
TH1D * PlotUtils::GetStackedSum(THStack *stk)
{
  const TList * ll = stk->GetHists();
  const TString tag = stk->GetName();
  TH1D * hout = (TH1D*)ll->At(0)->Clone(tag);
  hout->SetName(tag+"_sum");
  hout->SetTitle(tag);
  hout->SetDirectory(0);
  for(Int_t ii=1; ii<ll->GetEntries(); ii++){
    hout->Add((TH1D*)ll->At(ii));
  }

  hout->SetEntries(hout->Integral(0,10000));

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

int PlotUtils::GetColor(const int col)
{
  if(col>=1000){
    return col-1000;
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

void PlotUtils::gStyleSetup()
{
  gStyle->SetOptTitle(1);
  gStyle->SetStatY(0.87);
  gStyle->SetStatH(0.12);
  gStyle->SetStatX(0.83);
  gStyle->SetStatW(0.12);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetTitleX(0.55);
}


