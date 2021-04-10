#include "../include/PlotUtils.h"

void PlotUtils::DrawHist(TList *lout, const TString outdir)
{
  TCanvas * c1 = new TCanvas("c1", "", 1200, 800);
  PadSetup(c1);
  gStyleSetup();
  
  for(int ii=0; ii<lout->GetSize(); ii++){
    const TString tag = lout->At(ii)->GetName();
    TH1D * hh = (TH1D*)lout->FindObject(tag);
    hh->Draw("hist");
    c1->Print(outdir+"/"+tag+".png");
  } // End of for loop
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


