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
      // Check if the histogram name contians STK
      if(tag.Contains("STK") || tag.Contains("OVERLAY")){
        // Convert this tmp 2D histogram to a stack histogram
        THStack * stk = ConvertToStack(htmp,kMC);
        // Check if the stk exist 
        if(stk){
          // Add to lout for MC
          if(kMC) lout->Add(stk);   
          // Only need the sum of this stack histogram in data 
          else{
            // Get the sum of this hitogram
            TH1D * hsum = GetStackedSum(stk);
            lout->Add(hsum);
          }
        }
      } // End of STK tag

      // Check if the histogram name contians RES (resolution)
      else if(tag.Contains("RES") && kMC){
        // Column normalise each bin to easily see the maximum
        TH2D * hnor = NormalHist(htmp, 5, true);
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
        TH2D * hnor = NormalHist(htmp, 5, true);
        hnor->SetTitle(tag);
        lout->Add(hnor);
        if(tag.Contains("RawRecVSTruth")){
          getProfileFit(htmp);
          getSliceXDrawY(htmp);
        } 
      }
      // Do nothing (You can add more else if to process more tags)
      else {}
    } // End of if(htmp)

    // Stack histogram
    THStack * hstk = dynamic_cast<THStack *>(lout->At(ii));
    if(hstk){
      const TString tag = hstk->GetName();
      if(tag.Contains("stkTruth") && kMC == true){
        // Get the name of raw shower histogram
        TRegexp re("stk");
        TString tmp = tag;
        tmp(re) = "";

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
        // Draw 2D histogram from its name (MC only)
        if(!tag.Contains("proj") && (tag.Contains("RES") || tag.Contains("REG"))) {
          h2d->Draw("colz");
        }
      }
      // There is no h2d, we only have 1D histogram
      else {
        // Check if we have data overlay histogram
        if(holay){
          // Do the scaling for MC/Data
          if(plotscale!=1) hh->Scale(plotscale);
          hh->Draw("hist");
          DrawOverlay(holay);
          c1->Update();
        }
        // No data overlay histogram
        else{
          // Draw overlay for raw and corrected shower resolution
          if(!tag.Contains("RAW") && tag.Contains("proj")){
            // Get the name of raw shower histogram
            TRegexp re("_proj");
            TString name = tag;
            name(re) = "_RAW_proj";
            TH1D *holayRaw = (TH1D*)lout->FindObject(name);    
            TString fit = tag;
            fit(re) = "_FIT_proj";
            TH1D *holayFit = (TH1D*)lout->FindObject(fit);
            if(holayRaw){
              if(holayRaw->GetMaximum() > hh->GetMaximum()) hh->SetMaximum(holayRaw->GetMaximum()*1.2);
              else hh->SetMaximum(hh->GetMaximum()*1.2);
              hh->SetFillStyle(1);
              hh->SetLineColor(kBlue); 
              hh->SetLineWidth(3);
              hh->Draw("hist");
              c1->Update();
              holayRaw->SetFillStyle(1);
              holayRaw->SetLineColor(kRed);
              holayRaw->SetLineWidth(3);
	            holayRaw->SetStats(1);
              holayRaw->Draw("hist same"); 
              c1->Update();
              if(holayFit){
                holayFit->SetLineColor(kMagenta);
                holayFit->SetLineWidth(3);
                holayFit->Draw("hist same");
                c1->Update();
              }
            }
            else hh->Draw("hist");
          }
          else hh->Draw("hist");
        }
      }
    } // End of if(hh)

    ////=========================== stack histogram ===========================//
    else if (hstk) {
      // Get the data overlay histogram for stack
      holay = (TH1D*)overlayList->FindObject(tag+"_sum");
      hstk->Draw("hist");
      // Draw legend if needed
      if(tag.Contains("stkTruth")){
        const TString tag = hstk->GetName();
        // Get the name of raw shower histogram
        TRegexp re("stk");
        TString tmp = tag;
        tmp(re) = "";
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
        double lxoff = 0.5;
        if(tmp.Contains("Dalphat") || tmp.Contains("MomIniPi") || tmp.Contains("ThetaIniPi")) lxoff = 0.05;
        TLegend * lg = new TLegend(lxoff+0.15, 0.65, lxoff+0.35, 0.85);
        TString lheader("signal phase space");
        lg->SetHeader(lheader);
        lg->AddEntry(hh1p0n, "1p0n", "f");
        lg->AddEntry(hhNp0n, "Np0n", "f");
        lg->AddEntry(hh1pMn, "1pMn", "f");
        lg->AddEntry(hhNpMn, "NpMn", "f");
        lg->Draw("same");
      }
      if(holay){
        holay->Draw();//to generate statbox
        c1->Update(); 
        if(plotscale!=1) ScaleStack(hstk, plotscale); 
        hstk->SetMaximum(holay->GetMaximum()*1.2);
        hstk->Draw("hist");
        DrawOverlay(holay);
        c1->Update();
        
        const int overlayColor = kBlack;

        vector<TString> evtType;
        evtType.push_back("signal");
        evtType.push_back("background");
        evtType.push_back("non-#pi^{+} beam");
        evtType.push_back("data");

        vector<TString> htype;
        htype.push_back("f");
        htype.push_back("f");
        htype.push_back("f");
        htype.push_back("ple");

        int *cols=GetColorArray(4);
        cols[3]=overlayColor;
        const int mrks[]={1,1,1,6};
        TLegend * lg = 0x0;
        lg = DrawLegend(evtType, htype, tag, cols, mrks);
        lg->Draw("same");
  
      }
      // Spectial case
      else if(tag.Contains("OVERLAY")){
        TH1D * hsum = dynamic_cast<TH1D*> (hstk->GetStack()->Last());
        hstk->SetMaximum(hsum->GetMaximum()*1.2);
        hstk->Draw("nostack");
        hsum->SetLineColor(kBlack);
        hsum->SetFillStyle(0);
        hsum->Draw("same");
        auto* legend = new TLegend(0.65, 0.65, 0.85, 0.85);
        // Cheat Legend
        TH1D * E1 = 0x0;
        TH1D * E2 = 0x0;
        E1 = new TH1D("E1",  "", 20, -0.5, 19.5); 
        E2 = new TH1D("E2",  "", 20, -0.5, 19.5); 
        E1->SetFillStyle(3004);E2->SetFillStyle(3004);
        E1->SetFillColor(1509);E2->SetFillColor(1505);E1->SetLineColor(1509);E2->SetLineColor(1505);

        legend->AddEntry(E1, "SubLeading photon", "f");
        legend->AddEntry(E2, "Leading photon", "f");
        legend->AddEntry(hsum, "Photon Spectrum", "f");
        legend->Draw("same");
      }
      
    }
    else cout << "PlotUtils::DrawHist not found correct histogram!" << " name: " << tag << endl;
    c1->Print(outdir+"/"+tag+".png");
  } // End of for loop
}

THStack * PlotUtils::ConvertToStack(const TH2D * hh, const bool kMC)
{
  const TString tag = hh->GetName();
  const TString tit = hh->GetTitle();
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
  THStack * stk = new THStack(tag, tag);
  // Get color for each y value
  const int *col=GetColorArray(ny);

  //need to take into account of overflow in y
  for(int iy = y0; iy<=y1; iy++){
    const double toty = hh->Integral(0, 100000, iy, iy);
    if(toty<1E-12){
      continue;
    }

    TH1D * htmp = new TH1D(Form("%s%sy%d", typ.Data(), tag.Data(), iy), tit.Data(), nx, xmin, xmax);
    for(int ix=0; ix<=nx+1; ix++){
      const double ientry = hh->GetBinContent(ix, iy);
      newintegral += ientry;

      htmp->SetBinContent(ix, ientry);
    }

    const int icol = GetColor(col[iy-y0]);//need constant map between y and color
    htmp->SetFillColor(icol);
    htmp->SetLineColor(icol);
    htmp->SetMarkerSize(2);
    if(tag.Contains("OVERLAY")) htmp->SetFillStyle(3004);
    printf("PlotUtils::ConvertToStack %s adding y %f with color %d\n", tag.Data(), hh->GetYaxis()->GetBinCenter(iy), icol);
    stk->Add(htmp);
  }
  if(oldintegral!=newintegral){
    printf("PlotUtils::ConvertToStack integral not matched! %s old %f new %f\n", tag.Data(), oldintegral, newintegral); exit(1);
  }

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

TH2D * PlotUtils::NormalHist(const TH2D *hraw, const Double_t thres, const Bool_t kmax)
{
  TH2D *hh=(TH2D*)hraw->Clone(Form("%s_nor",hraw->GetName()));
  hh->Scale(0);

  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  Double_t hmax = -1e10;
  Double_t hmin = 1e10;
  Double_t nent = 0;
  for(Int_t ix=x0; ix<=x1; ix++){

    TH1D * sliceh = hraw->ProjectionY(Form("tmpnormalhist%sx%d", hh->GetName(), ix), ix, ix, "oe");
    const Double_t tot = sliceh->GetEntries();

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

  gStyle->SetOptTitle(1);
  gStyle->SetStatY(0.87);
  gStyle->SetStatH(0.12);
  gStyle->SetStatX(0.83);
  gStyle->SetStatW(0.12);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(1);
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

void PlotUtils::DrawOverlay(TH1D *holay)
{
  holay->SetMarkerStyle(8);
  holay->SetMarkerSize(1);
  holay->SetMarkerColor(kBlack);
  holay->SetLineColor(kBlack);
  holay->SetLineWidth(1);
  holay->Draw("same E");
}

void PlotUtils::getSliceXDrawY(TH2D * h2d)
{ 
  auto cc = new TCanvas("cc","cc",1600,1200);
  gStyle->SetOptStat(0);
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
    cc->Print("output/hSliceXDrawY.eps");
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

  TLegend * lg = new TLegend(0.6, 0.6, 0.88, 0.88);

  for(int ii=0; ii<nent; ii++){
    TH1D * hh=new TH1D(Form("h%d%s",ii,tag.Data()),"",1,0,1);
    const int col = GetColor(cols[ii]);
    hh->SetFillColor(col);
    hh->SetLineColor(col);
    hh->SetMarkerStyle(mkrs[ii]);
    hh->SetMarkerSize(3);
    hh->SetMarkerColor(col);
    lg->AddEntry(hh, entries[ii], htype[ii]);
  }

  lg->SetNColumns(ncol);

  return lg;
}
