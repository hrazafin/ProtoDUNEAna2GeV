#include "BetheBloch.h"
#include "util.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

int main(){

  SetProtoDUNEStyle();
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);


  // Create a BetheBloch object for muon
  BetheBloch bb(13);

  //cout<<bb.meandEdx(1000)<<" "<<bb.MPVdEdx(1000, 0.5)<<endl;
  //cout<<bb.KEAtLength(1000, 100)<<" "<<bb.KEAtLength(1000,500)<<endl;
  //https://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
  const int np = 29;
  double spline_KE[np] = {
    10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
    400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
    20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

  double spline_Range[np] = {
    9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1,
    1.063E2,  1.725E2, 2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3,
    2.297E3,  4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
    4.326E4,  5.768E4, 7.734E4, 1.060E5, 1.307E5};

  for (int i = 0; i<np; ++i) spline_Range[i] /= 1.396;

  double bbRange[np];
  double bbRangeSpline[np];

  for (int i = 0; i<np; ++i){
    //cout<<spline_KE[i]<<" "<<bb.RangeFromKE(spline_KE[i])<<endl;
    bbRange[i] = bb.RangeFromKE(spline_KE[i]);
    bbRangeSpline[i] = bb.RangeFromKESpline(spline_KE[i]);
    //cout<<spline_KE[i]<<" "<<spline_Range[i]<<" "<<bbRange[i]<<" "<<bbRangeSpline[i]<<endl;
  }

  TGraph *grmuKERangePDG = new TGraph(np, spline_Range, spline_KE);
  TGraph *grmuKERangeBB = new TGraph(np, bbRange, spline_KE);
  TGraph *grmuKERangeBBSpline = new TGraph(np, bbRangeSpline, spline_KE);

  TCanvas *c1 = new TCanvas("c1","c1");
  grmuKERangePDG->SetMarkerStyle(24);
  grmuKERangePDG->SetMarkerColor(2);
  grmuKERangePDG->SetLineColor(2);
  grmuKERangePDG->Draw("apc");
  grmuKERangePDG->SetTitle("Muon");
  grmuKERangePDG->GetXaxis()->SetTitle("Range (cm)");
  grmuKERangePDG->GetYaxis()->SetTitle("KE (MeV)");
  grmuKERangeBB->SetMarkerStyle(24);
  grmuKERangeBB->SetMarkerColor(4);
  grmuKERangeBB->SetLineColor(4);
  grmuKERangeBB->Draw("pc");
  grmuKERangeBBSpline->SetMarkerStyle(24);
  grmuKERangeBBSpline->SetMarkerColor(3);
  grmuKERangeBBSpline->SetLineColor(3);
  grmuKERangeBBSpline->Draw("pc");
  gPad->SetLogx();
  gPad->SetLogy();
  TLegend *leg1 = new TLegend(0.2,0.65,0.5,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(grmuKERangePDG, "PDG","lp");
  leg1->AddEntry(grmuKERangeBB, "Bethe Bloch","lp");
  leg1->AddEntry(grmuKERangeBBSpline, "Bethe Bloch spline","lp");
  leg1->Draw();
  c1->Print("muKERange.png");
  c1->Print("muKERange.pdf");

  gPad->SetLogx(1);
  gPad->SetLogy(0);
  const int np2 = 1000;
  double rr[np2];
  double KE[4][np2];
  double dEdx[4][np2];
  double pdg[4] = {13, 211, 321, 2212};
  char particle[4][100] = {"Muon", "Pion", "Kaon", "Proton"};
  TGraph *grrrdEdx[4];
  TGraph *grrrKE[4];
  for (int i = 0; i<4; ++i){
    if (bb.GetPdgCode() != pdg[i]){
      bb.SetPdgCode(pdg[i]);
    }
    for (int j = 0; j<np2; ++j){
      if (i==0){
        rr[j] = pow(10, log10(0.1)+j*log10(1000/0.1)/np2);
      }
      KE[i][j] = bb.KEFromRangeSpline(rr[j]);
      dEdx[i][j] = bb.meandEdx(KE[i][j]);
//      if (j==0) cout<<i<<" "<<rr[j]<<" "<<KE[i][j]<<" "<<dEdx[i][j]<<endl;
//      if (j==np2-1) cout<<i<<" "<<rr[j]<<" "<<KE[i][j]<<" "<<dEdx[i][j]<<endl;
//      if (i==3){
//        cout<<rr[j]<<" "<<KE[i][j]<<" "<<dEdx[i][j]<<endl;
//      }
    }
    grrrdEdx[i] = new TGraph(np2, rr, dEdx[i]);
    grrrKE[i] = new TGraph(np2, rr, KE[i]);
  }

  TLegend *leg2 = new TLegend(0.6,0.6,0.8,0.9);
  leg2->SetFillStyle(0);
  for (int i = 0; i<4; ++i){
    grrrdEdx[i]->SetLineColor(i+1);
    grrrdEdx[i]->SetLineWidth(2);
    grrrdEdx[i]->Draw(i==0?"ac":"c");
    if (i==0){
      grrrdEdx[i]->SetTitle("");
      grrrdEdx[i]->GetXaxis()->SetTitle("Residual range (cm)");
      grrrdEdx[i]->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
      grrrdEdx[i]->GetYaxis()->SetRangeUser(0,55);
    }
    leg2->AddEntry(grrrdEdx[i], particle[i], "l");
  }
  leg2->Draw();
  c1->Print("RRdEdx.png");
  c1->Print("RRdEdx.pdf");

  gPad->SetLogx(1);
  gPad->SetLogy(1);
  TLegend *leg3 = new TLegend(0.2,0.6,0.4,0.9);
  leg3->SetFillStyle(0);
  for (int i = 0; i<4; ++i){
    grrrKE[i]->SetLineColor(i+1);
    grrrKE[i]->SetLineWidth(2);
    grrrKE[i]->Draw(i==0?"ac":"c");
    if (i==0){
      grrrKE[i]->SetTitle("");
      grrrKE[i]->GetXaxis()->SetTitle("Range (cm)");
      grrrKE[i]->GetYaxis()->SetTitle("KE (MeV)");
    }
    leg3->AddEntry(grrrKE[i], particle[i], "l");
  }
  leg3->Draw();
  c1->Print("RRKE.png");
  c1->Print("RRKE.pdf");

  return 0;

}
