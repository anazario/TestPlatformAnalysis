#include <TLatex.h>
#include <TString.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>


inline void CMSmark(TString plotTitle){
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  // l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());                                                                                                         
  l.DrawLatex(0.51,0.91,plotTitle);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.11,0.91,"#bf{CMS} #it{Preliminary}");
}
/*
void SetPlotPar(TH1D* hist, TString xname, TString yname){

  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.06);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitle(xname);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.12);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitle(yname);
}

void Plot1D(TH1D* hist){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",600,800);

  can->SetLeftMargin(0.1);
  can->SetRightMargin(0.15);
  can->SetBottomMargin(0.15);
  can->Draw();
  can->cd();
  hist->Draw("");
}
*/


