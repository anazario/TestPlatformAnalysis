#include "Plot.h"

using namespace std;

void CMSmark(TString plotTitle){
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(0.51,0.91,plotTitle);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.11,0.91,"#bf{CMS} #it{Preliminary}");
}

void Plot1D(TH1D* hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog){

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(11111111);
  TCanvas can("can_"+name,"can_"+name,600,800);
  if(isLog)
    can.SetLogy();

  can.SetLeftMargin(0.12);
  can.SetRightMargin(0.3);
  can.SetBottomMargin(0.1);
  can.Update();
  
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->Draw();

  can.SaveAs(name+".pdf");
  can.Close();
}

void Plot2D(TH2D *hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  
  TCanvas can("can_"+name,"can_"+name,600,800);
  if(isLog)
    can.SetLogz();

  can.SetLeftMargin(0.12);
  can.SetRightMargin(0.2);
  can.SetBottomMargin(0.1);
  can.Update();
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetZaxis()->SetTitle("Events");
  hist->GetZaxis()->SetTitleOffset(1.5);

  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetZaxis()->CenterTitle();
  
  hist->Draw("colz");
  can.SaveAs(name+".pdf");
  can.Close();
}

void PlotGraph(TString name, TString xlabel, TString ylabel, int size, std::vector<double> xVec, std::vector<double> yVec){

  TGraph tg(size,&xVec[0],&yVec[0]);

  TCanvas can("can_"+name,"can_"+name,600,800);
  tg.SetLineWidth(2);
  tg.GetXaxis()->SetTitle(xlabel);
  tg.GetYaxis()->SetTitle(ylabel);
  tg.Draw("AL");
  can.SetLeftMargin(0.12);
  can.SaveAs(name+".pdf");

  can.Close();
}

TGraph* PlotScatter(TString name, TString xlabel, TString ylabel, std::vector<double> xVec, std::vector<double> yVec){

  int size = xVec.size();
  
  TGraph* tg = new TGraph(size,&xVec[0],&yVec[0]);

  //tg->SetMarkerStyle(8);
  tg->GetXaxis()->SetTitle(xlabel);
  tg->GetYaxis()->SetTitle(ylabel);
  tg->SetMarkerSize(0.8);
  //tg->Draw("AP");

  return tg;
}
