#include <TLatex.h>
#include <TString.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>


inline void CMSmark(TString plotTitle){
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

inline void Plot1D(TH1F hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog=false){

  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas can("can_"+name,"can_"+name,600,800);
  if(isLog)
    can.SetLogy();
  
  can.SetLeftMargin(0.12);
  can.SetRightMargin(0.1);
  can.SetBottomMargin(0.1);
  hist.GetXaxis()->SetTitle(xlabel);
  hist.GetYaxis()->SetTitle(ylabel);
  //can.cd();
  hist.SetLineColor(kBlack);
  hist.SetLineWidth(2);
  hist.Draw();
  can.SaveAs(name+".pdf");
  can.Close();
}



