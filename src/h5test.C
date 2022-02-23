#include <iostream>
#include <fstream>
#include <map>

#include "Scan.h"
#include "DataReaderH5.h"
#include "PulseList.h"
#include "PulseTools.h"
#include "Plot.h"
#include "UnbinnedFit.h"

#include <TROOT.h>
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
/*
#include "RooFit.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooLandau.h"
#include "RooKeysPdf.h"
#include "RooBinning.h"
*/
using namespace std;

int main(void){

  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int channel = 2;
  string path = "data/ScanCh2_1/hdf5/";
  
  Scan scan(path);

  int totalFiles = scan.GetTotalFiles();
  int count = 0;

  string dataLoc = scan.GetPathToData();

  //vector<double> mpvLXGlandau;
  //vector<double> mpvLXGgauss;
  
  vector<double> mpvBinned;
  vector<double> mpvLXG;
  vector<double> mpvL;
  vector<double> rateVec;
  
  scan.PrintInfo();

  TString name = ("Amplitude_ch"+to_string(channel)).c_str();
  //TH2D* rateVsMPV_landau = new TH2D("rateVsMPV_landau","rateVsMPV_landau", 20, 150.,160.,50,5.,30.);
  //TH2D* rateVsMPV_lxg = new TH2D("rateVsMPV_lxg","rateVsMPV_lxg", 20, 150.,160.,50,5.,30.);

  TCanvas cv("cv","cv",600,800);
  
  std::ofstream outfile ("test.txt");
  std::ofstream outfileL ("landau.txt");
  std::ofstream outfileLXG ("landauXgauss.txt");
  
  for(auto file : scan.GetFileList()){

    count++;
    fprintf(stdout, "\r  Processed files: %5d of %5d ", count, totalFiles);
    fflush(stdout);

    DataReaderH5 *data = new DataReaderH5((dataLoc+file).c_str());

    double x = data->GetXpos();
    double y = data->GetYpos();
    double rate = data->GetRate(channel);
    
    if(rate > 0.){
      rateVec.push_back(rate);
      
      Replace(file,".hdf5",""); 
      PulseList pulseList = data->GetPulseList(channel);
      vector<double> ampDist = pulseList.GetAmpDistribution();
      
      //binned fit
      TH1D* hist = GetHistFromVec(ampDist, name+"_"+file, 40, 0., 750.);
      TF1 *fitB = new TF1("fitB","landau",110.,450);
      hist->Fit("fitB", "Q", "", 110.,450.);
      double mpv = fitB->GetParameter(1);
      mpvBinned.push_back(mpv);
      delete hist;
      delete fitB;
      //unbinned fit
      UnbinnedFit fit(ampDist);
      fit.setRange("fitSignal_"+TString(file),125.,400.);
      //LXG fit
      fit.toLandauXgauss(155.,15.,0.,15.);
      mpv = fit.GetLandauMPV() + fit.GetGaussMean();
      outfileLXG << x << "\t" << y << "\t" << mpv << std::endl;
      mpvLXG.push_back(mpv);
      //mpvLXGlandau.push_back(fit.GetLandauMPV());
      //mpvLXGgauss.push_back(fit.GetGaussMean());
      //rateVsMPV_lxg->Fill(mpv,rate);
      //plot lxg fit
      /*
      TCanvas *can1 = fit.Plot(TString(file)+"_lxg");
      TString fitname = "MPV_"+file;
      can1->SaveAs("MPVplots/LandauXgauss/"+fitname+".pdf");
      can1->Close();
      delete can1;      
      */
      //Landau fit
      fit.toLandau(155.,15.);
      mpv = fit.GetLandauMPV();
      outfileL << x << "\t" << y << "\t" << mpv << std::endl;
      mpvL.push_back(mpv);
      //rateVsMPV_landau->Fill(mpv,rate);
      //plot landau fit
      /*
      TCanvas *can2 = fit.Plot(TString(file)+"_landau");
      fitname = "MPV_"+file;
      can2->SaveAs("MPVplots/Landau/"+fitname+".pdf");
      can2->Close();
      delete can2;
      */
    }
    else{
      outfileL << x << "\t" << y << "\t" << 0. << endl;
      outfileLXG << x << "\t" << y << "\t" << 0. << endl;
    }
    delete data;
  }
  cout << endl;
  outfile.close();

  name = "";
  TCanvas can("can_"+name,"can_"+name,600,800);

  TGraph *scatter1 = PlotScatter("scatter1", "MPV [mV]", "Rate [Hz]", mpvLXG, rateVec);
  TGraph *scatter2 = PlotScatter("scatter2", "MPV [mV]", "Rate [Hz]", mpvL, rateVec);
  TGraph *scatter3 = PlotScatter("scatter3", "MPV [mV]", "Rate [Hz]", mpvBinned, rateVec);
  //TGraph *scatter4 = PlotScatter("scatter4", "MPV [mV]", "Rate [Hz]", mpvLXGlandau, rateVec);
  //TGraph *scatter5 = PlotScatter("scatter5", "MPV [mV]", "Rate [Hz]", mpvLXGgauss, rateVec);

  //TF1* f1 = new TF1("f1","[0]/(x)+[1]",140.,170.);
  //scatter1->Fit(f1,"ROB","",155,170);
  //scatter1->Draw("AP");
  //f1->Draw("same");
  scatter1->SetMarkerStyle(8);
  scatter2->SetMarkerStyle(33);
  scatter2->SetMarkerColor(kGreen+2);
  //scatter2->Draw("sameP");
  scatter3->SetMarkerStyle(21);
  scatter3->SetMarkerColor(kBlue+1);
  //scatter3->Draw("sameP");
  //scatter4->SetMarkerStyle(29);
  //scatter4->SetMarkerColor(kRed+2);
  //scatter4->Draw("sameP");
  //scatter5->SetMarkerStyle(23);
  //scatter5->SetMarkerColor(kOrange);
  //scatter5->Draw("sameP");

  TMultiGraph *mg = new	TMultiGraph();

  mg->Add(scatter1);
  mg->Add(scatter2);
  mg->Add(scatter3);
  //mg->Add(scatter4);
  //mg->Add(scatter5);

  //mg->SetMaximum(60.);
  mg->GetXaxis()->SetTitle("MPV [mV]");
  mg->GetYaxis()->SetTitle("Rate [Hz]");
  mg->Draw("AP");

  TLegend leg(0.7,0.75,0.89,0.89);
  leg.AddEntry(scatter1,"UB LxG","p");
  leg.AddEntry(scatter2,"UB Landau","p");
  leg.AddEntry(scatter3,"B Landau","p");
  //leg.AddEntry(scatter4,"UB Landau MPV","p");
  leg.Draw("same");
  
  can.SetLeftMargin(0.12);
  name = "plots/RateVsMPV_lxg_ScanCh2_1";
  can.SaveAs(name+".pdf");

  can.Close();
  delete scatter1;
  
  //Plot2D(rateVsMPV_landau, "RateVsMPV_landau_ScanCh2_1", "MPV [mV]", "Rate [Hz]", false);
  //Plot2D(rateVsMPV_lxg, "RateVsMPV_lxg_ScanCh2_1", "MPV [mV]", "Rate [Hz]", false);
  //Plot1D(*hist, name, "Amplitude [mV]", "Events", true);
  //cout << allPulses.size() << endl;
  
  /*
  int x = 8;//x index of input file
  int y = 5;//y index of input file
  int s = 11;//index of sample
  int ch = 2;//select channel to analyze

  //create DataReaderH5 object and initialize with input file and channel
  DataReaderH5 data(("data/ScanCh2_1/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),ch);
  //print general scan info and file details
  data.PrintInfo();

  //get all data from DataReaderH5 object as a PulseList object
  PulseList pulseList = data.GetPulseList();

  vector<double> ampDist = pulseList.GetAmpDistribution();

  TString name = ("Amplitude_ch"+to_string(ch)).c_str();
  TH1D* hist = GetHistFromVec(ampDist, name, 52, 0., 780.);

  TF1 *fit = new TF1("fit","landau",0,700);
  hist->Fit("fit", "Q", "", 0,700);
  double mpv = fit->GetParameter(1);

  Plot1D(hist, name, "Amplitude [mV]", "Events", true);
  ////////////////////////////////////////////////////////////////////////////////////////////////
  Pulse *pulse = pulseList[s];

  vector<double> time;
  vector<double> waveForm;
  pulse->GetCenteredPulse(waveForm,time);

  int size = 300;
  double start = -10.;
  double end = 30.; 
  
  vector<double> interpolationTime = LinSpaceVec(start,end,size);
  vector<double> interpolatedPulse = pulse->GetInterpolation(size, start, end, pulse->GetMaxTime());
  vector<double> CDFpulse = pulse->GetCDF(size, start, end, pulse->GetMaxTime());

  PlotGraph("plots/pulse", "time [ns]", "Amplitude [mV]", waveForm.size(),time,waveForm);
  PlotGraph("plots/interpolation", "time [ns]", "Amplitude [mV]", interpolatedPulse.size(), interpolationTime, interpolatedPulse);
  PlotGraph("plots/cdf", "time [ns]", "Amplitude [mV]", CDFpulse.size(), interpolationTime, CDFpulse);
  */
  return 0;
}

