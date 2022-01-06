#include "DataReaderH5.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <fstream>
#include <map>

int main(void){

  int x = 8;//x index of input file
  int y = 5;//y index of input file
  int s = 1;//index of sample
  int ch = 1;//select channel to analyze

  //create DataReaderH5 object and initialize with input file and channel
  DataReaderH5 d(("data/Scan15/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),ch);
  //print general scan info and file details
  d.PrintInfo();

  //get all 1000 samples from DataReaderH5 object
  vector<vector<float>> samples = d.GetWaveForms();
  //get corresponding time arrays
  vector<vector<float>> time = d.GetTimeArr();

  //make plot with single pulse
  vector<float> xVec = time[s];
  vector<float> yVec = samples[s];

  TCanvas *cv = new TCanvas("cv","cv",600,800);
  TGraph *g = new TGraph(SAMPLE_SIZE,&xVec[0],&yVec[0]);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(0.5);
  g->GetXaxis()->SetTitle("time [ns]");
  g->GetYaxis()->SetTitle("Amplitude [mV]");

  gPad->SetGrid(1, 1); gPad->Update();
  g->Draw("AP");
  cv->SaveAs(("graph_pulse"+to_string(s)+".pdf").c_str());

  return 0;
}

