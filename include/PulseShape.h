#ifndef PULSESHAPE_H
#define PULSESHAPE_H

#define SAMPLE_SIZE 2002

#include <iostream>
#include <map>
#include <vector>
#include <math.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom.h"
#include "TEventList.h"
#include "TColor.h"
#include "Plot.h"
#include "TStyle.h"
#include "TLatex.h"
#include "PulseTools.h"



class PulseShape{

 public:

  PulseShape(vector<vector<float>> sampleCollection);

  virtual ~PulseShape();

  void PlotHists();
  void GetMeanErrors(vector<float>& mean, vector<float>& low_err, vector<float>& high_err, int CI);
  void GetErrors(vector<float>& low_err, vector<float>& high_err, float CI);
  void GetMaxT0andAmp0(float* sample, float& t0, float& amp0);
  float* GetMeanPulse(vector<float> t0, vector<float> max_amp, float step_size, std::function<float*(float*,float,float,float)> PulseType);
  static float* MakeInterpolation(float* sample, float t0, float amp_peak, float step_size);
  static float* GetPulseCDF(float* sample, float t0, float max_amp, float step_size);

  void PlotSinglePulse(int index);
  
  void PlotAllPulses();
  void PlotAllCDFs();
  void PlotPulseMean();
  void PlotCDFMean();
  void PlotPulseMeanErr();
  void PlotCDFMeanErr();
  void PlotRangeCDF(float xmin, float xmax, int total_slices);
  void PlotRangePulse(float xmin, float xmax, int total_slices);

 private:

  vector<vector<float>> samples_;
  int sampleSize_;
  vector<float> tPeak_;
  vector<float> maxAmp_;
  float stepSize_;

  vector<TH1F> histVec_;
  void GetMaxT0Amp0Vec();
  void PlotAll(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, 
	       std::function<float*(float*,float,float,float)> PulseType);
  void PlotMean(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, 
		std::function<float*(float*,float,float,float)> PulseType);
  void PlotMeanErr(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, TString opt, 
		   std::function<float*(float*,float,float,float)> PulseType);
  //void PlotRange(float xmin, float xmax, int total_slices, TString name, TString y_label, TString opt, 
  //std::function<float*(float*,float,float,float)> PulseType);
};

#endif
