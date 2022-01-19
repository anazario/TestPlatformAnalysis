#ifndef PULSEVARIATION_H_
#define PULSEVARIATION_H_

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TColor.h"
#include "TStyle.h"

#include "Plot.h"
#include "PulseTools.h"

class PulseVariation{

 public:
  
  PulseVariation(const vector<TH1F> histCollection, const vector<float> timeAxis);
  virtual ~PulseVariation();

  vector<float> GetMeanPulse();

  void GetErrors(vector<float>& lowErrors, vector<float>& highErrors, float CI);
  void PlotMeanErrors(const TString name, const TString ylabel);
  void PlotHistograms();

 private:

  int collectionSize_;
  vector<float> timeAxis_;
  vector<float> collectionMean_;
  vector<TH1F> histCollection_;
  
};

#endif
