#ifndef PULSEVARIATION_H_
#define PULSEVARIATION_H_

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TColor.h"
#include "TStyle.h"

#include "Plot.h"
#include "PulseList.h"
#include "PulseTools.h"

class PulseVariation{

 public:
  PulseVariation(const PulseList pulseList);  
  PulseVariation(const std::vector<TH1D> histCollection, const std::vector<double> timeAxis);
  virtual ~PulseVariation();

  std::vector<double> GetMeanPulse() const;

  void GetErrors(std::vector<double>& lowErrors, std::vector<double>& highErrors, double CI);
  void PlotMeanErrors(const TString name, const TString ylabel);
  void PlotHistograms();

 private:

  int collectionSize_;
  std::vector<double> timeAxis_;
  std::vector<double> collectionMean_;
  std::vector<TH1D> histCollection_;
  
};

#endif
