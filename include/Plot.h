#ifndef Plot_h
#define Plot_h

#include <iostream>
#include <vector>

#include <TLatex.h>
#include <TString.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>

using namespace std;

class Plot{

 public:

  Plot();
  virtual ~Plot();
  
  static void CMSmark(TString plotTitle);
  static void Plot1D(TH1F hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog=false);
  static void PlotGraph(TString name, TString xlabel, TString ylabel, int size, vector<float> xVec, vector<float> yVec);

};

#endif
