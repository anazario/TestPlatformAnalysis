#include <iostream>
#include <vector>

#include <TLatex.h>
#include <TString.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include "TLegend.h"
  
void CMSmark(TString plotTitle);
void Plot1D(TH1D* hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog);
void Plot2D(TH2D* hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog);
void PlotGraph(TString name, TString xlabel, TString ylabel, int size, std::vector<double> xVec, std::vector<double> yVec);
void PlotWaveform(const std::vector<double>& time, const std::vector<double>& waveform, const int waveform_id);

TGraph* PlotScatter(TString name, TString xlabel, TString ylabel, std::vector<double> xVec, std::vector<double> yVec);
void makeScatterPlot(const std::vector<double> &x_signal, const std::vector<double> &x_background,
                     const std::vector<double> &y_signal, const std::vector<double> &y_background,
                     std::string name, std::string xlabel, std::string ylabel);