#ifndef UnbinnedFit_h
#define UnbinnedFit_h

//standard libraries
#include <iostream>
#include <vector>

//Root libraries
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"

//Roofit libraries
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooLandau.h"

//Local libraries
#include "Plot.h"

class UnbinnedFit{
  
public:
  
  UnbinnedFit(const std::vector<double>& data, const double xMin = 0., const double xMax = 750.);
  virtual ~UnbinnedFit();

  double GetLandauMPV() const;
  double GetLandauSigma() const;

  double GetGaussMean() const;
  double GetGaussSigma() const;

  void setRange(const TString name, const double rangeMin, const double rangeMax);
  RooGaussian* toGauss(const double mean, const double sigma);
  RooLandau* toLandau(const double mean, const double sigma);
  RooFFTConvPdf* toLandauXgauss(const double meanL, const double sigmaL, const double meanG, const double sigmaG);

  RooDataSet* GetDataSet();
  RooRealVar* GetFitVar();

  TCanvas* Plot(const TString name, const int nBins = 40, const double xMin = 0., const double xMax = 750.);
  
private:

  double mpvLandau_ = -999.;
  double sigmaLandau_ = -999.;
  double meanGauss_ = -999.;
  double sigmaGauss_ = -999.;

  bool fitLandau_ = false;
  bool fitGauss_ = false;
  bool fitConv_ = false;  
  bool rangeOn_ = false;
  
  TString fitname_;
  TString rangeName_;
  
  RooDataSet* dataSet_ = nullptr;
  RooRealVar* var_ = nullptr;
  RooArgSet* varArgset_ = nullptr;
  //Landau fit
  RooRealVar* meanL_ = nullptr;
  RooRealVar* sigmaL_ = nullptr;
  RooLandau* landau_ = nullptr;
  //Gaussian fit                                                                                                                                          
  RooRealVar* meanG_ = nullptr;
  RooRealVar* sigmaG_ = nullptr;
  RooGaussian* gauss_ = nullptr;
  //LandauXgaussian convolution
  RooFFTConvPdf* lxg_ = nullptr;
  
  void SetupDataSet(const std::vector<double>& data, const double xMin, const double xMax);
  void SetLandauFit(const double mean, const double sigma);
  void SetGaussFit(const double mean, const double sigma);

  TCanvas* GetPlot(RooAbsPdf* fitPdf, const TString name, const int nBins, const double xMin, const double xMax);
};

#endif
