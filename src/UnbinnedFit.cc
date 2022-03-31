#include "UnbinnedFit.h"

UnbinnedFit::UnbinnedFit(const std::vector<double>& data, const double xMin, const double xMax){
  SetupDataSet(data, xMin, xMax);
}

UnbinnedFit::~UnbinnedFit(){
  if(dataSet_)
    delete dataSet_;
  if(var_)
    delete var_;
  if(varArgset_)
    delete varArgset_;
  if(meanL_)
    delete meanL_;
  if(sigmaL_)
    delete sigmaL_;
  if(landau_)
    delete landau_;
  if(meanG_)
    delete meanG_;
  if(sigmaG_)
    delete sigmaG_;
  if(gauss_)
    delete gauss_;
  if(lxg_)
    delete lxg_;
}

//public methods

double UnbinnedFit::GetLandauMPV() const{
  return mpvLandau_;
}

double UnbinnedFit::GetLandauSigma() const{
  return sigmaLandau_;
}

double UnbinnedFit::GetGaussMean() const{
  return meanGauss_;
}

double UnbinnedFit::GetGaussSigma() const{
  return sigmaGauss_;
}

RooDataSet* UnbinnedFit::GetDataSet(){
  return dataSet_;
}

RooRealVar* UnbinnedFit::GetFitVar(){
  return var_;
}

void UnbinnedFit::setRange(const TString name, const double rangeMin, const double rangeMax){
  rangeName_ = name;
  var_->setRange(rangeName_, rangeMin, rangeMax);
  rangeOn_ = true;
}

RooGaussian* UnbinnedFit::toGauss(const double mean, const double sigma){
  SetGaussFit(mean, sigma);

  if(rangeOn_)
    gauss_->fitTo(*dataSet_, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE),RooFit::Range(rangeName_));
  else
    gauss_->fitTo(*dataSet_, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE));

  meanGauss_ = meanG_->getValV();
  sigmaGauss_ = sigmaG_->getValV();

  fitLandau_ = false;
  fitGauss_ = true;
  fitConv_ = false;
  
  return gauss_;
}

RooLandau* UnbinnedFit::toLandau(const double mean, const double sigma){
  SetLandauFit(mean, sigma);

  if(rangeOn_)
    landau_->fitTo(*dataSet_, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE),RooFit::Range(rangeName_));
  else
    landau_->fitTo(*dataSet_, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE));

  mpvLandau_ = meanL_->getValV();
  sigmaLandau_ = sigmaL_->getValV();

  fitLandau_ = true;
  fitGauss_ = false;
  fitConv_ = false;

  return landau_;
}

RooFFTConvPdf* UnbinnedFit::toLandauXgauss(const double meanL, const double sigmaL, const double meanG, const double sigmaG){
  if(lxg_)
    delete lxg_;

  SetLandauFit(meanL, sigmaL);
  SetGaussFit(meanG, sigmaG);
  meanG_->setConstant(kTRUE);
  sigmaG_->setConstant(kFALSE);
  
  //LandauXGaussian convolution
  fitname_ = "lxg_";
  lxg_ = new RooFFTConvPdf(fitname_, fitname_, *var_, *landau_, *gauss_);

  if(rangeOn_)
    lxg_->fitTo(*dataSet_, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE),RooFit::Range(rangeName_));
  else
    lxg_->fitTo(*dataSet_, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE));

  mpvLandau_ = meanL_->getValV();
  sigmaLandau_ = sigmaL_->getValV();
  meanGauss_ = meanG_->getValV();
  sigmaGauss_ = sigmaG_->getValV();

  fitLandau_ = false;
  fitGauss_ = false;
  fitConv_ = true;
  
  return lxg_;
}

TCanvas* UnbinnedFit::Plot(const TString name, const int nBins, const double xMin, const double xMax){
  if(fitLandau_)
    return GetPlot(landau_, name, nBins, xMin, xMax);
  if(fitGauss_)
    return GetPlot(gauss_, name, nBins, xMin, xMax);
  if(fitConv_)
    return GetPlot(lxg_, name, nBins, xMin, xMax);
  else{
    std::cout << "No fit found for plotting! Returning nullptr" << std::endl;
    return nullptr;
  }
}

// private methods 
void UnbinnedFit::SetupDataSet(const std::vector<double>& data, const double xMin, const double xMax){

  var_ = new RooRealVar("amp","Signal Amplitude [mV]",xMin,xMax);
  varArgset_ = new RooArgSet(*var_);
  dataSet_ = new RooDataSet("dataset","dataset",*varArgset_);
  
  for(int s = 0; s < data.size(); s++){
    var_->setVal(data[s]);
    dataSet_->add(*varArgset_);
  }
}

void UnbinnedFit::SetLandauFit(const double mean, const double sigma){
  if(meanL_)
    delete meanL_;
  if(sigmaL_)
    delete sigmaL_;
  if(landau_)
    delete landau_;
  
  //Landau fit
  fitname_ = "meanL_";
  meanL_ = new RooRealVar(fitname_,"Landau MPV", mean, -1000000., 1000000.);
  fitname_ = "sigmaL_";
  sigmaL_ = new RooRealVar(fitname_,"Landau Width", sigma, 0., 10000.);
  fitname_ = "landau_";
  landau_ = new RooLandau(fitname_,fitname_,*var_,*meanL_,*sigmaL_);
}

void UnbinnedFit::SetGaussFit(const double mean, const double sigma){
  if(meanG_)
    delete meanL_;
  if(sigmaG_)
    delete sigmaL_;
  if(gauss_)
    delete landau_;
  
  //Gaussian fit
  fitname_ = "meanG_";
  meanG_ = new RooRealVar(fitname_,"Gauss Mean",mean);
  fitname_ = "sigmaG_";
  sigmaG_ = new RooRealVar(fitname_,"Gauss Width",sigma, 0., 10000.);
  fitname_ = "gauss_";
  gauss_ = new RooGaussian(fitname_,fitname_,*var_,*meanG_,*sigmaG_);
}

TCanvas* UnbinnedFit::GetPlot(RooAbsPdf* fitPdf, const TString name, const int nBins, const double xMin, const double xMax){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas("can_"+name,"can_"+name,600,800);

  //can->SetLogy();
  can->SetLeftMargin(0.15);
  can->SetGridx();
  can->SetGridy();

  RooPlot* rplot = var_->frame(xMin,xMax,nBins);
  dataSet_->plotOn(rplot);
  fitPdf->plotOn(rplot);
  fitPdf->plotOn(rplot,RooFit::Components(*fitPdf),RooFit::LineColor(kGreen));
  fitPdf->paramOn(rplot,RooFit::Format("NE", RooFit::AutoPrecision(2)), RooFit::Layout(0.65,0.89,0.89));
  rplot->getAttText()->SetTextSize(0.02);
  TString fitname_ = "fitSignal_"+name;
  rplot->SetTitle(fitname_);

  rplot->Draw();
  CMSmark("");

  return can;
}
