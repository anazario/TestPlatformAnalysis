#include "PulseVariation.h"

//Constructors
PulseVariation::PulseVariation(const type){}
PulseVariation::PulseVariation(const PulseList& pulseList, const type){

  collectionSize_ = pulseList.size();

  
}

PulseVariation::PulseVariation(const std::vector<TH1D> histCollection, const std::vector<double> timeAxis){

  collectionSize_ = timeAxis.size();
  timeAxis_ = timeAxis;
  histCollection_ = histCollection;
  collectionMean_ = GetMeanPulse();
}

//Destructor
PulseVariation::~PulseVariation(){}

//public methods
std::vector<double> PulseVariation::GetMeanPulse() const{

  int size = histCollection_.size();
  std::vector<double> meanPulse;
  
  for(int h = 0; h < size; h++)
    meanPulse.push_back(histCollection_[h].GetMean());

  return meanPulse;
}

void PulseVariation::GetErrors(std::vector<double>& lowErrors, std::vector<double>& highErrors, double CI){

  int size = histCollection_.size();
  double mean;
  double lowError  = -999.;
  double highError = -999.;

  for(int i = 0; i < size; i++){
    CalcInterval(histCollection_[i],CI,mean,lowError,highError);
    lowErrors.push_back(lowError);
    highErrors.push_back(highError);
  }
}

void PulseVariation::PlotMeanErrors(const TString name, const TString ylabel){

  std::vector<double> lowErr68;
  std::vector<double> highErr68;

  std::vector<double> lowErr95;
  std::vector<double> highErr95;

  std::vector<double> lowErr99;
  std::vector<double> highErr99;

  GetErrors(lowErr68,highErr68,0.68);
  GetErrors(lowErr95,highErr95,0.95);
  GetErrors(lowErr99,highErr99,0.99);

  TCanvas cv("cv","cv",600,800);
  TGraph mean(collectionSize_,&timeAxis_[0],&collectionMean_[0]);
  TGraphAsymmErrors tg68(collectionSize_,&timeAxis_[0],&collectionMean_[0],0,0,&lowErr68[0],&highErr68[0]);
  TGraphAsymmErrors tg95(collectionSize_,&timeAxis_[0],&collectionMean_[0],0,0,&lowErr95[0],&highErr95[0]);
  TGraphAsymmErrors tg99(collectionSize_,&timeAxis_[0],&collectionMean_[0],0,0,&lowErr99[0],&highErr99[0]);

  //Aesthetics
  int style = 3144;
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(kAvocado);
  gStyle->SetHatchesLineWidth(2);

  mean.SetLineWidth(3);
  mean.SetLineStyle(1);
  tg68.SetFillStyle(style);
  tg95.SetFillStyle(style);
  tg99.SetFillStyle(style);
  tg99.GetXaxis()->SetTitle("time (ns)");
  tg99.GetYaxis()->SetTitle(ylabel);
  tg99.GetYaxis()->SetTitleOffset(1.25);
  tg99.SetTitle("99% CL");
  tg95.SetTitle("95% CL");
  tg68.SetTitle("68% CL");
  mean.SetTitle("Mean Pulse");
  tg99.Draw("a4l PFC");
  tg95.Draw("samel4 PFC");
  tg68.Draw("samel4 PFC");
  mean.Draw("sameL");

  //Legend
  //TLegend leg(0.7,0.75,0.89,0.89);
  TLegend leg(0.11,0.75,0.3,0.89);
  /*
  if(opt == "pulse")
    leg = new TLegend(0.7,0.75,0.89,0.89);
  else
    leg = new TLegend(0.20,0.7,0.37,0.84);
  */
  leg.AddEntry(&mean,"Mean Pulse","l");
  leg.AddEntry(&tg68,"68% CL","f");
  leg.AddEntry(&tg95,"95% CL","f");
  leg.AddEntry(&tg99,"99% CL","f");
  leg.Draw("same");

  CMSmark("");
  gPad->SetGrid(1, 1); gPad->Update();

  cv.SaveAs(name+".pdf");
  cv.SaveAs(name+".gif");
  cv.Close();
}

void PulseVariation::PlotHistograms(){

  if(histCollection_.size() > 0){

    TCanvas cv("cv","cv",600,800);
    histCollection_[0].Draw();
    cv.SaveAs("plots/histCollection/pulse0.pdf");
    for(int i = 1; i < histCollection_.size(); i++){
      histCollection_[i].Draw();
      cv.SaveAs(Form("plots/histCollection/pulse%i.pdf",i));
    }
  }
  else{
    std::cout << "Histogram vector is empty!" << std::endl;
    return;
  }
}
