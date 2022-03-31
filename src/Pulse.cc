#include "Pulse.h"

//Constructor
Pulse::Pulse(const std::vector<double> pulse, const std::vector<double> time){

  sampleSize_ = time.size();
  timeLowEdge_ = time[0];
  timeHighEdge_ = time[sampleSize_-1];
  sampleRate_ = sampleSize_/(timeHighEdge_ - timeLowEdge_);

  pulse_ = pulse;
  time_ = time;
  FindMaximum();
}

//Destructor
Pulse::~Pulse(){}

//Public member functions
int Pulse::GetMaxIndex() const{
  return maxIndex_;
}

double Pulse::GetSampleRate() const{
  return sampleRate_;
}

double Pulse::GetMaxAmp() const{  
  return maxAmplitude_;
}

double Pulse::GetTimeLowEdge() const{
  return timeLowEdge_;
}

double Pulse::GetTimeHighEdge() const{
  return timeHighEdge_;
}

double Pulse::GetMaxTime() const{
  return maxTime_;
}

double Pulse::GetNoiseRMS() const{

  std::vector<double> noise;

  for(int i = 10; i < 710; i++)
    noise.push_back(pulse_[i]);
  
  return GetSampleRMS(noise);
}

std::vector<double> Pulse::GetPulse() const{
  return pulse_;
}

std::vector<double> Pulse::GetTime() const{
  return time_;
}

std::vector<double> Pulse::GetInterpolation(const int size, double start, double end, const double fixedTime){

  start += fixedTime;
  end += fixedTime;
  
  std::vector<double> interpolatedPulse = GetInterpolatedPulse(size, start, end, *this);
  NormalizeVec(interpolatedPulse);

  return interpolatedPulse;
}

std::vector<double> Pulse::GetCDF(const int size, double start, double end, const double fixedTime){

  double stepSize = (end - start)/double(size);

  std::vector<double> interpolatedPulse = GetInterpolation(size, start, end, fixedTime);
  std::vector<double> CDFpulse;

  double sum = 0;
  double integral = Integral(interpolatedPulse,stepSize);
  
  for(int i = 0; i < size; i++){

    if(interpolatedPulse[i] > 0){
      sum += interpolatedPulse[i]*stepSize/integral;
    }

    CDFpulse.push_back(sum);
  }
  return CDFpulse;
}

std::vector<double> Pulse::GetInterpolatedPulse(const int interpolationSize, const double pulseStart, const double pulseEnd,
						const Pulse& originalPulse){

  double sampleRate = originalPulse.GetSampleRate();

  std::vector<double> interpolatedPulse;
  std::vector<double> interpolationTime = LinSpaceVec(pulseStart,pulseEnd,interpolationSize);

  double timeDiff = originalPulse.GetTimeLowEdge();

  for(int i = 0; i < interpolationSize; i++){
    double amplitude = -GetInterpolatedPoint(originalPulse.GetPulse(), interpolationTime[i]-timeDiff, sampleRate);
    interpolatedPulse.push_back(amplitude);
  }

  return interpolatedPulse;
}

void Pulse::PlotCenterPulse(const std::string name){

  std::vector<double> shiftedTime;
  std::vector<double> centeredPulse;

  GetCenteredPulse(centeredPulse,shiftedTime);

  int size = centeredPulse.size();
  
  TCanvas can(TString("can_"+name),TString("can_"+name),600,800);
  can.SetLeftMargin(0.12);

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(11111111);
  TGraph pulseGraph(size,&shiftedTime[0],&centeredPulse[0]);

  pulseGraph.SetMarkerColor(kBlack);
  pulseGraph.SetMarkerStyle(8);
  pulseGraph.SetMarkerSize(.5);
  pulseGraph.GetXaxis()->SetTitle("time [ns]");
  pulseGraph.GetYaxis()->SetTitle("Amplitude [mV]");

  pulseGraph.Draw("AP");

  can.SaveAs(TString("plots/"+name+".pdf"));
  can.Close();
}

//Private Methods
void Pulse::GetCenteredPulse(std::vector<double> &centeredPulse, std::vector<double> &shiftedTime){

  int nLeft = -50;
  int nRight = 150;
  int size = nRight-nLeft+1;

  centeredPulse.clear();
  shiftedTime.clear();

  for(int i = nLeft; i <= nRight; i++){
    centeredPulse.push_back(-(pulse_[maxIndex_+i])/maxAmplitude_);
    shiftedTime.push_back(time_[maxIndex_+i]-maxTime_);
  }
}

void Pulse::FindMaximum(){

  int index = 999;
  double ampMin = 999.;
  double tMin = 999.;

  for (int s = 0 ; s < sampleSize_; s++){
    if (pulse_[s] < ampMin){
      index = s;
      ampMin = pulse_[s];
      tMin = time_[s];
    }
  }

  maxIndex_ = index;
  maxAmplitude_ = -ampMin;
  maxTime_ = tMin;  
}

