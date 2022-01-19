#include "Pulse.h"

//Public methods
//Constructor
Pulse::Pulse(const vector<float> pulse, const vector<float> time){

  sampleSize_ = time.size();
  sampleRate_ = sampleSize_/(time[sampleSize_-1]-time[0]);

  pulse_ = pulse;
  time_ = time;
  FindMaximum();
}

//Destructor
Pulse::~Pulse(){}

//Public member functions
float Pulse::GetMaxAmp() const{  
  return maxAmplitude_;
}

float Pulse::GetMaxTime() const{
  return maxTime_;
}

void Pulse::GetCenteredPulse(vector<float> &centeredPulse, vector<float> &shiftedTime){

  int nLeft = -50;
  int nRight = 250;
  int size = nRight-nLeft+1;
  
  centeredPulse.clear();
  shiftedTime.clear();

  for(int i = nLeft; i <= nRight; i++){
    centeredPulse.push_back(-(pulse_[maxIndex_+i])/maxAmplitude_);
    shiftedTime.push_back(time_[maxIndex_+i]-maxTime_);
  }
}

void Pulse::InterpolateTest(){

  float maxAmplitude = -999.;
  float maxTime = -999.;
  
  vector<float> shiftedTime;
  vector<float> centeredPulse;
  GetCenteredPulse(centeredPulse,shiftedTime);

  for(int i = 0; i < centeredPulse.size(); i++){
    if(maxAmplitude < centeredPulse[i]){
      maxAmplitude = centeredPulse[i];
      maxTime = shiftedTime[i];
    }
  }
  
  int sizeAroundPeak = 100;

  float initialTime = maxTime-shiftedTime[0]-2.;
  float finalTime = maxTime-shiftedTime[0]+5.;

  vector<float> timeAroundPeak = PulseTools::LinSpaceVec(initialTime,finalTime, sizeAroundPeak);

  maxAmplitude = -999.;
  maxTime = -999.;
  
  for(int i = 0; i < sizeAroundPeak; i++){
    float amplitude = PulseTools::GetInterpolatedPoint(centeredPulse, timeAroundPeak[i], sampleRate_);
    if(maxAmplitude < amplitude){
      maxAmplitude = amplitude;
      maxTime = timeAroundPeak[i];
    }
  }
  
  int interpolationSize = 500;
  float pulseStart = -10.;
  float pulseEnd = 50.;

  vector<float> interpolationTime = PulseTools::LinSpaceVec(pulseStart, pulseEnd, interpolationSize);
  vector<float> interpolatedPulse;
  
  for(int i = 0; i < interpolationSize; i++){
    float amplitude = (PulseTools::GetInterpolatedPoint(centeredPulse, interpolationTime[i]-pulseStart, sampleRate_))/maxAmplitude;
    interpolatedPulse.push_back(amplitude);
  }
  
  TCanvas can(TString("can"),TString("can"),600,800);
  can.SetLeftMargin(0.12);

  TGraph interpolationGraph(interpolationSize,&interpolationTime[0],&interpolatedPulse[0]);

  interpolationGraph.GetXaxis()->SetTitle("time [ns]");
  interpolationGraph.GetYaxis()->SetTitle("Amplitude[mV]");
  interpolationGraph.SetLineWidth(2);
  interpolationGraph.SetLineColor(kBlack);

  interpolationGraph.Draw("AL");

  can.SaveAs(TString("plots/Pulse_1_interpolated.pdf"));
  can.Close();
}

void Pulse::GetInterpolatedPulse(const int interpolationSize, vector<float> &interpolatedPulse, vector<float> &interpolationTime){

  interpolatedPulse.clear();
  interpolationTime.clear();

  /*  
  vector<float> shiftedTime;
  vector<float> centeredPulse;

  GetCenteredPulse(centeredPulse,shiftedTime);

  int size = centeredPulse.size();
  float ampMax = -999.;
  float minTime = shiftedTime[0];
  float maxTime = shiftedTime[size-1];
  interpolationTime = PulseTools::LinSpaceVec(minTime,maxTime,interpolationSize);

  for(int i = 0; i < interpolationSize; i++){
    interpolatedPulse.push_back(PulseTools::GetInterpolatedPoint(centeredPulse, (interpolationTime[i]-minTime), sampleRate_));
    if(ampMax < interpolatedPulse[i])
      ampMax = interpolatedPulse[i];
  }

  for(int i = 0; i < interpolationSize; i++)
    interpolatedPulse[i] /= ampMax;
  */
  
  //Find peak amplitude and time by using fine interpolation near the sample maximum
  int sizeAroundPeak = 100;

  float initialTime = maxTime_-time_[0]-2.;
  float finalTime   = maxTime_-time_[0]+5.;

  vector<float> timeAroundPeak = PulseTools::LinSpaceVec(initialTime,finalTime, sizeAroundPeak);

  int maxIndex;
  float maxAmplitude = -999.;
  float maxTime = -999.;

  for(int i = 0; i < sizeAroundPeak; i++){
    float amplitude = -PulseTools::GetInterpolatedPoint(pulse_, timeAroundPeak[i], sampleRate_);
    if(maxAmplitude < amplitude){
      maxAmplitude = amplitude;
      maxTime = timeAroundPeak[i];
    }
  }

  //make interpolation over 60 ns window around the maximum amplitude 
  float pulseStart = -10.;
  float pulseEnd = 50.;

  interpolationTime = PulseTools::LinSpaceVec(pulseStart, pulseEnd, interpolationSize);

  for(int i = 0; i < interpolationSize; i++){
    float amplitude = (-PulseTools::GetInterpolatedPoint(pulse_, interpolationTime[i]+maxTime, sampleRate_))/maxAmplitude;
    interpolatedPulse.push_back(amplitude);
  }
}

void Pulse::GetCDFpulse(const int interpolationSize, vector<float> &CDFpulse, vector<float> &CDFtime){

  CDFpulse.clear();
  CDFtime.clear();
  
  vector<float> interpolatedPulse;

  GetInterpolatedPulse(interpolationSize, interpolatedPulse, CDFtime);

  float stepSize = CDFtime[1] - CDFtime[0];

  float sum = 0;
  float integral = PulseTools::Integral(interpolatedPulse,stepSize);

  for(int i = 0; i < interpolationSize; i++){
    if(interpolatedPulse[i] > 0){
      sum += interpolatedPulse[i]*stepSize/integral;
    }

    CDFpulse.push_back(sum);
  }
}

void Pulse::PlotCenterPulse(const string name){

  vector<float> shiftedTime;
  vector<float> centeredPulse;

  GetCenteredPulse(centeredPulse,shiftedTime);
  
  int interpolationSize = 500;
  vector<float> interpolatedPulse;
  vector<float> interpolationTime;
  
  GetInterpolatedPulse(interpolationSize,interpolatedPulse,interpolationTime);
  //GetCDFpulse(interpolationSize,interpolatedPulse,interpolationTime);
  int size = centeredPulse.size();
  
  TCanvas can(TString("can_"+name),TString("can_"+name),600,800);
  can.SetLeftMargin(0.12);
  
  TGraph pulseGraph(size,&shiftedTime[0],&centeredPulse[0]);
  TGraph interpolationGraph(interpolationSize,&interpolationTime[0],&interpolatedPulse[0]);

  pulseGraph.SetMarkerColor(kBlack);
  pulseGraph.SetMarkerStyle(8);
  pulseGraph.SetMarkerSize(.5);

  interpolationGraph.GetXaxis()->SetTitle("time [ns]");
  interpolationGraph.GetYaxis()->SetTitle("Amplitude[mV]");
  interpolationGraph.SetLineWidth(2);
  interpolationGraph.SetLineColor(kGreen+2);
  interpolationGraph.SetTitle(TString(name));
  
  interpolationGraph.Draw("AL");
  pulseGraph.Draw("sameP");

  can.SaveAs(TString("plots/"+name+".pdf"));
  can.Close();
}

//Private Methods
void Pulse::FindMaximum(){

  int index = 999;
  float ampMin = 999.;
  float tMin = 999.;

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

