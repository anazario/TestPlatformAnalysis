#ifndef PULSE_H_
#define PULSE_H_

//#define SAMPLE_SIZE 2002

#include <iostream>
#include <vector>
#include <array>

#include "PulseTools.h"
#include "Plot.h"

using namespace std;

class Pulse{

 public:
  Pulse(const vector<float> pulse, const vector<float> time);
  virtual ~Pulse();
  
  float GetMaxAmp() const;
  float GetMaxTime() const;

  void InterpolateTest();
  void GetCenteredPulse(vector<float> &centeredPulse, vector<float> &shiftedTime);
  void GetInterpolatedPulse(const int interpolationSize, vector<float> &interpolatedPulse, vector<float> &interpolationTime);
  void GetCDFpulse(const int interpolationSize, vector<float> &interpolatedPulse, vector<float> &interpolationTime);
  void PlotCenterPulse(const string name);
  
 private:
  
  int maxIndex_;

  float sampleSize_;
  float sampleRate_;
  float maxAmplitude_;
  float maxTime_;

  vector<float> pulse_;
  vector<float> time_;

  void FindMaximum();
};

#endif
