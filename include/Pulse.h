#ifndef PULSE_H_
#define PULSE_H_

#include <iostream>
#include <vector>
#include <array>

#include "TH1.h"

#include "Interpolate.h"
#include "PulseTools.h"
#include "Plot.h"

class Pulse{

 public:
  Pulse(const std::vector<double> pulse, const std::vector<double> time);
  virtual ~Pulse();

  double GetSampleRate() const;
  double GetMaxAmp() const;
  double GetStartTime() const;
  double GetEndTime() const;
  double GetMaxTime() const;
  double GetNoiseRMS() const;

  std::vector<double> GetInterpolation(const int size, double start, double end, const double fixedTime = 0.);
  std::vector<double> GetCDF(const int size, double start, double end, const double fixedTime = 0.);
  std::vector<double> GetPulse() const;
  std::vector<double> GetTime() const;

  void GetCenteredPulse(std::vector<double> &centeredPulse, std::vector<double> &shiftedTime);
  void GetCDFpulse(const int interpolationSize, std::vector<double> &interpolatedPulse, std::vector<double> &interpolationTime);
  void PlotCenterPulse(const std::string name);
  
 private:
  
  int maxIndex_;

  double sampleSize_;
  double sampleRate_;
  double maxAmplitude_;
  double startTime_;
  double endTime_;
  double maxTime_;
  double noiseRMS_;
  
  std::vector<double> pulse_;
  std::vector<double> time_;
  std::vector<double> noise_;
  
  void FindMaximum();
};

#endif
