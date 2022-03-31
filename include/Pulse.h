#ifndef PULSE_H_
#define PULSE_H_

#include <iostream>
#include <vector>
#include <array>

#include "TH1.h"

#include "PulseTools.h"
#include "Plot.h"

class Pulse{

 public:
  Pulse(const std::vector<double> pulse, const std::vector<double> time);
  virtual ~Pulse();

  int GetMaxIndex() const;
  
  double GetSampleRate() const;
  double GetMaxAmp() const;
  double GetTimeLowEdge() const;
  double GetTimeHighEdge() const;
  double GetMaxTime() const;
  double GetNoiseRMS() const;

  std::vector<double> GetInterpolation(const int size, double start, double end, const double fixedTime = 0.);
  std::vector<double> GetCDF(const int size, double start, double end, const double fixedTime = 0.);
  std::vector<double> GetPulse() const;
  std::vector<double> GetTime() const;

  static std::vector<double> GetInterpolatedPulse(const int interpolationSize, const double pulseStart, const double pulseEnd,
                                                  const Pulse& originalPulse);
  
  void PlotCenterPulse(const std::string name);
  
 private:
  
  int maxIndex_;

  double sampleSize_;
  double sampleRate_;
  double maxAmplitude_;
  double timeLowEdge_;
  double timeHighEdge_;
  double maxTime_;
  double noiseRMS_;
  
  std::vector<double> pulse_;
  std::vector<double> time_;
  std::vector<double> noise_;

  void GetCenteredPulse(std::vector<double> &centeredPulse, std::vector<double> &shiftedTime);
  void FindMaximum();
};

#endif
