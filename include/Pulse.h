#ifndef Pulse_h
#define Pulse_h

#define SAMPLE_SIZE 2002

#include <iostream>
#include <vector>
#include <array>

using namespace std;

class Pulse{

 public:
  Pulse(const vector<float> pulse, const vector<float> time);
  virtual ~Pulse();
  
  float GetMaxAmp() const;
  float GetMaxTime() const;
  
 private:
  
  int maxIndex_;

  float maxAmplitude_;
  float maxTime_;

  vector<float> pulse_;
  vector<float> time_;
  
  void FindMaximum();
};

#endif
