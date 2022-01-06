#ifndef Pulse_h
#define Pulse_h

#include <iostream>
#include <vector>

using namespace std;

class Pulse{

 public:
  Pulse();
  virtual ~Pulse();
  
  float GetAmplitude(const float *pulse) const;

 private:
  float amplitude_;
  
};

#endif
