#include "Pulse.h"

//Public methods
//Constructor
Pulse::Pulse(const vector<float> pulse, const vector<float> time){

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

//Private Methods
void Pulse::FindMaximum(){

  int index = 999;
  float ampMin = 999.;
  float tMin = 999.;

  for (int s = 0 ; s < SAMPLE_SIZE; s++){
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
