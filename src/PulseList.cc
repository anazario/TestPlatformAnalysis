#include "PulseList.h"

PulseList::PulseList(){};

PulseList::PulseList(const std::vector<Pulse*>& pulseList){

  listSize_ = pulseList.size();
  
  for(auto pulse : pulseList)
    *this += pulse;
}

PulseList::~PulseList(){}

std::vector<double> PulseList::GetAmpDistribution() const{

  std::vector<double> ampDistribution;
  for(auto pulse : *this)
    ampDistribution.push_back(pulse->GetMaxAmp()*1e3);

  return ampDistribution;
}

PulseList& PulseList::operator += (Pulse* pulse){
  this->push_back(pulse);
  return *this;
}

PulseList& PulseList::operator += (const PulseList& pulseList){
  for(auto pulse : pulseList)
    this->push_back(pulse); 
  return *this;
}
