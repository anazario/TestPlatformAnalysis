#ifndef PulseList_h
#define PulseList_h

#include <iostream>
#include <vector>

#include "Pulse.h"

class PulseList : public std::vector<Pulse*>{

 public:
  PulseList();
  PulseList(const std::vector<Pulse*>& pulseList);
  virtual ~PulseList();

  std::vector<double> GetAmpDistribution() const;
  
  PulseList& operator += (Pulse* pulse);
  PulseList& operator += (const PulseList& pulseList);
  
 private:

  int listSize_;

};

#endif
