#ifndef PULSETOOLS_H
#define PULSETOOLS_H

#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom.h"
#include "TEventList.h"
#include "TColor.h"

using namespace std;

class PulseTools{

 public:
  
  PulseTools();

  virtual ~PulseTools();

  static int FindMinAbsolute(float* sample, int size);
  static int FindMaxAbsolute(float* sample, int size);
  static float* LinSpace(float xmin, float xmax, float step_size);
  static vector<float> LinSpaceVec(float xmin, float xmax, int size);
  static float Integral(float* sample, float step_size);
  static float Integral(vector<float> sample, float step_size);
  static vector<TH1F> MakeHistArr(int size, float xmin, float xmax, int nbins);
  static float* GetMeanArr(vector<TH1F> pulse_hist);
  static float InterpolateFunc(float* sample, int sample_size, float time);
  static float GetInterpolatedPoint(vector<float> pulse, float time, float frequency);
  static void CalcInterval(TH1F hist, float CI, float& mean, float& low, float& high);

};

#endif


    
