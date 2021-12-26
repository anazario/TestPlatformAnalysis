#ifndef DataReaderH5_h
#define DataReaderH5_h

#define N_TRIG 1000
#define SAMPLE_SIZE 2002

#include <iostream>
#include "H5Cpp.h"
#include "PulseTools.h"

using namespace std;
using namespace H5;

class DataReaderH5{

 public:
  DataReaderH5(string fileName, const int channel);
  virtual ~DataReaderH5();

  //get data from single channel given filename
  void GetChData(string fileName, const int channel);
  float GetXpos();
  float GetYpos();
  float GetRate();
  float** GetWaveForms();
  float** GetTimeArr();
  void PrintInfo();

private:

  //int nTrig_;
  int channel_;

  //h5 file
  H5File *file_;

  //h5 dataset objects
  DataSet *xDs_;
  DataSet *yDs_;
  DataSet *rateDs_;
  DataSet *samplesDs_;
  DataSet *horizOffsetDs_;
  DataSet *horizScaleDs_;
  DataSet *trigOffsetDs_;
  DataSet *trigTimeDs_;
  DataSet *vertOffsetDs_;
  DataSet *vertScaleDs_;

  //arrays for storing data from datasets
  float x_[1];
  float y_[1];
  float rate_[1];
  int samples_[N_TRIG][SAMPLE_SIZE];
  float horizOffset_[N_TRIG];
  float horizScale_[N_TRIG];
  float trigOffset_[N_TRIG];
  float trigTime_[N_TRIG];
  float vertOffset_[N_TRIG];
  float vertScale_[N_TRIG];

  //void GetScanInfo();
};

#endif
