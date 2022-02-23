#ifndef DataReaderH5_h
#define DataReaderH5_h

#define SAMPLE_SIZE 2002
#define N_TRIG 1000

#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include "H5Cpp.h"

#include "PulseList.h"
#include "PulseTools.h"

class DataReaderH5{

 public:
  DataReaderH5(const std::string filePath);
  virtual ~DataReaderH5();

  //Access data in data reader object
  int GetNtrig();
  double GetXpos();
  double GetYpos();
  double GetRate(const int channel);

  std::string GetFileName();
  
  PulseList GetPulseList(const int channel);
  
  //Print summary of loaded file
  void PrintInfo();

private:

  int nTrig_ = N_TRIG;

  //arrays for storing data from datasets
  double x_[1];//platform x position
  double y_[1];//platform y position
  double rate_[1];//trigger rate (different per channel)
  //waveform raw data
  int samples_[N_TRIG][SAMPLE_SIZE];
  double horizOffset_[N_TRIG];
  double horizScale_[N_TRIG];
  double trigOffset_[N_TRIG];
  double trigTime_[N_TRIG];
  double vertOffset_[N_TRIG];
  double vertScale_[N_TRIG];

  //file location
  std::string fileLoc_;
  //H5 file object
  H5::H5File file_;

  //private methods

  //get data from single channel given filename
  void GetFileData();
  void GetChData(const int channel);

  //waveform information
  std::vector<std::vector<double>> GetWaveForms();
  std::vector<std::vector<double>> GetTimeArr();

};

#endif
