#ifndef DataReaderH5_h
#define DataReaderH5_h

#define SAMPLE_SIZE 2002
#define N_TRIG 1000

#include <iostream>
#include <map>
#include <fstream>
#include <vector>
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

  //Access data in data reader object
  int GetNtrig();
  int GetXmaxIdx();
  int GetYmaxIdx();
  float GetXpos();
  float GetYpos();
  float GetRate();
  vector<vector<float>> GetWaveForms();
  vector<vector<float>> GetTimeArr();

  //Print summary of scan and specifics of loaded file
  void PrintInfo();

private:

  //data directory name
  string fileName_;
  string dirName_;
  string pathName_;
  string dataLoc_;
  
  //int nTrig_;
  int channel_;

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

  //map with general scan info
  map<string,float> scanInfo_;

  //variables from json file
  int nTrig_;
  int xMaxIdx_;
  int yMaxIdx_;
  float xMin_;
  float xMax_;
  float	yMin_;
  float	yMax_;
  
  //private methods
  void GetScanInfo();
  void SaveDirInfo(string inPath);
  void SplitFolder(string& fullPath, string& innerMostName);
  bool Replace(string& str, const string& from, const string& to);

};

#endif
