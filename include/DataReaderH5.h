#ifndef DataReaderH5_h
#define DataReaderH5_h

#define SAMPLE_SIZE 2002

#include <iostream>
#include <map>
#include <fstream>
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
  float** GetWaveForms();
  float** GetTimeArr();

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
  int **samples_;
  float *horizOffset_;
  float *horizScale_;
  float *trigOffset_;
  float *trigTime_;
  float *vertOffset_;
  float *vertScale_;

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
  void InitMembers();
  void GetScanInfo();
  void SaveDirInfo(string inPath);
  void SplitFolder(string& fullPath, string& innerMostName);
  bool Replace(string& str, const string& from, const string& to);

};

#endif
