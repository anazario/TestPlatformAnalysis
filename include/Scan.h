#ifndef Scan_h
#define Scan_h

#include <iostream>
#include <vector>

#include "PulseTools.h"
#include "DataReaderH5.h"

class Scan{
  
 public:

  Scan(const std::string filePath);
  virtual ~Scan();

  int GetNtrig() const;
  int GetXmaxIdx() const;
  int GetYmaxIdx() const;
  int GetTotalFiles() const;

  std::string GetPathToData() const;
  std::vector<std::string> GetFileList() const;
  
  void PrintInfo();

  static std::map<std::string,double> DictToMap(const std::string fileName);
  
 private:

  int nTrig_;
  int xMaxIdx_;
  int yMaxIdx_;
  double xMin_;
  double xMax_;
  double yMin_;
  double yMax_;

  std::string fileName_;
  std::string dirName_;
  std::string pathName_;
  std::string dataLoc_;

  std::vector<std::string> fileList_;
  
  std::map<std::string,double> scanInfo_;

  void GetScanInfo();
  void SaveDirInfo(std::string inPath);
};

#endif
