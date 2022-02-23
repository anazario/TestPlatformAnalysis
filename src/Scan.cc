#include "Scan.h"

Scan::Scan(const std::string filePath){
  SaveDirInfo(filePath);
  GetScanInfo();

  fileList_ = GetFileList();  
}

Scan::~Scan(){}

int Scan::GetNtrig() const{
  return nTrig_;
}
int Scan::GetXmaxIdx() const{
  return xMaxIdx_;
}
int Scan::GetYmaxIdx() const{
  return yMaxIdx_;
}
int Scan::GetTotalFiles() const{
  return xMaxIdx_*yMaxIdx_;
}

std::string Scan::GetPathToData() const{
  return (pathName_+"hdf5/").c_str();
}

std::vector<std::string> Scan::GetFileList() const{

  std::vector<std::string> fileList;

  for(int x = 0; x < xMaxIdx_; x++){
    for(int y = 0; y < yMaxIdx_; y++){
      fileList.push_back(("x"+std::to_string(x)+"_y"+std::to_string(y)+".hdf5").c_str());
    }
  }
  return fileList;
}

void Scan::PrintInfo(){
  std::cout << std::endl;
  std::cout << "General Details of " << dirName_ << ": " << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "\tNumber of Triggers: " << scanInfo_["nTrig"] << std::endl;
  std::cout << "\tMaximum index of x: " << scanInfo_["xMaxIdx"] << std::endl;
  std::cout << "\tMaximum index of y: " << scanInfo_["yMaxIdx"] << std::endl;
  std::cout << "\tx initial position: " << scanInfo_["xMin"] << " mm" << std::endl;
  std::cout << "\ty initial position: " << scanInfo_["yMin"] << " mm" << std::endl;
  std::cout << "\tx final position: " << scanInfo_["xMax"] << " mm" << std::endl;
  std::cout << "\ty final position: " << scanInfo_["yMax"] << " mm" << std::endl;
  std::cout << std::endl;
}

std::map<std::string,double> Scan::DictToMap(const std::string fileName){

  double value;
  std::fstream file;
  std::string  key;
  char c;

  std::map<std::string,double> dict;

  file.open(fileName.c_str());

  if(file >> c && c == '{'){
    while (file >> c)
      {
        if (c == '\"'){
          while (file >> c && c != '\"')
            key += c;
        }
        if(c == ':' && file >> value){
          dict.insert(std::pair<std::string,double>(key,value));
          key = "";
        }
      }
  }

  return dict;
}

void Scan::GetScanInfo(){
  scanInfo_ = DictToMap((dataLoc_+dirName_+"/"+dirName_+"_info.json").c_str());

  //save map to object member vars
  nTrig_   = int(scanInfo_["nTrig"]);
  xMaxIdx_ = int(scanInfo_["xMaxIdx"])+1;
  yMaxIdx_ = int(scanInfo_["yMaxIdx"])+1;
  xMin_    = scanInfo_["xMin"];
  xMax_    = scanInfo_["xMax"];
  yMin_    = scanInfo_["yMin"];
  yMax_    = scanInfo_["yMax"];
}

void Scan::SaveDirInfo(std::string inPath){

  pathName_ = inPath;
  SplitFolder(pathName_,fileName_);

  dataLoc_ = pathName_;
  Replace(dataLoc_,"/hdf5","");

  SplitFolder(dataLoc_,dirName_);
}
