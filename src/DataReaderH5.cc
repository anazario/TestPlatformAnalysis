#include "DataReaderH5.h"

using namespace std;
using namespace H5;

DataReaderH5::DataReaderH5(string filePath, const int channel){
  SaveDirInfo(filePath);
  GetScanInfo();
  GetChData(filePath,channel);
};

DataReaderH5::~DataReaderH5(){}

void DataReaderH5::GetChData(string filePath, const int channel){
  
  channel_ = channel;

  H5File file(filePath,H5F_ACC_RDONLY);

  //create datasets to access data in file
  DataSet xDs(file.openDataSet("x"));
  DataSet yDs(file.openDataSet("y"));
  DataSet rateDs(file.openDataSet(("ch"+to_string(channel_)+"_rate").c_str()));
  DataSet samplesDs(file.openDataSet(("ch"+to_string(channel_)+"_samples").c_str()));
  DataSet horizOffsetDs(file.openDataSet(("ch"+to_string(channel_)+"_horiz_offset").c_str()));
  DataSet horizScaleDs(file.openDataSet(("ch"+to_string(channel_)+"_horiz_scale").c_str()));
  DataSet trigOffsetDs(file.openDataSet(("ch"+to_string(channel_)+"_trig_offset").c_str()));
  DataSet trigTimeDs(file.openDataSet(("ch"+to_string(channel_)+"_trig_time").c_str()));
  DataSet vertOffsetDs(file.openDataSet(("ch"+to_string(channel_)+"_vert_offset").c_str()));
  DataSet vertScaleDs(file.openDataSet(("ch"+to_string(channel_)+"_vert_scale").c_str()));

  //read data from dataset
  xDs.read(x_,PredType::NATIVE_FLOAT);
  yDs.read(y_,PredType::NATIVE_FLOAT);
  rateDs.read(rate_,PredType::NATIVE_FLOAT);
  samplesDs.read(samples_,PredType::NATIVE_INT);
  horizOffsetDs.read(horizOffset_,PredType::NATIVE_FLOAT);
  horizScaleDs.read(horizScale_,PredType::NATIVE_FLOAT);
  trigOffsetDs.read(trigOffset_,PredType::NATIVE_FLOAT);
  trigTimeDs.read(trigTime_,PredType::NATIVE_FLOAT);
  vertOffsetDs.read(vertOffset_,PredType::NATIVE_FLOAT);
  vertScaleDs.read(vertScale_,PredType::NATIVE_FLOAT);

  file.close();
}

int DataReaderH5::GetNtrig(){
  return nTrig_;
}
int DataReaderH5::GetXmaxIdx(){
  return xMaxIdx_;
}
int DataReaderH5::GetYmaxIdx(){
  return yMaxIdx_;
}

float DataReaderH5::GetXpos(){
  return x_[0];
}

float DataReaderH5::GetYpos(){
  return y_[0];
}

float DataReaderH5::GetRate(){
  return rate_[0];
}

vector<vector<float>> DataReaderH5::GetWaveForms(){

  vector<vector<float>> waveForms;
  vector<float> tempVec;
  
  for(int i = 0; i < nTrig_; i++){
    for(int j = 0; j < SAMPLE_SIZE; j++){
      tempVec.push_back(vertOffset_[i]+samples_[i][j]*vertScale_[i]);
    }
    waveForms.push_back(tempVec);
    tempVec.clear();
  }

  return waveForms;
}

vector<vector<float>> DataReaderH5::GetTimeArr(){

  float tMin = -999.;
  float	tMax = -999.;
  vector<vector<float>> time;
  vector<float> t;
  
  for(int i = 0; i < nTrig_; i++){
    tMin = trigOffset_[i]*1e9;
    tMax = (horizScale_[i]*(SAMPLE_SIZE-1)-trigOffset_[i])*1e9;
    t = PulseTools::LinSpaceVec(tMin,tMax,SAMPLE_SIZE);
    time.push_back(t);
  }
  return time;
}

void DataReaderH5::PrintInfo(){
  cout << endl;
  cout << "General Details of " << dirName_ << ": " << endl;
  cout << "---------------------------------------" << endl;
  cout << "\tNumber of Triggers: " << scanInfo_["nTrig"] << endl;
  cout << "\tMaximum index of x: " << scanInfo_["xMaxIdx"] << endl;
  cout << "\tMaximum index of y: " << scanInfo_["yMaxIdx"] << endl;
  cout << "\tx initial position: " << scanInfo_["xMin"] << " mm" << endl;
  cout << "\ty initial position: " << scanInfo_["yMin"] << " mm" << endl;
  cout << "\tx final position: " << scanInfo_["xMax"] << " mm" << endl;
  cout << "\ty final position: " << scanInfo_["yMax"] << " mm" << endl;
  cout << endl;
  cout << "Details of loaded file: " << endl;
  cout << "---------------------------------------" << endl;
  cout << "\tFile name: " << fileName_ << endl;
  cout << "\tActive channel: " << channel_ << endl;
  cout << "\tLocation: x = " << x_[0] << " mm, y = " << y_[0] << " mm" << endl;
  cout << "\tTrigger rate: " << rate_[0] << " Hz" << endl;
  cout << endl;
}

//private methods
void DataReaderH5::GetScanInfo(){

  float value;
  fstream file;
  string  fileName,key;
  char c;

  fileName = (dataLoc_+dirName_+"/"+dirName_+"_info.json").c_str();
  file.open(fileName.c_str());

  if(file >> c && c == '{'){
    while (file >> c)
      {
        if (c == '\"'){
          while (file >> c && c != '\"')
            key += c;
        }
        if(c == ':' && file >> value){
          scanInfo_.insert(pair<string,float>(key,value));
          key = "";
        }
      }
  }
  //save map to object member vars
  nTrig_   = int(scanInfo_["nTrig"]);
  xMaxIdx_ = int(scanInfo_["xMaxIdx"]);
  yMaxIdx_ = int(scanInfo_["yMaxIdx"]);
  xMin_    = scanInfo_["xMin"];
  xMax_    = scanInfo_["xMax"];
  yMin_    = scanInfo_["yMin"];
  yMax_    = scanInfo_["yMax"];
}

void DataReaderH5::SaveDirInfo(string inPath){

  pathName_ = inPath;
  SplitFolder(pathName_,fileName_);
  
  dataLoc_ = pathName_;
  Replace(dataLoc_,"/hdf5","");

  SplitFolder(dataLoc_,dirName_);
}

void DataReaderH5::SplitFolder(string& fullPath, string& innerMostName){

  int length = fullPath.length();
  int slashIdx = fullPath.find_last_of('/');

  string path;
  
  if(slashIdx+1 == length){
    fullPath.replace(slashIdx,1,"");
    slashIdx = fullPath.find_last_of('/');
    length--;
  }

  for (int i = 0; i < slashIdx+1; i++)
    path += fullPath[i];
  for (int i = slashIdx+1; i < length; i++)
    innerMostName += fullPath[i];

  fullPath = path;
}

bool DataReaderH5::Replace(string& str, const string& from, const string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}
