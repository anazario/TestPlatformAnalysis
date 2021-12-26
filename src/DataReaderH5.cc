#include "DataReaderH5.h"

using namespace std;
using namespace H5;

DataReaderH5::DataReaderH5(string filePath, const int channel){
  SaveDirInfo(filePath);
  GetScanInfo();
  GetChData(filePath,channel);
};

DataReaderH5::~DataReaderH5(){
  delete xDs_;
  delete yDs_;
  delete rateDs_;
  delete samplesDs_;
  delete horizOffsetDs_;
  delete horizScaleDs_;
  delete trigOffsetDs_;
  delete trigTimeDs_;
  delete vertOffsetDs_;
  delete vertScaleDs_;
  file_->close();
  delete file_;

  delete samples_;
  delete horizOffset_;
  delete horizScale_;
  delete trigOffset_;
  delete trigTime_;
  delete vertOffset_;
  delete vertScale_;
}

void DataReaderH5::GetChData(string filePath, const int channel){
  
  channel_ = channel;

  InitMembers();
  
  file_ = new H5File(filePath,H5F_ACC_RDONLY);

  //create datasets to access data in file
  xDs_           = new DataSet(file_->openDataSet("x"));
  yDs_           = new DataSet(file_->openDataSet("y"));
  rateDs_        = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_rate").c_str()));
  samplesDs_     = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_samples").c_str()));
  horizOffsetDs_ = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_horiz_offset").c_str()));
  horizScaleDs_  = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_horiz_scale").c_str()));
  trigOffsetDs_  = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_trig_offset").c_str()));
  trigTimeDs_    = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_trig_time").c_str()));
  vertOffsetDs_  = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_vert_offset").c_str()));
  vertScaleDs_   = new DataSet(file_->openDataSet(("ch"+to_string(channel_)+"_vert_scale").c_str()));

  //read data from dataset
  xDs_->read(x_,PredType::NATIVE_FLOAT);
  yDs_->read(y_,PredType::NATIVE_FLOAT);
  rateDs_->read(rate_,PredType::NATIVE_FLOAT);
  samplesDs_->read(*samples_,PredType::NATIVE_INT);
  horizOffsetDs_->read(horizOffset_,PredType::NATIVE_FLOAT);
  horizScaleDs_->read(horizScale_,PredType::NATIVE_FLOAT);
  trigOffsetDs_->read(trigOffset_,PredType::NATIVE_FLOAT);
  trigTimeDs_->read(trigTime_,PredType::NATIVE_FLOAT);
  vertOffsetDs_->read(vertOffset_,PredType::NATIVE_FLOAT);
  vertScaleDs_->read(vertScale_,PredType::NATIVE_FLOAT);
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

float** DataReaderH5::GetWaveForms(){

  float **waveForms = new float*[nTrig_];
  for(int i = 0; i < nTrig_; i++){
    waveForms[i] = new float[SAMPLE_SIZE]; 
    for(int j = 0; j < SAMPLE_SIZE; j++){
      waveForms[i][j] = -vertOffset_[i]+samples_[i][j]*vertScale_[i];
    }
  }
  
  return waveForms;
}

float** DataReaderH5::GetTimeArr(){

  float tMin = -999.;
  float	tMax = -999.;
  float **time = new float*[nTrig_];

  for(int i = 0; i < nTrig_; i++){
    time[i] = new float[SAMPLE_SIZE];
    tMin = trigOffset_[i]*1e9;
    tMax = (horizScale_[i]*(SAMPLE_SIZE-1)-trigOffset_[i])*1e9;
    time[i] = PulseTools::LinSpace(tMin,tMax,1./SAMPLE_SIZE);
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

void DataReaderH5::InitMembers(){
  samples_ = new int*[nTrig_];
  for (int i = 0; i < nTrig_; i++)
    samples_[i] = new int[SAMPLE_SIZE];
  
  horizOffset_ = new float[nTrig_];
  horizScale_ = new float[nTrig_];
  trigOffset_ = new float[nTrig_];
  trigTime_ = new float[nTrig_];
  vertOffset_ = new float[nTrig_];
  vertScale_ = new float[nTrig_];
}

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
