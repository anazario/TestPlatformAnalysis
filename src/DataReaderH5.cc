#include "DataReaderH5.h"

DataReaderH5::DataReaderH5(const std::string filePath){
  fileLoc_ = filePath;
  H5::H5File file(fileLoc_,H5F_ACC_RDONLY);
  file_ = file;
  GetFileData();
}

DataReaderH5::~DataReaderH5(){
  file_.close();
}

//public methods

int DataReaderH5::GetNtrig(){
  return nTrig_;
}

double DataReaderH5::GetXpos(){
  return x_[0];
}

double DataReaderH5::GetYpos(){
  return y_[0];
}

double DataReaderH5::GetRate(const int channel){
  H5::DataSet rateDs(file_.openDataSet(("ch"+std::to_string(channel)+"_rate").c_str()));
  rateDs.read(rate_,H5::PredType::NATIVE_DOUBLE);
  return rate_[0]; 
}

std::string DataReaderH5::GetFileName(){
  std::string fileName;
  SplitFolder(fileLoc_,fileName);
  return GetFileName();  
}

PulseList DataReaderH5::GetPulseList(const int channel){
  GetChData(channel);
  std::vector<std::vector<double>> waveForms = GetWaveForms();
  std::vector<std::vector<double>> times = GetTimeArr();

  PulseList pulseList;

  for(int n = 0; n < waveForms.size(); n++){
    Pulse* pulse = new Pulse(waveForms[n],times[n]);
    pulseList += pulse;
  }
  
  return pulseList;
}

void DataReaderH5::PrintInfo(){
  std::cout << std::endl;
  std::cout << "Details of loaded file: " << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "\tFile Location: " << fileLoc_ << std::endl;
  std::cout << "\tLocation: x = " << x_[0] << " mm, y = " << y_[0] << " mm" << std::endl;
  std::cout << "\tTrigger rate: " << rate_[0] << " Hz" << std::endl;
  std::cout << std::endl;
}

//private methods

void DataReaderH5::GetFileData(){
  H5::DataSet xDs(file_.openDataSet("x"));
  H5::DataSet yDs(file_.openDataSet("y"));
  
  xDs.read(x_,H5::PredType::NATIVE_DOUBLE);
  yDs.read(y_,H5::PredType::NATIVE_DOUBLE);
}

void DataReaderH5::GetChData(const int channel){
  using namespace H5;

  //H5File file(fileLoc_,H5F_ACC_RDONLY);

  //create datasets to access data in file
  DataSet samplesDs(file_.openDataSet(("ch"+std::to_string(channel)+"_samples").c_str()));
  DataSet horizOffsetDs(file_.openDataSet(("ch"+std::to_string(channel)+"_horiz_offset").c_str()));
  DataSet horizScaleDs(file_.openDataSet(("ch"+std::to_string(channel)+"_horiz_scale").c_str()));
  DataSet trigOffsetDs(file_.openDataSet(("ch"+std::to_string(channel)+"_trig_offset").c_str()));
  DataSet trigTimeDs(file_.openDataSet(("ch"+std::to_string(channel)+"_trig_time").c_str()));
  DataSet vertOffsetDs(file_.openDataSet(("ch"+std::to_string(channel)+"_vert_offset").c_str()));
  DataSet vertScaleDs(file_.openDataSet(("ch"+std::to_string(channel)+"_vert_scale").c_str()));

  //read data from dataset
  samplesDs.read(samples_,PredType::NATIVE_INT);
  horizOffsetDs.read(horizOffset_,PredType::NATIVE_DOUBLE);
  horizScaleDs.read(horizScale_,PredType::NATIVE_DOUBLE);
  trigOffsetDs.read(trigOffset_,PredType::NATIVE_DOUBLE);
  trigTimeDs.read(trigTime_,PredType::NATIVE_DOUBLE);
  vertOffsetDs.read(vertOffset_,PredType::NATIVE_DOUBLE);
  vertScaleDs.read(vertScale_,PredType::NATIVE_DOUBLE);
}

std::vector<std::vector<double>> DataReaderH5::GetWaveForms(){

  std::vector<std::vector<double>> waveForms;

  for(int i = 0; i < nTrig_; i++){
    std::vector<double> pulse;
    for(int j = 0; j < SAMPLE_SIZE; j++){
      pulse.push_back(-vertOffset_[i]+samples_[i][j]*vertScale_[i]);
    }
    waveForms.push_back(pulse);
  }

  return waveForms;
}

std::vector<std::vector<double>> DataReaderH5::GetTimeArr(){

  double tMin = -999.;
  double tMax = -999.;
  std::vector<std::vector<double>> time;

  for(int i = 0; i < nTrig_; i++){
    tMin = trigOffset_[i]*1e9;
    tMax = (horizScale_[i]*(SAMPLE_SIZE-1)-trigOffset_[i])*1e9;
    time.push_back(LinSpaceVec(tMin,tMax,SAMPLE_SIZE));
  }
  return time;
}

