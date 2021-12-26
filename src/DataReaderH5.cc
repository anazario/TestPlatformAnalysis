#include "DataReaderH5.h"


using namespace std;
using namespace H5;

DataReaderH5::DataReaderH5(string fileName, const int channel){
  GetChData(fileName,channel);
  //GetScanInfo();
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
}

void DataReaderH5::GetChData(string fileName, const int channel){

  channel_ = channel;
  
  file_ = new H5File(fileName,H5F_ACC_RDONLY);

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
  samplesDs_->read(samples_,PredType::NATIVE_INT);
  horizOffsetDs_->read(horizOffset_,PredType::NATIVE_FLOAT);
  horizScaleDs_->read(horizScale_,PredType::NATIVE_FLOAT);
  trigOffsetDs_->read(trigOffset_,PredType::NATIVE_FLOAT);
  trigTimeDs_->read(trigTime_,PredType::NATIVE_FLOAT);
  vertOffsetDs_->read(vertOffset_,PredType::NATIVE_FLOAT);
  vertScaleDs_->read(vertScale_,PredType::NATIVE_FLOAT);
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

  float **waveForms = new float*[N_TRIG];

  for(int i = 0; i < N_TRIG; i++){
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
  float **time = new float*[N_TRIG];

  for(int i = 0; i < N_TRIG; i++){
    time[i] = new float[SAMPLE_SIZE];
    tMin = trigOffset_[i]*1e9;
    tMax = (horizScale_[i]*(SAMPLE_SIZE-1)-trigOffset_[i])*1e9;
    time[i] = PulseTools::LinSpace(tMin,tMax,1./SAMPLE_SIZE);
  }
  
  return time;
}

void DataReaderH5::PrintInfo(){
  cout << "x: " << x_[0] << ", y: " << y_[0] << endl;
  cout << "rate: " << rate_[0] << endl;
}
