#include "DataReaderH5.h"

using namespace std;
using namespace H5;

DataReaderH5::DataReaderH5(){
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
  xDs_->read(x_,PredType::NATIVE_INT);
  yDs_->read(y_,PredType::NATIVE_INT);
  rateDs_->read(rate_,PredType::NATIVE_FLOAT);
  samplesDs_->read(samples_,PredType::NATIVE_INT);
  horizOffsetDs_->read(horizOffset_,PredType::NATIVE_FLOAT);
  horizScaleDs_->read(horizScale_,PredType::NATIVE_FLOAT);
  trigOffsetDs_->read(trigOffset_,PredType::NATIVE_FLOAT);
  trigTimeDs_->read(trigTime_,PredType::NATIVE_FLOAT);
  vertOffsetDs_->read(vertOffset_,PredType::NATIVE_FLOAT);
  vertScaleDs_->read(vertScale_,PredType::NATIVE_FLOAT);

}

void DataReaderH5::PrintInfo(){
  cout << "x: " << x_[0] << ", y: " << y_[0] << endl;
  cout << "rate: " << rate_[0] << endl;
}
