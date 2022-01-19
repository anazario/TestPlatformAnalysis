#include "PulseTools.h"

using namespace std;

PulseTools::PulseTools(){};

PulseTools::~PulseTools(){};

/***************************************************/
//Name: FindMinAbsolute
//Input: 1.Array of floats (float*) and 2.Total size 
//of the array. 
//Output: index of the array element with the lowest
//float value.
/***************************************************/
int PulseTools::FindMinAbsolute(float* sample, int size){
  float xmin = 999;
  
  int index_min = 0;
  
  for (int i = 0 ; i < size; i++){
    if (sample[i]<xmin){// && sample[i+1] < 0.5*sample[i]){
      xmin = sample[i]; //minimum
      index_min = i; //index number of minimum
    }
  }
  //cout << "min value: " << xmin << " at index : " << index_min << endl;
  return index_min;
}

/***************************************************/
//Name: FindMaxAbsolute                                                                                                                    
//Input: 1.Array of floats (float*) and 2.Total size                                                                                       
//of the array.                                                                                                                            
//Output: index of the array element with the highest                                                                                       
//float value.                                                                                                                             
/***************************************************/
int PulseTools::FindMaxAbsolute(float* sample, int size){
  float xmax = sample[20];
  
  int index_max = 0;
  
  for (int i = 20 ; i < size-20; i++){
    if (sample[i]>xmax){
      xmax = sample[i]; //maximum
      index_max = i; //index number of maximum
    }
  }
  return index_max;
}

/**************************************************************/
//Name: LinSpace
//Input: 1.Step size (float), 2.Minimum x-range value (float) 
//and 3.Maximum x-range value (float).
//Output: Array of evenly spaced elements between the values of 
//xmin and xmax (args 2. and 3.) with a separation of step_size 
//(arg 1.).
/**************************************************************/
float* PulseTools::LinSpace(float xmin, float xmax, float step_size){
  
  const int nElements = (xmax - xmin)/step_size + 1;
  float* time_arr = new float[nElements];
  
  for(int i = 0; i < nElements; i++){
    time_arr[i] = round((xmin+float(i)*step_size)*100.)/100.;
  }
  return time_arr;
}

vector<float> PulseTools::LinSpaceVec(float xmin, float xmax, int size){

  const float stepSize = (xmax - xmin)/size;

  vector <float> time_vec;
  
  for(int i = 0; i < size; i++){
    time_vec.push_back(round((xmin+float(i)*stepSize)*100.)/100.);
  }
  return time_vec;
}

float PulseTools::Integral(float* sample, float step_size){
  
  int total = int(2./step_size)+1;
  float sum = 0;
  
  for(int i = 0; i < total; i++){
    cout << sample[i] << endl;
    if(sample[i] > 0){
      sum += sample[i]*step_size;
      //cout << sum << endl;
    }
  }
  return sum;
}

float PulseTools::Integral(vector<float> sample, float step_size){

  int total = sample.size();
  float sum = 0;

  for(int i = 0; i < total; i++){
    if(sample[i] > 0)
      sum += sample[i]*step_size;
  }
  return sum;
}

/***************************************************************************/
//Name: MakeHistArr
//Input: 1.Total size of vector to be produced (int), 2.Minimum x-range value 
//of 1D histograms added to vector, 3.Maximum x-range of 1D histograms and 
//4.Total amount of bins in each histogram
//Output: Vector with a total of 'size' (arg 1.) elements. All TH1F objects 
//in the vector have the same range in x, between xmin and xmax (arg 2. and 
//arg 3), and the same amount of bins, given by nbins (arg 4.).
/***************************************************************************/
vector<TH1F> PulseTools::MakeHistArr(int size, float xmin, float xmax, int nbins){
  TString name;
  
  vector<TH1F> pulse_hist;
  TH1F temp_hist;
  
  for(int i = 0; i < size; i++){
    name.Form("pulse_hist%i",i);
    temp_hist = TH1F(name,name,nbins,xmin,xmax);
    pulse_hist.push_back(temp_hist);
  }
  return pulse_hist;
}

/************************************************/
//Name: GetMeanArr
//Input: 1.Vector of TH1F objects
//Output: Array of the corresponding mean value of 
//each of the TH1F objects in the input vector.
/************************************************/
float* PulseTools::GetMeanArr(vector<TH1F> pulse_hist){
  
  int arr_size = pulse_hist.size();
  float* mean_arr = new float[arr_size];
  
  for (int i = 0; i < arr_size; i++){
    mean_arr[i] = float(pulse_hist[i].GetMean());
    //cout << mean_arr[i] << endl;
  }
  return mean_arr;
}

/**********************************************************************************/
//Name: InterpolateFunc
//Input: 1.Pulse sample to be interpolated (float*), 2.Total elements in the sample 
//(int) and 3.Point in time to be interpolated from sample (float). 
//Output: Amplitude value (float) of the interpolated sample at time t (arg 3.).
//Uses the Nyquist-Shannon sampling theorem.
/**********************************************************************************/
float PulseTools::InterpolateFunc(float* sample, int sample_size, float time){
  
  float frequency = 10.5e9;//10./(2e-9);
  float omega_c = 2.*M_PI*frequency;
  float sampling_time = 1./frequency;
  
  time *= (1e-9);
  
  float sum = 0;
  
  for(int i = 0; i < sample_size; i++){
    sum += sample[i]*sin(omega_c*(time - i*sampling_time))/(omega_c*(time - i*sampling_time));
  }
  return (isnan(sum))?0:sum;
}

float PulseTools::GetInterpolatedPoint(vector<float> pulse, float time, float frequency){

  int size = pulse.size();
  double bandLimit = frequency/2;//half the sampling frequency (GHz)

  double sum = 0.;
  for(int i = 0; i < size; i++){
    double argument = M_PI*(2*bandLimit*time - i);
    double sinc = pulse[i]*sin(argument)/argument;
    if(isnan(sinc)){
      sum += 0;
      cout << i	<< endl;
      cout << time << endl;
      cout << argument << endl;
    }
    else{
      sum += sinc;
    }
  }

  return (isnan(sum))?0:sum;
}

/********************************************************************************/
//Name: CalcInterval
//Input: 1.Filled histogram (TH1F), 2.Confidence interval percentage (values 
//between 0 - 1.0) around mean, 3.Empty variable for mean value of histogram 
//(float&), 4.Empty variable for errors to the left side of the mean (float&) and 
//5.Empty variable for errors on the right side of the mean (float&).
//Output: The input mean, low and up (args 2,3 and 4) are filled and passed by
//reference.
/********************************************************************************/
void PulseTools::CalcInterval(TH1F hist, float CI, float& mean, float& low, float& high){
  if(CI < 0. || CI >= 1.)
    return;
  
  hist.Scale(1./hist.Integral());
  int Nbin = hist.GetNbinsX();
  
  mean = hist.GetMean();
  int imean = hist.GetXaxis()->FindBin(mean);
  
  float prob = hist.GetBinContent(imean);
  
  float probLeft, probRight;
  int binIdxLeft  = imean-1;
  int binIdxRight = imean+1;
  int lastBinIdxLeft  = binIdxLeft;
  int lastBinIdxRight = binIdxRight;

  while(prob < CI){
    
    if(binIdxLeft >= 1)
      probLeft = hist.GetBinContent(binIdxLeft);
    else
      probLeft = 0;
    
    if(binIdxRight <= Nbin)
      probRight = hist.GetBinContent(binIdxRight);
    else
      probRight = 0;
    
    if(binIdxLeft <= 1 || binIdxRight >= Nbin){
      cout << "Breaking with total probability: " << prob << " (CI = " << CI << ")" << endl;
      break;
    }
   
    if(probLeft == 0. && probRight == 0.){
      binIdxLeft--;
      binIdxRight++;
      continue;
    }  
    
    if(probLeft > probRight){
      prob += probLeft;
      lastBinIdxLeft = binIdxLeft;
      binIdxLeft--;
    }
    else {
      prob += probRight;
      lastBinIdxRight = binIdxRight;
      binIdxRight++;
    }
  }
  
  low = mean - hist.GetXaxis()->GetBinCenter(lastBinIdxLeft);
  high  = hist.GetXaxis()->GetBinCenter(lastBinIdxRight) - mean;
}



