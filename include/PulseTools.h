#ifndef PULSETOOLS_H
#define PULSETOOLS_H

#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom.h"
#include "TEventList.h"
#include "TColor.h"

using namespace std;

class PulseTools{

 public:
  
  PulseTools() {};

  virtual ~PulseTools(){ };

  /***************************************************/
  //Name: FindMinAbsolute
  //Input: 1.Array of floats (float*) and 2.Total size 
  //of the array. 
  //Output: index of the array element with the lowest
  //float value.
  /***************************************************/
  static int FindMinAbsolute(float* sample, int size){
    float xmin = 999;

    int index_min = 0;

    for (int i = 100 ; i < size-2; i++){
      if (sample[i]<xmin && sample[i+1] < 0.5*sample[i]){
        xmin = sample[i]; //minimum
        index_min = i; //index number of minimum
      }
    }
    return index_min;
  }

  /***************************************************/
  //Name: FindMaxAbsolute                                                                                                                    
  //Input: 1.Array of floats (float*) and 2.Total size                                                                                       
  //of the array.                                                                                                                            
  //Output: index of the array element with the highest                                                                                       
  //float value.                                                                                                                             
  /***************************************************/
  static int FindMaxAbsolute(float* sample, int size){
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
  static float* LinSpace(float xmin, float xmax, float step_size){

    const int nElements = (xmax - xmin)/step_size + 1;
    float* time_arr = new float[nElements];

    for(int i = 0; i < nElements; i++){
      time_arr[i] = round((xmin+float(i)*step_size)*100.)/100.;
    }
    return time_arr;
  }

  static float Integral(float* sample, float step_size){

    int total = int(2./step_size)+1;
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
  static vector<TH1F> MakeHistArr(int size, float xmin, float xmax, int nbins){
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
  static float* GetMeanArr(vector<TH1F> pulse_hist){

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
  //Uses the Shannon-Nyquist sampling theorem.
  /**********************************************************************************/
  static float InterpolateFunc(float* sample, int sample_size, float time){

    float frequency = 5e9;
    float omega_c = 2*M_PI*frequency;
    float sampling_time = 1/(2*frequency);

    time *= (1e-9);

    float sum = 0;

    for(int i = 0; i < sample_size; i++){
      sum += sample[i]*sin(omega_c*(time - i*sampling_time))/(omega_c*(time - i*sampling_time));
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
  static void CalcInterval(TH1F hist, float CI, float& mean, float& low, float& high){
    if(CI < 0. || CI >= 1.)
      return;

    hist.Scale(1./hist.Integral());
    int Nbin = hist.GetNbinsX();

    mean = hist.GetMean();
    int imean = hist.GetXaxis()->FindBin(mean);

    float prob = hist.GetBinContent(imean);
    
    float pl, pr;
    int cl = imean-1;
    int cr = imean+1;
    int cl_non0 = cl;
    int cr_non0 = cr;
    while(prob < CI){
      if(cl >= 1)
        pl = hist.GetBinContent(cl);
      else
        pl = 0;

      if(cr <= Nbin)
        pr = hist.GetBinContent(cr);
      else
        pr = 0;

      if(cl <= 1 || cr >= Nbin){
        break;
      }

      if(pl == 0. && pr == 0.){
	cl--;
	cr++;
	continue;
      }  

      if(pl > pr){
        prob += pl;
	cl_non0 = cl;
        cl--;
      } else {
        prob += pr;
	cr_non0 = cr;
        cr++;
      }
    }
    low = mean - hist.GetXaxis()->GetBinCenter(cl_non0);
    high  = hist.GetXaxis()->GetBinCenter(cr_non0) - mean;
  }
};

#endif
