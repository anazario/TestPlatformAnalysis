#include "PulseTools.h"

using namespace std;

int FindMaxIndex(const std::vector<double>& sample){

  int maxIndex = -999;
  double maxSample = -999;
  
  for(int s = 0; s < sample.size(); s++){
    if(sample[s] > maxSample){
      maxSample = sample[s];
      maxIndex = s;
    } 
  }
  return maxIndex;
}

double FindMax(const std::vector<double>& sample){
  int maxIndex = FindMaxIndex(sample);
  return sample[maxIndex];
}

double GetSampleRMS(const std::vector<double>& sample){

  TH1D sampleHist("sample","sample",40,-0.1,0.1);

  for(int i = 0; i < sample.size(); i++)
    sampleHist.Fill(sample[i]);

  return sampleHist.GetRMS();
}

double Integral(const vector<double>& sample, const double step_size){

  int total = sample.size();
  double sum = 0;

  for(int i = 0; i < total; i++){
    if(sample[i] > 0)
      sum += sample[i]*step_size;
  }
  return sum;
}

double GetInterpolatedPoint(const vector<double>& pulse, const double time, const double frequency){

  int size = pulse.size();

  double sum = 0.;
  for(int i = 0; i < size; i++){
    double argument = M_PI*(frequency*time - double(i));
    double sinc = -999.;

    if(fabs(argument) < 1e-8){
      sinc = 1.;
    }
    else
      sinc = sin(argument)/argument;

    sum += pulse[i]*sinc;
  }
  
  return (isnan(sum))?0:sum;
}

vector<double> LinSpaceVec(const double xmin, const double xmax, const int size){

  const double stepSize = (xmax - xmin)/size;

  vector<double> time_vec;
  
  for(int i = 0; i < size; i++){
    time_vec.push_back(round((xmin+double(i)*stepSize)*100.)/100.);
  }
  return time_vec;
}

vector<double> GetInterpolatedPulse(const int interpolationSize, const double pulseStart, const double pulseEnd,
				    const Pulse& originalPulse){

  double sampleRate = originalPulse.GetSampleRate();

  vector<double> interpolatedPulse;
  vector<double> interpolationTime = LinSpaceVec(pulseStart,pulseEnd,interpolationSize);

  double timeDiff = originalPulse.GetStartTime();

  for(int i = 0; i < interpolationSize; i++){    
    double amplitude = -GetInterpolatedPoint(originalPulse.GetPulse(), interpolationTime[i]-timeDiff, sampleRate);
    interpolatedPulse.push_back(amplitude);
  }

  return interpolatedPulse;
}

TH1D* GetHistFromVec(const vector<double>& inputVector, const TString name, const int bins, const double xInitial, const double xFinal){

  TH1D* distribution = new TH1D(name,name,bins,xInitial,xFinal);

  for(int i = 0; i < inputVector.size(); i++)
    distribution->Fill(inputVector[i]);

  return distribution;
  
}

vector<TH1D> MakeHistArr(const int size, const double xmin, const double xmax, const int nbins){
  TString name;
  
  vector<TH1D> pulse_hist;
  TH1D temp_hist;
  
  for(int i = 0; i < size; i++){
    name.Form("pulse_hist%i",i);
    temp_hist = TH1D(name,name,nbins,xmin,xmax);
    pulse_hist.push_back(temp_hist);
  }
  return pulse_hist;
}

void CalcInterval(TH1D hist, const double CI, double& mean, double& low, double& high){
  if(CI < 0. || CI >= 1.)
    return;
  
  hist.Scale(1./hist.Integral());
  int Nbin = hist.GetNbinsX();
  
  mean = hist.GetMean();
  int imean = hist.GetXaxis()->FindBin(mean);
  
  double prob = hist.GetBinContent(imean);
  
  double probLeft, probRight;
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

void NormalizeVec(std::vector<double>& sample){
  double max = FindMax(sample);

  for(int s = 0; s < sample.size(); s++)
    sample[s] = sample[s]/max;
}

void SplitFolder(std::string& fullPath, std::string& innerMostName){

  int length = fullPath.length();
  int slashIdx = fullPath.find_last_of('/');

  std::string path;

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

bool Replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}
