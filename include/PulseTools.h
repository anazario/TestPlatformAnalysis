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

#include "Pulse.h"
#include "Plot.h"

int FindMaxIndex(const std::vector<double>& sample);

double FindMax(const std::vector<double>& sample);
double GetSampleRMS(const std::vector<double>& sample);
double Integral(const std::vector<double>& sample, const double step_size);
double GetInterpolatedPoint(const std::vector<double>& pulse, const double time, double frequency);

std::vector<double> LinSpaceVec(const double xmin, const double xmax, const int size);
std::vector<double> GetInterpolatedPulse(const int interpolationSize, const double pulseStart, const double pulseEnd,
					 const Pulse& originalPulse);

TH1D* GetHistFromVec(const std::vector<double>& inputVector, const TString name, const int bins, const double xInitial, const double xFinal);
std::vector<TH1D> MakeHistArr(const int size, const double xmin, const double xmax, const int nbins);

void NormalizeVec(std::vector<double>& sample);
void CalcInterval(TH1D hist, const double CI, double& mean, double& low, double& high);

void SplitFolder(std::string& fullPath, std::string& innerMostName);
bool Replace(std::string& str, const std::string& from, const std::string& to);

    
