#ifndef PULSETOOLS_H
#define PULSETOOLS_H

#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <complex>
#include <functional>
#include <random>
#include <limits>
#include <armadillo>
#include <fftw3.h>
#include <chrono>

#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom.h"
#include "TEventList.h"
#include "TColor.h"

//#include "Pulse.h"
#include "NelderMead.h"
#include "BSpline.h"
#include "Plot.h"

int FindMaxIndex(const std::vector<double>& sample);

double FindMax(const std::vector<double>& sample);
double GetSampleRMS(const std::vector<double>& sample);
double Integral(const std::vector<double>& sample, double step_size);
double GetInterpolatedPoint(const std::vector<double>& pulse, double time, double frequency);
double Mean(const std::vector<double>& inputVec);
void PedestalSubtraction(std::vector<double>& signal, int subset_size);
double CalculateCFDTime(const std::vector<double>& sample, const std::vector<double>& time, double cfd_fraction);

std::vector<double> Diff(const std::vector<double>& inputVec);
std::vector<double> LinSpaceVec(double xmin, double xmax, int size);
std::vector<double> ArangeVec(double start, double stop, double step);
  
TH1D* GetHistFromVec(const std::vector<double>& inputVector, const TString& name, int bins, double xInitial, double xFinal);
std::vector<TH1D> MakeHistArr(int size, double xmin, double xmax, int nbins);

void NormalizeVec(std::vector<double>& sample);
void CalcInterval(TH1D hist, double CI, double& mean, double& low, double& high);

void SplitFolder(std::string& fullPath, std::string& innerMostName);
bool Replace(std::string& str, const std::string& from, const std::string& to);

void sort_non_decreasing(std::vector<double>& v);
std::vector<int> find_matching_indices(const std::vector<double>& data, const std::vector<double>& values_to_check);
void copy_values(const std::vector<int>& indices, const std::vector<double>& source, std::vector<double>& dest);
bool checkVector(const std::vector<double>& input_vec);
  
std::vector<std::vector<double>> CreateSimplex(const std::vector<double>& input_coords, const std::vector<int>& fixed_idx);
std::vector<std::vector<double>> CreateSimplexRand(const std::vector<double>& input_coords, double std_dev = 0.1);
std::vector<double> MinimizeKnots(const std::vector<double> &knot_vector, const std::vector<double> &samples,
                                  int max_iter = 100, double stddev = 1e-3, bool isRand = false);
double goldenSectionSearch(const std::function<double(double)> &f, double& a, double& b, double tol = 1e-6);
double BrentsMethod(const std::function<double(double)> &func, double a, double b, double tol = 1e-8, int max_iter = 1e3);

void FillDensityMatrix(std::vector<std::vector<double>> &density_matrix, const std::vector<double> &sample);
std::vector<double> GetEigenvectorAtIndex(const std::vector<std::vector<double>>& matrix, int index = 0);

double SumGauss(const std::vector<double>& inputVector, double stddev);
std::function<double(std::vector<double>)> BSplineErr(const std::vector<double>& samples, double stddev);
std::function<double(double)> WindowFitErr(const std::function<std::vector<double>(const std::vector<double>&)> &interp_function,
                                           std::vector<double> knot_vector,
                                           double time_begin, double time_end,
                                           double window_size=25,
                                           double step_interp=1e-3);

double WindowFitErr(double t,
                    std::vector<double> knot_vector,
                    std::vector<double> &samples,
                    double time_begin, double time_end,
                    double window_size=25,
                    double step_interp=1e-3);


std::vector<double> sincInterp(std::vector<double> x, std::vector<double> y, std::vector<double> t);

void ShannonNyquistInterp(const std::vector<double>& xdata, const std::vector<double>& ydata,
                          const std::vector<double>& xnew, std::vector<double>& ynew);
std::function<std::vector<double>(const std::vector<double>&)> ShannonNyquistInterp(
    const std::vector<double>& xdata,
    const std::vector<double>& ydata);

std::function<std::vector<double>(const std::vector<double>&)> ShannonNyquistInterpFFT(
    const std::vector<double>& xdata,
    const std::vector<double>& ydata);

void ShannonNyquistInterpFFT(const std::vector<double>& xdata, const std::vector<double>& ydata,
                             const std::vector<double>& xnew, std::vector<double>& ynew);

std::function<std::vector<double>(const std::vector<double>&)> LinearInterp(
    const std::vector<double>& xdata,
    const std::vector<double>& ydata);

std::function<std::vector<double>(const std::vector<double>&)> CubicInterpolation(
    const std::vector<double>& xdata,
    const std::vector<double>& ydata);

void computeSplineCoefficients(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& a,
                               std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);

std::function<std::vector<double>(const std::vector<double>&)> cubicSplineInterpolate(const std::vector<double>& x, const std::vector<double>& y);

#endif
