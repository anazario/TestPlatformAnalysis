#ifndef BSpline_h
#define BSpline_h

#include <vector>
#include <armadillo>
#include <TGraph.h>

#include "PulseTools.h"

class BSpline{

 public:

  explicit BSpline(const std::vector<double> &knot_vector, int order=3, double precision=1e-4);

  virtual ~BSpline();

  double GetError() const;
  std::vector<double> SplineLS(const std::vector<double> &samples);
  std::vector<double> GetFitCoefficients();
  
private:

  int n_knots_;
  int order_;
  double precision_;
  double ls_error_;

  std::vector<double> knot_vector_;

  arma::dvec fit_coeff_;

  int FindInterval(double value) const;
  arma::dvec BasisValuesAt(double t) const;
  arma::dmat CollocationMatrix(int size=100) const;
  static arma::dvec CoeffLeastSquare(const arma::dvec &samples, const arma::dmat &collocation_matrix) ;
  arma::dvec& SplineFunction(arma::dvec &spline, const arma::dvec &samples, const arma::dvec &coeff, const arma::dmat &collocation_matrix);
  
};

void PlotSplineFit(const std::vector<double> &waveform,
		   const std::vector<double> &spline,
		   const std::vector<double> &knot_vector,
		   const std::string &name);

#endif
