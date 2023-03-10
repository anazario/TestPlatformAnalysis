#ifndef NELDERMEAD_H
#define NELDERMEAD_H

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "PulseTools.h"

class NelderMead {
public:
  NelderMead(int ndim, std::vector<std::vector<double>> init_points, std::vector<int> fixed_idx={}, double reflection_coef=1.0,
	     double expansion_coef=2.0, double contraction_coef=0.5, double shrinkage_coef=0.5, double tol=1e-8);

  virtual ~NelderMead();
  
  std::vector<double> optimize(std::function<double(std::vector<double>)> &func, int max_it = 500);

private:
  int ndim_;
  int iter_;
  std::vector<int> order_;
  std::vector<int> fixed_idx_;
  std::vector<double> f_values_;
  std::vector<std::vector<double>> simplex_;
  std::vector<std::vector<double>> init_points_;

  std::function<double(std::vector<double>)> *objective_function_ = nullptr;
  
  double reflection_coef_;
  double expansion_coef_;
  double contraction_coef_;
  double shrinkage_coef_;
  double tol_;

  void initialize_simplex();
  void evaluate();
  void order_points();
  void shrink_simplex();
  //void enforce_unique(std::vector<double>& row);
  void update_simplex();
  bool check_convergence();
  std::vector<double> calculate_centroid();
  std::vector<double> reflect_point();
  std::vector<double> expand_point();
  std::vector<double> contract_point();
};

#endif
