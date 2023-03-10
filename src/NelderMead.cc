#include "NelderMead.h"

#include <utility>

NelderMead::NelderMead(int ndim, std::vector<std::vector<double>> init_points, std::vector<int> fixed_idx, double reflection_coef,
		       double expansion_coef, double contraction_coef, double shrinkage_coef, double tol)
  : ndim_(ndim),
    order_(ndim + 1),
    simplex_(ndim + 1, std::vector<double>(ndim)),
    f_values_(ndim + 1),
    init_points_(std::move(init_points)),
    fixed_idx_(std::move(fixed_idx)),
    reflection_coef_(reflection_coef),
    expansion_coef_(expansion_coef),
    contraction_coef_(contraction_coef),
    shrinkage_coef_(shrinkage_coef),
    tol_(tol),
    iter_(0.){

  initialize_simplex(); 
}

NelderMead::~NelderMead() = default;

std::vector<double> NelderMead::optimize(std::function<double(std::vector<double>)> &func, const int max_iter) {
    objective_function_ = &func;
    evaluate();
    order_points();

    if(iter_ != 0)
        iter_ = 0;

    while (!check_convergence() && iter_ < max_iter) {
        std::vector<double> centroid = calculate_centroid();
        std::vector<double> reflected_point = reflect_point();
        double f_reflected = (*objective_function_)(reflected_point);
    
    if (f_reflected < f_values_[order_[0]]) {
        std::vector<double> expanded_point = expand_point();
        double f_expanded = (*objective_function_)(expanded_point);

        if (f_expanded < f_reflected) {
            simplex_[order_[ndim_]] = expanded_point;
            f_values_[order_[ndim_]] = f_expanded;
        }
        else {
            simplex_[order_[ndim_]] = reflected_point;
            f_values_[order_[ndim_]] = f_reflected;
        }
    }
    else {
        if (f_reflected < f_values_[order_[ndim_-1]]) {
            simplex_[order_[ndim_]] = reflected_point;
            f_values_[order_[ndim_]] = f_reflected;
        }
        else {
            std::vector<double> contracted_point = contract_point();
            double f_contracted = (*objective_function_)(contracted_point);
            if (f_contracted < f_values_[order_[ndim_]]) {
                simplex_[order_[ndim_]] = contracted_point;
                f_values_[order_[ndim_]] = f_contracted;
            }
            else {
                shrink_simplex();
            }
        }
    }
    update_simplex();
    iter_++;
  }
  std::cout << "Finished in " << iter_ << " iterations." << std::endl;

  /*
  for(auto s : samples)
    std::cout << s << " ";
  std::cout << std::endl;
  */
  return simplex_[order_[0]];
}

void NelderMead::initialize_simplex() {
  //simplex_[0] = init_points_[0];
  for (int i = 0; i < ndim_ + 1; i++) {
    simplex_[i] = init_points_[i];
    if (simplex_[i].size() != ndim_) {
      throw std::invalid_argument("Invalid initialization points");
    }
  }
}

void NelderMead::evaluate(){
  for (int i = 0; i < ndim_ + 1; i++){
    f_values_[i] = (*objective_function_)(simplex_[i]);
  }
}

void NelderMead::order_points() {
  // Create a vector of indices and sort them by the corresponding function values
  std::vector<int> indices(simplex_.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](int i, int j) { return f_values_[i] < f_values_[j]; });
  
  // Reorder simplex and f_values_ according to the sorted indices
  std::vector<std::vector<double>> sorted_simplex(simplex_.size(), std::vector<double>(ndim_));
  std::vector<double> sorted_f_values(f_values_.size());
  for (int i = 0; i < indices.size(); ++i) {
    int index = indices[i];
    sorted_simplex[i] = simplex_[index];
    sort_non_decreasing(sorted_simplex[i]);
    sorted_f_values[i] = f_values_[index];
  }
  simplex_ = sorted_simplex;
  f_values_ = sorted_f_values;
  order_ = indices;
}

std::vector<double> NelderMead::calculate_centroid() {
  // Compute the centroid of all the vertices except the worst one
  std::vector<double> centroid(ndim_, 0.0);
  for (int i = 0; i < simplex_.size() - 1; ++i) {
    for (int j = 0; j < ndim_; ++j) 
      centroid[j] += simplex_[i][j];
  } 
  for (int j = 0; j < ndim_; ++j) 
    centroid[j] /= static_cast<double>(ndim_);
  
  return centroid;
}

std::vector<double> NelderMead::reflect_point(){//const std::vector<double> &centroid) {
  // Compute the reflection of the worst point with respect to the centroid
  std::vector<double> centroid = calculate_centroid();
  std::vector<double> reflected_point(ndim_);
  
  for (int j = 0; j < ndim_; ++j) 
        reflected_point[j] = centroid[j] + reflection_coef_ * (centroid[j] - simplex_.back()[j]);

  sort_non_decreasing(reflected_point);
  
  return reflected_point;
}

std::vector<double> NelderMead::expand_point(){//const std::vector<double> &centroid, const std::vector<double> &reflected_point) {
  // Compute the expansion of the reflected point
  std::vector<double> centroid = calculate_centroid();
  std::vector<double> expanded_point(ndim_);
  for (int j = 0; j < ndim_; ++j) 
    expanded_point[j] = centroid[j] + expansion_coef_ * (reflect_point()[j] - centroid[j]);

  return expanded_point;
}

std::vector<double> NelderMead::contract_point() {
  // Compute the contraction of the worst point with respect to the centroid
  std::vector<double> centroid = calculate_centroid();
  std::vector<double> contracted_point(ndim_);
  for (int j = 0; j < ndim_; ++j) 
    contracted_point[j] = centroid[j] + contraction_coef_ * (simplex_.back()[j] - centroid[j]);

  return contracted_point;
}

void NelderMead::shrink_simplex() {
  for (int i = 1; i < ndim_ + 1; ++i) {
    for (int j = 0; j < ndim_; ++j) 
      simplex_[i][j] = simplex_[0][j] + shrinkage_coef_ * (simplex_[i][j] - simplex_[0][j]);    
  }
  evaluate();
}

/*void NelderMead::enforce_unique(std::vector<double>& row) {
  int n = row.size();
  for (int i = 0; i < n-1; i++) {
    if (std::find(fixed_idx_.begin(), fixed_idx_.end(), i) != fixed_idx_.end()) {
      continue; // skip if this index is fixed
    }
    if (std::abs(row[i] - row[i-1]) < 1e-3) {
      if (i == n-2 || std::find(fixed_idx_.begin(), fixed_idx_.end(), i+2) != fixed_idx_.end()) {
	row[i+1] += (row[i] - row[i-1])/2.0; // average with previous element
      } else {
	row[i] += (row[i+1] - row[i-1])*3/2.0; // average with next element
      }
    }
  }
}*/

void NelderMead::update_simplex() {
  for (int i = 0; i < ndim_ + 1; ++i) {
    copy_values(fixed_idx_, init_points_[0], simplex_[i]);
    //enforce_unique(simplex_[i]);
  }
  order_points();
  /*
  std::cout << "Simplex " << iter_ << ": " << std::endl;
  for(int i = 0; i < f_values_.size(); i++){
    for(int j = 0; j < simplex_[i].size(); j++)
      std::cout << simplex_[i][j] << " ";
    std::cout << f_values_[i] << std::endl;
  }
  std::cout << std::endl;
  */
}

bool NelderMead::check_convergence() {
  double max_diff = 0.0;
  for (int i = 1; i < ndim_ + 1; i++) {
    double diff = std::abs(f_values_[order_[i]] - f_values_[order_[0]]);
    if (diff > max_diff) {
      max_diff = diff;
    }
  }
  //std::cout << "Convergence value: " << max_diff << std::endl;
  return (max_diff < tol_);
}
