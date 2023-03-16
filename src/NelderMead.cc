#include "NelderMead.h"

NelderMead::NelderMead(int n_dim, std::vector<std::vector<double>> init_simplex, std::vector<int> fixed_idx,
                       double reflection_coefficient, double expansion_coefficient, double contraction_coefficient,
                       double shrinkage_coefficient, double tol)
        : n_dim_(n_dim),
          order_(n_dim + 1),
          simplex_(n_dim + 1, std::vector<double>(n_dim)),
          f_values_(n_dim + 1),
          init_points_(std::move(init_simplex)),
          fixed_idx_(std::move(fixed_idx)),
          reflection_coefficient_(reflection_coefficient),
          expansion_coefficient_(expansion_coefficient),
          contraction_coefficient_(contraction_coefficient),
          shrinkage_coefficient_(shrinkage_coefficient),
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
                simplex_[order_[n_dim_]] = expanded_point;
                f_values_[order_[n_dim_]] = f_expanded;
            }
            else {
                simplex_[order_[n_dim_]] = reflected_point;
                f_values_[order_[n_dim_]] = f_reflected;
            }
        }
        else {
            if (f_reflected < f_values_[order_[n_dim_ - 1]]) {
                simplex_[order_[n_dim_]] = reflected_point;
                f_values_[order_[n_dim_]] = f_reflected;
            }
            else {
                std::vector<double> contracted_point = contract_point();
                double f_contracted = (*objective_function_)(contracted_point);
                if (f_contracted < f_values_[order_[n_dim_]]) {
                    simplex_[order_[n_dim_]] = contracted_point;
                    f_values_[order_[n_dim_]] = f_contracted;
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

    return simplex_[order_[0]];
}

void NelderMead::initialize_simplex() {
    for (int i = 0; i < n_dim_ + 1; i++) {
        simplex_[i] = init_points_[i];
        if (simplex_[i].size() != n_dim_) {
            throw std::invalid_argument("Invalid initialization points");
        }
    }
}

void NelderMead::evaluate(){
    for (int i = 0; i < n_dim_ + 1; i++){
        f_values_[i] = (*objective_function_)(simplex_[i]);
    }
}

void NelderMead::order_points() {
    // Create a vector of indices and sort them by the corresponding function values
    std::vector<int> indices(simplex_.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { return f_values_[i] < f_values_[j]; });

    // Reorder simplex and f_values_ according to the sorted indices
    std::vector<std::vector<double>> sorted_simplex(simplex_.size(), std::vector<double>(n_dim_));
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
    std::vector<double> centroid(n_dim_, 0.0);
    for (int i = 0; i < simplex_.size() - 1; ++i) {
        for (int j = 0; j < n_dim_; ++j)
            centroid[j] += simplex_[i][j];
    }
    for (int j = 0; j < n_dim_; ++j)
        centroid[j] /= static_cast<double>(n_dim_);

    return centroid;
}

std::vector<double> NelderMead::reflect_point(){
    // Compute the reflection of the worst point with respect to the centroid
    std::vector<double> centroid = calculate_centroid();
    std::vector<double> reflected_point(n_dim_);

    for (int j = 0; j < n_dim_; ++j)
        reflected_point[j] = centroid[j] + reflection_coefficient_ * (centroid[j] - simplex_.back()[j]);

    sort_non_decreasing(reflected_point);

    return reflected_point;
}

std::vector<double> NelderMead::expand_point(){
    // Compute the expansion of the reflected point
    std::vector<double> centroid = calculate_centroid();
    std::vector<double> expanded_point(n_dim_);
    for (int j = 0; j < n_dim_; ++j)
        expanded_point[j] = centroid[j] + expansion_coefficient_ * (reflect_point()[j] - centroid[j]);

    return expanded_point;
}

std::vector<double> NelderMead::contract_point() {
    // Compute the contraction of the worst point with respect to the centroid
    std::vector<double> centroid = calculate_centroid();
    std::vector<double> contracted_point(n_dim_);
    for (int j = 0; j < n_dim_; ++j)
        contracted_point[j] = centroid[j] + contraction_coefficient_ * (simplex_.back()[j] - centroid[j]);

    return contracted_point;
}

void NelderMead::shrink_simplex() {
    for (int i = 1; i < n_dim_ + 1; ++i) {
        for (int j = 0; j < n_dim_; ++j)
            simplex_[i][j] = simplex_[0][j] + shrinkage_coefficient_ * (simplex_[i][j] - simplex_[0][j]);
    }
    evaluate();
}

void NelderMead::update_simplex() {
    for (int i = 0; i < n_dim_ + 1; ++i)
        copy_values(fixed_idx_, init_points_[0], simplex_[i]);
    order_points();
}

bool NelderMead::check_convergence() {
    double max_diff = 0.0;
    for (int i = 1; i < n_dim_ + 1; i++) {
        double diff = std::abs(f_values_[order_[i]] - f_values_[order_[0]]);
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    return (max_diff < tol_);
}
