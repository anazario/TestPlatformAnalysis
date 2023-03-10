#include "BSpline.h"

using namespace arma;

BSpline::BSpline(const std::vector<double> &knot_vector, const int order, const double precision):
  knot_vector_(knot_vector),
  n_knots_(int(knot_vector.size())),
  precision_(precision),
  order_(order),
  ls_error_(-999.){

  if(order_ > n_knots_-2 || order_ < 0)
    throw std::invalid_argument("Polynomial order has to be a positive integer no larger than the number of knots minus 2!");
  
}

BSpline::~BSpline() = default;

//public
double BSpline::GetError() const{
  return ls_error_;
}

std::vector<double> BSpline::SplineLS(const std::vector<double> &samples){

  int sample_size = int(samples.size());

  dmat collocation_matrix = CollocationMatrix(sample_size);
  dvec fit_coeff = CoeffLeastSquare(samples, collocation_matrix);

  dvec spline;
  SplineFunction(spline, samples, fit_coeff, collocation_matrix);

  return {spline.begin(), spline.end()};
}

//private
int BSpline::FindInterval(const double value) const{

  if(value < knot_vector_[0] || value > knot_vector_[n_knots_-1])
    throw std::invalid_argument("Value must be within the domain of the knot vector!");

  int interval = -1;
  
  for(int i = 0; i < n_knots_-1; i++)
    if (value >= knot_vector_[i] && value < knot_vector_[i+1])
      interval = i;

  return interval;
}

dvec BSpline::BasisValuesAt(const double t) const{

  auto omega = [](double x, double ti, double tik){
    if(ti == tik)
      return 0.;
    else
      return (x-ti)/(tik-ti);
  };

  dmat basis = zeros<dmat>(n_knots_-1, n_knots_-1); 
  basis(0, FindInterval(t)) = 1.;

  for(uword k = 1; k < order_+1; k++){
    drowvec B = basis.row(k-1);
    for(uword i = 0; i < n_knots_-k-1; i++)
      basis(k,i) = omega(t, knot_vector_[i], knot_vector_[i+k])*B(i) + (1. - omega(t, knot_vector_[i+1], knot_vector_[i+k+1]))*B(i+1);
  }

  dvec basis_slice = basis(order_, span(0, n_knots_-order_-2)).t();

  return basis_slice;
}

dmat BSpline::CollocationMatrix(const int size) const{

  vec domain = linspace(knot_vector_[0], knot_vector_[n_knots_-1]-precision_, size);

  dmat collocation_matrix = zeros<dmat>(n_knots_-order_-1, size);
  
  for(uword i = 0; i < domain.size(); i++)
    collocation_matrix.col(i) = BasisValuesAt(domain(i));

  return collocation_matrix;
}

dvec BSpline::CoeffLeastSquare(const dvec &samples, const dmat &collocation_matrix) {
  
  //return (inv_sympd(collocation_matrix*collocation_matrix.t())*collocation_matrix)*samples;
  return solve(collocation_matrix*collocation_matrix.t(), collocation_matrix)*samples;
}

dvec& BSpline::SplineFunction(dvec &spline, const dvec &samples, const dvec &coeff, const dmat &collocation_matrix){

  int sample_size = int(samples.size());
  spline = zeros<dvec>(sample_size);

  for(int i = 0; i < sample_size; i++)
    spline(i) = dot(coeff, collocation_matrix.col(i));

  dvec diff = spline-conv_to<dvec>::from(samples);
  ls_error_ = dot(diff, diff);
  
  return spline;
}


void PlotSplineFit(const std::vector<double> &waveform,
		   const std::vector<double> &spline,
		   const std::vector<double> &knot_vector,
		   const std::string &name){

  std::vector<int> colors = {kBlack, kOrange+1};

  // Create TGraph objects for each data set
  const int data_size = int(waveform.size());
  std::vector<double> x_data = LinSpaceVec(knot_vector.front(), knot_vector.back(), data_size);
  std::vector<double> knot_values(knot_vector.size(), 0.);
  ShannonNyquistInterp(x_data, waveform, knot_vector, knot_values);

  std::vector<TGraph*> graphs(3, nullptr);
  graphs[0] = new TGraph(data_size, &x_data[0], &waveform[0]); 
  graphs[1] = new TGraph(data_size, &x_data[0], &spline[0]);
  auto g_knots = new TGraph(int(knot_vector.size()), &knot_vector[0], &knot_values[0]);
  
  // Create a TMultiGraph object to hold the TGraphs
  auto multi_graph = new TMultiGraph();
  for(int i = 0; i < graphs.size()-1; i++){
    graphs[i]->SetLineWidth(3);
    graphs[i]->SetLineColor(short(colors[i]));
    graphs[i]->SetMarkerSize(0);
    multi_graph->Add(graphs[i]);
  }

  g_knots->SetMarkerStyle(20);
  g_knots->SetMarkerColor(kBlack);
  g_knots->SetLineWidth(0);
  //multi_graph->Add(graphs[2]);

  // Create a canvas to draw the TMultiGraph
  auto canvas = new TCanvas("canvas", "canvas");

  // create the legend
  auto legend = new TLegend(0.6, 0.6, 0.8, 0.8);
  legend->AddEntry(graphs[0], "Interpolated Data", "l");
  legend->AddEntry(graphs[1], "Spline", "l");
  legend->AddEntry(g_knots, "Knots", "p");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  
  // Draw the TMultiGraph with proper axis labels and title
  multi_graph->SetTitle("");
  multi_graph->GetXaxis()->SetTitle("Domain");
  multi_graph->GetXaxis()->CenterTitle();
  multi_graph->GetYaxis()->SetTitle("Amplitude");
  multi_graph->GetYaxis()->CenterTitle();
  multi_graph->Draw("AL");
  g_knots->Draw("Psame");

  legend->Draw("same");
  canvas->Update();
  // Save the canvas to a PDF file with the given name
  std::string filename = name + ".pdf";
  canvas->SetGrid();
  canvas->SaveAs(filename.c_str());
  
  // Clean up the TGraphs, TMultiGraph, and TCanvas
  for (TGraph* graph : graphs) {
    delete graph;
  }
  delete g_knots;
  delete multi_graph;
  delete canvas;
}
