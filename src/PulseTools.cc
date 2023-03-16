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

    for(double i : sample)
        sampleHist.Fill(i);

    return sampleHist.GetRMS();
}

double Integral(const vector<double>& sample, const double step_size){

    int total = int(sample.size());
    double sum = 0;

    for(int i = 0; i < total; i++){
        if(sample[i] > 0)
            sum += sample[i]*step_size;
    }
    return sum;
}

double GetInterpolatedPoint(const vector<double>& pulse, const double time, const double frequency){

    int size = int(pulse.size());

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

double Mean(const vector<double>& inputVec){

    double sum = 0;

    for(double i : inputVec)
        sum += i;

    return sum/double(inputVec.size());
}

void PedestalSubtraction(std::vector<double>& signal, const int subset_size) {
    // Calculate the mean of the subset of the signal
    double pedestal_mean = std::accumulate(signal.begin(), signal.begin() + subset_size, 0.0) / subset_size;
    // Subtract the mean from the entire signal
    std::transform(signal.begin(), signal.end(), signal.begin(), [&](double value) { return value - pedestal_mean; });
}

double CalculateCFDTime(const std::vector<double>& sample, const std::vector<double>& time, const double cfd_fraction) {
    // Find the maximum sample value
    double max_sample = *std::max_element(sample.begin(), sample.end());
    double min_sample = *std::min_element(sample.begin(), sample.end());
    // Calculate the CFD threshold value
    double threshold = cfd_fraction * max_sample;

    if(threshold < min_sample)
        return std::numeric_limits<double>::quiet_NaN();

    // Find the time at which the CFD signal first crosses the threshold value
    for (size_t i = 1; i < sample.size(); i++) {
        if (sample[i] > threshold) {
            // Linearly interpolate the crossing time using the two adjacent samples
            double slope = (sample[i] - sample[i-1]) / (time[i] - time[i-1]);
            double crossing_time = time[i-1] + (threshold - sample[i-1]) / slope;
            return crossing_time;
        }
    }
    // If the CFD threshold is not crossed, return NaN
    return std::numeric_limits<double>::quiet_NaN();
}

vector<double> Diff(const vector<double>& inputVec){

    vector<double> diffVec;
    for(int i = 1; i < inputVec.size(); i++)
        diffVec.push_back(inputVec[i]-inputVec[i-1]);

    return diffVec;
}

vector<double> LinSpaceVec(const double xmin, const double xmax, const int vec_size){

    double h = (xmax - xmin) / static_cast<double>(vec_size-1);
    std::vector<double> xs(vec_size);
    //std::vector<double>::iterator x;
    //double val;
    for(int i = 0; i < vec_size; i++){
        xs[i] = xmin + i*h;
        //for (x = xs.begin(), val = xmin; x != xs.end(); ++x, val += h) {
        //*x = val;
    }
    return xs;
}

std::vector<double> ArangeVec(const double start, const double stop, const double step) {

    std::vector<double> vec;
    double n = ceil((stop - start) / step) + 1;
    if (n <= 1) {
        vec.push_back(start);
        return vec;
    }
    vec.reserve(n);
    double val = start;
    vec.push_back(val);
    for (int i = 1; i < n; i++) {
        val += step;
        vec.push_back(val);
        if (i == n - 1 && val != stop) {
            vec.pop_back();
            vec.push_back(stop);
        }
    }
    return vec;
}


TH1D* GetHistFromVec(const vector<double>& inputVector, const TString& name, const int bins, const double xInitial, const double xFinal){

    auto distribution = new TH1D(name,name,bins,xInitial,xFinal);

    for(double i : inputVector)
        distribution->Fill(i);

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

    for(double &s : sample)
        s = s/max;
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

void sort_non_decreasing(std::vector<double>& v) {
    std::sort(v.begin(), v.end(),
              [](double x, double y) {
                  return x < y || (x == y && &x < &y);
              });
}

std::vector<int> FindMatchingIndices(const std::vector<double>& data, const std::vector<double>& values_to_check) {
    std::vector<int> matching_indices;

    // Loop through each value to check
    for (int i = 0; i < values_to_check.size(); i++) {
        // Loop through each data value to find a match
        for (int j = 0; j < data.size(); j++) {
            if (values_to_check[i] == data[j]) {
                matching_indices.push_back(j);
                //break;  // Found a match, move on to the next value to check
            }
        }
    }
    return matching_indices;
}

void copy_values(const std::vector<int>& indices, const std::vector<double>& source, std::vector<double>& dest) {
    // Loop through each index in the indices vector
    for (int i = 0; i < indices.size(); i++) {
        int index = indices[i];
        // Check if the index is within bounds of the source vector
        if (index < source.size()) {
            // Check if the destination vector has enough space to hold the new value
            if (i >= dest.size()) {
                // Throw an exception if the destination vector is not large enough
                throw std::out_of_range("Destination vector is not large enough to hold all values");
            }
            // Copy the value from the source vector to the destination vector
            dest[index] = source[index];
        }
        else
            throw std::out_of_range("Index out of range in source vector");
    }
}

bool checkVector(const std::vector<double>& input_vec) {
    for (const auto& val : input_vec) {
        if (std::isnan(val) || std::isinf(val)) {
            return false;  // vector contains NaN or Inf
        }
    }
    return true;  // vector does not contain NaN or Inf
}


std::vector<std::vector<double>> CreateSimplex(const std::vector<double>& input_coords, const std::vector<int>& fixed_idx) {
    const int N = input_coords.size();
    const double nonzdelt = 0.05;
    const double zdelt = 0.00025;
    std::vector<std::vector<double>> sim(N + 1, std::vector<double>(N));

    sim[0] = input_coords;
    for (int k = 0; k < N; ++k) {
        std::vector<double> y = input_coords;
        if (std::find(fixed_idx.begin(), fixed_idx.end(), k) == fixed_idx.end()) {
            if (y[k] != 0) {
                y[k] = (1 + nonzdelt) * y[k];
            } else {
                y[k] = zdelt;
            }
        }
        sim[k + 1] = y;
    }
    return sim;
}


// Function to create simplex from input coordinates
std::vector<std::vector<double>> CreateSimplexRand(const std::vector<double>& input_coords, const double std_dev){
    // Seed random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0, std_dev);

    // Determine length of input coordinates vector
    int ndim = input_coords.size();

    // Initialize simplex matrix
    std::vector<std::vector<double>> simplex(ndim+1, std::vector<double>(ndim));

    // Create variations of input coordinates and add to simplex
    simplex[0] = input_coords;
    for (int i = 1; i < ndim+1; i++) {
        for (int j = 0; j < ndim; j++) {
            if (input_coords[j] == 0. || input_coords[j] == 1.) {
                simplex[i][j] = input_coords[j];
            } else {
                simplex[i][j] = input_coords[j] + dist(gen);
                while (simplex[i][j] < 0. || simplex[i][j] > 1.) {
                    simplex[i][j] = input_coords[j] + dist(gen);
                }
            }
        }
    }

    // Sort simplex in non-decreasing order
    for (int i = 0; i < ndim+1; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = j+1; k < ndim; k++) {
                if (simplex[i][k] < simplex[i][j]) {
                    std::swap(simplex[i][j], simplex[i][k]);
                }
            }
        }
    }

    return simplex;
}

void FillDensityMatrix(std::vector<std::vector<double>>& density_matrix, const std::vector<double>& sample) {
    // Convert the sample to an arma::dvec
    arma::dvec sample_vec = arma::conv_to<arma::dvec>::from(sample);

    // Normalize the sample by the square root of its inner product
    double norm_factor = std::sqrt(arma::dot(sample_vec, sample_vec));
    arma::dvec normalized_vec = sample_vec / norm_factor;

    //int cnt = 0.;
    // Compute the outer product of the normalized sample and add the values to density_matrix
    for (arma::uword i = 0; i < normalized_vec.n_elem; i++) {
        for (arma::uword j = 0; j < normalized_vec.n_elem; j++) {
            density_matrix[i][j] += normalized_vec(i) * normalized_vec(j);
        }
    }
}

std::vector<double> GetEigenvectorAtIndex(const std::vector<std::vector<double>>& matrix, const int index) {
    // Convert matrix to Armadillo matrix
    arma::mat A(matrix.size(), matrix[0].size());
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            A(i, j) = matrix[i][j];
        }
    }
    // Divide matrix by its trace
    double trace = arma::trace(A);
    A /= trace;

    // Get eigenvalues and eigenvectors using Schur decomposition
    arma::mat Q, T;
    arma::schur(Q, T, A);

    // Extract eigenvector at given index
    arma::vec eigenvector = Q.col(index);

    // Convert eigenvector back to std::vector<double> and return
    std::vector<double> eigenvector_stdvec(eigenvector.begin(), eigenvector.end());

    return eigenvector_stdvec;
}

double AreaOutsideWindow(const std::vector<double>& amplitudes, const std::vector<double>& times, double t1, double t2) {

    double area = 0.0;
    double dt = times[1] - times[0]; // assume equal time step
    double t_edge_min = times.back();
    double t_edge_max = times.front();
    int i_edge_min = -1.;
    int i_edge_max = -1.;
    int n = amplitudes.size();
    for (int i = 0; i < n; i++) {
        if (times[i] < t1 || times[i] > t2) {
            area += amplitudes[i];
        }
        else{
            if(t_edge_min >= times[i]) {
                t_edge_min = times[i];
                i_edge_min = i;
            }
            if(t_edge_max < times[i]) {
                t_edge_max = times[i];
                i_edge_max = i;
            }
        }
    }
    //Add the points at the edge of the time window (shannon-nyquist interpolation)
    area += GetInterpolatedPoint(amplitudes, t1, 1/dt);
    area += GetInterpolatedPoint(amplitudes, t2, 1/dt);
    return area;
}

/*std::function<double(double)> WindowFitErr(const std::function<std::vector<double>(const std::vector<double>&)> &interp_function,
        //const std::tuple<double, double> time_edges,
                                           const std::vector<double> knot_vector,
                                           const double time_begin, const double time_end,
                                           const double window_size,
                                           const double step_interp) {

    return [&interp_function, knot_vector, time_begin, time_end, window_size, step_interp](const double t) {
        //std::cout << "t: " << t << endl;
        //std::cout << "time begin: " << std::get<0>(time_edges) << ", time end: " << std::get<1>(time_edges) << " and t: " << t << endl;
        if(t < time_begin || t > time_end-window_size)
            throw std::runtime_error("input time is out of range for given interpolation vector.");

        auto start_time = std::chrono::high_resolution_clock::now();
        //Get spline error
        std::vector<double> window_seg = ArangeVec(t, t + window_size, step_interp);
        std::vector<double> window_interp = interp_function(ArangeVec(t, t + window_size, step_interp));
        //auto end_time = std::chrono::high_resolution_clock::now();
        BSpline fit(knot_vector);
        fit.SplineLS(window_interp);
        double error = fit.GetError();
        //auto end_time = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        //std::cout << "Duration: " << duration.count()*1e-6 << " seconds" << std::endl;
        //interpolation before and after the time window used to calculate the penalty
        std::vector<double> pre_interp = interp_function(ArangeVec(time_begin+step_interp, t, step_interp));
        std::vector<double> post_interp = interp_function(ArangeVec(t+window_size+step_interp, time_end, step_interp));

        const double penalty = std::accumulate(post_interp.begin(), post_interp.end(), std::accumulate(pre_interp.begin(), pre_interp.end(), 0.));

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        std::cout << "Duration: " << duration.count()*1e-6 << " seconds" << std::endl;

        return error + penalty;
    };
}*/

std::function<double(double)> WindowFitErr(const std::vector<double> &knot_vector,
                                           const std::vector<double> &samples,
                                           const std::vector<double> &time,
                                           const double window_size,
                                           const double step_interp) {

    //std::vector<double> const time = LinSpaceVec(time_begin, time_end, samples.size());
    //auto fit = new BSpline(knot_vector);

    return [&samples, &knot_vector, &time, window_size, step_interp](const double t) {

        if(t < time.front() || t > time.back() - window_size)
            throw std::runtime_error("input time is out of range for given interpolation vector.");

        //std::vector<double> time = LinSpaceVec(time_begin, time_end, samples.size());

        std::vector<double> window_seg = ArangeVec(t, t + window_size, step_interp);
        std::vector<double> window_interp;// = interp_function(ArangeVec(t, t + window_size, step_interp));
        ShannonNyquistInterp(time, samples, window_seg, window_interp);

        BSpline fit(knot_vector);
        fit.SplineLS(window_interp);
        double error = fit.GetError();
        const double penalty = AreaOutsideWindow(samples, time, t, t + window_size);
        return error + penalty;
    };
}

double WindowFitErr(const double t,
                    const std::vector<double> knot_vector,
                    const std::vector<double> &samples,
                    const double time_begin, const double time_end,
                    const double window_size,
                    const double step_interp){

    if(t < time_begin || t > time_end-window_size)
        throw std::runtime_error("input time is out of range for given interpolation vector.");

    //auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<double> time_vec = LinSpaceVec(time_begin, time_end, samples.size());

    std::vector<double> window_seg = ArangeVec(t, t + window_size, step_interp);
    std::vector<double> window_interp;// = interp_function(ArangeVec(t, t + window_size, step_interp));
    ShannonNyquistInterp(time_vec, samples, window_seg, window_interp);

    //auto start_time = std::chrono::high_resolution_clock::now();

    BSpline fit(knot_vector);
    fit.SplineLS(window_interp);
    double error = fit.GetError();

    //auto end_time = std::chrono::high_resolution_clock::now();

    /*std::vector<double> pre_interp;
    ShannonNyquistInterp(time_vec, samples, ArangeVec(time_begin+step_interp, t, step_interp), pre_interp);
    std::vector<double> post_interp;
    ShannonNyquistInterp(time_vec, samples, ArangeVec(t+window_size+step_interp, time_end, step_interp), post_interp);
*/

    const double penalty = AreaOutsideWindow(samples, time_vec, t, t + window_size);//std::accumulate(post_interp.begin(), post_interp.end(), std::accumulate(pre_interp.begin(), pre_interp.end(), 0.));

    //auto end_time = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    //std::cout << "Duration: " << duration.count()*1e-6 << " seconds" << std::endl;

    return error + penalty;
}

double KnotSeparationPenalty(const std::vector<double>& knot_vector, const double stddev) {
    double sum = 0.0;
    double norm = std::sqrt(2*M_PI*pow(stddev, 2));

    for (int i = 0; i < knot_vector.size() - 1; i++) {
        double diff = knot_vector[i + 1] - knot_vector[i];
        sum += exp(-pow(diff/(2*stddev), 2))/norm;
    }
    return sum;
}

std::function<double(std::vector<double>)> BSplineErr(const std::vector<double>& samples, const double stddev) {
    return [&samples, stddev](const std::vector<double>& knot_vector) {
        BSpline fit(knot_vector);
        fit.SplineLS(samples);
        double error = fit.GetError();
        return error + KnotSeparationPenalty(knot_vector, stddev);
    };
}

std::vector<double> MinimizeKnots(const std::vector<double> &knot_vector, const std::vector<double> &samples,
                                  const int max_iter, const double stddev, const bool isRand){

    int n_dim = knot_vector.size();
    auto error_func = BSplineErr(samples, stddev);

    std::vector<int> fixed_idx = FindMatchingIndices(knot_vector, {knot_vector.front(), knot_vector.back()});

    std::vector<std::vector<double>> init_points;
    if(isRand)
        init_points = CreateSimplexRand(knot_vector, stddev);
    else
        init_points = CreateSimplex(knot_vector, fixed_idx);

    NelderMead fit(n_dim, init_points, fixed_idx);
    std::vector<double> optimized_knots = fit.optimize(error_func, max_iter);

    return optimized_knots;
}

double goldenSectionSearch(const std::function<double(double)> &f, double& a, double& b, const double tol) {
    const double phi = (1.0 + sqrt(5.0)) / 2.0; // golden ratio
    double c = b - (b - a) / phi; // first intermediate point
    double d = a + (b - a) / phi; // second intermediate point
    double fc = f(c); // function value at first intermediate point
    double fd = f(d); // function value at second intermediate point

    while (std::abs(b - a) > tol) {
        if (fc < fd) { // update bounds and intermediate points
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) / phi;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) / phi;
            fd = f(d);
        }
    }
    return (a + b) / 2.0; // return midpoint of final bracket as estimate
}


double BrentsMethod(const std::function<double(double)> &func, double a, double b, const double tol, const int max_iter) {

    if(b < a)
        std::swap(a,b);

    const double eps = std::numeric_limits<double>::epsilon(); // machine epsilon
    const double golden_ratio = (1 + std::sqrt(5)) / 2;
    const double rho = 1 - (1/golden_ratio);

    double c = b - rho * (b - a);
    double d = a + rho * (b - a);

    // Set initial values
    double x = 0.5 * (d - c);
    double w = x;
    double v = x;
    double fx = func(x);
    double fw = fx;
    double fv = fx;
    double e = 0.0;
    double u;

    // Main loop
    for (int i = 0; i <= max_iter; ++i) {
        //std::cout << x << std::endl;
        double xm = 0.5 * (a + b);
        double tol1 = tol * std::abs(x) + eps;
        double tol2 = 2.0 * tol1;

        if (std::abs(x - xm) <= tol2 - 0.5 * (b - a)) {
            // Convergence
            //std::cout << "Converged in " << i+1 << " iterations." << std::endl;
            return x;
        }

        if (fabs(e) > tol1) { // check if we can use inverse quadratic interpolation
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) {
                p = -p;
            }
            q = fabs(q);
            double e_temp = e;
            e = d;
            if (fabs(p) >= fabs(0.5 * q * e_temp) || p <= q * (a - x) || p >= q * (b - x)) {
                // use bisection instead
                e = (x >= xm) ? a - x : b - x;
                d = rho*e;
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = std::copysign(tol1, xm - x);
            }
        }
        else{
            e = (x >= xm) ? a - x : b - x;
            d = rho*e;
        }

        u = (abs(d) >= tol1 ? x + d : x + std::copysign(tol1, d));
        double fu = func(u);

        if (fu <= fx) {
            if (u >= x) {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    // Maximum number of iterations exceeded
    throw std::runtime_error("Brent's method failed to converge");
}

void ShannonNyquistInterp(const std::vector<double>& xdata, const std::vector<double>& ydata,
                          const std::vector<double>& xnew, std::vector<double>& ynew) {

    const int N = xdata.size();
    const double dx = xdata.back() - xdata.front();
    const double freq = static_cast<double>(N-1)/dx;

    //std::vector<double> y(x.size());
    if(ynew.size() == 0. || ynew.size() != xnew.size())
        ynew.assign(xnew.size(), 0.);

    for (size_t i = 0; i < xnew.size(); i++) {
        ynew[i] = GetInterpolatedPoint(ydata, xnew[i] - xdata.front(), freq);

    }
}

// Function to compute the sinc interpolation of a signal using FFTW3
std::vector<double> sincInterp(std::vector<double> x, std::vector<double> y, std::vector<double> t){

    // Determine the size of the input and output vectors
    const int N = x.size();
    const int M = t.size();

    // Compute the FFT of the input signal x
    fftw_complex* X = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_r2c_1d(N, x.data(), X, FFTW_ESTIMATE);
    fftw_execute(p);

    // Compute the frequency vector
    std::vector<double> f(N);
    for (int i = 0; i < N; i++) {
        if (i < N / 2) {
            f[i] = i / (N * (y[1] - y[0]));
        } else {
            f[i] = (i - N) / (N * (y[1] - y[0]));
        }
    }

    // Compute the interpolation factor
    std::vector<double> h(N);
    for (int i = 0; i < N; i++) {
        if (std::abs(f[i]) > 1) {
            h[i] = 0;
        } else {
            h[i] = std::sin(M_PI * f[i]) / (M_PI * f[i]);
        }
    }

    // Compute the inverse FFT of the product of the FFT of the input signal and the interpolation factor
    fftw_complex* Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for (int i = 0; i < N; i++) {
        Y[i][0] = X[i][0] * h[i];
        Y[i][1] = X[i][1] * h[i];
    }
    fftw_plan q = fftw_plan_dft_c2r_1d(N, Y, t.data(), FFTW_ESTIMATE);
    fftw_execute(q);

    // Normalize the output signal
    const double scale = (y[1] - y[0]) / (t[1] - t[0]);
    for (int i = 0; i < M; i++) {
        t[i] = t[i] * scale;
    }

    // Free the memory used by the FFTW3 plans and arrays
    fftw_destroy_plan(p);
    fftw_free(X);
    fftw_destroy_plan(q);
    fftw_free(Y);

    return t;
}

std::function<std::vector<double>(const std::vector<double>&)> ShannonNyquistInterp(
        const std::vector<double>& xdata,
        const std::vector<double>& ydata) {

    const int N = xdata.size();
    const double dx = xdata.back() - xdata.front();
    const double freq = (N-1)/dx;

    return [&xdata, &ydata, N, freq](const std::vector<double>& x)->std::vector<double>{
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i)
            y[i] = GetInterpolatedPoint(ydata, x[i]-xdata.front(), freq);
        return y;
    };
}

std::function<std::vector<double>(const std::vector<double>&)> ShannonNyquistInterpFFT(
        const std::vector<double>& x,
        const std::vector<double>& y) {

    auto sinc = [](double x){
        if(fabs(x) < 1e-8)
            return 1.;
        else
            return sin(x) / x;
    };

    int n = x.size();
    int n_interp = n * 10;
    double dx = (x.back() - x.front()) / (n - 1);
    double x0 = x.front();
    double nyquist = M_PI / dx;
    vector<double> y_interp(n_interp);
    vector<double> freq(n_interp);
    for (int i = 0; i < n_interp; i++) {
        double xf = x0 + i * dx / 10;
        freq[i] = xf <= nyquist ? xf : -1;
    }
    fftw_complex *in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n_interp));
    fftw_complex *out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n_interp));
    fftw_plan plan = fftw_plan_dft_1d(n_interp, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (int i = 0; i < n_interp; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) {
            s += y[j] * sinc(M_PI * (freq[i] - x[j]) / dx);
        }
        in[i][0] = s / n;
        in[i][1] = 0.0;
    }
    fftw_execute(plan);
    for (int i = 0; i < n_interp; i++) {
        y_interp[i] = out[i][0];
    }
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return [=](const vector<double>& x_values) -> std::vector<double>{
        int n_values = x_values.size();
        vector<double> y_values(n_values);
        for (int i = 0; i < n_values; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                s += y[j] * sinc(M_PI * (x_values[i] - x[j]) / dx);
            }
            y_values[i] = s;
        }
        return y_values;
    };
}

void ShannonNyquistInterpFFT(const std::vector<double>& xdata, const std::vector<double>& ydata,
                             const std::vector<double>& xnew, std::vector<double>& ynew) {

    // Calculate the number of data points and new points
    const size_t N = xdata.size();
    const size_t M = xnew.size();

    ynew.assign(M, 0.);

    // Calculate the sampling frequency and the Nyquist frequency
    const double Fs = 1.0 / (xdata[1] - xdata[0]);
    const double Nyq = Fs / 2.0;

    // Create an FFT plan for the input and output arrays
    fftw_plan input_plan = fftw_plan_r2r_1d(N, const_cast<double*>(ydata.data()), const_cast<double*>(ydata.data()), FFTW_R2HC, FFTW_ESTIMATE);
    fftw_plan output_plan = fftw_plan_r2r_1d(M, const_cast<double*>(ynew.data()), const_cast<double*>(ynew.data()), FFTW_HC2R, FFTW_ESTIMATE);

    // Allocate memory for the FFT input and output arrays
    double* input_array = fftw_alloc_real(N);
    double* output_array = fftw_alloc_real(M);

    // Copy the input data into the FFT input array
    std::copy(ydata.begin(), ydata.end(), input_array);

    // Execute the FFT on the input data
    fftw_execute(input_plan);

    // Calculate the frequency axis for the FFT output
    std::vector<double> freq_axis(int(N / 2) + 1);
    const double df = Fs / N;
    for (size_t i = 0; i < freq_axis.size(); i++){
        freq_axis[i] = i * df;
    }

    // Interpolate the FFT output at the new frequencies
    for (size_t i = 0; i < M; i++){
        // Calculate the corresponding frequency index for the new frequency
        const double f = xnew[i] * Nyq;
        const size_t index = std::floor(f / df + 0.5);
        std::cout << f << " " << index << std::endl;
        // Check that the index is within the FFT frequency range
        if (index >= freq_axis.size()){
            ynew[i] = 0.0;
        }
        else{
            // Interpolate the FFT output at the new frequency
            const double A = std::sqrt(std::pow(input_array[0], 2) / 2.0);
            const double B = input_array[index] / std::sqrt(N);
            const double C = std::sqrt(std::pow(input_array[N - index], 2) / 2.0);
            ynew[i] = 2.0 * (A + B + C);
            std::cout << ynew[i] << std::endl;
        }
    }

    // Execute the inverse FFT on the interpolated FFT output
    fftw_execute(output_plan);

    // Normalize the output of the inverse FFT
    const double norm_factor = 1.0 / M;
    for (size_t i = 0; i < M; i++){
        ynew[i] *= norm_factor;
    }

    // Free memory and destroy the FFT plans
    fftw_destroy_plan(input_plan);
    fftw_destroy_plan(output_plan);
    fftw_free(input_array);
    fftw_free(output_array);
}


/*
void ShannonNyquistInterpFFT(const std::vector<double>& xdata, const std::vector<double>& ydata,
                             const std::vector<double>& xnew, std::vector<double>& ynew) {

    // Determine the number of data points
    const int n = xdata.size();

    // Calculate the spacing between data points
    const double dx = xdata[1] - xdata[0];

    // Create a FFTW3 plan for the input data
    fftw_plan plan = fftw_plan_r2r_1d(n, &ydata[0], &ynew[0], FFTW_R2HC, FFTW_ESTIMATE);

    // Create temporary arrays for the FFT and inverse FFT
    double *fft_data = new double[n];
    double *ifft_data = new double[n];

    // Initialize the output array to zero
    ynew.assign(xnew.size(), 0.0);

    // Loop over the output data points
    for (int i = 0; i < xnew.size(); i++) {

        // Calculate the FFT of the input data
        for (int j = 0; j < n; j++) {
            fft_data[j] = ydata[j] * exp(-2.0 * M_PI * xdata[j] * xnew[i] * I);
        }

        fftw_execute(plan, fft_data, ifft_data);

        // Scale the result and store it in the output array
        ynew[i] = creal(ifft_data[0]) / n;
    }

    // Clean up the FFTW3 plan and temporary arrays
    fftw_destroy_plan(plan);
    delete[] fft_data;
    delete[] ifft_data;
}
*/

/*
void ShannonNyquistInterpFFT(const std::vector<double>& xdata,
                             const std::vector<double>& ydata,
			     const std::vector<double>& x_interp,
                             std::vector<double>& interpolation){

  auto sinc = [](double x){
    if(fabs(x) < 1e-8)
      return 1.;
    else
      return sin(x) / x;
  };

  int n = x.size();
  int n_interp = n * 10;
  double dx = (x.back() - x.front()) / (n - 1);
  double x0 = x.front();
  double nyquist = M_PI / dx;
  vector<double> y_interp(n_interp);
  vector<double> freq(n_interp);
  for (int i = 0; i < n_interp; i++) {
    double xf = x0 + i * dx / 10;
    freq[i] = xf <= nyquist ? xf : -1;
  }
  fftw_complex *in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n_interp));
  fftw_complex *out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n_interp));
  fftw_plan plan = fftw_plan_dft_1d(n_interp, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  for (int i = 0; i < n_interp; i++) {
    double s = 0.0;
    for (int j = 0; j < n; j++) {
      s += y[j] * sinc(M_PI * (freq[i] - x[j]) / dx);
    }
    in[i][0] = s / n;
    in[i][1] = 0.0;
  }
  fftw_execute(plan);
  for (int i = 0; i < n_interp; i++) {
    y_interp[i] = out[i][0];
  }
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  int n_values = x_interp.size();
  interpolation.clear();
  //vector<double> y_values(n_values);
  for (int i = 0; i < n_values; i++) {
    double s = 0.0;
    for (int j = 0; j < n; j++) {
      s += y[j] * sinc(M_PI * (x_values[i] - x[j]) / dx);
    }
    interpolation.push_back(s);
  }
}
*/
std::function<std::vector<double>(const std::vector<double>&)> LinearInterp(
        const std::vector<double>& xdata,
        const std::vector<double>& ydata) {

    const int N = xdata.size();
    const double dx = xdata.back() - xdata.front();
    const double dt = dx / (N - 1);

    return [xdata, ydata, N, dx, dt](const std::vector<double>& x)->std::vector<double>{
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            double t = (x[i] - xdata.front()) / dx;
            int k = static_cast<int>(t * (N - 1));
            double f = t * (N - 1) - k;
            y[i] = (1 - f) * ydata[k] + f * ydata[k + 1];
        }
        return y;
    };
}

std::function<std::vector<double>(const std::vector<double>&)> CubicInterpolation(
        const std::vector<double>& xdata,
        const std::vector<double>& ydata){

    // Check that the input arrays have the same size
    if (xdata.size() != ydata.size())
        throw std::invalid_argument("x and y arrays must have the same size");


    // Compute the slopes
    std::vector<double> h(xdata.size() - 1), b(xdata.size() - 1), u(xdata.size() - 1), v(xdata.size() - 1);
    for (size_t i = 0; i < xdata.size() - 1; i++){
        h[i] = xdata[i+1] - xdata[i];
        b[i] = (ydata[i+1] - ydata[i]) / h[i];
    }

    // Compute the u and v vectors
    u[0] = 2 * (h[0] + h[1]);
    v[0] = 6 * (b[1] - b[0]);
    for (size_t i = 1; i < xdata.size() - 2; i++){
        u[i] = 2 * (h[i] + h[i-1]) - (h[i-1] * h[i-1]) / u[i-1];
        v[i] = 6 * (b[i] - b[i-1]) - (h[i-1] * v[i-1]) / u[i-1];
    }

    // Compute the coefficients of the cubic polynomials
    std::vector<double> z(xdata.size(), 0);
    z[xdata.size()-1] = 0;
    for (int i = xdata.size() - 3; i >= 0; i--)
        z[i+1] = (v[i] - h[i] * z[i+2]) / u[i];

    z[0] = 0;

    // Define the interpolation function
    return [xdata, ydata, z, h](const std::vector<double>& xs) -> std::vector<double> {
        std::vector<double> ys(xs.size(), 0);
        for (size_t i = 0; i < xs.size(); i++){
            // Find the index of the largest x value less than xs[i]
            int j = std::distance(xdata.begin(), std::lower_bound(xdata.begin(), xdata.end(), xs[i])) - 1;

            // Compute the coefficients of the cubic polynomial
            double A = z[j+1] - z[j];
            double B = -h[j] * (z[j+1] + 2*z[j]) + 3*(ydata[j+1] - ydata[j]);
            double C = h[j] * z[j];
            double D = ydata[j];

            // Compute the interpolated value
            double dx = xs[i] - xdata[j];
            ys[i] = A*dx*dx*dx + B*dx*dx + C*dx + D;
            ys[0] = ydata[0];
            ys[ys.size()-1] = ydata[ydata.size()-1];
        }
        return ys;
    };
}


// Compute the coefficients of the cubic polynomials for cubic spline interpolation
void computeSplineCoefficients(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& a,
                               std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
    // Initialize the coefficients as vectors of zeros
    a = std::vector<double>(x.size(), 0.0);
    b = std::vector<double>(x.size(), 0.0);
    c = std::vector<double>(x.size(), 0.0);
    d = std::vector<double>(x.size(), 0.0);

    // Compute the h_i and delta_i values
    std::vector<double> h(x.size() - 1);
    std::vector<double> delta(x.size() - 1);
    for (int i = 0; i < x.size() - 1; i++) {
        h[i] = x[i + 1] - x[i];
        delta[i] = (y[i + 1] - y[i]) / h[i];
    }

    // Compute the tridiagonal system of equations
    std::vector<double> alpha(x.size() - 1);
    std::vector<double> beta(x.size() - 1);
    alpha[0] = 0.0;
    beta[0] = 0.0;
    for (int i = 1; i < x.size() - 1; i++) {
        alpha[i] = h[i] / (h[i - 1] + h[i] - alpha[i - 1] * h[i - 1]);
        beta[i] = (delta[i] - beta[i - 1] * h[i - 1]) / (h[i - 1] + h[i] - alpha[i - 1] * h[i - 1]);
    }

    // Solve for the coefficients c_i
    c[x.size() - 1] = 0.0;
    for (int i = x.size() - 2; i >= 0; i--) {
        c[i] = alpha[i] * c[i + 1] + beta[i];
    }

    // Compute the remaining coefficients a_i, b_i, and d_i
    for (int i = 0; i < x.size() - 1; i++) {
        a[i] = y[i];
        b[i] = delta[i] - c[i] * h[i];
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    }
}

// Perform cubic spline interpolation on the given x and y vectors and return a lambda function that can interpolate new x values
std::function<std::vector<double>(const std::vector<double>&)> cubicSplineInterpolate(const std::vector<double>& x, const std::vector<double>& y) {
    // Compute the coefficients of the cubic polynomials for cubic spline interpolation
    std::vector<double> a, b, c, d;
    computeSplineCoefficients(x, y, a, b, c, d);

    // Return a lambda function that performs the interpolation
    return [a, b, c, d, x](const std::vector<double>& xi) -> std::vector<double> {
        std::vector<double> yi(xi.size());

        for (int i = 0; i < xi.size(); i++) {
            int j = 0;
            while (j < x.size() - 1 && x[j+1] < xi[i]) {
                j++;
            }

            double dx = xi[i] - x[j];
            yi[i] = a[j] + b[j] * dx + c[j] * dx * dx + d[j] * dx * dx * dx;
        }

        return yi;
    };
}
