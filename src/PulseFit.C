#include <vector>
#include <sstream>
#include <algorithm>
#include <functional>
#include <unistd.h>

#include <TChain.h>
#include <TKey.h>
#include <iostream>

#include "RootInterface.h"
#include "Plot.h"
#include "PulseTools.h"

using namespace std;

int main(int argc, char* argv[]) {

    std::string filename;
    int num_entries = 1000;
    int interpolation_size = 100;
    double time_window = 30;//ns
    double cfd_fraction = -1.0;

    int c;
    while ((c = getopt(argc, argv, "i:n:s:t:c:")) != -1) {
        switch (c) {
            case 'i':
                filename = optarg;
                break;
            case 'n':
                num_entries = std::stoi(optarg);
                break;
            case 's':
                interpolation_size = std::stoi(optarg);
                break;
            case 't':
                time_window = std::stod(optarg);
                break;
            case 'c':
                cfd_fraction = std::stod(optarg);
                break;
            case '?':
                std::cerr << "Unknown option: -" << optopt << std::endl;
                return 1;
            default:
                std::cerr << "Unexpected error" << std::endl;
                return 1;
        }
    }

    auto tch = new TChain("data");
    tch->Add(filename.c_str());

    auto data = new RootInterface(tch, 1);

    int sample_size = data->GetSampleSize();
    std::vector<double> waveform(sample_size);
    std::vector<double> time;
    std::vector<double> interpolated_waveform(interpolation_size, 0.);
    std::vector<double> interpolation_time;
    //std::vector<double> knot_vector = {0, 0.0566526, 0.0706713, 0.118964, 0.13563, 0.257123, 0.507515, 1, 1, 1, 1 };
    std::vector<double> knot_vector = {0., 0.05, 0.085, 0.1, 0.14, 0.2, 0.4, 1., 1., 1., 1.};
    //std::vector<double> knot_vector = {0, 0.0507, 0.0862, 0.09, 0.142, 0.203, 0.405, 1, 1, 1, 1};
    //std::vector<double> knot_vector = {0, 0.0089733, 0.0609578, 0.101097, 0.153922, 0.305948, 0.490758, 1, 1, 1, 1};
    std::vector<double> spline;

    std::vector<std::vector<double> > density_matrix(interpolation_size, std::vector<double>(interpolation_size, 0.));

    auto start_time = std::chrono::high_resolution_clock::now();
    int count = 0;
    for(int e = 0; e < tch->GetEntries(); e++){

        data->fChain->GetEntry(e);

        //double horiz_offset = data->horiz_offset;
        double horiz_scale = data->horiz_scale;
        double trig_offset = data->trig_offset;
        //double trig_time = data->trig_time;
        double vert_offset = data->vert_offset;
        double vert_scale = data->vert_scale;

        for(int i = 0; i < sample_size; i++)
            waveform[i] = -(-vert_offset+double(data->samples[i])*vert_scale);

        double time_begin = trig_offset*1e9;
        double time_end = (horiz_scale*(sample_size-1)-trig_offset)*1e9;
        time = LinSpaceVec(time_begin, time_end, sample_size);
        auto max_elem = std::max_element(waveform.begin(), waveform.end());
        auto min_elem = std::min_element(waveform.begin(), waveform.end());
        int max_index = std::distance(waveform.begin(), max_elem);
        int min_index = std::distance(waveform.begin(), min_elem);
        double max_time = time[max_index];

        PedestalSubtraction(waveform, int(waveform.size()/3));

        //skip bad waveforms, noise and saturated samples
        if(waveform[max_index] > 0.65 ||  waveform[max_index] < 0.1 || fabs(waveform[min_index]) > fabs(waveform[max_index]))
            continue;

        //define interpolation function to use
        auto interp_function = ShannonNyquistInterp(time, waveform);

        //make interpolation for cfd time
        double cfd_time = -1.;
        if(cfd_fraction > 0.05) {
            interpolation_time = LinSpaceVec(max_time - 4., max_time + 10., interpolation_size);
            interpolated_waveform = interp_function(interpolation_time);
            cfd_time = CalculateCFDTime(interpolated_waveform, interpolation_time, cfd_fraction);
        }

        //time fit over window
        /*auto time_func = WindowFitErr(knot_vector, waveform, time,
                                      time_window, 0.15);
        double fit_time = BrentsMethod(time_func, time_begin, time_end-time_window);*/

        //make interpolation window at CFD time or time from fit
        if(!isnan(cfd_time) && cfd_fraction > 0.)
            interpolation_time = LinSpaceVec(cfd_time - 2., cfd_time + time_window - 2., interpolation_size);
        else if(cfd_fraction < 0.) {
            auto time_func = WindowFitErr(knot_vector, waveform, time,
                                          time_window, 0.15);
            double fit_time = BrentsMethod(time_func, time_begin, time_end - time_window);
            interpolation_time = LinSpaceVec(fit_time, fit_time + time_window, interpolation_size);
        }
        else
            continue;

        bool noNans = !std::any_of(interpolation_time.begin(), interpolation_time.end(),
                                   [](double d) { return std::isnan(d); });

        if(noNans){
            interpolated_waveform = interp_function(interpolation_time);
            FillDensityMatrix(density_matrix, interpolated_waveform);
            count++;
        }

        //print progress to command line
        if (count % 100 == 0) {
            fprintf(stdout, "\r  Processed events: %8d of %8d ", count, num_entries);
        }
        fflush(stdout);
        if(count == num_entries)
            break;
    }

    cout << "\nProcessed " << count << " waveforms\n" << endl;

    vector<double> lead_eigenvector = GetEigenvectorAtIndex(density_matrix);

    //if eigenvector is upside down, flip it
    auto max_elem = std::max_element(lead_eigenvector.begin(), lead_eigenvector.end());
    auto min_elem = std::min_element(lead_eigenvector.begin(), lead_eigenvector.end());
    if(fabs(*min_elem) > fabs(*max_elem))
        std::transform(lead_eigenvector.begin(), lead_eigenvector.end(), lead_eigenvector.begin(),
                       [](double d) { return -1.0 * d; });

    //subtract first element from entire sample
    std::transform(lead_eigenvector.begin(), lead_eigenvector.end(), lead_eigenvector.begin(),
                   [&lead_eigenvector](double d) { return d - (lead_eigenvector[0]); });

    knot_vector = MinimizeKnots(knot_vector, lead_eigenvector, 1e3);

    BSpline bsp(knot_vector);
    spline = bsp.SplineLS(lead_eigenvector);

    cout << "Minimized knot vector: ";
    for(auto elem : knot_vector)
        cout << elem << " ";
    cout << endl;

    PlotSplineFit(lead_eigenvector, spline, knot_vector, "pulse");

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Duration: " << duration.count()*1e-6 << " seconds" << std::endl;

    return 0;
}
