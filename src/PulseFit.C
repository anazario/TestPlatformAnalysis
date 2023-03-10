#include <vector>
#include <sstream>
#include <algorithm>
#include <functional>

#include <TChain.h>
#include <TKey.h>
#include <iostream>

#include "Pulse.h"
#include "RootInterface.h"
#include "Plot.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: print_branches <filename.root> <max_entries>" << std::endl;
        return 1;
    }

    const char* filename = argv[1];
    const int nEntries = std::atoi(argv[2]);

    auto tch = new TChain("data");
    tch->Add(filename);

    auto data = new RootInterface(tch, 1);

    int interp_size = 200;
    int sample_size = data->GetSampleSize();
    std::vector<double> waveform(sample_size);
    std::vector<double> time;
    std::vector<double> interpolated_waveform(interp_size, 0.);
    std::vector<double> interpolation_time;
    //std::vector<double> knot_vector = {0, 0.0566526, 0.0706713, 0.118964, 0.13563, 0.257123, 0.507515, 1, 1, 1, 1 };
    //std::vector<double> knot_vector = {0., 0.05, 0.085, 0.1, 0.14, 0.2, 0.4, 1., 1., 1., 1.};
    std::vector<double> knot_vector = {0, 0.0507, 0.0862, 0.09, 0.142, 0.203, 0.405, 1, 1, 1, 1};
    std::vector<double> spline;

    std::vector<std::vector<double> > density_matrix(interp_size, std::vector<double>(interp_size, 0.));

    for(int e = 0; e < nEntries; e++){

        if (e % 100 == 0) {
            fprintf(stdout, "\r  Processed events: %8d of %8d ", e, nEntries);
        }
        fflush(stdout);

        data->fChain->GetEntry(e);

        double horiz_offset = data->horiz_offset;
        double horiz_scale = data->horiz_scale;
        double trig_offset = data->trig_offset;
        double trig_time = data->trig_time;
        double vert_offset = data->vert_offset;
        double vert_scale = data->vert_scale;

        for(int i = 0; i < sample_size; i++)
            waveform[i] = -(-vert_offset+double(data->samples[i])*vert_scale);

        double time_begin = trig_offset*1e9;
        double time_end = (horiz_scale*(sample_size-1)-trig_offset)*1e9;
        time = LinSpaceVec(time_begin, time_end, sample_size);
        auto max_elem = std::max_element(waveform.begin(), waveform.end());
        double max_time = time[std::distance(waveform.begin(), max_elem)];

        PedestalSubtraction(waveform, int(waveform.size()/3));

        auto interp_function = ShannonNyquistInterp(time, waveform);

        //make interpolation for cfd time
        interpolation_time = LinSpaceVec(max_time-10., max_time+10., interp_size);
        //ShannonNyquistInterp(time, waveform, interpolation_time, interpolated_waveform);
        interpolated_waveform = interp_function(interpolation_time);
        double cfd_time = CalculateCFDTime(interpolated_waveform, interpolation_time, 0.4);

        //make interpolation window at CFD time
        if(!isnan(cfd_time)) {
            interpolation_time = LinSpaceVec(cfd_time - 2., cfd_time + 23., interp_size);
            //ShannonNyquistInterp(time, waveform, interpolation_time, interpolated_waveform);
            interpolated_waveform = interp_function(interpolation_time);
            FillDensityMatrix(density_matrix, interpolated_waveform);
        }

        double time_window = 25;//ns
        auto time_func = WindowFitErr(interp_function, knot_vector, time_begin, time_end, time_window, 0.01);
        std::vector<double> time_test = LinSpaceVec(time_begin, time_end-time_window, 100);
    }

    std::cout << "Finished!" << std::endl;
    std::vector<double> lead_eigvector = GetEigenvectorAtIndex(density_matrix);
    std::cout << "here" << std::endl;
    knot_vector = MinimizeKnots(knot_vector, lead_eigvector, 1e3);

    BSpline bsp(knot_vector);
    spline = bsp.SplineLS(lead_eigvector);

    for(auto elem : knot_vector)
        std::cout << elem << " ";
    std::cout << std::endl;

    PlotSplineFit(lead_eigvector, spline, knot_vector, "pulse");

    knot_vector = {0, 0.0506818, 0.0861591, 0.09, 0.141909, 0.202727, 0.405455, 1, 1, 1, 1};

    double stddev = 1e-3;
    auto obj_func = BSplineErr(lead_eigvector, stddev);

    double val = obj_func(knot_vector);
    double temp = -999;

    for(int i = 0; i < 1e3; i++){
        temp = obj_func(knot_vector);
        if(temp > val || temp < val)
            std::cout << temp << std::endl;
    }

    return 0;
}
