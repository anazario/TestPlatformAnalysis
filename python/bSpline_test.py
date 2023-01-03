import numpy as np
import matplotlib.pyplot as plt
import sys
import ROOT

from pulseTools import *
import pulse as pls
import BSpline as bsp

if __name__ == '__main__':
    import h5py
    import optparse
    from scipy import optimize
    from scipy.linalg import *
    from scipy import interpolate

    usage = "usage: %prog [-n]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default=None)
    parser.add_option("--ch", type=int, dest="channel",
                      help="channel to process (between 1 and 8)", default=1)
    parser.add_option("-o", type=int, dest="order",
                      help="polynomial order for spline basis", default=3)
    (options, args) = parser.parse_args()

    order = options.order
    channel = options.channel
    input_file = options.infile
    
    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    time_data = f[f'ch{channel}_time']
    fixed_time_array = f[f'ch{channel}_fixed_time_array']
    density_matrix = f[f'ch{channel}_density_matrix']
    eig_value_array, eig_vector_array = schur(density_matrix)

    lead_eigvector = eig_vector_array[:,0]
    eig_value_array = np.array([val[idx] for idx, val in enumerate(eig_value_array)])
    if np.max(lead_eigvector) < -np.min(lead_eigvector):
        eig_vector_array = -1*eig_vector_array.T
    else:
        eig_vector_array = eig_vector_array.T

    lead_eigvector = eig_vector_array[0]

    pulse = pls.Pulse(lead_eigvector, 0., 1.)
    max_pos = pulse.get_time_at_max()
    
    #knot_vector = list(np.linspace(0., 1., 8))
    #knot_vector = [0., 0., 0., 0., 0.14285714285714285, 0.2857142857142857, 0.42857142857142855, 0.5714285714285714, 0.7142857142857142, 0.8571428571428571, 1.0, 1.0, 1.0, 1.0]
    #knot_vector = [0., 0., 0., 0., 0.05, 0.085, 0.1, 0.14, 0.2, 0.4, 1., 1., 1., 1.]
    knot_vector = [0., 0.05, 0.085, 0.1, 0.14, 0.2, 0.4, 1., 1., 1., 1.]
    #knot_vector = [0., 0.05, 0.085, 0.1, 0.2, 0.4, 1., 1., 1., 1.]
    #knot_vector = [0., 0.05546708, 0.08395173, 0.16736349, 1., 1., 1., 1.]
    #knot_vector = [0., 0.05, 0.085, 0.1, 0.14, 0.2, 0.4, 1.]
    #knot_vector_test = [0, 1, 1, 3, 4, 6, 6, 6]
    #knot_vector_test = [0, 0, 1, 2, 3,3]
    #knot_vector = [0., 0.07, 0.1, 0.15, 0.3, 0.5, 0.75, 1.]
    #knot_size = len(knot_vector)
    #knot_vector = [0., 0., 0., 0., 0.04634193, 0.0820082, 0.09746398, 0.12442732, 0.19711181, 0.68400499, 1., 1., 1., 1.]
    '''
    for k in range(3):
        knot_vector.append(knot_vector[-1])
    print(knot_vector)
    
    bounds = []
    for idx, knot in enumerate(knot_vector):
        if knot == knot_vector[0] or knot == knot_vector[-1]:
            bounds.append((knot, knot))
        else:
            #bounds.append((knot_vector[0]+1e-5, knot_vector[idx+1]-1e-5))
            bounds.append((knot_vector[0], knot_vector[-1]))
            
    print(bounds)
    '''    
    #plot_basis('test_knot_vector', knot_vector_test, order=2)

    spline_size = 200
    domain = np.linspace(0., 1., spline_size)

    c_list = [[] for i in range(len(knot_vector)-4)]
    err_list = []
    
    for trig, waveform in enumerate(samples):

        if trig%1000 == 0:
            sys.stdout.write(f"\rProcessing waveform {trig}")
            sys.stdout.flush()
            
        if np.any(waveform):# and fixed_time_array[trig] > -1.0:
            pulse = pls.Pulse(waveform, time_data[trig][0], time_data[trig][-1])
            max_amplitude = pulse.get_max_amplitude()
            
            interpolated_function = interpolate.interp1d(x=np.linspace(0, 1, len(waveform)), y=waveform, kind='cubic')
            interpolated_waveform = interpolated_function(domain)
            fit = bsp.BSpline(interpolated_waveform, knot_vector, precision=1e-8)
            error = fit.get_error()
            fit_coeff = fit.get_fit_coeff()

            err_list.append(error)
            
            for i in range(len(fit_coeff)):
                c_list[i].append(fit_coeff[i]/max_amplitude)

            if error > 0.05 and fit_coeff[0] > 1.:
                fit.plot_fit(f'pulse_plots/fit_pulse{trig}')
            '''
            else:
                fit.plot_fit(f'pulse_plots/fit_good_pulse{trig}')
            '''
                
    for idx, ls in enumerate(c_list):
        plot_hist(ls, bins=100, range_tuple=(np.min(ls)-1e-3, np.max(ls)+1e-4), name=f'coeff{idx}_refch1_hist_optimized', xlabel='Coefficient Value', ylabel='Entries', title='')
        plot_scatter(err_list, ls, name=f'coeff{idx}_vs_error_refch1_optimized', xlabel='Least Squares Error', ylabel='Coefficient Value')
    plot_hist(err_list, bins=100, range_tuple=(np.min(err_list)-1e-3, np.max(err_list)+1e-4), name=f'fit_error_refch1_hist_optimized', xlabel='Error', ylabel='Entries', title='')
   
    #plot_basis('test_LGAD', knot_vector, 3, 1000, 0.001)
    
    '''
    bounds = []
    for idx, knot in enumerate(knot_vector):
        if knot == knot_vector[0] or knot == knot_vector[-1]:
            bounds.append((knot, knot))
        else:
            #bounds.append((knot_vector[0]+1e-5, knot_vector[idx+1]-1e-5))
            bounds.append((knot_vector[0], knot_vector[-1]))

    print(bounds)
    
    pulse = pls.Pulse(samples[2], 0,1)#fixed_time_array[0], fixed_time_array[-1])
    #pulse.save_plot('test_pulse')
    #waveform = pulse.get_waveform()
    waveform = lead_eigvector

    spline_size = 200
    domain = np.linspace(0,1, spline_size)
    interpolated_function = interpolate.interp1d(x=np.linspace(0, 1, len(lead_eigvector)), y=waveform, kind='cubic')
    interpolated_waveform = interpolated_function(domain)
    interpolated_waveform_at_knots = interpolated_function(knot_vector)
    #test_coeff = [0.3,0.38,0.12,0.07]

    res = optimize.minimize(bspline_err, knot_vector, args=interpolated_waveform, method='Nelder-Mead', bounds=bounds, tol=1e-8, options={'maxfev':1e5})
    knot_vector = res['x']
    print(res)
    #print(f'optimized knot vector: {knot_vector}')
    #spline, error = bspline(interpolated_waveform, knot_vector)

    knot_vector.sort()
    fit = bsp.BSpline(interpolated_waveform, knot_vector)
    spline = fit.spline_function()
    error = fit.get_error()
    fit_coeff = fit.get_fit_coeff()
    fit.plot_basis('test_LGAD', 1000)
    fit.plot_fit('test_fit')
    
    print(f'fit coefficients: {fit_coeff}')
    #0., 0.05, 0.085, 0.1, 0.14, 0.2, 0.4, 1.

    '''
    #error, spline = spline_function(interpolated_waveform, knot_vector)
    '''
    print(f'Error: {error}')
    
    plt.xlabel('Domain')
    plt.ylabel('Amplitude')
    plt.plot(knot_vector, interpolated_waveform_at_knots, 'o', color='black', label='Interpolation at Knots')
    plt.plot(domain, interpolated_waveform, label='Interpolated data')
    plt.plot(domain, spline, label='B-Spline Fit')
    plt.grid()
    plt.text(0.75,0.25,f'Error: {round(error,5)}')
    legend=plt.legend(loc="upper right", fancybox=True)
    frame = legend.get_frame()
    frame.set_edgecolor('black')
    frame.set_alpha(1)
    plt.show()
    '''
    
