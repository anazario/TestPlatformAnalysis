import sys
import h5py
import math
import numpy as np
from scipy.linalg import *
import matplotlib.pyplot as plt

import pulse as pls
from InvKLJfilter import *
from pulseTools import *
#from Andres_cern import cfd

if __name__ == '__main__':

    channel = 1
    pulse_index = 0
    smooth=False
    kernel_sigma=1
    kernel_width=5*kernel_sigma
    
    input_file = "transformed_ch1_coincidence_test9_CFD40_51samples.h5"

    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    time_data = f[f'ch{channel}_time']
    num_trig = len(f['ch1_samples'])
    num_samples  = len(f['ch1_samples'][0])
    baseline_n = int(num_samples*0.3)

    fixed_time_array = f[f'ch{channel}_fixed_time_array']
    density_matrix = f[f'ch{channel}_density_matrix']
    eig_value_array, eig_vector_array = schur(density_matrix)

    #modify eigendecomposition results for ease of use
    eig_value_array = np.array([val[idx] for idx, val in enumerate(eig_value_array)])
    lead_eigvector = eig_vector_array[:,0]    
    if np.max(lead_eigvector) < -np.min(lead_eigvector):
        eig_vector_array = -1*eig_vector_array.T
    else:
        eig_vector_array = eig_vector_array.T
    
    pulse = pls.Pulse(eig_vector_array[1], 0, num_samples)
    pulse.save_plot('next_to_lead_eigvector')
    
    max_trig = 50000

    time_array = []
    coefficient_array = []
    max_amp_array = []
    
    for trig, waveform in enumerate(samples[:max_trig]):
        if np.any(waveform):
            max_amp_array.append(np.max(waveform))
            norm_waveform = waveform/math.sqrt(np.inner(waveform, waveform))
            coefficient = np.inner(eig_vector_array[0], norm_waveform)
            coefficient_array.append(coefficient)
            time_array.append(fixed_time_array[trig])
            if np.max(waveform) < 0.3 and fixed_time_array[trig] > 0.:
                print(coefficient)
                pulse = pls.Pulse(waveform, time_data[trig][0], time_data[trig][-1])
                pulse.save_plot(f'pulse_plots/pulse_{trig}')

    plot_scatter(time_array, coefficient_array, 'coeff_vs_CFDtime_CFD40_51samples', 'CFD Time w.r.t Trigger [ns]',
                 'Leading Eigenvector Coefficient Value', f'Channel {channel} (CFD 40%)')
    plot_scatter(coefficient_array, max_amp_array, 'maxAmp_vs_coeff_CFD40_51samples', 'Leading Eigenvector Coefficient Value',
                 'Maximum Amplitude [V]', f'Channel {channel} (CFD 40%)')
    plot_scatter(time_array, max_amp_array, 'maxAmp_vs_CFDtime_CFD40_51samples', 'CFD Time w.r.t Trigger [ns]',
                 'Maximum Amplitude [V]', f'Channel {channel} (CFD 40%)')
    plot_line(np.linspace(0, num_samples, num_samples), np.gradient(eig_vector_array[0]), 'derivative_of_lead_eigvector_CFD40_51samples',
              'Sample Number', 'Derivative of Lead Eigenvector', f'Channel {channel} (CFD 40%)')

    y_list = [np.gradient(eig_vector_array[0]), eig_vector_array[1]]
    leg_labels = ['Derivative of Lead Eigenvector', 'Second Eigenvector']

    plot_line(np.linspace(0, num_samples, num_samples), y_list, 'derivative_of_lead_eigvector_vs_NLEigvector_CFD40_51samples',
              'Sample Number', 'No Units', f'Channel {channel} (CFD 40%)', leg_labels)

    plot_hist(data=max_amp_array, bins=28, range_tuple=(0.,0.7), name='max_amplitude_dist_CFD40_51samples', xlabel='Maximum Amplitude [V]',
              ylabel='Entries', title=f'Channel {channel} (CFD 40%)')
    
    #for idx, vec in enumerate(eig_vector_array):
    #plot_line(np.linspace(0, num_samples, num_samples), vec, f'pulse_plots/eigvector{idx}', 'Sample Number', 'Eigenvector Value', f'Channel {channel} (CFD 40%)')
    
    #pulse = pls.Pulse(eig_vector_array[:,0], 0, num_samples)
    #pulse.save_plot('lead_eigvector')
    '''
    time = np.linspace(0, num_samples, num_samples)
    eig_vec_idx = 0
    #print(test_waveform[50])
    fig, ax = plt.subplots()
    #ax.set_ylabel('Leading Eigenvector Coefficient Value')
    #ax.set_xlabel('CFD Time w.r.t Trigger [ns]')
    #ax.set_title(f'Channel {channel} (CFD 40%)')
    #plt.scatter(cfd_list, c_list)
    #plt.hist(cfd_list, bins=7, range=(-0.2,0.5))
    #plt.plot(time, eig_vector_array[:,0])
    #plt.plot(time, test_waveform, 'o')
    x = fixed_time_array[:max_trig]
    y = coefficient_array
    plt.scatter(x,y)
    plt.savefig(f'test_coeff_vd_CFDtime.pdf')
    '''
    '''
    waveform = samples[:][0]
    
    time_begin = f[f'ch{channel}_trig_offset'][0]*1e9
    time_end = f[f'ch{channel}_horiz_scale'][0]*(num_samples-1)*1e9+time_begin

    pulse = pls.Pulse(waveform, time_begin, time_end)
    max_time, max_amplitude = pulse.get_interpolated_max(50)
    pedestal = np.mean(waveform[:int(len(waveform)/3)])
    interpolation_time, interpolated_pulse = pulse.get_interpolated_waveform(interpolation_size, -10, 10, max_time)
    mod_waveform = np.array(interpolated_pulse)-pedestal
    norm_waveform = mod_waveform/math.sqrt(np.inner(mod_waveform, mod_waveform))
    
    reco_waveform = np.zeros(interpolation_size)
    for idx, vector in enumerate(eig_vector_array):
        c = np.inner(vector, norm_waveform)
        print(c)
        temp_vec = c*vector
        reco_waveform += temp_vec

    fig, ax = plt.subplots()
    plt.plot(time, reco_waveform)
    plt.savefig('test_reco_waveform.pdf')
    '''
    '''
    for idx, vector in enumerate(eig_vector_array.T):
        normalization = np.inner(vector, vector)
        #print(normalization)
        eig_value_x_vector = np.matmul(density_matrix, vector)
        print(np.inner(vector, eig_value_x_vector))
    '''
