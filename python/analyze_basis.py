import h5py
import math
import numpy as np
from scipy.linalg import *
import matplotlib.pyplot as plt
import ROOT

import pulse as pls
from InvKLJfilter import *
from pulseTools import *
#from Andres_cern import cfd

if __name__ == '__main__':
    import optparse
    import sys

    usage = "usage: %prog [-n]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default=None)
    parser.add_option("--ch", type=int, dest="channel",
                      help="channel to process (between 1 and 8)", default=1)
    parser.add_option("--vector", type=int, dest="vector",
                      help="choose eigenvector to calculate coefficient value (starts at 1)", default=1)
    parser.add_option("--trig", type=int, dest="max_trig",
                      help="triggers to loop over (set to maximum trigger if it exceeds it)", default=None)
    parser.add_option("--cfd", type=float, dest="cfd",
                      help="fraction of maximum for CFD time estimation", default=None)
    
    (options, args) = parser.parse_args()

    channel = options.channel
    cfd_fraction = options.cfd
    if cfd_fraction is None:
        tag='maxAmp'
    else:
        tag=f'CFD{int(100*cfd_fraction)}'
    
    if options.infile is None:
        sys.exit("Please specify an input file (-i [input_file_name])!")
    if channel < 1 or channel > 8:
        sys.exit(f'{options.channel} is not a valid channel!')
    
    pulse_index = 0
    smooth=False
    kernel_sigma=1
    kernel_width=5*kernel_sigma
    
    input_file = options.infile

    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    time_data = f[f'ch{channel}_time']
    num_trig = len(f[f'ch{channel}_samples'])
    num_samples  = len(f[f'ch{channel}_samples'][0])
    
    #set max trigger to loop through
    max_trig = options.max_trig
    if max_trig is None or max_trig > num_trig:
        max_trig = num_trig
    
    fixed_time_array = f[f'ch{channel}_fixed_time_array']
    density_matrix = f[f'ch{channel}_density_matrix']
    eig_value_array, eig_vector_array = schur(density_matrix)

    print(eig_value_array[0], eig_value_array[1])
    
    #modify eigendecomposition results for ease of use
    eig_value_array = np.array([val[idx] for idx, val in enumerate(eig_value_array)])
    lead_eigvector = eig_vector_array[:,0]
    if np.max(lead_eigvector) < -np.min(lead_eigvector):
        eig_vector_array = -1*eig_vector_array.T
    else:
        eig_vector_array = eig_vector_array.T

    h2f_leadcoeff_vs_time = ROOT.TH2F('h2f_c_vs_t', 'h2f_c_vs_t', 32, -1.75, 0.25, 60, 0.85, 1.00)
    h2f_maxamp_vs_time = ROOT.TH2F('h2f_a_vs_t', 'h2f_a_vs_t', 90, -1.5, 0.75, 75, 0., 0.75)
    
    time_array = []
    coefficient_array = []
    max_amp_array = []

    eig_vector_idx = options.vector-1
    if eig_vector_idx > num_samples:
        eig_vector_idx = num_samples-1
    
    for trig, waveform in enumerate(samples[:max_trig]):
        if np.any(waveform) and fixed_time_array[trig] > -1.0:
            norm_waveform = waveform/math.sqrt(np.inner(waveform, waveform))
            coefficient = np.inner(eig_vector_array[eig_vector_idx], norm_waveform)
            max_amp_array.append(np.max(waveform))
            coefficient_array.append(coefficient)
            time_array.append(fixed_time_array[trig])

            h2f_leadcoeff_vs_time.Fill(fixed_time_array[trig], coefficient)
            h2f_maxamp_vs_time.Fill(fixed_time_array[trig], np.max(waveform))
            
            if np.max(waveform) < 0.1 and fixed_time_array[trig] < -1.25:
                print(coefficient)
                pulse = pls.Pulse(waveform, time_data[trig][0], time_data[trig][-1])
                pulse.save_plot(f'pulse_plots/pulse_{trig}')

    plot_scatter(time_array, coefficient_array, f'coeff{eig_vector_idx+1}_vs_CFDtime_{tag}_{num_samples}samples', 'CFD Time w.r.t Trigger [ns]',
                 f'Eigenvector {eig_vector_idx+1} Coefficient Value', f'Channel {channel}')
    plot_scatter(coefficient_array, max_amp_array, f'maxAmp_vs_coeff{eig_vector_idx+1}_{tag}_{num_samples}samples', f'Eigenvector {eig_vector_idx+1} Coefficient Value',
                 'Maximum Amplitude [V]', f'Channel {channel}')
    plot_scatter(time_array, max_amp_array, f'maxAmp_vs_time_{tag}_{num_samples}samples', 'CFD Time w.r.t Trigger [ns]',
                 'Maximum Amplitude [V]', f'Channel {channel}', fit=2)
    plot_scatter(max_amp_array, time_array, f'time_vs_maxAmp_{tag}_{num_samples}samples', 'Maximum Amplitude [V]',
                 'CFD Time w.r.t Trigger [ns]', f'Channel {channel}', fit=6)

    
    #plot_line(np.linspace(0, num_samples, num_samples), np.gradient(eig_vector_array[0]), f'derivative_of_lead_eigvector_{tag}_{num_samples}samples',
    #'Sample Number', 'Derivative of Lead Eigenvector', f'Channel {channel}')

    #y_list = [np.gradient(eig_vector_array[0]), eig_vector_array[1]]
    #leg_labels = ['Derivative of Lead Eigenvector', 'Second Eigenvector']

    #plot_line(np.linspace(0, num_samples, num_samples), y_list, f'derivative_of_lead_eigvector_vs_NLEigvector_{tag}_{num_samples}samples',
    #'Sample Number', 'No Units', f'Channel {channel}', leg_labels)

    y_list = eig_vector_array[:,1]
    plot_line(np.linspace(0, num_samples, num_samples), y_list, f'NL_eigvector_{tag}_{num_samples}samples',                                       
    'Sample Number', 'Eigvector Value', f'Channel {channel}')
    
    plot_hist(data=max_amp_array, bins=28, range_tuple=(0.,0.7), name=f'max_amplitude_dist_{tag}_{num_samples}samples', xlabel='Maximum Amplitude [V]',
              ylabel='Entries', title=f'Channel {channel}')
    plot_hist(data=time_array, bins=16, range_tuple=(-1.75, 0.25), name=f'time_dist_{tag}_{num_samples}samples', xlabel='Time [ns]',
              ylabel='Entries', title=f'Channel {channel}')

    plot_2D_hist(h2f_leadcoeff_vs_time, name='coefficient_vs_trigtime', x_label=f'CFD Time ({int(100*cfd_fraction)}%) w.r.t Trigger [ns]',
                 y_label='Leading Coefficient Value', setLogz=True)
    plot_2D_hist(h2f_maxamp_vs_time, name='maxamp_vs_trigtime', x_label=f'CFD Time ({int(100*cfd_fraction)}%) w.r.t Trigger [ns]',
                 y_label='Maximum Amplitude [V]', setLogz=True)
    pulse = pls.Pulse(samples[:][0], time_data[0][0], time_data[0][-1])
    pulse.save_plot(f'pulse_plots/pulse_0', title=f'Channel {channel}: Pulse 0', isDot=True)
