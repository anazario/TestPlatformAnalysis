import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import schur

import pulse as pls
from InvKLJfilter import *

if __name__ == '__main__':

    channel = 1
    pulse_index = 0
    smooth=False
    kernel_sigma=1
    kernel_width=5*kernel_sigma

    input_file = "source_coincidence_test9.h5"
    
    f = h5py.File(input_file,'r')

    num_trig = len(f['ch1_samples'])
    num_samples  = len(f['ch1_samples'][0])
    #time_data = f[f'ch{channel}_time']
    #time = np.linspace(time_begin, time_end, num_samples)
    baseline = int(num_samples/3)
    
    samples = f[f'ch{channel}_samples']
    #horiz_scale = f[f'ch{channel}_horiz_scale']
    #trig_offset = f[f'ch{channel}_trig_offset']

    interpolation_size=51
    for idx in range(1):

        time_begin = f[f'ch{channel}_trig_offset'][idx]*1e9
        time_end = f[f'ch{channel}_horiz_scale'][idx]*(num_samples-1)*1e9+time_begin
        time = np.linspace(time_begin, time_end, num_samples)
        waveform = samples[:][idx]
        '''
        t1 = time[63]
        t2 = time[64]
        v1 = -waveform[63]
        v2 = -waveform[64]
        fraction = np.max(-waveform)*0.4
        print(fraction/0.4)
        print(t1+((t2 - t1)*(fraction-v1)/(v2 - v1)))
        '''
        #print(-waveform[63], -waveform[64])
        #print(time[63], time[64])
        
        pulse = pls.Pulse(waveform, time_begin, time_end)
        pulse.subtract_pedestal()
        pulse.flip_sign()
        pulse.save_plot(f'pulse_plots/pulseOG_{idx}', title=f'Channel {channel}: Pulse {idx}', isDot=True)
        
        max_time, max_amplitude = pulse.get_interpolated_max(50)
        cfd_time = pulse.get_cfd_time(fraction=0.4, baseline=baseline)
        transformation_time, transformed_waveform = pulse.get_interpolated_waveform(size=interpolation_size, time_window=2, fixed_time=cfd_time)

        time_point = [cfd_time]
        amp_point = [max_amplitude*0.4]
        fig, ax = plt.subplots()
        plt.plot(transformation_time, transformed_waveform, 'o')
        plt.plot(time_point, amp_point, 'o', color='black')
        plt.savefig('test_cfd_point.pdf')
        #print(max_amplitude, max_time)
        #print(np.max(transformed_waveform), transformation_time[np.argmax(transformed_waveform)])
        
        print(cfd_time, transformed_waveform[10])
        #print(transformed_waveform)
        #print(transformation_time)
        
        #transformed_waveform/=np.max(transformed_waveform)
        pulse = pls.Pulse(transformed_waveform, transformation_time[0], transformation_time[-1])
        pulse.save_plot(f'pulse_plots/pulse_{idx}', title=f'Channel {channel}: Pulse {idx}', isDot=True)

    #transformed dataset
    input_file = "transformed_ch1_coincidence_test9_CFD40_51samples.h5"
    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    fixed_time_array = f[f'ch{channel}_fixed_time_array']
    density_matrix = f[f'ch{channel}_density_matrix']
    eig_value_array, eig_vector_array = schur(density_matrix)
    '''
    for idx, waveform in enumerate(samples[:]):
        if any(waveform):
            print(np.correlate(eig_vector_array[0],waveform), fixed_time_array[idx])
    '''
    '''
    pulse = pls.Pulse(samples[:][pulse_index],time_begin,time_end)

    max_time, max_amplitude = pulse.get_interpolated_max(50)

    #print(max_amplitude)
    
    interpolation_size = 100
    xaxis = np.linspace(max_time-10, max_time+10, interpolation_size)
    interpolation_time, interpolated_pulse = pulse.get_interpolated_waveform(interpolation_size, max_time-10, max_time+10)

    print(pulse.get_interpolated_max(50)) 
    
    fig, ax = plt.subplots()
    ax.set_ylabel('Amplitude [V]')
    ax.set_xlabel('Time [ns]')
    ax.set_title(f'Channel {channel}: Pulse {pulse_index}')
    plt.plot(xaxis, interpolated_pulse)
    #plt.plot(peak_time, fit_res, color = 'red')
    
    plt.savefig('test.pdf')

    #pulse_list = pls.PulseList([pulse])
    #pulse_list.test_function()
    '''
    '''    
    eig_values, eig_vectors = GetBasisInvKLJ(samples)
    #print(eig_values)
    
    cfd_list = []
    c_list = []
    
    for trig, waveform in enumerate(samples[:1000]):
    
        waveform = np.mean(waveform[:int(num_samples/3)]) - np.array(waveform)
        cfd_val = cfd(waveform, 0.1, base_line=baseline_n, fraction=0.4, smooth=smooth, kernel_width=kernel_width, kernel_sigma=kernel_sigma)

        if cfd_val is not None:
            converted_cfd_time = (cfd_val*horiz_scale[trig]+trig_offset[trig])*1e9
            coefficient = np.inner(eig_vectors[0], waveform)
            cfd_list.append(converted_cfd_time)
            c_list.append(coefficient)

            if converted_cfd_time < 0 and coefficient < 0:
                time_begin = f[f'ch{channel}_trig_offset'][trig]*1e9
                time_end = f[f'ch{channel}_horiz_scale'][trig]*(num_samples-1)*1e9+time_begin
                pulse = pls.Pulse(-1*samples[:][trig], time_begin, time_end)
                pulse.save_plot(f'Channel_{channel}_Trigger_{trig}')
                print(f'index: {trig}, Max Amplitude: {-np.min(samples[:][trig])}')
                
    fig, ax = plt.subplots()
    ax.set_ylabel('Leading Eigenvector Coefficient Value')
    ax.set_xlabel('CFD Time w.r.t Trigger [ns]')
    ax.set_title(f'Channel {channel} (CFD 40%)')
    plt.scatter(cfd_list, c_list)
    #plt.hist(cfd_list, bins=7, range=(-0.2,0.5))
    plt.savefig('test.pdf')
    '''
    '''
    #reco_waveform = np.zeros(num_samples)

    for idx, vector in enumerate(eig_vectors):
        c = np.inner(vector, waveform)
        temp_vec = c*vector
        reco_waveform += temp_vec    

    fig, ax = plt.subplots()
    ax.set_ylabel('Amplitude [V]')
    ax.set_xlabel('Time [ns]')
    ax.set_title(f'Channel {channel}: Pulse {pulse_index}')
    plt.plot(time, reco_waveform)

    plt.savefig('test.pdf')
    '''
