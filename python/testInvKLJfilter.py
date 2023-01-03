import h5py
import math
import numpy as np
import matplotlib.pyplot as plt

from InvKLJfilter import *

def PlotLeadingEigVectors(input_file, max_idx = 5):
    colors = ['black','red', 'blue', 'green', 'orange']
    f = h5py.File(input_file,'r')
    channels = [1,2,3,4]

    xaxis = np.linspace(0,126,127)
    for channel in channels:
        fig, ax = plt.subplots()
        ax.set_xlabel('Nth Sample')
        ax.set_ylabel('Eigenvector Sample Value')
        ax.set_title(f'Channel {channel}')
    
        samples = f[f'ch{channel}_samples']
        num_samples = len(f[f'ch{channel}_samples'][0])
        
        eig_values, eig_vectors = GetBasisInvKLJ(samples)
        for idx in range(max_idx):
            temp_eig_vector = eig_vectors[:,idx]
            plt.plot(xaxis, temp_eig_vector, color=colors[idx], label=f'λ = {round(eig_values[idx],4)}')    

        plt.legend()    
        plt.savefig(f"channel{channel}_{max_idx}_leading_eigenvectors.pdf") 
        plt.clf()

def PlotWaveform(input_file, channel, pulse_index, useFilter = False):

    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    num_samples = len(f[f'ch{channel}_samples'][pulse_index])

    time_begin = f[f'ch{channel}_trig_offset'][pulse_index]*1e9
    time_end = f[f'ch{channel}_horiz_scale'][pulse_index]*(num_samples-1)*1e9+time_begin
    time = np.linspace(time_begin, time_end, num_samples)

    print(time_end)
    print(time_begin)
    print((time_end-time_begin)/num_samples)
    
    waveform = samples[:][pulse_index]
    waveform = np.mean(waveform[:int(num_samples/3)]) - waveform

    if(useFilter):
        filter_invKLJ = KLJ_filter(samples, 5)
        waveform = np.matmul(filter_invKLJ, waveform)
    
    fig, ax = plt.subplots()
    ax.set_ylabel('Amplitude [V]')
    ax.set_xlabel('Time [ns]')
    ax.set_title(f'Channel {channel}: Pulse {pulse_index}')
    plt.plot(time, waveform)

    label = 'unfiltered'
    if useFilter:
        label = 'filtered'

    plt.savefig(f'channel{channel}_pulse{pulse_index}_{label}.pdf')

if __name__ == '__main__':

    input_file = "source_coincidence_test9.h5"

    PlotWaveform(input_file, 1, 0, False)
    #PlotLeadingEigVectors(input_file, 1)
    '''
    channel = 1

    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    eig_values, eig_vectors = GetBasisInvKLJ(samples)
    first_derivative = np.gradient(eig_vectors[:,0])
    first_derivative /= math.sqrt(np.inner(first_derivative, first_derivative))#np.max(first_derivative)
    eig_vector2 = eig_vectors[:,1]

    xaxis = np.linspace(0,126,127)
    
    fig, ax = plt.subplots()
    ax.set_xlabel('Nth Sample')
    ax.set_ylabel('Eigenvector Sample Value')
    ax.set_title(f'Channel {channel}')
    plt.plot(xaxis, first_derivative, color = 'blue', label = 'derivative of 1st eigenvector')
    plt.plot(xaxis, eig_vector2, color = 'red', label = 'second eigenvector')
    plt.legend()
    plt.savefig(f'channel{channel}_2ndEigenvector.pdf')
    '''
    '''
    max_idx = 3
    colors = ['black','red', 'blue', 'green', 'orange']
    f = h5py.File(input_file,'r')

    num_trig = len(f['ch1_samples'])
    num_samples  = len(f['ch1_samples'][0])
    
    channels = [2,3,4]
    for channel in channels:
        samples = []
        for trig in range(num_trig):
            waveform = f[f'ch{channel}_samples'][trig]
            max_amplitude = np.max(waveform)
            if max_amplitude > 0.1 and max_amplitude < 0.65:
                samples.append(waveform)
        eig_values, eig_vectors = GetBasisInvKLJ(samples)

        xaxis = np.linspace(0,126,127)
        fig, ax = plt.subplots()
        ax.set_xlabel('Nth Sample')
        ax.set_ylabel('Eigenvector Sample Value')
        ax.set_title(f'Channel {channel}')
        for idx in range(max_idx):
            temp_eig_vector = eig_vectors[:,idx]
            plt.plot(xaxis, temp_eig_vector, color=colors[idx], label=f'λ = {round(eig_values[idx],4)}')
            
        plt.legend()
        plt.savefig(f"channel{channel}_{max_idx}_leading_eigenvectors_ampCut.pdf")
        plt.clf()
    '''
    '''
    noise = {}
    signal = {}
    
    for trig in range(num_trig):
        max_ch = -1
        temp_amplitude = -1
        for channel in channels:
            waveform = f[f'ch{channel}_samples'][trig]
            max_amplitude = np.max(waveform)
            if amplitude > 0.1 and amplitude < 0.65:
                if max_ch == -1 or (max_ch > 0 and max_amplitude > temp_amplitude):
                    max_ch = channel
                    temp_amplitude = max_amplitude
        signal[max_ch]
    '''        
    

