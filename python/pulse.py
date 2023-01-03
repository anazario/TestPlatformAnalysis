import numpy as np
import math
import matplotlib.pyplot as plt
import sys

from pulseTools import *

class Pulse:

    def __init__(self, waveform, time_start, time_end, flipSign=False, isNormal=False):
        self.waveform = waveform
        self.num_samples = len(waveform)
        self.time_start = time_start
        self.time_end = time_end
        self.sampling_frequency = self.num_samples/(time_end-time_start)
        self.noise_rms = np.std(waveform[:int(0.3*self.num_samples)])
        self.sign_inverted = False

        if flipSign:
            self.flip_sign()
        if isNormal:
            self.normalize()

    def get_waveform(self):
        return self.waveform
            
    def normalize(self):
        self.waveform /= np.max(self.waveform)

    def subtract_pedestal(self, baseline=3):
        waveform = self.waveform
        noise_rms = np.mean(waveform[:int(len(waveform)/baseline)])
        self.waveform -= noise_rms
        
    def flip_sign(self):
        self.waveform *= -1

    def get_num_samples(self):
        return len(self.waveform)

    def get_time_array(self):
        return np.linspace(self.time_start, self.time_end, self.num_samples)

    def get_time_at_max(self):
        waveform = self.waveform
        time_arr = self.get_time_array()
        return(time_arr[np.argmax(waveform)])

    def get_max_amplitude(self):
        waveform = self.waveform
        return np.max(waveform)

    def save_plot(self, name='pulse', title='', isDot=False):
        fig, ax = plt.subplots()
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Amplitude [mV]')
        ax.set_title(title)
        if isDot:
            plt.plot(self.get_time_array(), self.waveform, 'o')
        else:
            plt.plot(self.get_time_array(), self.waveform)
        plt.savefig(name+".pdf")
        plt.clf()
        print(f'Saved plot: {name}.pdf')
        
    def get_cfd_time(self, fraction=0.3, threshold=None, baseline=None, smooth=False):
        if threshold is None:
            threshold=(self.noise_rms)*6
        if baseline is None:
            baseline=int(0.3*self.num_samples)
        time = self.get_time_array()
        return cfd(self.waveform, time, threshold, baseline, fraction, smooth=smooth)


    def get_interpolated_waveform(self, size, time_window, fixed_time=None, t_below=2, t_above=8):

        if fixed_time == None:
            fixed_time = self.get_time_at_max()
        
        waveform = self.waveform
        sample_rate = self.sampling_frequency

        step_size = time_window/size
        interpolation_time = np.linspace(fixed_time-step_size*int(t_below*(size-1)/10), fixed_time+step_size*int(t_above*(size-1)/10), size)

        interpolated_pulse = []
        time_diff = self.time_start

        for idx in range(size):
            amplitude = get_interpolated_point(waveform, interpolation_time[idx]-time_diff, sample_rate)
            interpolated_pulse.append(amplitude)
            
        return interpolation_time, interpolated_pulse

    def get_interpolated_max(self, size):

        approx_max_time = self.get_time_at_max()

        interpolation_time, interpolated_waveform = self.get_interpolated_waveform(size, 2, approx_max_time, t_below=5, t_above=5)

        max_idx = np.argmax(interpolated_waveform)
        max_amplitude = interpolated_waveform[max_idx]
        max_time = interpolation_time[max_idx]
        
        return max_time, max_amplitude        
    

class PulseList(list):

    def __init__(self, iterable=None):
        super(PulseList, self).__init__()
        if iterable:
            for item in iterable:
                self.append(item)

    def append(self, item):
        if isinstance(item, Pulse):
            super(PulseList, self).append(item)
        else:
            raise ValueError('Only Pulse objects are allowed!')

    def insert(self, index, item):
        if isinstance(item, Pulse):
            super(PulseList, self).insert(index, item)
        else:
            raise ValueError('Only Pulse objects are allowed!')

    def __add__(self, item):
        if isinstance(item, Pulse):
            super(PulseList, self).__add__(item)
        else:
            raise ValueError('Only Pulse objects are allowed!')

    def __iadd__(self, item):
        if isinstance(item, Pulse):
            super(PulseList, self).__iadd__(item)
        else:
            raise ValueError('Only Pulse objects are allowed!')
                
'''    
def get_pulse(h5_file, channel, pulse_index):

    import h5py
    
    file_h5 = h5py.File(h5_file, 'r')

    samples      = file_h5[f'ch{channel}_samples']
    horiz_offset = file_h5[f'ch{channel}_horiz_offset']
    horiz_scale  = file_h5[f'ch{channel}_horiz_scale']
    trig_time    = file_h5[f'ch{channel}_trig_time']
    trig_offset  = file_h5[f'ch{channel}_trig_offset']

    n_samples   = (len(samples))
    sample_size = (len(samples[0]))

    if (pulse_index > n_samples):
        sys.exit("Specified index is out of bounds!")

    time_begin = trig_offset[pulse_index]*1e9
    time_end   = (horiz_scale[pulse_index]*(sample_size-1)-trig_offset[pulse_index])*1e9
    time       = np.linspace(time_begin, time_end, sample_size)

    pulse = samples[pulse_index]

    return pulse, time
'''
'''
import h5py

f = h5py.File('output.h5', 'r')
samples = f['ch1_samples']
time_begin = f['ch1_trig_offset'][0]*1e9
num_samples = len(samples[0])
time_end = f['ch1_horiz_scale'][0]*(num_samples-1)*1e9+time_begin


pulse_ch1 = Pulse(samples[0]*1e3, time_begin, time_end)
print(pulse_ch1.get_cfd_index(0.4))
#pulse_ch1.save_plot("test_plot", True, True)
'''
