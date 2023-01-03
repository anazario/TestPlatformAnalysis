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

class H5Interface:

    def __init__(self, h5_file, channel):
        self.h5_file = h5_file
        self.channel = channel
        self.samples = h5_file[f'ch{channel}_samples']
        #import waveform conversion factors if 
        if self.samples.dtype == 'int8':
            self.isConverted = False
            self.vert_offset = infile_h5[f'ch{channel}_vert_offset']
            self.vert_scale = infile_h5[f'ch{channel}_vert_scale']
            
        self.trig_offset = h5_file[f'ch{channel}_trig_offset']
        self.horiz_offset = h5_file[f'ch{channel}_horiz_offset']
        self.horiz_scale = h5_file[f'ch{channel}_horiz_scale']

    def change_channel(self, channel):
        self.channel = channel
        self.samples = self.h5_file[f'ch{channel}_samples']
        self.trig_offset = self.h5_file[f'ch{channel}_trig_offset']
        self.vert_offset = self.h5_file[f'ch{channel}_vert_offset']
        self.vert_scale = self.h5_file[f'ch{channel}_vert_scale']
        self.horiz_offset = self.h5_file[f'ch{channel}_horiz_offset']
        self.horiz_scale = self.h5_file[f'ch{channel}_horiz_scale']

    def get_current_channel(self):
        return self.channel
        
    def get_samples(self):
        return self.samples

    def get_trig_offset(self):
        return self.trig_offset

    def	get_vert_offset(self):
        return self.vert_offset

    def get_vert_scale(self):
        return self.vert_scale

    def	get_horiz_offset(self):
        return self.horiz_offset

    def	get_horiz_scale(self):
        return self.horiz_scale
        
    def get_waveform(self, trig):
        waveform = -self.vert_offset[trig]+self.samples[trig]*self.vert_scale[trig]
        return waveform

    def get_time_interval(self, trig):
        num_samples = len(self.samples[0])
        time_begin = self.trig_offset[trig]*1e9
        time_end = self.horiz_scale[trig]*(num_samples-1)*1e9+time_begin
        return time_begin, time_end

     def get_pulse(self, trig):
        waveform = self.get_waveform(trig)
        time_begin, time_end = self.get_time_interval(trig)
        pulse = pls.Pulse(waveform, time_begin, time_end)
        return pulse
    
if __name__ == '__main__':
    import optparse
    import sys

    usage = "usage: %prog [-n]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default=None)
    parser.add_option("--trig", type=int, dest="trigger_ch",
                      help="trigger channel (between 1 and 8)", default=1)
    parser.add_option("--events", type=int, dest="events",
                      help="amount of events to loop over (set to maximum amount if it exceeds it)", default=None)
    parser.add_option("--cfd", type=float, dest="cfd",
                      help="fraction of maximum for CFD time estimation", default=None)
    parser.add_option("--thresh", type=float, dest="threshold",
                      help="sample amplitude threshold value", default=0.7)
    parser.add_option("--baseline", type=float, dest="baseline",
                      help="sample fraction before waveform corresponding to noise", default=3)
    (options, args) = parser.parse_args()

    trigger_ch = options.trigger_ch
    events = options.events
    cfd_fraction = options.cfd
    
    if cfd_fraction is not None and (cfd_fraction < 1e-5 or cfd_fraction > 1.0):
        sys.exit("CFD fraction must be a number between 0 and 1!")
    
    if options.infile is None:
        sys.exit("Please specify an input file (-i [input_file_name])!")
    if trigger_ch < 1 or trigger_ch > 8:
        sys.exit(f'{options.channel} is not a valid channel!')
    
    pulse_index = 0
    smooth=False
    kernel_sigma=1
    kernel_width=5*kernel_sigma
    
    input_file = options.infile

    f = h5py.File(input_file,'r')

    h5 = H5Interface(f, trigger_ch)
    
    samples = h5.get_samples()#f[f'ch{trigger_ch}_samples']
    time_begin_arr = h5.get_trig_offset()
    time_end_arr = h5.get_horiz_scale()
    num_trig = len(samples)
    num_samples  = len(samples[0])
    vert_offset = h5.get_vert_offset()
    vert_scale = h5.get_vert_scale()
    baseline = int(num_samples/options.baseline)
    
    #set max trigger to loop through
    max_trig = options.events
    if max_trig is None or max_trig > num_trig:
        max_trig = num_trig

    ref_channels = [i for i in range(1,9) if i is not trigger_ch]
    coincidences = np.zeros(len(ref_channels)+1)

    for trig, waveform in enumerate(samples[:max_trig]):
        if trig%1000 == 0:
            sys.stdout.write(f"\rProcessing waveform {trig}")
            sys.stdout.flush()

        waveform = -vert_offset[trig]+waveform*vert_scale[trig]
        time_begin = time_begin_arr[trig]*1e9
        time_end = time_end_arr[trig]*(num_samples-1)*1e9+time_begin
        
        #Create pulse object from each waveform in dataset
        pulse = pls.Pulse(waveform, time_begin, time_end)
        pulse.subtract_pedestal()
        pulse.flip_sign()

        if pulse.get_max_amplitude() < options.threshold:
            cfd_time = pulse.get_cfd_time(fraction=cfd_fraction, baseline=baseline)#, smooth=True)                                                                            
            if cfd_time is not None:# and abs(cfd_time) < 5:
                for channel in ref_channels:
                    h5.change_channel(channel)
                    ref_waveform = h5.get_waveform(trig)
                    time_begin, time_end = h5.get_time_interval(trig)
                    pulse = pls.Pulse(ref_waveform, time_begin, time_end)
                    pulse.subtract_pedestal()
                    pulse.flip_sign()
                    max_amplitude = pulse.get_max_amplitude()
                    
                    ref_cfd_time = pulse.get_cfd_time(fraction=cfd_fraction, baseline=baseline)
                    if ref_cfd_time is not None and max_amplitude > 0.1:
                        #pulse.save_plot(name=f'plots/test_pulse{trig}_ch{channel}', title='', isDot=True)
                        coincidences[channel-1] += 1

    print(f'Coincidences in channels {ref_channels} were {coincidences}')
