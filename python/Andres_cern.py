import h5py
import numpy as np
import matplotlib.pyplot as plt

import optparse
import sys
import os
from InvKLJfilter import KLJ_filter 

from iminuit import Minuit
from probfit import UnbinnedLH, gaussian

from scipy import signal
from scipy import interpolate
from scipy.ndimage import filters

def cfd(samples, threshold=0.02, base_line=200, fraction=0.3, hysteresis=0.001, smooth=True, kernel_width=10, kernel_sigma=5):
    if samples.max()-samples.min() < threshold:
        return None
    if smooth:
        kernel = signal.gaussian(kernel_width, kernel_sigma)
        samples = filters.convolve(samples, kernel)
    samples -= np.mean(samples[:base_line])
    if abs(np.max(samples)) < abs(np.min(samples)):
        maximum = abs(samples.min())
        samples /= -maximum
    else:
        maximum = samples.max()
        samples /= maximum
    
    hysteresis = hysteresis/maximum
    above = False
    lockedForHysteresis = False
    numberOfPeaks = 0
    for i,v in enumerate(samples):
        if not above and not lockedForHysteresis and v>fraction:
            firstCrossingIndex = i
            above = True
            lockedForHysteresis = True
        if above and lockedForHysteresis and v>fraction+hysteresis:
            lockedForHysteresis = False
        if above and not lockedForHysteresis and v<fraction:
            numberOfPeaks += 1
            above = False

    if numberOfPeaks==1:
        t1 = firstCrossingIndex
        t2 = firstCrossingIndex-1
        v1 = samples[firstCrossingIndex]
        v2 = samples[firstCrossingIndex-1]
        return (t1 + (t2 - t1)*(fraction-v1)/(v2 - v1))
    else:
        return None


# In[3]:


if __name__ == '__main__':
    import optparse
    import sys
    import pulse

    usage = "usage: %prog [-n]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default="output.h5")
    parser.add_option("--ref", type=int, dest="ref_ch",
                      help="channel that is the reference sensor", default=1)
    parser.add_option("--ch", type=int, dest="channels",
                      help="amount of channels saved in input file", default=1)
    parser.add_option("--cfd", type=float, dest="cfd",
                      help="percentage for CFD time estimation", default=0.3)
    parser.add_option("--thresh", type=float, dest="threshold",
                      help="sample amplitude threshold value", default=0.025)
    parser.add_option("--baseline", type=float, dest="baseline",
                      help="baseline time in nanoseconds", default=20)
    parser.add_option("--plot", action="store_true", dest="isPlot",
                      help="if called, save the time difference plot", default=False)
    parser.add_option("--filter", action="store_true", dest="useFilter",
                      help="if called, use inverse KLJ filter on the input dataset", default=False)
    (options, args) = parser.parse_args()

    if(options.channels < 1 or options.channels > 8):
        sys.exit(f'{options.channels} is not a valid amount of channels!')

    if(options.cfd < 0.01 or options.cfd > 1.0):
        sys.exit("CFD percentage must be a number between 0 and 1!")

    channels = [1,2,3,4]
    baseline_fraction = 0.3
    ref_ch = 1
    cfd_fraction = 0.4
    infile = options.infile
    amplitudeLimits = {1:(0.1,0.65), 2:(0.1,0.65), 3:(0.1,0.65), 4:(0.1,0.65)}
    
    smooth=False
    kernel_sigma=1
    kernel_width=5*kernel_sigma

    f = h5py.File(infile,'r')

    num_triggers = len(f['ch1_samples'])
    #print(num_triggers)
    # num_triggers = 10000
    num_samples = len(f['ch1_samples'][0])
    time_begin = f['ch1_trig_offset'][0]*1e9
    time_end = f['ch1_horiz_scale'][0]*(num_samples-1)*1e9+time_begin
    time = np.linspace(time_begin, time_end, num_samples)
    sampling_freq = num_samples/(time_end-time_begin)
    
    amplitudes = {}
    amplitudesSelected = {}
    cfd_times = {}
    noise = {}
    riseTime = {}
    
    for channel in channels:
        temp_cfd = []
        temp_amp = []
        temp_amp_selected = []
        temp_noise = []
        temp_riseTime = []
        samples = f[f'ch{channel}_samples']
        horiz_scale = f[f'ch{channel}_horiz_scale']
        trig_offset = f[f'ch{channel}_trig_offset']

        if options.useFilter:
            inv_KLJ_filter = KLJ_filter(samples, 15)
        
        for trig in range(num_triggers):
            waveform = samples[trig]
            if options.useFilter and channel == ref_ch:
                waveform = np.matmul(inv_KLJ_filter, waveform)
            baseline_n = int(len(waveform)*baseline_fraction)
            waveform = np.mean(waveform[:baseline_n]) - waveform
            temp_noise.append(np.std(waveform[:baseline_n]))
        
            amp = np.max(waveform)
            temp_amp.append(amp)
        
            if amp>amplitudeLimits[channel][0] and amp<amplitudeLimits[channel][1]:
                cfd_val = cfd(waveform, threshold=amplitudeLimits[channel][0], base_line=baseline_n, fraction=cfd_fraction,
                              smooth=smooth, kernel_width=kernel_width, kernel_sigma=kernel_sigma)
                if cfd_val is not None:
                    temp_cfd.append((cfd_val*horiz_scale[trig]+trig_offset[trig])*1e9)
                    temp_amp_selected.append(amp)
                    #                 print(cfd_val)
                else:
                    temp_cfd.append(-900)
                    # Rise time:
                    cfd_rt20 = cfd(waveform, threshold=amplitudeLimits[channel][0], base_line=baseline_n, fraction=0.2, smooth=False)
                    cfd_rt80 = cfd(waveform, threshold=amplitudeLimits[channel][0], base_line=baseline_n, fraction=0.8, smooth=False)
                    if cfd_rt20 is not None and cfd_rt80 is not None:
                        temp_riseTime.append( (cfd_rt80-cfd_rt20)*horiz_scale[trig]*1e9 )
                    else:
                        temp_riseTime.append(-900)
            
            else:
                temp_cfd.append(-900)
                temp_riseTime.append(-900)
        
        cfd_times[channel] = np.array(temp_cfd)
        amplitudes[channel] = np.array(temp_amp)
        amplitudesSelected[channel] = np.array(temp_amp_selected)
        noise[channel] = np.array(temp_noise)
        riseTime[channel] = np.array(temp_riseTime)
    
    fig, axes = plt.subplots(2, 2, figsize=(20,20))
    axes = axes.flatten()
    for ch,ax in zip(cfd_times.keys(), axes):
        filteredCfdTimes = cfd_times[ch][cfd_times[ch]>-900]
        ax.hist(filteredCfdTimes, bins=30)
        ax.set_title(f'CFD time of ch {ch}')
        ax.set_xlabel('t (ns)')
    plt.savefig(f'CFD_time_CDF{int(cfd_fraction*100)}_{os.path.splitext(infile)[0]}.pdf')

    fig, axes = plt.subplots(2, 2, figsize=(20,20))
    axes = axes.flatten()
    for ch,ax in zip(channels, axes):
        ax.hist(amplitudes[ch], bins=30, range=(0,1), label='All')
        bins = ax.hist(amplitudesSelected[ch], bins=30, range=(0,1), label='Selected', alpha=0.6, hatch='/')    
        ax.set_title(f'Amplitude of ch {ch}')
        ax.set_xlabel('A (V)')
        ax.legend()
        ax.set_ylim(0, 1.5*np.max(bins[0]))
    plt.savefig(f'amplitude_{os.path.splitext(infile)[0]}.pdf')

    fig, axes = plt.subplots(2, 2, figsize=(20,20))
    axes = axes.flatten()
    for ch,ax in zip(channels, axes):
        ax.hist(1e3*noise[ch], bins=30)  
        ax.set_title(f'Noise of ch {ch}: {1e3*np.mean(noise[ch]):.3f} mV')
        ax.set_xlabel('noise (mV)')
    plt.savefig(f'noise_{os.path.splitext(infile)[0]}.pdf')


    fig, axes = plt.subplots(2, 2, figsize=(20,20))
    axes = axes.flatten()
    for ch,ax in zip(channels, axes):
        rtFiltered = riseTime[ch][riseTime[ch]>-500]
        ax.hist(rtFiltered, bins=30)  
        ax.set_title(f'Rise Time 20-80 of ch {ch}: {np.mean(rtFiltered):.3f} ns')
        ax.set_xlabel('rt20-80 (ns)')
    plt.savefig(f'risetime_{os.path.splitext(infile)[0]}.pdf')

    time_diff = {}

    for channel in channels:
        temp = []
        if(channel is not ref_ch):
            for t,r,rt in zip(cfd_times[channel], cfd_times[ref_ch], riseTime[channel]):
                if abs(t)<4 and t>0 and -10<r<10:
                    temp.append(t-r)
            time_diff[channel] = temp

    fig, axes = plt.subplots(2, 2, figsize=(20,20))
    axes = axes.flatten()
    for ch,d in time_diff.items():
        tmpArray = np.array(time_diff[ch])
        axes[ch-1].hist(tmpArray, bins=20)
        axes[ch-1].set_title(f'Mean: {np.mean(tmpArray):.3f} ns   std: {np.std(tmpArray):.3f} ns')
    plt.savefig(f'time_diff_noFit_CDF{int(cfd_fraction*100)}_{os.path.splitext(infile)[0]}.pdf')

    for ch,d in time_diff.items():
        print(f'Channel {ch}:')

        unbinned_likelihood = UnbinnedLH(gaussian, np.array(d))
        minuit = Minuit(unbinned_likelihood, mean=0.1, sigma=1.1)
        minuit.migrad()
        parameters = minuit.values
        errors = minuit.errors
        mean = parameters['mean']
        sigma = parameters['sigma']
        mean_err = errors['mean']
        sigma_err = errors['sigma']
        print(f'Mean: {mean:.3} +/- {mean_err:.3} ns')
        print(f'Standard Dev: {sigma:.3} +/- {sigma_err:.3} ns')
        
        n_bin = 33
        fig, ax = plt.subplots(figsize=(20,10))
        ax.set_title(f'Time Difference Distribution  (CFD {int(cfd_fraction*100)}%)')
        ax.set_xlabel('Time Difference [ns]')
        ax.set_ylabel('Triggers')
        unbinned_likelihood.draw(minuit, bins=n_bin, show_errbars=None, no_plot=False)
        plt.hist(d, n_bin)
        plt.savefig(f'ch{ch}_time_diff_CDF{int(cfd_fraction*100)}_{os.path.splitext(infile)[0]}.pdf')
