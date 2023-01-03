import numpy as np
import h5py

def cfd(voltage, threshold=0.02, base_line=30, fraction=0.3, hysteresis=0.001, chooseFirst=True):
    samples = voltage.astype(float)
    #print(base_line*sample_rate)
    # If peak2peak < threshold there is no peak for sure...
    if samples.max()-samples.min() < threshold:
        return None
    # Work only with positive and normalized signals 
    samples -= np.mean(samples[:base_line])
    #print(np.mean(samples[:int(base_line*sample_rate)]))
    if abs(samples.max()) < abs(samples.min()):
        maximum = abs(samples.min())
        samples *= -1/maximum
    else:
        maximum = samples.max()
        samples *= 1/maximum
    #print(f'Maximum: {maximum} V')
    #print(f'Time at maximum: {np.mean(time[np.where(samples == samples.max())])}')
    
    #threshold = threshold/maximum
    hysteresis = hysteresis/maximum
    
    above = False
    lockedForHysteresis = False
    numberOfPeaks = 0
    for i,v in enumerate(samples):
        #print("%d: %.3f"%(i,v))
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
        #print((fraction-v1))
        return (t1 + (t2 - t1)*(fraction-v1)/(v2 - v1))
    else:
        return None
    
def get_pulse(h5_file, channel, pulse_index):

    file_h5 = h5py.File(h5_file, 'r')

    samples      = file_h5[f'ch{channel}_samples']
    horiz_offset = file_h5[f'ch{channel}_horiz_offset']
    horiz_scale  = file_h5[f'ch{channel}_horiz_scale']
    vert_offset  = file_h5[f'ch{channel}_vert_offset']
    vert_scale   = file_h5[f'ch{channel}_vert_scale']
    trig_time    = file_h5[f'ch{channel}_trig_time']
    trig_offset  = file_h5[f'ch{channel}_trig_offset']

    n_samples   = (len(samples))
    sample_size = (len(samples[0]))

    if (pulse_index > n_samples):
        sys.exit("Specified index is out of bounds!")

    time_begin = trig_offset[pulse_index]*1e9
    time_end   = (horiz_scale[pulse_index]*(sample_size-1)-trig_offset[pulse_index])*1e9
    time       = np.linspace(time_begin, time_end, sample_size)

    pulse = -vert_offset[pulse_index]+samples[pulse_index,:]*vert_scale[pulse_index]
    
    return pulse, time
'''
pulse, time = get_pulse('x8_y5.hdf5', 1, 11)

print(cfd(time, pulse))

plt.plot(time, pulse*-1)
plt.show()
'''
