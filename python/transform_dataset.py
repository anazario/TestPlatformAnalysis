import h5py
import math
import numpy as np
import pulse as pls

from scipy import interpolate
from scipy.linalg import *

if __name__ == '__main__':
    import optparse
    import sys

    usage = "usage: %prog [-i <input file>] [-o <output file>] [--output_size <size of interpolation>] [--ch <channel to process>]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default="source_coincidence_test9.h5")
    parser.add_option("-o", type=str, dest="outfile",
                      help="output file name", default="output_test.h5")
    parser.add_option("--output_size", type=int, dest="interpolation_size",
                      help="size of output waveforms in transformed dataset", default=101)
    parser.add_option("--ch", type=int, dest="channel",
                      help="channel to process (between 1 and 8)", default=1)
    parser.add_option("--cfd", type=float, dest="cfd",
                      help="fraction of maximum for CFD time estimation", default=None)
    parser.add_option("--thresh", type=float, dest="threshold",
                      help="sample amplitude threshold value", default=0.7)
    parser.add_option("--baseline", type=float, dest="baseline",
                      help="sample fraction before waveform corresponding to noise", default=3)
    (options, args) = parser.parse_args()

    channel = options.channel
    interpolation_size = options.interpolation_size
    cfd_fraction = options.cfd
    
    if channel < 1 or channel > 8:
        sys.exit(f'{options.channel} is not a valid channel!')

    if cfd_fraction is not None and (cfd_fraction < 1e-5 or cfd_fraction > 1.0):
        sys.exit("CFD fraction must be a number between 0 and 1!")

    if cfd_fraction is None:
        file_label = 'maxVal'
    else:
        file_label = f'CFD{int(cfd_fraction*10)}'
        
    input_file = options.infile
    output_file = options.outfile+'.h5'
    
    infile_h5 = h5py.File(input_file, 'r')
    outfile_h5 = h5py.File(output_file, 'w')

    samples = infile_h5[f'ch{channel}_samples']
    time_begin_arr = infile_h5[f'ch{channel}_trig_offset']
    time_end_arr = infile_h5[f'ch{channel}_horiz_scale']
    num_trig = len(infile_h5[f'ch{channel}_samples'])
    num_samples  = len(infile_h5[f'ch{channel}_samples'][0])
    baseline = int(num_samples/options.baseline)
    
    outfile_h5.create_dataset(f'ch{channel}_samples', (num_trig, interpolation_size), dtype='f8')
    outfile_h5.create_dataset(f'ch{channel}_time', (num_trig, interpolation_size), dtype='f8')
    
    density_matrix = np.zeros((interpolation_size, interpolation_size))
    fixed_time_array = np.zeros(num_trig)

    #If sample input dataset is not converted, get vertical conversion factors
    isConverted = True
    vert_offset = None
    vert_scale = None

    if samples.dtype ==	'int8':
        isConverted = False
        vert_offset = infile_h5[f'ch{channel}_vert_offset']
        vert_scale = infile_h5[f'ch{channel}_vert_scale']
        
    for trig, waveform in enumerate(samples[:]):

        if trig%1000 == 0:
            sys.stdout.write(f"\rProcessing waveform {trig}")
            sys.stdout.flush()

        time_begin = time_begin_arr[trig]*1e9
        time_end = time_end_arr[trig]*(num_samples-1)*1e9+time_begin

        if not isConverted:
            waveform = -vert_offset[trig]+waveform*vert_scale[trig]
        
        #Create pulse object from each waveform in dataset 
        pulse = pls.Pulse(waveform, time_begin, time_end)
        pulse.subtract_pedestal()
        pulse.flip_sign()

        #create empty buffer for interpolated waveform in case CFD fails
        transformed_waveform = np.zeros(interpolation_size)
        transformation_time = np.zeros(interpolation_size)

        if pulse.get_max_amplitude() < options.threshold:

            if cfd_fraction is None:
                max_time, max_amplitude = pulse.get_interpolated_max(50)
                if abs(max_time) < 5:
                    transformation_time, transformed_waveform = pulse.get_interpolated_waveform(size=interpolation_size, time_window=20, fixed_time=max_time)
                    fixed_time_array[trig] = max_time

            else:
                cfd_time = pulse.get_cfd_time(fraction=cfd_fraction, baseline=baseline)#, smooth=True)
                if cfd_time is not None:# and abs(cfd_time) < 5:
                    transformation_time = np.linspace(cfd_time-2, cfd_time+20, interpolation_size)
                    transformation_function = interpolate.interp1d(x=pulse.get_time_array(), y=pulse.get_waveform(), kind='cubic')
                    transformed_waveform = transformation_function(transformation_time)
                    #transformation_time, transformed_waveform = pulse.get_interpolated_waveform(size=interpolation_size, time_window=20, fixed_time=cfd_time)
                    fixed_time_array[trig] = cfd_time
                
        outfile_h5[f'ch{channel}_samples'][trig] = transformed_waveform
        outfile_h5[f'ch{channel}_time'][trig] = transformation_time

        if np.any(transformed_waveform):
            norm_pulse = np.array(transformed_waveform)/math.sqrt(np.inner(transformed_waveform, transformed_waveform))
            #pulse = pls.Pulse(norm_pulse, transformation_time[0], transformation_time[-1])
            #pulse.save_plot(f'pulse_plots/pulse_{trig}')
            density_matrix += np.outer(norm_pulse, norm_pulse)
        
    density_matrix /= density_matrix.trace()
    print(density_matrix.size)

    eig_value_array, eig_vector_array = schur(density_matrix)
    pulse = pls.Pulse(eig_vector_array[:,0], 0, num_samples)
    pulse.save_plot(name=options.outfile+'_lead_eigvector', isDot=True)

    outfile_h5.create_dataset(f'ch{channel}_fixed_time_array', data=fixed_time_array)
    outfile_h5.create_dataset(f'ch{channel}_density_matrix', data=density_matrix)
    outfile_h5.close()
    infile_h5.close()
