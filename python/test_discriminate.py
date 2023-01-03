import h5py
import numpy as np
import pulse as pls
import DensityMatrix as dm

from pulseTools import *

if __name__ == '__main__':
    import optparse
    import sys

    usage = "usage: %prog [-n]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default=None)
    parser.add_option("--bkgch", type=int, dest="bkg_channel",
                      help="background channel (between 1 and 8)", default=2)
    parser.add_option("--sigch", type=int, dest="sig_channel",
                      help="signal channel (between 1 and 8)", default=3)
    parser.add_option("--thresh", type=float, dest="threshold",
                      help="sample amplitude threshold value", default=0.1)
    parser.add_option("--output_size", type=int, dest="interpolation_size",
                      help="size of output waveforms in transformed dataset", default=101)
    (options, args) = parser.parse_args()

    bkg_channel = options.bkg_channel
    sig_channel = options.sig_channel
    interpolation_size = options.interpolation_size
    cfd_fraction = 0.4

    output_file = 'test_discrimination1.h5'
    
    if options.infile is None:
        sys.exit("Please specify an input file (-i [input_file_name])!")
        
    if bkg_channel == sig_channel:
        sys.exit("background and signal channels must be different!")

    if bkg_channel < 1 or bkg_channel > 8:
        sys.exit("background channel out of bounds! Please choose a number between 1 and 8")

    if sig_channel < 1 or sig_channel > 8:
        sys.exit("signal channel out of bounds! Please choose a number between 1 and 8")

    input_file = options.infile

    infile_h5 = h5py.File(input_file,'r')
    outfile_h5 = h5py.File(output_file, 'w')

    sig_samples = infile_h5[f'ch{sig_channel}_samples']
    sig_time_begin_arr = infile_h5[f'ch{sig_channel}_trig_offset']
    sig_time_end_arr = infile_h5[f'ch{sig_channel}_horiz_scale']

    bkg_samples = infile_h5[f'ch{bkg_channel}_samples']
    bkg_time_begin_arr = infile_h5[f'ch{bkg_channel}_trig_offset']
    bkg_time_end_arr = infile_h5[f'ch{bkg_channel}_horiz_scale']
    
    num_trig = len(infile_h5[f'ch{sig_channel}_samples'])
    num_samples  = len(infile_h5[f'ch{sig_channel}_samples'][0])
    baseline = int(num_samples/3)

    isConverted = True
    vert_offset = None
    vert_scale = None

    if sig_samples.dtype == 'int8':
        isConverted = False
        sig_vert_offset = infile_h5[f'ch{sig_channel}_vert_offset']
        sig_vert_scale = infile_h5[f'ch{sig_channel}_vert_scale']
        bkg_vert_offset = infile_h5[f'ch{bkg_channel}_vert_offset']
        bkg_vert_scale = infile_h5[f'ch{bkg_channel}_vert_scale']

    sig_cnt = 0
    bkg_cnt = 0

    sig_density_matrix = np.zeros((interpolation_size, interpolation_size))
    bkg_density_matrix = np.zeros((num_samples, num_samples))
    #bkg_density_matrix = np.zeros((interpolation_size, interpolation_size))
    
    proc_sig_list = []
    time_sig_list = []
    
    proc_bkg_list = []
    time_bkg_list = []

    for trig, waveform in enumerate(sig_samples[:]):

        if trig%1000 == 0:
            sys.stdout.write(f"\rProcessing waveform {trig}")
            sys.stdout.flush()

        sig_time_begin = sig_time_begin_arr[trig]*1e9
        sig_time_end = sig_time_end_arr[trig]*(num_samples-1)*1e9+sig_time_begin

        bkg_time_begin = bkg_time_begin_arr[trig]*1e9
        bkg_time_end = bkg_time_end_arr[trig]*(num_samples-1)*1e9+bkg_time_begin

        bkg_waveform = bkg_samples[trig]
        
        if not isConverted:
            waveform = -sig_vert_offset[trig]+waveform*sig_vert_scale[trig]
            bkg_waveform = -bkg_vert_offset[trig]+bkg_waveform*bkg_vert_scale[trig]

        sig_pulse = pls.Pulse(waveform, sig_time_begin, sig_time_end)
        bkg_pulse = pls.Pulse(bkg_waveform, bkg_time_begin, bkg_time_end)

        sig_pulse.subtract_pedestal()
        sig_pulse.flip_sign()
        sig_max_amplitude =  sig_pulse.get_max_amplitude()

        bkg_pulse.subtract_pedestal()
        bkg_pulse.flip_sign()

        proc_sig = np.zeros(interpolation_size)
        time_sig = np.zeros(interpolation_size)
        
        proc_bkg = np.zeros(num_samples)
        time_bkg = np.zeros(num_samples)

        if sig_max_amplitude > options.threshold and sig_max_amplitude < 0.65:
            
            time_sig, proc_sig = cfd_aligned_waveform(sig_pulse, interpolation_size, cfd_fraction, baseline)
            #time_bkg, proc_bkg = cfd_aligned_waveform(bkg_pulse, interpolation_size, cfd_fraction, baseline)
            #time, proc = cfd_aligned_waveform(bkg_pulse, interpolation_size, cfd_fraction, baseline)
            
            time_bkg = np.linspace(bkg_time_begin, bkg_time_end, num_samples)
            proc_bkg = bkg_waveform
            
            if np.any(proc_sig):
                proc_sig_list.append(proc_sig)
                time_sig_list.append(time_sig)
                sig_cnt +=1
                
            if np.any(proc_bkg) and np.max(proc_bkg) < options.threshold:
                proc_bkg_list.append(proc_bkg)
                time_bkg_list.append(time_bkg)
                bkg_cnt +=1

        if np.any(proc_sig):
            norm_pulse = np.array(proc_sig)/math.sqrt(np.inner(proc_sig, proc_sig))
            sig_density_matrix += np.outer(norm_pulse, norm_pulse)

        if np.any(proc_bkg):
            norm_pulse = np.array(proc_bkg)/math.sqrt(np.inner(proc_bkg, proc_bkg))
            bkg_density_matrix += np.outer(norm_pulse, norm_pulse)

    sig_density_matrix /= sig_density_matrix.trace()
    bkg_density_matrix /= bkg_density_matrix.trace()
    
    print('\n')
    outfile_h5.create_dataset(f'signal_samples_ch{sig_channel}', data=proc_sig_list)
    outfile_h5.create_dataset(f'signal_time_ch{sig_channel}', data=time_sig_list)
    outfile_h5.create_dataset(f'sig_density_matrix_ch{sig_channel}', data=sig_density_matrix)

    outfile_h5.create_dataset(f'background_samples_ch{bkg_channel}', data=proc_bkg_list)
    outfile_h5.create_dataset(f'background_time_ch{bkg_channel}', data=time_bkg_list)
    outfile_h5.create_dataset(f'bkg_density_matrix_ch{bkg_channel}', data=bkg_density_matrix)

    den_sig = dm.DensityMatrix(sig_density_matrix)
    den_sig.plot_eig_vector(0, 'sig')

    den_bkg = dm.DensityMatrix(bkg_density_matrix)
    den_bkg.plot_eig_vector(0, 'bkg')

    print(den_bkg.get_eig_values())
    
    print(f'\nSignal: {sig_cnt}, Background: {bkg_cnt}')
        
