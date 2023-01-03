import h5py

if __name__ == '__main__':
    import optparse
    import sys

    usage = "usage: %prog [-i <input file>] [-o <output file>] "
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default=None)
    parser.add_option("-o", type=str, dest="outfile",
                      help="output file name", default=None)
    (options, args) = parser.parse_args()
    
    infile = options.infile
    outfile = options.outfile

    if infile is None:
        print(usage)
        sys.exit('No input file specified!')

    if outfile is None:
        outfile = f'converted_{infile}'
    
    with h5py.File(outfile,'w') as f_dest:
        with h5py.File(infile,'r') as f_src:

            for channel in range(1,8):
                num_trig     = len(f_src[f'ch{channel}_samples'])
                num_samples  = len(f_src[f'ch{channel}_samples'][0])
                samples      = f_src[f'ch{channel}_samples']
                vert_offset  = f_src[f'ch{channel}_vert_offset']
                vert_scale   = f_src[f'ch{channel}_vert_scale']

                converted_samples=[]
                for trig, waveform in enumerate(samples):
                    converted_samples.append(-vert_offset[trig]+waveform*vert_scale[trig])

                f_dest.create_dataset(f'ch{channel}_samples', data=converted_samples)
                f_dest.create_dataset(f'ch{channel}_horiz_offset', data=f_src[f'ch{channel}_horiz_offset'])
                f_dest.create_dataset(f'ch{channel}_scale', data=f_src[f'ch{channel}_horiz_scale'])
                f_dest.create_dataset(f'ch{channel}_trig_time', data=f_src[f'ch{channel}_trig_time'])

            
