import numpy as np
import matplotlib.pyplot as plt
import sys

import pulse as pls

def find_interval(value, knot_vector):

    n_knots = len(knot_vector)
    if value < knot_vector[0] or value > knot_vector[-1]:
        print(value)
        sys.exit("Value must be within the domain of the knot vector!")

    for i in range(n_knots-1):
        if value >= knot_vector[i] and value < knot_vector[i+1]:
            return i

def get_bSpline_value(value, knot_vector, order):

    n_knots = len(knot_vector)
    if order > n_knots-2:
        order = n_knots-2
        print(f'Warning: Order reduced to maximum possible!')
    if order < 0:
        sys.exit("The polynomial order has to be a positive integer!")

    basis = np.zeros((n_knots-1,n_knots-1))
    basis[0][find_interval(value, knot_vector)] = 1
    
    for k in range(1, order+1):
        N = basis[k-1]

        for i in range(n_knots-k-1):
            basis[k][i] = (value-knot_vector[i])*N[i]/(knot_vector[i+k]-knot_vector[i]) + (knot_vector[i+k+1]-value)*N[i+1]/(knot_vector[i+k+1]-knot_vector[i+1])
    
    return basis[order][0:n_knots-order-1]

def get_bSpline_array(knot_vector, order=3, size=100, precision=0.001):

    n_knots = len(knot_vector)
    domain = np.linspace(0., 1.-precision, size)
    bSplineMatrix = np.zeros((n_knots-order-1, len(domain)))

    for i,x in enumerate(domain):
        bSplineMatrix[:,i] = get_bSpline_value(x, knot_vector, order)

    return bSplineMatrix

def fit_spline(bSplineMatrix, samples):
    aTa = np.matmul(bSplineMatrix,bSplineMatrix.T)
    aTa_inv = np.linalg.inv(aTa)
    aTa_inv_aT = np.matmul(aTa_inv,bSplineMatrix)
    test_coeff = np.matmul(aTa_inv_aT, samples)
    
    return test_coeff

def spline_function(samples, knot_vector, order=3, size=100, precision=0.001):
    size = len(samples)
    n_knots = len(knot_vector)
    domain = np.linspace(0, 1.-precision, size)
    bSplineMatrix = get_bSpline_array(knot_vector, order, size)
    fit_coeff = fit_spline(bSplineMatrix, samples)
    error = 0.
    
    #create spline function
    spline = np.zeros(size)
    for i in range(size):
        total = 0.
        for idx, basis_val in enumerate(bSplineMatrix.T[i]):
            total+=basis_val*fit_coeff[idx]
        spline[i] = total
        diff = spline[i]-samples[i]
        error += diff*diff

    return error, spline
        
def plot_basis(name, knot_vector, order=3, size=100, precision=0.001):

    n_knots = len(knot_vector)
    domain = np.linspace(0., 1.-precision, size)
    bSplineMatrix = get_bSpline_array(knot_vector, order, size)

    for idx, arr in enumerate(bSplineMatrix):
        plt.plot(domain, arr, label=f'B({idx})')

    y_knots = np.zeros(n_knots)
    plt.plot(knot_vector, y_knots, 'o', color='black')

    plt.xlabel('Domain')
    plt.ylabel('Spline Basis Value')
    plt.grid()
    legend=plt.legend(loc="upper right", title=f'k={order}', fontsize='small', fancybox=True)
    frame = legend.get_frame()#sets up for color, edge, and transparency
    frame.set_edgecolor('black') #edge color of legend
    frame.set_alpha(1) #deals with transparency
    plt.savefig(f'{name}_basis_order{order}_size{int(size)}.pdf')
    print(f'Saved plot: {name}_basis_order{order}_size{int(size)}.pdf')
    plt.clf()

if __name__ == '__main__':
    import h5py
    import optparse
    from scipy.linalg import *
    from scipy import interpolate

    usage = "usage: %prog [-n]"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-i", type=str, dest="infile",
                      help="input file name", default=None)
    parser.add_option("--ch", type=int, dest="channel",
                      help="channel to process (between 1 and 8)", default=1)
    parser.add_option("-o", type=int, dest="order",
                      help="polynomial order for spline basis", default=3)
    (options, args) = parser.parse_args()

    order = options.order
    channel = options.channel
    input_file = options.infile
    
    f = h5py.File(input_file,'r')

    samples = f[f'ch{channel}_samples']
    fixed_time_array = f[f'ch{channel}_fixed_time_array']
    density_matrix = f[f'ch{channel}_density_matrix']
    eig_value_array, eig_vector_array = schur(density_matrix)

    lead_eigvector = eig_vector_array[:,0]
    eig_value_array = np.array([val[idx] for idx, val in enumerate(eig_value_array)])
    if np.max(lead_eigvector) < -np.min(lead_eigvector):
        eig_vector_array = -1*eig_vector_array.T
    else:
        eig_vector_array = eig_vector_array.T

    lead_eigvector = eig_vector_array[0]

    knot_vector = [0., 0.05, 0.085, 0.1, 0.14, 0.25, 0.4, 1.]
    #knot_vector = [0., 0.07, 0.1, 0.15, 0.3, 0.5, 0.75, 1.]
    #knot_size = len(knot_vector)
    '''
    for trig, waveform in enumerate(samples):

        if trig%1000 == 0:
            sys.stdout.write(f"\rProcessing waveform {trig}")
            sys.stdout.flush()

        if np.any(waveform) and fixed_time_array[trig] > -1.0:
            norm_waveform = waveform/math.sqrt(np.inner(waveform, waveform))
    '''
    #plot_basis('test_LGAD', knot_vector, 3, 1000, 0.001)

    pulse = pls.Pulse(samples[0], 0,1)#fixed_time_array[0], fixed_time_array[-1])
    pulse.save_plot('test_pulse')
    #waveform = pulse.get_waveform()
    waveform = lead_eigvector
    
    spline_size = 200
    domain = np.linspace(0,1, spline_size)
    interpolated_function = interpolate.interp1d(x=np.linspace(0, 1, len(lead_eigvector)), y=waveform, kind='cubic')
    interpolated_waveform = interpolated_function(domain)
    interpolated_waveform_at_knots = interpolated_function(knot_vector)
    #test_coeff = [0.3,0.38,0.12,0.07]

    error, spline = spline_function(interpolated_waveform, knot_vector)
    print(f'Error: {error}')
    
    plt.xlabel('Domain')
    plt.ylabel('Amplitude')
    plt.plot(knot_vector, interpolated_waveform_at_knots, 'o', color='black', label='Interpolation at Knots')
    plt.plot(domain, interpolated_waveform, label='Interpolated data')
    plt.plot(domain, spline, label='B-Spline Fit')
    plt.grid()
    plt.text(0.75,0.25,f'Error: {round(error,5)}')
    legend=plt.legend(loc="upper right", fancybox=True)
    frame = legend.get_frame()
    frame.set_edgecolor('black')
    frame.set_alpha(1)
    plt.show()
    
    
