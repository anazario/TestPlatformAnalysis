import numpy as np
import matplotlib.pyplot as plt
import sys

class BSpline:

    def __init__(self, samples, knot_vector, order=3, precision=1e-8, clampRight=False, clamp=False):
        
        self.samples = samples
        self.knot_vector = list(knot_vector)
        self.order = order
        self.precision = precision

        self.clamp = clamp
        if clamp:
            self.og_knot_vector	= knot_vector
            self.__clamp_ends()

        self.clampRight = clampRight 
        if clampRight and not clamp:
            self.og_knot_vector = knot_vector
            self.__clamp_right()
        
        self.spline = self.__spline_function()

    def __clamp_ends(self):
        if self.clamp:
            for i in range(self.order):
                self.knot_vector.insert(0, self.knot_vector[0])
                self.knot_vector.insert(-1, self.knot_vector[-1])

    def __clamp_right(self):
        if self.clampRight:
            for i in range(self.order):
                self.knot_vector.insert(-1, self.knot_vector[-1])
                
    def __find_interval(self, value):

        knot_vector = self.knot_vector
        n_knots = len(knot_vector)
    
        if value < knot_vector[0] or value > knot_vector[-1]:
            sys.exit("Value must be within the domain of the knot vector!")

        for i in range(n_knots-1):
            if value >= knot_vector[i] and value < knot_vector[i+1]:
                return i

    def __get_basis_value(self, value):

        order = self.order
        knot_vector = self.knot_vector
        n_knots = len(knot_vector)
        
        if order > n_knots-2:
            order = n_knots-2
            print(f'Warning: Order reduced to maximum possible!')
        if order < 0:
            sys.exit("The polynomial order has to be a positive integer!")

        def omega(x, ti, tik):
            if ti == tik:
                return 0.
            else:
                return (x-ti)/(tik-ti)
            
        basis = np.zeros((n_knots-1,n_knots-1))
        basis[0][self.__find_interval(value)] = 1

        for k in range(1, order+1):
            N = basis[k-1]

            for i in range(n_knots-k-1):
                basis[k][i] = omega(value, knot_vector[i], knot_vector[i+k])*N[i] + (1. - omega(value, knot_vector[i+1], knot_vector[i+k+1]))*N[i+1]
                    
        return basis[order][0:n_knots-order-1]
    
    def __get_collocation_matrix(self, size=100, precision=1e-8):

        order = self.order
        n_knots = len(self.knot_vector)
        domain = np.linspace(self.knot_vector[0], self.knot_vector[-1]-precision, size)
        collocationMatrix = np.zeros((n_knots-order-1, size))

        for i,x in enumerate(domain):
            collocationMatrix[:,i] = self.__get_basis_value(x)

        return collocationMatrix

    def __fit_spline(self, collocationMatrix):

        aTa = np.matmul(collocationMatrix, collocationMatrix.T)
        aTa_inv = np.linalg.inv(aTa)
        aTa_inv_aT = np.matmul(aTa_inv, collocationMatrix)
        fit_coeff = np.matmul(aTa_inv_aT, self.samples)

        return fit_coeff

    def __spline_function(self):

        collocationMatrix = self.__get_collocation_matrix(size=len(self.samples), precision=self.precision)

        #self.plot_basis('test_basis')
        self.fit_coeff = self.__fit_spline(collocationMatrix)

        #construct spline function using fit coefficients
        spline = [np.matmul(self.fit_coeff, basis_slice) for basis_slice in collocationMatrix.T]

        #calculate error
        diff = np.subtract(spline, self.samples)
        self.error = 0.5*np.matmul(diff, diff)
        
        return spline

    def get_collocation_matrix(self, size=100, precision=1e-8):
        return self.__get_collocation_matrix(size=size, precision=precision)
    
    def spline_function(self):
        return self.spline

    def get_fit_coeff(self):
        return self.fit_coeff

    def get_error(self):
        #print(f'knot vector: {self.knot_vector}, error: {self.error}')
        return self.error

    def plot_fit(self, name):
        from scipy import interpolate

        knot_vector = self.knot_vector
        size = len(self.samples)
        domain = np.linspace(knot_vector[0], knot_vector[-1], size)
        interpolated_function = interpolate.interp1d(x=domain, y=self.samples, kind='cubic')
        interpolated_waveform_at_knots = interpolated_function(knot_vector)
        
        plt.xlabel('Domain')
        plt.ylabel('Amplitude')
        plt.plot(knot_vector, interpolated_waveform_at_knots, 'o', color='black', label='Interpolation at Knots')
        plt.plot(domain, self.samples, label='Interpolated data')
        plt.plot(domain, self.spline, label='B-Spline Fit')
        plt.grid()
        plt.text(0.75,0.25,f'Error: {round(self.error,5)}')
        legend=plt.legend(loc="upper right", fancybox=True)
        frame = legend.get_frame()
        frame.set_edgecolor('black')
        frame.set_alpha(1)
        plt.savefig(name+'.pdf')
        print(f'Saved plot: {name}.pdf')
        plt.clf()
    
    def plot_basis(self, name, size=100, precision=1e-8):

        order = self.order
        knot_vector = self.knot_vector
        n_knots = len(knot_vector)
        domain = np.linspace(knot_vector[0], knot_vector[-1]-precision, size)
        collocationMatrix = self.__get_collocation_matrix(size)

        for idx, arr in enumerate(collocationMatrix):
            plt.plot(domain, arr, label=f'B({idx})')

        y_knots = np.zeros(n_knots)
        plt.plot(knot_vector, y_knots, 'o', color='black')

        plt.xlabel('Domain')
        plt.ylabel('Spline Basis Value')
        plt.grid()
        legend=plt.legend(loc="upper right", title=f'k={order}', fontsize='small', fancybox=True)
        frame = legend.get_frame()
        frame.set_edgecolor('black')
        frame.set_alpha(1)
        plt.savefig(f'{name}_basis_order{order}_size{int(size)}.pdf')
        print(f'Saved plot: {name}_basis_order{order}_size{int(size)}.pdf')
        plt.clf()
        
