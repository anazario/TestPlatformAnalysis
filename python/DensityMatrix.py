import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import schur
from pulseTools import *

class DensityMatrix:

    def __init__(self, density_matrix):
        self.eig_value_array, self.eig_vector_array = schur(density_matrix)
        self.eig_value_array = np.array([val[idx] for idx, val in enumerate(self.eig_value_array)])
        self.size = self.eig_value_array.size
        
        lead_eigvector = self.eig_vector_array[:,0]
        
        if np.max(lead_eigvector) < -np.min(lead_eigvector):
            self.eig_vector_array = -1*self.eig_vector_array.T
        else:
            self.eig_vector_array = self.eig_vector_array.T
        
    def get_eig_values(self):
        return self.eig_value_array
        
    def get_eig_vectors(self):
        return self.eig_vector_array

    def plot_eig_vector(self, index, name):
        domain=np.linspace(0., 1., self.size)
        name += f'_eig_vector_index{index}'
        plot_line(x_array=domain, y_array=self.eig_vector_array[index], name=name, xlabel='Domain', ylabel='Eigenvector Value', title='', leg_labels=None)
        
