import math
import numpy as np

def GetBasisInvKLJ(sample_dataset, base_line = -1):

    num_triggers = len(sample_dataset)
    num_samples  = len(sample_dataset[0])

    if base_line == -1:
        base_line = num_samples

    density_matrix = np.zeros((base_line,base_line))

    for waveform in sample_dataset:
        waveform = np.mean(waveform[:int(num_samples/3)]) - waveform
        norm_waveform = waveform/math.sqrt(np.inner(waveform, waveform))
        density_matrix += np.outer(norm_waveform, norm_waveform)

    density_matrix /= density_matrix.trace()

    eig_value_array, eig_vector_array = np.linalg.eig(density_matrix)

    return eig_value_array, eig_vector_array

def KLJ_filter(sample_dataset, eig_begin = 5, base_line = -1):

    num_samples = len(sample_dataset[0])

    if base_line == -1:
        base_line = num_samples
    
    eig_value_array, eig_vector_array = GetBasisInvKLJ(sample_dataset)
    
    I = np.identity(base_line)
    noise_matrix = np.zeros((base_line,base_line))

    for idx in range(eig_begin, base_line):
        noise_vec = eig_vector_array[:,idx]
        noise_matrix += np.outer(noise_vec, noise_vec)

    filter_invKLJ = I - noise_matrix

    return filter_invKLJ

