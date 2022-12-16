from scipy import signal
from scipy.ndimage import filters
import numpy as np
#algorithms for one (1) waveform
class algos:
    def __init__(self,samples,verbose=False):
        if len(samples.shape) > 1:
            print("This class is for analysis of one waveform.")
        self.samples = samples
        self.verbose = verbose

    def cfd(self, threshold=0.02, base_line=200, fraction=0.3, hysteresis=0.001, smooth=True, kernel_width=10, kernel_sigma=5):
        if len(self.samples.shape) > 1:
            if self.verbose:
                print("This class is for analysis of one waveform.")
            return None
        if self.samples.max()-self.samples.min() < threshold:
            if self.verbose:
                print("Samples did not meet threshold.")
            return None
        if smooth: #smoothing not working - filter size doesn't match?
            kernel = signal.gaussian(kernel_width, kernel_sigma)
            self.samples = filters.convolve(self.samples, kernel)
        self.samples -= np.mean(self.samples[:base_line])
        if abs(np.max(self.samples)) < abs(np.min(self.samples)):
            maximum = abs(self.samples.min())
            self.samples /= -maximum
        else:
            maximum = self.samples.max()
            self.samples /= maximum
        
        hysteresis = hysteresis/maximum
        above = False
        lockedForHysteresis = False
        numberOfPeaks = 0
        for i,v in enumerate(self.samples):
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
            v1 = self.samples[firstCrossingIndex]
            v2 = self.samples[firstCrossingIndex-1]
            return (t1 + (t2 - t1)*(fraction-v1)/(v2 - v1))
        else:
            if self.verbose:
                print("More or less than one peak. Number of peaks =",numberOfPeaks)
            return None
