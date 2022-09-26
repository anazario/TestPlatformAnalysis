import numpy as np
import h5py
from algorithms import algos

#class to prepare hdf5 files for python analysis
class DataPrep:
	def __init__(self,infile):
		self.file = h5py.File(infile,'r')
		self.num_triggers = len(self.file['ch1_samples'])
		self.num_samples = len(self.file['ch1_samples'][0])
		time_begin = self.file['ch1_trig_offset'][0]*1e9
		time_end = self.file['ch1_horiz_scale'][0]*(self.num_samples-1)*1e9+time_begin
		self.time = np.linspace(time_begin, time_end, self.num_samples)
		self.sampling_freq = self.num_samples/(time_end-time_begin)
		self.channels = [1,2,3,4]
		self.ref_ch = 1e9
		self.amplitudeLimits = {1:(0.1,0.65), 2:(0.1,0.65), 3:(0.1,0.65), 4:(0.1,0.65)}
		self.smooth = False
		self.kernel_sigma = 1
		self.kernel_width = 5*self.kernel_sigma

		self.baseline_fraction = 0.3
		self.cfd_fraction = 0.4

		self.samples = []
		self.horiz_scale = []
		self.trig_offset = []

		self.cfd_times = []
		self.rise_times = []

		self.waveforms = []

		print("CFD fraction:",self.cfd_fraction,"baseline fraction:",self.baseline_fraction)

	def setCFDfraction(self,frac):
		self.cfd_fraction = frac

	def setBaselineFraction(self,frac):
		self.baseline_fraction = frac


	def getSamplingFrequency(self):
		return self.sampling_freq

	def prepAndCalculate(self,channel):
		#prep samples per channel
		self.samples = self.file[f'ch{channel}_samples']
		self.horiz_scale = self.file[f'ch{channel}_horiz_scale']
		self.trig_offset = self.file[f'ch{channel}_trig_offset']
		for trig in range(self.num_triggers):
			waveform = self.samples[trig]
			baseline_n = int(len(waveform)*self.baseline_fraction)
			waveform = np.mean(waveform[:baseline_n]) - waveform
			amp = np.max(waveform)
			a = algos(waveform)

			#do CFD calculation
			if amp>self.amplitudeLimits[channel][0] and amp<self.amplitudeLimits[channel][1]:
				cfd_val = a.cfd(threshold=self.amplitudeLimits[channel][0],base_line = baseline_n, fraction=self.cfd_fraction, smooth=self.smooth, kernel_width=self.kernel_width, kernel_sigma=self.kernel_sigma)
				if cfd_val is not None:
					self.cfd_times.append((cfd_val*self.horiz_scale[trig]+self.trig_offset[trig])*1e9)
				else:
					self.cfd_times.append(-900)

			#do risetime calculation
			cfd_rt20 = a.cfd(threshold=self.amplitudeLimits[channel][0], base_line=baseline_n, fraction=0.2, smooth=False)
			cfd_rt80 = a.cfd(threshold=self.amplitudeLimits[channel][0], base_line=baseline_n, fraction=0.8, smooth=False)
			if cfd_rt20 is not None and cfd_rt80 is not None:
				self.rise_times.append( (cfd_rt80-cfd_rt20)*self.horiz_scale[trig]*1e9 )

			#append prepped waveform to class object
			self.waveforms.append(waveform)
	



#get array of CFD times (1 per trigger) in ns
	def getCFDtimes(self):
		if self.cfd_times != []:
			return self.cfd_times


	def getRiseTimes(self):
		if self.rise_times != []:
			return self.rise_times

	def getWaveforms(self):
		if self.waveforms != []:
			return self.waveforms

#returns time axis
	def getTimeAxis(self):
		if self.time != []:
			return self.time
            




