import numpy as np
class DataPrep:
	def __init__(self,dataset):
		self.dataset = np.load(dataset, allow_pickle=True)
		x, y = self.dataset['X'], self.dataset['y']
		self.x_ns = x*156
		self.y_ns = y*156
		
	def getSamples(self):
		return self.x_ns, self.y_ns