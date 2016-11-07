import numpy as np
import matplotlib.pyplot as plt
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

class Plotter():
	def __init__(self, savefile):
		self.savefile = savefile	# If True, saves the plots to a file

	def read_data(self, filename_open):
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				data.append(data_set)
			i += 1
		filename.close()

		N = len(data)
		N_half = int(N/2.0)
		self.L = int(data[0][0])
		self.MC_cycles = np.zeros(N_half)
		self.Mean_E_T1 = np.zeros(N_half)
		self.Mean_E_T24 = np.zeros(N_half)
		for i in range(0,N_half):
			self.MC_cycles[i] = float(data[i][1])
			self.Mean_E_T1[i] = float(data[i][4])
			self.Mean_E_T24[i] = float(data[i+N_half][4])

	def plot_state(self):
		self.read_data("4c.txt")
		plt.plot(np.log(self.MC_cycles), self.Mean_E_T1)
		plt.show()

solver = Plotter(False)
solver.plot_state()