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
		self.T = np.array([float(data[0][2]), float(data[-1][2])])
		self.MC_cycles = np.zeros(N_half)
		self.E_mean_T1 = np.zeros(N_half)
		self.E_mean_T24 = np.zeros(N_half)
		self.M_mean_T1 = np.zeros(N_half)
		self.M_mean_T24 = np.zeros(N_half)

		for i in range(0,N_half):
			self.MC_cycles[i] = float(data[i][1])
			self.E_mean_T1[i] = float(data[i][4])
			self.E_mean_T24[i] = float(data[i+N_half][4])
			self.M_mean_T1[i] = float(data[i][6])
			self.M_mean_T24[i] = float(data[i+N_half][6])
			
	def plot_state(self):
		self.read_data("4c.txt")
		plt.plot(np.log(self.MC_cycles), self.E_mean_T1, 'b-')
		plt.hold("on")
		plt.plot(np.log(self.MC_cycles), self.M_mean_T1, 'r-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E\rangle$')
		plt.show()

solver = Plotter(False)
solver.plot_state()