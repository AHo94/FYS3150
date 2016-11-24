import numpy as np
import matplotlib.pyplot as plt
import os
# Change directory based based on work place, i.e home or at UiO computer.
file_directory = '../build-Project5_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'
#file_directory = '../build-Project5_cpp_program-Desktop_Qt_5_7_0_GCC_64bit-Debug'

class Plotter():
	def __init__(self, savefile):
		self.savefile = savefile	# If True, saves the plots to a file

	def read_data(self, filename_open):
		""" Reads data, specific for task 4c """
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
		self.Energy = np.zeros(N)
		self.alpha = np.zeros(N)
		for j in range(0, N):
			self.Energy[j] = float(data[j][1])
			self.alpha[j] = float(data[j][0])

	def plot_energy_alpha(self):
		self.read_data("Energy_alpha.txt")
		plt.plot(self.alpha, self.Energy)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the energy as a functio no $\alpha$ with $\omega$ = 1')
		plt.show()	

## Comment out the functions to plot what you want
solver = Plotter(False)
solver.plot_energy_alpha()
