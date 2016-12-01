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
			data_set = line.split()
			if i != 0:
				data.append(data_set)
			i += 1
		filename.close()

		N = len(data)
		self.omega = float(data[0][-1])
		self.Energy = np.zeros(N)
		self.alpha = np.zeros(N)
		self.Energy_Exact = 3*self.omega*np.ones(N)
		self.MeanDistance = np.zeros(N)
		for j in range(0, N):
			self.alpha[j] = float(data[j][0])
			self.Energy[j] = float(data[j][1])
			self.MeanDistance[j] = float(data[j][2])


	def plot_energy_alpha(self):
		self.read_data("Energy_Alpha_Mdistance.txt")
		fig1 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.hold("on")
		plt.plot(self.alpha, self.Energy_Exact)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.g. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig2 = plt.figure()
		plt.plot(self.alpha, self.MeanDistance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$r_{12}$')
		plt.title(r'Plot of the mean distance $r_{12}$ as a function of $\alpha$ with $\omega$ = %.g' %(self.omega))
		if self.savefile == True:
			fig1.savefig('../Plots/Energy_alpha_plot.pdf')
			fig2.savefig('../Plots/MeanDistance_alpha_plot.pdf')
		else:
			plt.show()

## Comment out the functions to plot what you want
solver = Plotter(False)
solver.plot_energy_alpha()
