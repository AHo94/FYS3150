import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project2/build-Project2_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

class Proj2_plot_script:
	def __init__(self, save_figure):
		self.n = 0
		self.save_fig = save_figure		# Saves as a figure if True
		# The omega values of the filename
		self.omega_value = np.array(['0.01', '0.50', '1.00', '5.00'])
		self.omega_filenames = []
		# For loops that creates filenames that is used for datafiles
		for omega_points in self.omega_value:
			self.omega_filenames.append('Eigenvector_data_omega_'+omega_points+'.txt')

	def read_data(self, filename_open):
		""" Saves the datapoints and returns the values in arrays """
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			n_value = line.split()
			if i >= 2:
				data.append(n_value)
			i += 1
		filename.close()
		self.n = int(data[0][0])
		rho_max = float(data[0][1])
		omega = float(data[0][2])
		self.rho = np.zeros(self.n+2)
		self.rho[-1] = rho_max
		eig = np.zeros(self.n+2)
		for j in range(0, self.n, 1):
			self.rho[j+1] = float(data[j+1][0])
			eig[j+1] = float(data[j+1][1])

		return eig

	def plot_data(self):
		""" Plots the data """
		fig1 = plt.figure()
		plt.hold("on")
		for filename in self.omega_filenames:
			eigenvec = self.read_data(filename)
			plt.plot(self.rho, eigenvec)

		fig_filename = 'Plot_groundstate_n'+str(self.n)+'.pdf'
		if self.save_fig:
			fig1.savefig('../Plots/'+fig_filename)
		else:
			plt.xlabel(r'$\rho$')
			plt.ylabel('$|\psi|^2$')
			plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$, ground state, n=%g' %(self.n))
			plt.legend(['$\omega=0.01$','$\omega=0.5$','$\omega=1$', '$\omega=5$'])
			plt.hold("off")
			plt.show()


solve = Proj2_plot_script(True)
solve.plot_data()