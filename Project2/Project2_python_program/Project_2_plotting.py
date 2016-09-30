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
		self.omega = float(data[0][2])
		self.rho = np.zeros(self.n+2)
		self.rho[-1] = rho_max
		eig = np.zeros(self.n+2)
		self.V = np.zeros(self.n+2)
		self.E = np.zeros(self.n+2)
		for j in range(0, self.n, 1):
			self.rho[j+1] = float(data[j+1][0])
			eig[j+1] = float(data[j+1][1])
			self.V[j+1] = float(data[j+1][2])
			self.E[j+1] = float(data[j+1][3])
		# Normalizing
		new_eig = np.sum(eig)*(rho_max)/(self.n+1)
		factor = 1.0/new_eig
		eig *= factor

		# Placing boundary to the potential and E (to prevent them for being 0)
		self.V[0] = self.omega**2*self.rho[0]**2 + 1.0/self.rho[0]
		self.V[-1] = self.omega**2*self.rho[-1]**2 + 1.0/self.rho[-1]
		self.E[0] = self.E[1]
		self.E[-1] = self.E[-2]
		return eig

	def plot_data(self):
		""" Plots the data """
		fig1 = plt.figure()
		plt.hold("on")
		for filename in self.omega_filenames:
			eigenvec = self.read_data(filename)
			plt.plot(self.rho, eigenvec)

		plt.xlabel(r'$\rho$')
		plt.ylabel('$|\psi|^2$')
		plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$, ground state, n=%g' %(self.n))
		plt.legend([r'$\omega=0.01, \rho_{max} = 40$',r'$\omega=0.5, \rho_{max}=5$',r'$\omega=1, \rho_{max}=4$',r'$\omega=5, \rho_{max}=2$'])
		plt.hold("off")
		# Specified filename for a set number of meshpoints n
		fig_filename = 'Plot_groundstate_n'+str(self.n)+'.pdf'

		# Plotting example of tunneling effect
		fig2, ax1 = plt.subplots()
		Tunneling_eigvec = self.read_data(self.omega_filenames[2])
		ax1.plot(self.rho, Tunneling_eigvec)
		ax1.set_xlabel(r'$\rho$')
		ax1.set_ylabel('$|\psi|^2$')

		ax2 = ax1.twinx()
		ax2.plot(self.rho, self.V, 'r-')
		plt.hold("on")
		ax2.plot(self.rho, self.E, 'g-')
		ax2.set_ylabel('$E$')
		plt.hold("off")
		plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$. n=%g, $\omega$ = %f' %(self.n, self.omega))
		ax2.legend(['Potential', 'Energy'])

		tunneling_filename = 'Plot_groundstate_tunneling_effect.pdf'

		if self.save_fig:
			fig1.savefig('../Plots/'+fig_filename)
			fig2.savefig('../Plots/'+tunneling_filename)
		else:
			plt.show()


solve = Proj2_plot_script(True)
solve.plot_data()