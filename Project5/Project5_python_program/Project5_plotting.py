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
		""" Reads data and saves values to variables and arrays"""
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
		self.Variance = np.zeros(N)
		self.Energy_Exact = 3*self.omega*np.ones(N)
		self.MeanDistance = np.zeros(N)
		for j in range(0, N):
			self.alpha[j] = float(data[j][0])
			self.Energy[j] = float(data[j][1])
			self.Variance[j] = float(data[j][2])
			self.MeanDistance[j] = float(data[j][3])


	def plot_energy_alpha(self):
		self.read_data("Energy_Alpha_Mdistance_omega_0.01.txt")
		fig1 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.hold("on")
		plt.plot(self.alpha, self.Energy_Exact)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig2 = plt.figure()
		plt.plot(self.alpha, self.MeanDistance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$r_{12}$')
		plt.title(r'Plot of the mean distance $r_{12}$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		fig3 = plt.figure()
		plt.plot(self.alpha, self.Variance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\sigma_E$')
		plt.title(r'Plot of the variance $\sigma_E$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		self.read_data("Energy_Alpha_Mdistance_omega_0.50.txt")
		fig4 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.hold("on")
		plt.plot(self.alpha, self.Energy_Exact)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig5 = plt.figure()
		plt.plot(self.alpha, self.MeanDistance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$r_{12}$')
		plt.title(r'Plot of the mean distance $r_{12}$ as a function of $\alpha$ with $\omega$ = %.2f.' %(self.omega))


		fig6 = plt.figure()
		plt.plot(self.alpha, self.Variance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\sigma_E$')
		plt.title(r'Plot of the variance $\sigma_E$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		self.read_data("Energy_Alpha_Mdistance_omega_1.00.txt")
		fig7 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.hold("on")
		plt.plot(self.alpha, self.Energy_Exact)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig8 = plt.figure()
		plt.plot(self.alpha, self.MeanDistance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$r_{12}$')
		plt.title(r'Plot of the mean distance $r_{12}$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		fig9 = plt.figure()
		plt.plot(self.alpha, self.Variance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\sigma_E$')
		plt.title(r'Plot of the variance $\sigma_E$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))
		if self.savefile == True:
			fig1.savefig('../Plots/Energy_alpha_plot_omega001.pdf')
			fig2.savefig('../Plots/MeanDistance_alpha_plot_omega001.pdf')
			fig3.savefig('../Plots/Variance_alpha_plot_omega001.pdf')

			fig4.savefig('../Plots/Energy_alpha_plot.pdf_omega05')
			fig5.savefig('../Plots/MeanDistance_alpha_plot_omega05.pdf')
			fig6.savefig('../Plots/Variance_alpha_plot_omega05.pdf')

			fig7.savefig('../Plots/Energy_alpha_plot.pdf_omega1')
			fig8.savefig('../Plots/MeanDistance_alpha_plot_omega1.pdf')
			fig9.savefig('../Plots/Variance_alpha_plot_omega1.pdf')
		else:
			plt.show()

	def testplot(self):
		self.read_data("Testing.txt")
		fig1 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.hold("on")
		plt.plot(self.alpha, self.Energy_Exact)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig2 = plt.figure()
		plt.plot(self.alpha, self.MeanDistance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$r_{12}$')
		plt.title(r'Plot of the mean distance $r_{12}$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))
		plt.show()
## Comment out the functions to plot what you want
solver = Plotter(False)
solver.plot_energy_alpha()
#solver.testplot()