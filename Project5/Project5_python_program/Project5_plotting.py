import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
# Change directory based based on work place, i.e home or at UiO computer.
file_directory = '../build-Project5_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Release'
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

		self.N = len(data)
		self.omega = float(data[0][-2])
		self.Energy = np.zeros(self.N)
		self.alpha = np.zeros(self.N)
		self.beta = np.zeros(self.N)
		self.Variance = np.zeros(self.N)
		self.MeanDistance = np.zeros(self.N)
		self.AcceptedConfigs = np.zeros(self.N)
		for j in range(0, self.N):
			self.alpha[j] = float(data[j][0])
			self.beta[j] = float(data[j][1])
			self.Energy[j] = float(data[j][2])
			self.Variance[j] = float(data[j][3])
			self.MeanDistance[j] = float(data[j][4])
			self.AcceptedConfigs[j] = float(data[j][-1])

	def read_data_virial(self, filename_open):
		""" Reads data and saves to variables. Specific for Virial problem """
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			data_set = line.split()
			if i != 0:
				data.append(data_set)
			i += 1
		filename.close()

		self.N = len(data)
		self.alpha = float(data[0][0])
		self.beta = float(data[0][1])
		self.omega = np.zeros(self.N)

		self.KineticExpect = np.zeros(self.N)
		self.PotentialExpect = np.zeros(self.N)

		for j in range(0, self.N):
			self.omega[j] = float(data[j][-2])
			self.KineticExpect[j] = float(data[j][2])
			self.PotentialExpect[j] = float(data[j][4])

	def stability_test(self):
		self.read_data("Stability_test_MC_100.txt")
		fig1 = plt.figure()
		ax = plt.subplot(111)
		ax.plot(self.alpha, self.Energy, label='$N_{MC}$ = 100')
		plt.hold("on")
		self.read_data("Stability_test_MC_1000.txt")
		ax.plot(self.alpha, self.Energy, label='$N_{MC}$ = 1000')
		self.read_data("Stability_test_MC_10000.txt")
		ax.plot(self.alpha, self.Energy, label='$N_{MC}$ = 10000')
		self.read_data("Stability_test_MC_100000.txt")
		ax.plot(self.alpha, self.Energy, label='$N_{MC}$ = 100000')
		self.read_data("Stability_test_MC_1000000.txt")
		ax.plot(self.alpha, self.Energy, label='$N_{MC}$ = 1000000')
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		ax.legend(loc='upper right', bbox_to_anchor=(1,1), ncol=1, fancybox=True)
		plt.title(r'Mean energy as a function of $\alpha$ with $\omega = 1$.'+' \n Running for different number of MC cycles.')
		if self.savefile == True:
			fig1.savefig('../Plots/Stability_check.pdf')
		else:
			plt.show()
	def plot_energy_alpha(self):
		" Function used to plot the energy and variance as a function of alpha for different omegas "
		# Plots for omega = 0.01
		self.read_data("Energy_Alpha_Mdistance_omega_0.01.txt")
		fig1 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig3 = plt.figure()
		plt.plot(self.alpha, self.Variance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\sigma_E$')
		plt.title(r'Plot of the variance $\sigma_E$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		### Plots for omega = 0.5
		self.read_data("Energy_Alpha_Mdistance_omega_0.50.txt")
		fig4 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig6 = plt.figure()
		plt.plot(self.alpha, self.Variance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\sigma_E$')
		plt.title(r'Plot of the variance $\sigma_E$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		### Plots for omega = 1
		self.read_data("Energy_Alpha_Mdistance_omega_1.00.txt")
		fig7 = plt.figure()
		plt.plot(self.alpha, self.Energy)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\langle H \rangle$')
		plt.title(r'Plot of the $\langle H \rangle$ as a function of $\alpha$ with $\omega$ = %.2f. Using $\psi_{T_1}$' %(self.omega))
		plt.legend(['Numerical value','Exact'])

		fig9 = plt.figure()
		plt.plot(self.alpha, self.Variance)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'$\sigma_E$')
		plt.title(r'Plot of the variance $\sigma_E$ as a function of $\alpha$ with $\omega$ = %.2f' %(self.omega))

		### Plots accepted configs
		fig10 = plt.figure()
		ax1 = plt.subplot(111)
		self.read_data("Energy_Alpha_Mdistance_omega_0.01.txt")
		ax1.plot(self.alpha, self.AcceptedConfigs, label=r'$\omega = %.2f$' %self.omega)
		plt.hold("on")
		self.read_data("Energy_Alpha_Mdistance_omega_0.50.txt")
		ax1.plot(self.alpha, self.AcceptedConfigs, label=r'$\omega = %.2f$' %self.omega)
		self.read_data("Energy_Alpha_Mdistance_omega_1.00.txt")
		ax1.plot(self.alpha, self.AcceptedConfigs, label=r'$\omega = %.2f$' %self.omega)
		plt.xlabel(r'$\alpha$')
		plt.ylabel(r'Accepted configurations (in %)')
		plt.title(r'Number of accepted configurations (in %) as a function of $\alpha$')
		ax1.legend(loc='upper right', bbox_to_anchor=(1,1), ncol=1, fancybox=True)
		
		if self.savefile == True:
			fig1.savefig('../Plots/Energy_alpha_plot_omega001.pdf')
			fig3.savefig('../Plots/Variance_alpha_plot_omega001.pdf')

			fig4.savefig('../Plots/Energy_alpha_plot_omega05.pdf')
			fig6.savefig('../Plots/Variance_alpha_plot_omega05.pdf')

			fig7.savefig('../Plots/Energy_alpha_plot_omega1.pdf')
			fig9.savefig('../Plots/Variance_alpha_plot_omega1.pdf')

			fig10.savefig('../Plots/AcceptedConfigs.pdf')
		else:
			plt.show()
		

	def Find_Optimal_AlphaBeta(self):
		fig1 = plt.figure()
		self.read_data("Optimal_AlphaBeta_omega_1.00.txt")
		ax = fig1.add_subplot(111, projection='3d')
		x = np.linspace(np.amin(self.alpha), np.amax(self.alpha), self.N)
		y = np.linspace(np.amin(self.beta), np.amax(self.beta), self.N)
		Alpha,Beta = np.meshgrid(x,y)
		ax.plot_surface(Alpha, Beta, self.Energy)
		ax.set_xlabel(r'$\alpha$')
		ax.set_ylabel(r'$\beta$')
		ax.set_zlabel(r'$\langle E \rangle$')

		# Find optimal alpha and beta for the lowest energy
		EnergyMinimum = np.amin(self.Energy)
		MinIndex = np.where(self.Energy == EnergyMinimum)[0][0]
		MinAlpha = self.alpha[MinIndex]
		MinBeta = self.beta[MinIndex]
		StdEnergy = np.sqrt(self.Variance[MinIndex])
		print 'For omega = ', self.omega
		print 'Lowest energy = ', EnergyMinimum
		print 'Optimal Alpha value = ', MinAlpha
		print 'Optimal Beta value = ', MinBeta
		print 'Standard deviation = ', StdEnergy
		print ''

		self.read_data("Optimal_AlphaBeta_omega_0.50.txt")
		EnergyMinimum = np.amin(self.Energy)
		MinIndex = np.where(self.Energy == EnergyMinimum)[0][0]
		MinAlpha = self.alpha[MinIndex]
		MinBeta = self.beta[MinIndex]
		StdEnergy = np.sqrt(self.Variance[MinIndex])
		print 'For omega = ', self.omega
		print 'Lowest energy = ', EnergyMinimum
		print 'Optimal Alpha value = ', MinAlpha
		print 'Optimal Beta value = ', MinBeta
		print 'Standard deviation = ', StdEnergy
		print ''

		self.read_data("Optimal_AlphaBeta_omega_0.01.txt")
		EnergyMinimum = np.amin(self.Energy)
		MinIndex = np.where(self.Energy == EnergyMinimum)[0][0]
		MinAlpha = self.alpha[MinIndex]
		MinBeta = self.beta[MinIndex]
		StdEnergy = np.sqrt(self.Variance[MinIndex])
		print 'For omega = ', self.omega
		print 'Lowest energy = ', EnergyMinimum
		print 'Optimal Alpha value = ', MinAlpha
		print 'Optimal Beta value = ', MinBeta
		print 'Standard deviation = ', StdEnergy

		if self.savefile == True:
			fig1.savefig('../Plots/OptimalAlphaBeta_3DPlot.pdf')
		else:
			plt.show()

	def Virial_plotting(self):
		""" Function made to plot for the virial problem """
		self.read_data_virial("Virial_data_V0.txt")
		fig1 = plt.figure()
		ax = plt.subplot(111)
		ax.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$\alpha = %.2f, \beta = %.2f$'%(self.alpha, self.beta))
		plt.hold("on")
		self.read_data_virial("Virial_data_V1.txt")
		ax.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$\alpha = %.2f, \beta = %.2f$'%(self.alpha, self.beta))
		self.read_data_virial("Virial_data_V2.txt")
		ax.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$\alpha = %.2f, \beta = %.2f$'%(self.alpha, self.beta))
		plt.xlabel(r'$\omega$')
		plt.ylabel(r'$\langle T \rangle/\langle V \rangle$')
		plt.title(r'Plot of $\langle T \rangle/\langle V \rangle$ as a function of $\omega$')
		ax.legend(loc='upper center', bbox_to_anchor=(0.6,0.5), ncol=1, fancybox=True)

		fig2 = plt.figure()
		ax1 = plt.subplot(111)
		self.read_data_virial("Virial_NoCoulomb_data_V0.txt")
		ax1.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$\alpha = %.2f, \beta = %.2f$'%(self.alpha, self.beta))
		self.read_data_virial("Virial_NoCoulomb_data_V1.txt")
		ax1.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$\alpha = %.2f, \beta = %.2f$'%(self.alpha, self.beta))
		self.read_data_virial("Virial_NoCoulomb_data_V2.txt")
		ax1.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$\alpha = %.2f, \beta = %.2f$'%(self.alpha, self.beta))
		plt.xlabel(r'$\omega$')
		plt.ylabel(r'$\langle T \rangle/\langle V \rangle$')
		plt.title(r'Plot of $\langle T \rangle/\langle V \rangle$ as a function of $\omega$. No Coulomb interaction')
		ax1.legend(loc='upper center', bbox_to_anchor=(0.5,0.5), ncol=1, fancybox=True)

		fig3 = plt.figure()
		ax2 = plt.subplot(111)
		plt.hold("on")
		self.read_data_virial("Virial_stability_MC_10000.txt")
		ax2.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$N_{MC} = 10000$')
		self.read_data_virial("Virial_stability_MC_100000.txt")
		ax2.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$N_{MC} = 100000$')
		self.read_data_virial("Virial_stability_MC_1000000.txt")
		ax2.plot(self.omega, self.KineticExpect/self.PotentialExpect, label=r'$N_{MC} = 1000000$')
		plt.xlabel(r'$\omega$')
		plt.ylabel(r'$\langle T \rangle/\langle V \rangle$')
		plt.title(r'Plot of $\langle T \rangle/\langle V \rangle$ as a function of $\omega$'+'\n Running for different Monte Carlo cycles.')
		ax2.legend(loc='lower right', bbox_to_anchor=(1,0), ncol=1, fancybox=True)		
		
		if self.savefile == True:
			fig1.savefig('../Plots/Virial_Plot.pdf')
			fig2.savefig('../plots/Virial_Plot_NoCoulombInt.pdf')
			fig3.savefig('../Plots/Virial_stability_test.pdf')
		else:
			plt.show()
## Comment out the functions to plot what you want
solver = Plotter(True)
#solver.stability_test()
#solver.plot_energy_alpha()
#solver.Find_Optimal_AlphaBeta()
solver.Virial_plotting()