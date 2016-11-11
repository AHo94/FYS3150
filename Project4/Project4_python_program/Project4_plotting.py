import numpy as np
import matplotlib.pyplot as plt
import os
# Change directory based based on work place, i.e home or at UiO computer.
#file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'
file_directory = '/uio/hume/student-u70/aleh/H2016/FYS3150/FYS3150_Projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_GCC_64bit-Debug'
class Plotter():
	def __init__(self, savefile):
		self.savefile = savefile	# If True, saves the plots to a file

	def read_data_4c(self, filename_E, filename_M, Temp_check):
		""" Reads data, specific for task 4c """
		filename = open(os.path.join(file_directory, filename_E), 'r')
		i = 0
		data_E = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				data_E.append(data_set)
			i += 1
		filename.close()

		filename = open(os.path.join(file_directory, filename_M), 'r')
		i = 0
		data_M = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				data_M.append(data_set)
			i += 1
		filename.close()
		N = len(data_E)-1
		self.MC_max = int(float(data_E[0][0]))
		self.L = float(data_E[0][1])
		self.MC_cycles = np.linspace(1, self.MC_max, N)
		if Temp_check == 1:
			self.T1 = float(data_E[0][2])
			self.E_expectation_1 = np.zeros(N)
			self.M_expectation_1 = np.zeros(N)
			self.E_counter_1 = np.zeros(N)

			for i in range(0, N):
				self.E_expectation_1[i] = float(data_E[i+1][0])
				self.M_expectation_1[i] = float(data_M[i+1][0])
				self.E_counter_1[i] = float(data_E[i+1][1])
		
		elif Temp_check == 2:
			self.T2 = float(data_E[0][2])
			self.E_expectation_2 = np.zeros(N)
			self.M_expectation_2 = np.zeros(N)
			self.E_counter_2 = np.zeros(N)

			for i in range(0, N):
				self.E_expectation_2[i] = float(data_E[i+1][0])
				self.M_expectation_2[i] = float(data_M[i+1][0])
				self.E_counter_2[i] = float(data_E[i+1][1])
		else:
			raise ValueError('Temp_check not set properly')

	
	def read_data_4d(self, filename_E1, filename_E2):
		""" Reads data, specific for 4d """
		filename = open(os.path.join(file_directory, filename_E1), 'r')
		i = 0
		data_E1 = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				data_E1.append(data_set)
			i += 1
		filename.close()

		filename = open(os.path.join(file_directory, filename_E2), 'r')
		i = 0
		data_E2 = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				data_E2.append(data_set)
			i += 1
		filename.close()
		N = len(data_E1)-1
		self.E_values_T1 = np.zeros(N)
		self.E_values_T24 = np.zeros(N)

		for i in range(0, N):
			self.E_values_T1[i] = float(data_E1[i+1][0])
			self.E_values_T24[i] = float(data_E2[i+1][0])
		

	def plot_state(self):
		""" Function that plots all plots in task 4c """
		# Plots the expecation values for T = 1
		self.read_data_4c("Mean_E_T1.00.txt", "Mean_M_T1.00.txt", 1)
		fig1 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T1, self.MC_max))
		fig2 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |M| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T1, self.MC_max))

		# Plots the expecation values for T = 2.4
		self.read_data_4c("Mean_E_T2.40.txt", "Mean_M_T2.40.txt", 2)
		fig3 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T2, self.MC_max))
		fig4 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |M| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T2, self.MC_max))

		# Plots the expecation values for T = 1, now with all states pointing up
		self.read_data_4c("Mean_E_AllUpState_T1.00.txt", "Mean_M_AllUpState_T1.00.txt", 1)
		fig1 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up' \
							 %(self.T1, self.MC_max))
		fig2 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |M| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T1, self.MC_max))

		# Plots the expecation values for T = 2.4, now with all states pointing up
		self.read_data_4c("Mean_E_AllUpState_T2.40.txt", "Mean_M_AllUpState_T2.40.txt", 2)
		fig3 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T2, self.MC_max))
		fig4 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |M| \rangle$')
		plt.title('Plot of magnetization as a function of MC Cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T2, self.MC_max))



		# Plots the accepted configurations as a function of MC cycles
		fig5 = plt.figure()
		plt.semilogy(self.MC_cycles, self.E_counter_1, 'b-')
		plt.hold("on")
		plt.semilogy(self.MC_cycles, self.E_counter_2, 'r-')
		plt.title('Number of accepted configurations as a function of Monte Carlo cycles')
		plt.xlabel('$N_{MC}$')
		plt.ylabel('Accepted spin flips')
		plt.legend(['T=1.0','T=2.40'])

		# Plots accepted configurations as a function of temperature
		fig6 = plt.figure()
		plt.semilogy(np.array([self.T1, self.T2]), np.array([self.E_counter_1[-1], self.E_counter_2[-1]]))
		plt.xlabel('Temperature - T')
		plt.ylabel('Accepted spin flips - Logscale')
		plt.title('Number of accepted spin configurations as a function of temperature')
		if self.savefile == True:
			fig1.savefig('../Plots/Energy_stability_T1.pdf')
			fig2.savefig('../Plots/Magnetization_stability_T1.pdf')
			fig3.savefig('../Plots/Energy_stability_T24.pdf')
			fig4.savefig('../Plots/Magnetization_stability_T24.pdf')
			fig5.savefig('../Plots/Accepted_configurations_wrt_MC_cycles.pdf')
			fig6.savefig('../Plots/Accepted_configs_wrt_temp.pdf')
		else:
			plt.show()
	
		
	def plot_probability(self):
		""" Function that plots the probability distribution """
		self.read_data_4d("Mean_E_T1.00.txt", "Mean_E_T2.40.txt")
		fig1 = plt.figure()
		plt.hist(self.E_values_T1, bins=100)
		plt.hold("on")
		plt.hist(self.E_values_T24, bins=100)
		plt.xlabel('Energies')
		plt.ylabel('Number of times energy is calculated')
		plt.legend(['T = 1.0','T = 2.4'])
		if self.savefile == True:
			fig1.savefig('../Plots/Probability_distribution_T1.pdf')
			fig2.savefig('../Plots/Probability_distribution_T2.pdf')
		else:
			plt.show()

solver = Plotter(False)
solver.plot_state()
#solver.plot_probability()
