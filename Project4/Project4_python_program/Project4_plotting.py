import numpy as np
import matplotlib.pyplot as plt
import os
#file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'
file_directory = '/uio/hume/student-u70/aleh/H2016/FYS3150/FYS3150_Projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_GCC_64bit-Debug'
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
		self.E_expectation_T1 = np.zeros(N_half)
		self.E_expectation_T24 = np.zeros(N_half)
		self.M_expectation_T1 = np.zeros(N_half)
		self.M_expectation_T24 = np.zeros(N_half)

		for i in range(0,N_half):
			self.MC_cycles[i] = float(data[i][1])
			self.E_expectation_T1[i] = float(data[i][4])
			self.E_expectation_T24[i] = float(data[i+N_half][4])
			self.M_expectation_T1[i] = float(data[i][6])
			self.M_expectation_T24[i] = float(data[i+N_half][6])

	def read_data_4c(self, filename_E, filename_M, Temp_check):
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
		MC_max = int(float(data_E[0][0]))
		self.L = float(data_E[0][1])
		self.MC_cycles = np.linspace(1, MC_max, N)
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

	def plot_state(self):
		# Plots the expecation values for T = 1
		self.read_data_4c("Mean_E_T1.00.txt", "Mean_M_T1.00.txt", 1)
		plt.plot(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of monte carlo cycles. T = %.2f' %(self.T1))
		plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |M| \rangle$')
		plt.title('Plot of magnetization as a function of monte carlo cycles. T = %.2f' %(self.T1))
		# Plots the expecation values for T = 2.4
		self.read_data_4c("Mean_E_T2.40.txt", "Mean_M_T2.40.txt", 2)
		plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of monte carlo cycles. T = %.2f' %(self.T2))
		plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |M| \rangle$')
		plt.title('Plot of magnetization as a function of monte carlo cycles. T = %.2f' %(self.T2))
		# Plots the accepted configurations as a function of MC cycles
		plt.figure()
		plt.plot(self.MC_cycles, self.E_counter_1, 'b-')
		plt.hold("on")
		plt.plot(self.MC_cycles, self.E_counter_2, 'r-')
		plt.title('Number of accepted spin flips as a function of Monte Carlo cycles')
		plt.xlabel('$N_{MC}$')
		plt.ylabel('Accepted spin flips')
		plt.legend(['T=1.0','T=2.40'])
		# Plots accepted configurations as a function of temperature
		plt.figure()
		plt.plot(np.array([self.T1, self.T2]), np.array([self.E_counter_1[-1], self.E_counter_2[-1]]))
		plt.show()

		
	def plot_probability(self):
		self.read_data_4c("Mean_E_T1.00.txt", "Mean_M_T1.00.txt")
		
		weights = np.ones_like(self.E_mean)/len(self.E_mean)
		plt.hist(self.E_mean, weights=weights, bins = 100)
		plt.show()

solver = Plotter(False)
solver.plot_state()
#solver.plot_probability()
