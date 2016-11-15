import numpy as np
import matplotlib.pyplot as plt
import os
# Change directory based based on work place, i.e home or at UiO computer.
#file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'
#file_directory = '/uio/hume/student-u70/aleh/H2016/FYS3150/FYS3150_Projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_GCC_64bit-Debug'
file_directory = '../build-Project4_cpp_program-Desktop_Qt_5_7_0_GCC_64bit-Debug'

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

			for j in range(0, N):
				self.E_expectation_1[j] = float(data_E[j+1][0])
				self.M_expectation_1[j] = float(data_M[j+1][0])
				self.E_counter_1[j] = float(data_E[j+1][1])
		
		elif Temp_check == 2:
			self.T2 = float(data_E[0][2])
			self.E_expectation_2 = np.zeros(N)
			self.M_expectation_2 = np.zeros(N)
			self.E_counter_2 = np.zeros(N)

			for j in range(0, N):
				self.E_expectation_2[j] = float(data_E[j+1][0])
				self.M_expectation_2[j] = float(data_M[j+1][0])
				self.E_counter_2[j] = float(data_E[j+1][1])
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
		self.E_variance_T1 = float(data_E1[0][3])
		self.E_variance_T24 = float(data_E2[0][3])
		for j in range(0, N):
			self.E_values_T1[j] = float(data_E1[j+1][0])
			self.E_values_T24[j] = float(data_E2[j+1][0])

	def read_data_parallellization(self, fileL40, fileL60, fileL100, fileL140):
		""" Function specifically made to read data for the parallelization part """
		filename = open(os.path.join(file_directory, fileL40), 'r')
		i = 0
		dataL40 = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				dataL40.append(data_set)
			i += 1
		filename.close()
		filename = open(os.path.join(file_directory, fileL60), 'r')
		i = 0
		dataL60 = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				dataL60.append(data_set)
			i += 1
		filename.close()
		filename = open(os.path.join(file_directory, fileL100), 'r')
		i = 0
		dataL100 = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				dataL100.append(data_set)
			i += 1
		filename.close()
		filename = open(os.path.join(file_directory, fileL140), 'r')
		i = 0
		dataL140 = []
		for line in filename:
			if i != 0:
				data_set = line.split()
				dataL140.append(data_set)
			i += 1
		filename.close()

		N = len(dataL40)
		self.MC_max_parallell = int(dataL40[0][0])
		self.T_parallell = np.zeros(N)

		self.E_L40 = np.zeros(N)
		self.E_L60 = np.zeros(N)
		self.E_L100 = np.zeros(N)
		self.E_L140 = np.zeros(N)

		self.M_L40 = np.zeros(N)
		self.M_L60 = np.zeros(N)
		self.M_L100 = np.zeros(N)
		self.M_L140 = np.zeros(N)

		self.C_vL40 = np.zeros(N)
		self.C_vL60 = np.zeros(N)
		self.C_vL100 = np.zeros(N)
		self.C_vL140 = np.zeros(N)

		self.ChiL40 = np.zeros(N)
		self.ChiL60 = np.zeros(N)
		self.ChiL100 = np.zeros(N)
		self.ChiL140 = np.zeros(N)
		for j in range(0, N):
			self.T_parallell[j] = float(dataL40[j][1]) 

			self.E_L40[j] = float(dataL40[j][3])
			self.E_L60[j] = float(dataL60[j][3])
			self.E_L100[j] = float(dataL100[j][3])
			self.E_L140[j] = float(dataL140[j][3])

			self.M_L40[j] = float(dataL40[j][4])
			self.M_L60[j] = float(dataL60[j][4])
			self.M_L100[j] = float(dataL100[j][4])
			self.M_L140[j] = float(dataL140[j][4])

			self.C_vL40[j] = float(dataL40[j][5])
			self.C_vL60[j] = float(dataL60[j][5])
			self.C_vL100[j] = float(dataL100[j][5])
			self.C_vL140[j] = float(dataL140[j][5])

			self.ChiL40[j] = float(dataL40[j][6])
			self.ChiL60[j] = float(dataL60[j][6])
			self.ChiL100[j] = float(dataL100[j][6])
			self.ChiL140[j] = float(dataL140[j][6])

	def plot_state(self):
		""" Function that plots all plots in task 4c """
		# Plots the expecation values for T = 1
		self.read_data_4c("E_expect_T1.00.txt", "M_expect_T1.00.txt", 1)
		fig1 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T1, self.MC_max))
		fig2 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T1, self.MC_max))

		# Plots the expecation values for T = 2.4
		self.read_data_4c("E_expect_T2.40.txt", "M_expect_T2.40.txt", 2)
		fig3 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T2, self.MC_max))
		fig4 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T2, self.MC_max))

		# Plots the expecation values for T = 1, now with all states pointing up
		self.read_data_4c("E_expect_AllUpState_T1.00.txt", "M_expect_AllUpState_T1.00.txt", 1)
		fig5 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up' \
							 %(self.T1, self.MC_max))
		fig6 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T1, self.MC_max))

		# Plots the expecation values for T = 2.4, now with all states pointing up
		self.read_data_4c("E_expect_AllUpState_T2.40.txt", "M_expect_AllUpState_T2.40.txt", 2)
		fig7 = plt.figure()
		plt.plot(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T2, self.MC_max))
		fig8 = plt.figure()
		plt.plot(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC Cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T2, self.MC_max))

		# Plots the accepted configurations as a function of MC cycles
		fig9 = plt.figure()
		plt.semilogy(self.MC_cycles, self.E_counter_1, 'b-')
		plt.hold("on")
		plt.semilogy(self.MC_cycles, self.E_counter_2, 'r-')
		plt.title('Number of accepted configurations as a function of Monte Carlo cycles')
		plt.xlabel('$N_{MC}$')
		plt.ylabel('Accepted spin flips')
		plt.legend(['T=1.0','T=2.40'])

		# Plots accepted configurations as a function of temperature
		fig10 = plt.figure()
		plt.semilogy(np.array([self.T1, self.T2]), np.array([self.E_counter_1[-1], self.E_counter_2[-1]]))
		plt.xlabel('Temperature - T')
		plt.ylabel('Accepted spin flips - Logscale')
		plt.title('Number of accepted spin configurations as a function of temperature')
		if self.savefile == True:
			fig1.savefig('../Plots/Energy_stability_T1.pdf')
			fig2.savefig('../Plots/Magnetization_stability_T1.pdf')
			fig3.savefig('../Plots/Energy_stability_T24.pdf')
			fig4.savefig('../Plots/Magnetization_stability_T24.pdf')

			fig5.savefig('../Plots/Energy_stability_UpInitSpin_T1.pdf')
			fig6.savefig('../Plots/Magnetization_stability_UpInitSpin_T1.pdf')
			fig7.savefig('../Plots/Energy_stability_UpInitSpin_T24.pdf')
			fig8.savefig('../Plots/Magnetization_stability_UpInitSpin_T24.pdf')
			
			fig9.savefig('../Plots/Accepted_configurations_wrt_MC_cycles.pdf')
			fig10.savefig('../Plots/Accepted_configs_wrt_temp.pdf')
		else:
			plt.show()
	
		
	def plot_probability(self):
		""" Function that plots the probability distribution """
		self.read_data_4d("E_expect_T1.00.txt", "E_expect_T2.40.txt")
		fig1 = plt.figure()
		plt.hist(self.E_values_T1, bins=100)
		plt.hold("on")
		plt.hist(self.E_values_T24, bins=100)
		plt.xlabel('Energies')
		plt.ylabel('Number of times energy is calculated')
		plt.legend(['T = 1.0','T = 2.4'])
		plt.title('Probability distribution of the energies')
		
		# Splits the probability distribution in two separate plots
		fig2 = plt.figure()
		plt.hist(self.E_values_T1, bins=100)
		plt.xlabel('Energies')
		plt.ylabel('Number of times energy is calculated')
		plt.title('Probability distribution of the energy, for $T = 1.0$')

		fig3 = plt.figure()
		plt.hist(self.E_values_T24, bins=100)
		plt.xlabel('Energies')
		plt.ylabel('Number of times energy is calculated')
		plt.title('Probability distribution of the energy, for $T=2.4$')

		print 'Computed variance = ', self.E_variance_T1, ', for T = 1.0'
		print 'Computed variance = ', self.E_variance_T24, ', for T = 2.4'
		print 'Numpy variance = ', np.var(self.E_values_T1), ',for T = 1.0'
		print 'Numpy variance = ', np.var(self.E_values_T24), ',for T = 2.4'
		
		if self.savefile == True:
			fig1.savefig('../Plots/Probability_distribution_merged.pdf')
			fig2.savefig('../Plots/Probability_distribution_T1.pdf')
			fig3.savefig('../Plots/Probability_distribution_T2.pdf')
		else:
			plt.show()

	def plot_state_logarithmic(self):
		""" Plots the stability of E and M, with logarithmic x axis """
		self.read_data_4c("E_expect_T1.00.txt", "M_expect_T1.00.txt", 1)
		fig1 = plt.figure()
		plt.semilogx(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T1, self.MC_max))
		fig2 = plt.figure()
		plt.semilogx(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T1, self.MC_max))

		# Plots the expecation values for T = 2.4
		self.read_data_4c("E_expect_T2.40.txt", "M_expect_T2.40.txt", 2)
		fig3 = plt.figure()
		plt.semilogx(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T2, self.MC_max))
		fig4 = plt.figure()
		plt.semilogx(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Abritary initial state' \
							 %(self.T2, self.MC_max))

		# Plots the expecation values for T = 1, now with all states pointing up
		self.read_data_4c("E_expect_AllUpState_T1.00.txt", "M_expect_AllUpState_T1.00.txt", 1)
		fig5 = plt.figure()
		plt.semilogx(self.MC_cycles, self.E_expectation_1, 'b-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up' \
							 %(self.T1, self.MC_max))
		fig6 = plt.figure()
		plt.semilogx(self.MC_cycles, self.M_expectation_1, 'r-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T1, self.MC_max))

		# Plots the expecation values for T = 2.4, now with all states pointing up
		self.read_data_4c("E_expect_AllUpState_T2.40.txt", "M_expect_AllUpState_T2.40.txt", 2)
		fig7 = plt.figure()
		plt.semilogx(self.MC_cycles, self.E_expectation_2, 'b-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  E \rangle$')
		plt.title('Plot of the energies as a function of MC cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T2, self.MC_max))
		fig8 = plt.figure()
		plt.semilogx(self.MC_cycles, self.M_expectation_2, 'r-')
		plt.xlabel('$\log(N_{MC})$')
		plt.ylabel(r'$\langle  |\mathcal{M}| \rangle$')
		plt.title('Plot of magnetization as a function of MC Cycles. T = %.2f, $N_{mc}$ = %.g. \n Initial state, all spins up'\
							 %(self.T2, self.MC_max))

		if self.savefile == True:
			fig1.savefig('../Plots/Energy_stability_logarithmic_T1.pdf')
			fig2.savefig('../Plots/Magnetization_stability_logarithmic_T1.pdf')
			fig3.savefig('../Plots/Energy_stability_logarithmic_T24.pdf')
			fig4.savefig('../Plots/Magnetization_stability_logarithmic_T24.pdf')

			fig5.savefig('../Plots/Energy_stability_UpInitSpin_logarithmic_T1.pdf')
			fig6.savefig('../Plots/Magnetization_stability_UpInitSpin_logarithmic_T1.pdf')
			fig7.savefig('../Plots/Energy_stability_UpInitSpin_logarithmic_T24.pdf')
			fig8.savefig('../Plots/Magnetization_stability_UpInitSpin_logarithmic_T24.pdf')
			
		else:
			plt.show()

	def plot_parallellization(self, compute_TC):
		""" Function made specifically for the parallelization part of plotting """
		self.read_data_parallellization("4e_data_L40.txt","4e_data_L60.txt","4e_data_L100.txt","4e_data_L140.txt")
		# Plots energy
		fig1 = plt.figure()
		ax = plt.subplot(111)
		line, = ax.plot(self.T_parallell, self.E_L40, label='$L=40$')
		plt.hold("on")
		line, = ax.plot(self.T_parallell, self.E_L60, label='$L=60$')
		line, = ax.plot(self.T_parallell, self.E_L100, label='$L=100$')
		line, = ax.plot(self.T_parallell, self.E_L140, label='$L=140$')
		plt.xlabel('$T$')
		plt.ylabel(r'$\langle E \rangle/L^2$')
		plt.title(r'Plot of $\langle E \rangle/L^2$ as a function of $T$. $N_{mc} = %.g$' %(self.MC_max_parallell))
		ax.legend(loc='upper center', bbox_to_anchor=(0.5,1.00), ncol=2, fancybox=True)
		
		# Plots magnetization
		fig2 = plt.figure()
		ax = plt.subplot(111)
		line, = ax.plot(self.T_parallell, self.M_L40, label='$L=40$')
		plt.hold("on")
		line, = ax.plot(self.T_parallell, self.M_L60, label='$L=60$')
		line, = ax.plot(self.T_parallell, self.M_L100, label='$L=100$')
		line, = ax.plot(self.T_parallell, self.M_L140, label='$L=140$')
		plt.xlabel('$T$')
		plt.ylabel(r'$\langle |\mathcal{M}| \rangle/L^2$')
		plt.title(r'Plot of $\langle |\mathcal{M}| \rangle/L^2$ as a function of $T$. $N_{mc} = %.g$' %(self.MC_max_parallell))
		ax.legend(loc='upper center', bbox_to_anchor=(0.3,0.2), ncol=2, fancybox=True)
		
		# Plots non-interpolated heat capacity
		fig3 = plt.figure()
		ax = plt.subplot(111)
		line, = ax.plot(self.T_parallell, self.C_vL40, label='$L=40$')
		plt.hold("on")
		line, = ax.plot(self.T_parallell, self.C_vL60, label='$L=60$')
		line, = ax.plot(self.T_parallell, self.C_vL100, label='$L=100$')
		line, = ax.plot(self.T_parallell, self.C_vL140, label='$L=140$')
		plt.xlabel('$T$')
		plt.ylabel(r'$C_V/L^2$')
		plt.title(r'Plot of $C_V/L^2$ as a function of $T$. $N_{mc} = %.g$' %(self.MC_max_parallell))
		ax.legend(loc='upper center', bbox_to_anchor=(0.4,1.00), ncol=2, fancybox=True)
		
		# Plots non-interpolated suceptibility
		fig4 = plt.figure()
		ax = plt.subplot(111)
		line, = ax.plot(self.T_parallell, self.ChiL40, label='$L=40$')
		plt.hold("on")
		line, = ax.plot(self.T_parallell, self.ChiL60, label='$L=60$')
		line, = ax.plot(self.T_parallell, self.ChiL100, label='$L=100$')
		line, = ax.plot(self.T_parallell, self.ChiL140, label='$L=140$')
		plt.xlabel('$T$')
		plt.ylabel(r'$\chi/L^2$')
		plt.title(r'Plot of $\chi/L^2$ as a function of $T$. $N_{mc} = %.g$' %(self.MC_max_parallell))
		ax.legend(loc='upper center', bbox_to_anchor=(0.3,1.00), ncol=2, fancybox=True)
		
		# Plots interpolated heat capacity
		fig5 = plt.figure()
		ax = plt.subplot(111)
		# Calculate interpolation of the heat capacity
		xvals = np.linspace(2.1, 2.35, 70)
		C_v_interp_L40 = np.interp(xvals, self.T_parallell, self.C_vL40)
		C_v_interp_L60 = np.interp(xvals, self.T_parallell, self.C_vL60)
		C_v_interp_L100 = np.interp(xvals, self.T_parallell, self.C_vL100)
		C_v_interp_L140 = np.interp(xvals, self.T_parallell, self.C_vL140)

		line, = ax.plot(self.T_parallell, self.C_vL40, 'o', label='$L=40$')
		plt.hold("on")
		line, = ax.plot(self.T_parallell, self.C_vL60, 'o', label='$L=60$')
		line, = ax.plot(self.T_parallell, self.C_vL100, 'o', label='$L=100$')
		line, = ax.plot(self.T_parallell, self.C_vL140, 'o', label='$L=140$')
		# Plots interpolation
		line, = ax.plot(xvals, C_v_interp_L40, '-x', label='$L=40$')
		line, = ax.plot(xvals, C_v_interp_L60, '-x', label='$L=60$')
		line, = ax.plot(xvals, C_v_interp_L100, '-x', label='$L=100$')
		line, = ax.plot(xvals, C_v_interp_L140, '-x', label='$L=140$')
		plt.xlabel('$T$')
		plt.ylabel(r'$C_V/L^2$')
		plt.title(r'Plot of $C_V/L^2$ as a function of $T$. Interpolated points. $N_{mc} = %.g$' %(self.MC_max_parallell))
		ax.legend(loc='upper center', bbox_to_anchor=(0.3,1.00), ncol=2, fancybox=True)
		
		# Plots interpolated heat capacity
		fig6 = plt.figure()
		ax = plt.subplot(111)
		# Calculate interpolation of the susceptibility
		Chi_interp_L40 = np.interp(xvals, self.T_parallell, self.ChiL40)
		Chi_interp_L60 = np.interp(xvals, self.T_parallell, self.ChiL60)
		Chi_interp_L100 = np.interp(xvals, self.T_parallell, self.ChiL100)
		Chi_interp_L140 = np.interp(xvals, self.T_parallell, self.ChiL140)

		line, = ax.plot(self.T_parallell, self.ChiL40, 'o', label='$L=40$')
		plt.hold("on")
		line, = ax.plot(self.T_parallell, self.ChiL60, 'o', label='$L=60$')
		line, = ax.plot(self.T_parallell, self.ChiL100, 'o', label='$L=100$')
		line, = ax.plot(self.T_parallell, self.ChiL140, 'o', label='$L=140$')
		# Plots interpolation
		line, = ax.plot(xvals, Chi_interp_L40, '-x', label='$L=40$')
		line, = ax.plot(xvals, Chi_interp_L60, '-x', label='$L=60$')
		line, = ax.plot(xvals, Chi_interp_L100, '-x', label='$L=100$')
		line, = ax.plot(xvals, Chi_interp_L140, '-x', label='$L=140$')
		plt.xlabel('$T$')
		plt.ylabel(r'$\chi/L^2$')
		plt.title(r'Plot of $\chi/L^2$ as a function of $T$. Interpolated points. $N_{mc} = %.g$' %(self.MC_max_parallell))
		ax.legend(loc='upper center', bbox_to_anchor=(0.3,1.00), ncol=2, fancybox=True)

		if compute_TC == True:
			# Calculate the critical temperature for L -> infinity
			T_C = [2.28874, 2.28119, 2.27788, 2.7788]
			L = [40, 60, 100, 140]
			T_C_inf = 0
			norm = 0
			for i in range(0, 3):
				for j in range(i+1, 4):
					norm += 1
					T_C_inf += T_C[i] - (T_C[i]-T_C[j])/(1 - (L[i]/L[j]))

			T_C_inf = T_C_inf/(norm)
			print 'Computed critical temperature for L -> infinity:'
			print 'T_C(L = infty) = ', T_C_inf
		else:
			# Plots or saves the files.
			if self.savefile == True:
				fig1.savefig('../Plots/Energy_parallellization.pdf')
				fig2.savefig('../Plots/Magnetization_parallellization.pdf')
				fig3.savefig('../Plots/Heat_capacity_parallellization.pdf')
				fig4.savefig('../Plots/Susceptibility_parallellization.pdf')			
				fig5.savefig('../Plots/Heat_capacity_parallellization_interpolated.pdf')
				fig6.savefig('../Plots/Susceptibility_parallellization_interpolated.pdf')
			else:
				plt.show()

## Comment out the functions to plot what you want
solver = Plotter(True)
#solver.plot_state()
#solver.plot_probability()
#solver.plot_state_logarithmic()
solver.plot_parallellization(False)
