import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project2/build-Project2_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

# Defines lists that contains the number of points used
# Used for filename convention
n = 0
omega_value = np.array(['0.01', '0.50', '1.00', '5.00'])
omega_filenames = []
# For loops that creates filenames that is used for datafiles
for omega_points in omega_value:
	omega_filenames.append('Eigenvector_data_omega_'+omega_points+'.txt')

def save_and_plot(filename_open):
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
	n = int(data[0][0])
	rho_max = float(data[0][1])
	omega = float(data[0][2])
	rho = np.zeros(n+2)
	rho[-1] = rho_max
	eig = np.zeros(n+2)
	for j in range(0, n, 1):
		rho[j+1] = float(data[j+1][0])
		eig[j+1] = float(data[j+1][1])
	return rho, eig

rho, o1_eig1, omega1 = save_and_plot(omega_filenames[0])
rho, o2_eig1, omega2 = save_and_plot(omega_filenames[1])
rho, o3_eig1, omega3 = save_and_plot(omega_filenames[2])
rho, o4_eig1, omega4 = save_and_plot(omega_filenames[3])

plotting = False

fig1 = plt.figure()
plt.plot(rho, o1_eig1, 'b-')
plt.hold("on")
plt.plot(rho, o2_eig1, 'r-')
plt.plot(rho, o3_eig1, 'g-')
plt.plot(rho, o4_eig1, 'k-')
plt.xlabel(r'$\rho$')
plt.ylabel('$|\psi|^2$')
plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$, ground state')
plt.legend(['$\omega=0.01$','$\omega=0.5$','$\omega=1$', '$\omega=5$'])
plt.hold("off")
"""
fig2 = plt.figure()
plt.plot(rho, o1_eig2, 'b-')
plt.hold("on")
plt.plot(rho, o2_eig2, 'r-')
plt.plot(rho, o3_eig2, 'g-')
plt.plot(rho, o4_eig2, 'k-')
plt.xlabel(r'$\rho$')
plt.ylabel('$|\psi|^2$')
plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$, first excited state')
plt.legend(['$\omega=0.01$','$\omega=0.5$','$\omega=1$', '$\omega=5$'])
plt.hold("off")

fig3 = plt.figure()
plt.plot(rho, o1_eig3, 'b-')
plt.hold("on")
plt.plot(rho, o2_eig3, 'r-')
plt.plot(rho, o3_eig3, 'g-')
plt.plot(rho, o4_eig3, 'k-')
plt.xlabel(r'$\rho$')
plt.ylabel('$|\psi|^2$')
plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$, second excited state')
plt.legend(['$\omega=0.01$','$\omega=0.5$','$\omega=1$', '$\omega=5$'])
plt.hold("off")
"""
"""
if plotting:
	plt.show()
else:
	fig1.savefig('../Plots/Plot_groundstate.pdf')
	#fig2.savefig('../Plots/Plot_first_excitedstate.pdf')
	#fig3.savefig('../Plots/Plot_second_excitedstate.pdf')
"""

class Proj2_plot_script:
	def __init__(self):
		self.n = 0

	def read_data(self, fielane_open):
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
		rho = np.zeros(n+2)
		rho[-1] = rho_max
		eig1 = np.zeros(n+2)
		eig2 = np.zeros(n+2)
		eig3 = np.zeros(n+2)
		for j in range(0, n, 1):
			rho[j+1] = float(data[j+1][0])
			eig1[j+1] = float(data[j+1][1])