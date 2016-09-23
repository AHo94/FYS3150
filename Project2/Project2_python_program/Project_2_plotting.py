import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project2/build-Project2_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

# Defines lists that contains the number of points used
# Used for filename convention
omega_value = np.array([0.01, 0.50, 1.00, 5.00])
omega_filenames = []
# For loops that creates filenames that is used for datafiles
for omega_points in omega_value:
	omega_filenames.append('Eigenvector_data_omega_'+str(omega_points)+'.txt')



def save_and_plot(filename_open):
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
	eig1 = np.zeros(n+2)
	eig2 = np.zeros(n+2)
	eig3 = np.zeros(n+2)
	for j in range(0, n, 1):
		rho[j+1] = float(data[j+1][0])
		eig1[j+1] = float(data[j+1][1])
		eig2[j+1] = float(data[j+1][2])
		eig3[j+1] = float(data[j+1][3])

	plt.plot(rho, eig1, 'b-')
	plt.xlabel(r'$\rho$')
	plt.ylabel('$|\psi|^2$')
	plt.title(r'Plot of $|\psi|^2$ as a function of $\rho$, with $\omega$ = %.g' %(omega))

save_and_plot(omega_filenames[0])
"""
fig_err = plt.figure()
relative_error(filename_error)
fig_err.savefig('../Data_plots/Relative_error.png')
"""
plt.show()
