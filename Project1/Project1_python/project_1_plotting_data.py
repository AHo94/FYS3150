import matplotlib.pyplot as plt
import numpy as np
import os

file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project1/build-Project1_cpp-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

# Defines lists that contains the number of points used
# Used for filename convention
general_points = np.array([10, 100, 1000])
simplified_points = np.array([10, 100, 1000, 10000, 100000, 1000000])
error_points = np.array([10, 100, 1000, 10000, 100000, 1000000, 10000000])
filename_general = []
filename_simplified = []
filename_error = []

# For loops that creates filenames that is used for datafiles
for general_n_pts in general_points:
	filename_general.append('General_data_n'+str(general_n_pts)+'.txt')

for simplified_n_pts in simplified_points:
	filename_simplified.append('Simplified_data_n'+str(simplified_n_pts)+'.txt')

for error_n_pts in error_points:
	filename_error.append('Error_data_n'+str(error_n_pts)+'.txt')

def save_and_plot(filename_open, simplified_plot):
	""" 
	Function that reads the datafile, saves the values of n, x and v
	Also solves the exact solution u(x) for the given values of x
	Plots both v and u as a function of x
	"""
	filename = open(os.path.join(file_directory, filename_open), 'r')
	i = 0
	data = []
	for line in filename:
		n_value = line.split()
		if i != 0:
			data.append(n_value)
		i += 1
	filename.close()
	n = int(data[0][0])
	x = np.zeros(n)
	v = np.zeros(n)
	if n <= 1000:
		maxpoints = n
	elif n > 1000:
		maxpoints = 1000

	for j in range(0,maxpoints,1):
		x[j] = float(data[j+1][0])
		v[j] = float(data[j+1][1])

	u_exact = 1-(1-np.exp(-10))*x - np.exp(-10*x)

	plt.plot(x,v,'b-')
	plt.hold("on")
	plt.plot(x,u_exact,'r-')
	plt.xlabel('x')
	plt.ylabel('u')
	plt.legend(['Numerical solution', 'Exact solution'])
	if simplified_plot == True:
		plt.title('Plot of $v$ and $u$, special algorithm, , n = %g' %(n))
	else:
		plt.title('Plot of $v$ and $u$, general algorithm, , n = %g' %(n))

def relative_error(error_list):
	""" Function that solves and plots the relative error """
	#n_values = zeros(len(error_list))
	h_values = np.zeros(len(error_list))
	rel_error = np.zeros(len(error_list))
	h_count = 0
	for file_error in error_list:
		filename = open(os.path.join(file_directory, file_error), 'r')
		i = 0
		data = []
		for line in filename:
			n_value = line.split()
			if i != 0:
				data.append(n_value)
			i += 1
		filename.close()
		n = int(data[0][0])
		v = np.zeros(n)
		x = np.zeros(n)
		h = 1.0/(n-1)

		if n <= 1000:
			maxpoints = n
		elif n > 1000:
			maxpoints = 1000

		for j in range(0,maxpoints,1):
			x[j] = float(data[j+1][0])
			v[j] = float(data[j+1][1])

		u_exact = 1-(1-np.exp(-10))*x - np.exp(-10*x)

		h_values[h_count] = np.log(h)/np.log(10)
		rel_error[h_count] = np.max(np.log(np.abs((v[1:999]-u_exact[1:999])/u_exact[1:999]))/np.log(10))
		h_count += 1

	
	plt.figure()
	plt.plot(h_values, rel_error, 'g-o')
	plt.xlabel('$log_{10}(h)$')
	plt.ylabel('Relative error, in $log_{10}$ scale')
	plt.title('Relative error')
# Plots the figures
"""
for file_gen in filename_general:
	plt.figure()
	save_and_plot(file_gen, False)

for file_simpl in filename_simplified:
	plt.figure()
	save_and_plot(file_simpl, True)
"""
relative_error(filename_error)	
plt.show()
