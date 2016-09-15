import matplotlib.pyplot as plt
import numpy as np
import os

file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project1/build-Project1_cpp-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

general_points = np.array([10, 100, 1000])
simplified_points = np.array([10, 100, 1000, 10000, 100000, 1000000])
filename_general = []
filename_simplified = []
for general_n_pts in general_points:
	filename_general.append('General_data_n'+str(general_n_pts)+'.txt')

for simplified_n_pts in simplified_points:
	filename_simplified.append('Simplified_data_n'+str(simplified_n_pts)+'.txt')

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
	step = int(data[0][1])
	x = np.zeros(n)
	v = np.zeros(n)
	if step <= 1:
		maxpoints = n
	elif step > 1:
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

def relative_error(filename):
	""" Function that solves and plots the relative error """
	i = 0
	data = []
	for line in filename:
		n_value = line.split()
		if i != 1:
			data.append(n_value)
		i += 1
	filename.close()
	n = int(data[0][-1])
	x = np.zeros(n)
	v = np.zeros(n)

	for j in range(0,n,1):
		x[j] = float(data[j+1][0])
		v[j] = float(data[j+1][1])

	u_exact = 1-(1-np.exp(-10))*x - np.exp(-10*x)

for file_gen in filename_general:
	plt.figure()
	save_and_plot(file_gen, False)

for file_simpl in filename_simplified:
	plt.figure()
	save_and_plot(file_simpl, True)

plt.show()
"""
save_and_plot(fn10)
plt.figure()
save_and_plot(fn100)
plt.figure()
save_and_plot(fn1000)	
plt.figure()
save_and_plot(f_simplified)
plt.figure()
#relative_error(f_rel_error)
plt.show()
"""