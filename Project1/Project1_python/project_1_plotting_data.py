import matplotlib.pyplot as plt
import numpy as np
import os

file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project1/build-Project1_cpp-Desktop_Qt_5_7_0_MinGW_32bit-Debug'
filename_n10 = 'Project1_data_n10.txt'
filename_n100 = 'Project1_data_n100.txt'
filename_n1000 = 'Project1_data_n1000.txt'
filename_simplified = 'Project1c_data_simplified.txt'
filename_relative_error = 'Project1d_relative_error.txt'
fn10 = open(os.path.join(file_directory, filename_n10), 'r')
fn100 = open(os.path.join(file_directory, filename_n100), 'r')
fn1000 = open(os.path.join(file_directory, filename_n1000), 'r')
f_simplified = open(os.path.join(file_directory, filename_simplified), 'r')
f_rel_error = open(os.path.join(file_directory, filename_relative_error), 'r')

def save_and_plot(filename):
	""" 
	Function that reads the datafile, saves the values of n, x and v
	Also solves the exact solution u(x) for the given values of x
	Plots both v and u as a function of x
	"""
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

	plt.plot(x,v,'b-')
	plt.hold("on")
	plt.plot(x,u_exact,'r-')
	plt.xlabel('x')
	plt.ylabel('u')
	plt.legend(['Numerical solution', 'Exact solution'])
	if filename == f_simplified:
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

save_and_plot(fn10)
plt.figure()
save_and_plot(fn100)
plt.figure()
save_and_plot(fn1000)	
plt.figure()
save_and_plot(f_simplified)
plt.figure()
relative_error(f_rel_error)
plt.show()