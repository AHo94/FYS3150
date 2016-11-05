import numpy as np
import matplotlib.pyplot as plt
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project4/build-Project4_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

class Plotter():
	def __init__(self, savefile):
		self.savefile = savefile	# If True, saves the plots to a file

	def read_data(self, filename_open):
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			data_set = line.split()
			data.append(data_set)
		filename.close()

		N = len(data)
		self.L = data[0][0]
		self.MC_cycles = data[0][1]
		self.T = np.zeros(N)
		self.C_v = np.zeros(N)
		self.Chi = np.zeros(N)
		self.MC_range = np.linspace(0,MC_cycles, MC_cycles)
		for i in range(0,N):
			self.T[i] = data[i][2]
			self.C_v = data[i][4]
			self.Chi = data[i][6]


	def plot_state(self):
		self.read_data("4c.txt")
