import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project3/build-Project3_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

class Plotter():
	def __init__(self):
		self.filename = "Celestial_positions.txt"
		self.read_data(self.filename)

	def read_data(self, filename_open):
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			data_set = line.split()
			data.append(data_set)
		filename.close()

		self.N = len(data)
		NumberofObjects = (len(data[0])/3.0)
		self.Sun_pos = np.zeros((3,self.N))
		self.Earth_pos = np.zeros((3,self.N))
		#if NumberofObjects <= 3:
		self.Jupiter_pos = np.zeros((3, self.N))

		for i in range(0,3):
			for j in range(0,self.N):
				self.Sun_pos[i][j] = data[j][i]
				self.Earth_pos[i][j] = data[j][i+3]
				self.Jupiter_pos[i][j] = data[j][i+6]



	def plotting(self):
		plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
		plt.hold("on")
		plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
		plt.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:])
		plt.xlabel('x - [AU]')
		plt.ylabel('y - [AU]')
		plt.legend(['Earth', 'Sun', 'Jupiter'])
		plt.show()

	def animation(self):
		# Animates 
		a = 3

solve = Plotter()
solve.plotting()