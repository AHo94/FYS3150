import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project3/build-Project3_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

class Plotter():
	def __init__(self, savefile):
		self.filename = "Celestial_positions.txt"
		self.read_data(self.filename)
		self.savefile = savefile

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
		self.Jupiter_pos = np.zeros((3, self.N))
		self.Mercury_pos = np.zeros((3, self.N))
		self.Mars_pos = np.zeros((3, self.N))
		self.Saturn_pos = np.zeros((3, self.N))
		self.Uranus_pos = np.zeros((3, self.N))
		self.Neptune_pos = np.zeros((3, self.N))
		self.Pluto_pos = np.zeros((3,self.N))

		for i in range(0,3):
			for j in range(0,self.N):
				self.Sun_pos[i][j] = data[j][i]
				self.Earth_pos[i][j] = data[j][i+3]
				self.Jupiter_pos[i][j] = data[j][i+6]
				self.Mercury_pos[i][j] = data[j][i+9]
				self.Mars_pos[i][j] = data[j][i+12]
				self.Saturn_pos[i][j] = data[j][i+15]
				self.Uranus_pos[i][j] = data[j][i+18]
				self.Neptune_pos[i][j] = data[j][i+21]
				self.Pluto_pos[i][j] = data[j][i+24]

	def Earth_Jupiter_test(self):	
		fig = plt.figure()
		if self.savefile:
			fig.savefig('../Plots/'+'Earth_Jupiter_test.pdf')
		else:		
			plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
			plt.hold("on")
			plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
			plt.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:])
			plt.xlabel('x - [AU]')
			plt.ylabel('y - [AU]')
			plt.legend(['Earth', 'Sun', 'Jupiter'])
			plt.title('Orbits for a system with the Sun, Earth and Jupiter.')
			plt.show()

	def plotting_2D(self):
		fig = plt.figure()
		if self.savefile:
			fig.savefig('../Plots/'+'All_planets_plot.pdf')
		else:
			plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
			plt.hold("on")
			plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
			plt.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:])
			plt.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:])
			plt.plot(self.Mars_pos[0][:], self.Mars_pos[1][:])
			plt.plot(self.Saturn_pos[0][:], self.Saturn_pos[1][:])
			plt.plot(self.Uranus_pos[0][:], self.Uranus_pos[1][:])
			plt.plot(self.Neptune_pos[0][:], self.Neptune_pos[1][:])
			plt.plot(self.Pluto_pos[0][:], self.Pluto_pos[1][:])
			plt.xlabel('x - [AU]')
			plt.ylabel('y - [AU]')
			plt.legend(['Earth', 'Sun', 'Jupiter', 'Mercury', 'Mars', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
			plt.title('Orbits of the planets. Timescale in years.')
			plt.show()

	def animation(self):
		# Animates 
		a = 3

solve = Plotter(False)
solve.plotting_2D()