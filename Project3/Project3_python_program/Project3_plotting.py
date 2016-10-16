import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
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
			if i == 0:
				self.Nsteps = float(data_set[0])
				self.dt = float(data_set[1])
				i += 1
			else:
				data.append(data_set)

		filename.close()

		self.N = len(data)
		self.Nstep = int(self.Nsteps)

		self.Sun_pos = np.zeros((3,self.N))
		self.Earth_pos = np.zeros((3,self.N))
		self.Jupiter_pos = np.zeros((3, self.N))
		self.Mercury_pos = np.zeros((3, self.N))
		self.Venus_pos = np.zeros((3, self.N))
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
				self.Venus_pos[i][j] = data[j][i+12]
				self.Mars_pos[i][j] = data[j][i+15]
				self.Saturn_pos[i][j] = data[j][i+18]
				self.Uranus_pos[i][j] = data[j][i+21]
				self.Neptune_pos[i][j] = data[j][i+24]
				self.Pluto_pos[i][j] = data[j][i+27]

	def read_data_ESJ_sys(self, filename_open, Add_jupiter):
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			data_set = line.split()
			if i == 0:
				self.Nsteps = float(data_set[0])
				self.dt = float(data_set[1])
				i += 1
			else:
				data.append(data_set)
		filename.close()
		self.N = len(data)
		self.Nstep = int(self.Nsteps)
		self.Sun_pos = np.zeros((3,self.N))
		self.Earth_pos = np.zeros((3,self.N))

		if Add_jupiter == True:
			self.Jupiter_pos = np.zeros((3, self.N))
			for i in range(0,3):
				for j in range(0,self.N):
					self.Sun_pos[i][j] = data[j][i]
					self.Earth_pos[i][j] = data[j][i+3]
					self.Jupiter_pos[i][j] = data[j][i+6]
		else:
			for i in range(0,3):
				for j in range(0,self.N):
					self.Sun_pos[i][j] = data[j][i]
					self.Earth_pos[i][j] = data[j][i+3]


	def read_data_mercury_GR(self, filename_open):
		filename = open(os.path.join(file_directory, filename_open), 'r')
		i = 0
		data = []
		for line in filename:
			data_set = line.split()
			if i == 0:
				self.Nsteps = float(data_set[0])
				self.dt = float(data_set[1])
				i += 1
			else:
				data.append(data_set)
		filename.close()

		self.N = len(data)
		self.Nsteps = int(self.Nsteps)

		self.Sun_pos_GR = np.zeros((3,self.N))
		self.Mercury_pos_GR = np.zeros((3, self.N))
		for i in range(0,3):
			for j in range(0,self.N):
				self.Sun_pos_GR[i][j] = data[j][i]
				self.Mercury_pos_GR[i][j] = data[j][i+3]

	def Earth_Sun_sys(self):
		method = ['Euler_method', 'Euler_Cromer_method', 'Verlet_method']
		filename_method = ["Earth_Sun_sys_euler.txt", "Earth_Sun_sys_eulercromer.txt", "Earth_Sun_sys_verlet.txt"]
		for i in range(0, len(method)):
			fig = plt.figure()
			self.read_data_ESJ_sys(filename_method[i], False)
			plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
			plt.hold("on")
			plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
			plt.xlabel('X - [AU]')
			plt.ylabel('Y - [AU]')
			plt.legend(['Earth', 'Sun'])
			plt.title('Earth and Sun system, '+method[i]+' N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

			if self.savefile:
				fig.savefig('../Plots/Earth_sun_'+method[i]+'.pdf')
			else:		
				plt.show()

	def ESJ_System(self):	
		fig = plt.figure()
		self.read_data_ESJ_sys("ESJ_sys.txt", True)
		plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
		plt.hold("on")
		plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
		plt.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:])
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		plt.legend(['Earth', 'Sun', 'Jupiter'])
		plt.title('Earth, Jupiter and Sun system. N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		if self.savefile:
			fig.savefig('../Plots/Earth_Sun_Jupiter.pdf')
		else:		
			plt.show()

	def plotting_2D(self):
		fig = plt.figure()
		plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], 'b--')
		plt.hold("on")
		plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], 'r--')
		plt.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], 'g-')
		plt.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:], 'k--')
		plt.plot(self.Venus_pos[0][:], self.Venus_pos[1][:], 'm--')
		plt.plot(self.Mars_pos[0][:], self.Mars_pos[1][:], 'c--')
		plt.plot(self.Saturn_pos[0][:], self.Saturn_pos[1][:], 'b-')
		plt.plot(self.Uranus_pos[0][:], self.Uranus_pos[1][:], 'r-')
		plt.plot(self.Neptune_pos[0][:], self.Neptune_pos[1][:], 'm-')
		plt.plot(self.Pluto_pos[0][:], self.Pluto_pos[1][:], 'k-')
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		plt.legend(['Earth', 'Sun', 'Jupiter', 'Mercury', 'Venus', 'Mars', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
		plt.title('Orbits of the planets. Timescale in years.')

		if self.savefile:
			fig.savefig('../Plots/'+'All_planets_2D_plot.pdf')
		else:
			plt.show()

	def plotting_3D(self):		
		fig1 = plt.figure()
		ax = fig1.gca(projection='3d')
		ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], self.Earth_pos[2][:], 'b--')
		plt.hold("on")
		ax.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:], 'r--')
		ax.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], self.Jupiter_pos[2][:], 'g-')
		ax.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:], self.Mercury_pos[2][:], 'k--')
		ax.plot(self.Venus_pos[0][:], self.Venus_pos[1][:], self.Venus_pos[2][:], 'm--')
		ax.plot(self.Mars_pos[0][:], self.Mars_pos[1][:], self.Mars_pos[2][:], 'c--')
		ax.plot(self.Saturn_pos[0][:], self.Saturn_pos[1][:], self.Saturn_pos[2][:], 'b-')
		ax.plot(self.Uranus_pos[0][:], self.Uranus_pos[1][:], self.Uranus_pos[2][:], 'r-')
		ax.plot(self.Neptune_pos[0][:], self.Neptune_pos[1][:], self.Neptune_pos[2][:], 'm-')
		ax.plot(self.Pluto_pos[0][:], self.Pluto_pos[1][:], self.Pluto_pos[2][:], 'k-')
		ax.set_xlabel('X - [AU]')
		ax.set_ylabel('Y - [AU]')
		ax.set_zlabel('Z - [AU]')
		ax.legend(['Earth', 'Sun', 'Jupiter', 'Mercury', 'Venus', 'Mars', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
		ax.set_title('Orbits of all planets')
		fig2 = plt.figure()
		ax2 = fig2.gca(projection='3d')
		ax2.plot(self.Earth_pos[0][1:100], self.Earth_pos[1][1:100], self.Earth_pos[2][1:100], 'b-')
		plt.hold("on")
		ax2.plot(self.Sun_pos[0][1:100], self.Sun_pos[1][1:100], self.Sun_pos[2][1:100], 'r-')
		ax2.plot(self.Mercury_pos[0][1:100], self.Mercury_pos[1][1:100], self.Mercury_pos[2][1:100], 'k-')
		ax2.plot(self.Venus_pos[0][1:100], self.Venus_pos[1][1:100], self.Venus_pos[2][1:100], 'm-')
		ax2.plot(self.Mars_pos[0][1:100], self.Mars_pos[1][1:100], self.Mars_pos[2][1:100], 'c-')
		ax2.set_xlabel('X - [AU]')
		ax2.set_ylabel('Y - [AU]')
		ax2.set_zlabel('Z - [AU]')
		ax2.legend(['Earth', 'Sun', 'Mercury', 'Venus', 'Mars'])
		ax2.set_title('Plot of the planets inside the asteroid belt')


		ax.scatter(self.Jupiter_pos[0][0], self.Jupiter_pos[1][0], self.Jupiter_pos[2][0], color='green', s=200)
		ax.scatter(self.Saturn_pos[0][0], self.Saturn_pos[1][0], self.Saturn_pos[2][0], color='green', s=150)
		ax.scatter(self.Neptune_pos[0][0], self.Neptune_pos[1][0], self.Neptune_pos[2][0], color='green', s=90)
		ax.scatter(self.Uranus_pos[0][0], self.Uranus_pos[1][0], self.Uranus_pos[2][0], color='green', s=100)
		ax.scatter(self.Pluto_pos[0][0], self.Pluto_pos[1][0], self.Pluto_pos[2][0], color='green', s=5)
		#ax.add_artist(circle_Jupiter)
		if self.savefile:
			fig1.savefig('../Plots/'+'All_planets_3D_plot.pdf')
		else:
			plt.show()
	def plot_mercury_GR(self):
		self.read_data_mercury_GR("Mercury_GR.txt")
		fig1 = plt.figure()
		ax = fig1.gca(projection='3d')
		plt.hold("on")
		ax.plot(self.Sun_pos_GR[0][:], self.Sun_pos_GR[1][:], self.Sun_pos_GR[2][:], 'r-')
		ax.plot(self.Mercury_pos_GR[0][:], self.Mercury_pos_GR[1][:], self.Mercury_pos_GR[2][:], 'b-')
		ax.set_xlabel('X - [AU]')
		ax.set_ylabel('Y - [AU]')
		ax.set_zlabel('Z - [AU]')
		ax.legend(['Sun', 'Mercury'])
		ax.set_title('Orbital path of Mercury around the Sun')
		ax.scatter(self.Sun_pos_GR[0][0], self.Sun_pos_GR[1][0], self.Sun_pos_GR[2][0], color='red', s=200)
		ax.scatter(self.Mercury_pos_GR[0][0], self.Mercury_pos_GR[1][0], self.Mercury_pos_GR[2][0], color='blue', s=30)
		fig2 = plt.figure()
		plt.plot(self.Sun_pos_GR[0][:], self.Sun_pos_GR[1][:], 'r-')
		plt.hold("on")
		plt.plot(self.Mercury_pos_GR[0][:], self.Mercury_pos_GR[1][:], 'b-')
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		plt.legend(['Sun', 'Mercury'])
		plt.show()



	def animate(self):
		# Animates 
		fig = plt.figure()
		#ax = fig.gca(projection='3d')
		ax = fig.add_subplot(111, projection='3d')
		ax.set_xlabel('X - [AU]')
		ax.set_ylabel('Y - [AU]')
		ax.set_zlabel('Z - [AU]')

		theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		Jx = 0.1*np.cos(theta)*np.sin(phi)+self.Jupiter_pos[0][0]
		Jy = 0.1*np.sin(theta)*np.sin(phi)+self.Jupiter_pos[1][0]
		Jz = 0.1*np.cos(phi)+self.Jupiter_pos[2][0]
		Jsphere =ax.plot_wireframe(Jx,Jy,Jz)
		Sx = 0.6*np.cos(theta)*np.sin(phi)+self.Sun_pos[0][0]
		Sy = 0.6*np.sin(theta)*np.sin(phi)+self.Sun_pos[1][0]
		Sz = 0.6*np.cos(phi)+self.Sun_pos[2][0]
		#Ssphere = ax.plot_wireframe(Sx, Sy, Sz, color='r')
		#plt.show()
		#jupiter = ax.scatter(self.Jupiter_pos[0][0], self.Jupiter_pos[1][0], self.Jupiter_pos[2][0], s=100, animated=True)
		
		def init():
			return

		def update(i):
			#jupiter = ax.scatter(self.Jupiter_pos[0][i], self.Jupiter_pos[1][i], self.Jupiter_pos[2][i], s=100, animated=True)
			#jupiter._offsets3D(self.Jupiter_pos[0][i], self.Jupiter_pos[1][i], self.Jupiter_pos[2][i])
			#jupiter.set_data(self.Jupiter_pos[0][i], self.Jupiter_pos[1][i], self.Jupiter_pos[2][i])
			Jx = 0.1*np.cos(theta)*np.sin(phi)+self.Jupiter_pos[0][i]
			Jy = 0.1*np.sin(theta)*np.sin(phi)+self.Jupiter_pos[1][i]
			Jz = 0.1*np.cos(phi)+self.Jupiter_pos[2][i]
			Jsphere = ax.plot_wireframe(Jx, Jy, Jz, animated=True)
			return Jsphere,


		anim = animation.FuncAnimation(fig, update, frames=300, interval=20, blit=True)
		
		plt.show()
solve = Plotter(False)
#solve.Earth_Sun_sys()
solve.ESJ_System()
#solve.plotting_3D()
#solve.animate()
#solve.plot_mercury_GR()