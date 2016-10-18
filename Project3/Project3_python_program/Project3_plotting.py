import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os
file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/Project3/build-Project3_cpp_program-Desktop_Qt_5_7_0_MinGW_32bit-Debug'

class Plotter():
	def __init__(self, savefile):
		self.savefile = savefile	# If True, saves the plots to a file

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

	def read_data_escape_vel(self, filename_open):
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
		self.Sun_pos_escape_vel = np.zeros((3,self.N))
		self.Earth_pos_escape_vel = np.zeros((3,self.N))
		for i in range(0,3):
			for j in range(0,self.N):
				self.Sun_pos_escape_vel[i][j] = data[j][i]
				self.Earth_pos_escape_vel[i][j] = data[j][i+3]

	def Earth_Sun_sys(self):
		""" 
		Function that plots for the Earth - Sun system.
		This function plots for two different values of dt.
		"""
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
				fig.savefig('../Plots/Earth_Sun_'+method[i]+'.pdf')
			else:		
				plt.show()

		larger_dt_method = ['Euler_method', 'Verlet_method']
		larger_dt_filename_method = ["Earth_Sun_sys_euler_largerdt.txt", "Earth_Sun_sys_verlet_largerdt.txt"]	
		for i in range(0, len(larger_dt_method)):
			fig = plt.figure()
			self.read_data_ESJ_sys(larger_dt_filename_method[i], False)
			plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
			plt.hold("on")
			plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
			plt.xlabel('X - [AU]')
			plt.ylabel('Y - [AU]')
			plt.legend(['Earth', 'Sun'])
			plt.title('Earth and Sun system, '+larger_dt_method[i]+' N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

			if self.savefile:
				fig.savefig('../Plots/Earth_Sun_'+larger_dt_method[i]+'_larger_dt'+'.pdf')
			else:		
				plt.show()

	def ESJ_System(self):
		""" Function that plots the system consisting of Earth, Sun and Jupiter. """
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

		fig2 = plt.figure()
		plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
		plt.hold("on")
		plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:])
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		plt.legend(['Earth', 'Sun'])
		plt.title('Plot of Earth and Sun in the Earth-Sun-Jupiter system.\n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		fig3 = plt.figure()
		ax = fig3.gca(projection='3d')
		ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], self.Earth_pos[2][:])
		plt.hold("on")
		ax.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:])
		ax.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], self.Jupiter_pos[2][:])
		ax.set_xlabel('X - [AU]')
		ax.set_ylabel('Y - [AU]')
		ax.set_zlabel('Z - [AU]')
		ax.legend(['Earth', 'Sun', 'Jupiter'])
		ax.set_title('3D Orbits of the Earth-Jupiter-Sun system. N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		if self.savefile:
			fig.savefig('../Plots/Earth_Sun_Jupiter.pdf')
			fig2.savefig('../Plots/ESJ_EarthandSun.pdf')
			fig3.savefig('../Plots/Earth_Sun_Jupiter_3D.pdf')
		else:		
			plt.show()

	def plot_escape_velocity(self):
		""" 
		Function that plots the escape velocity of a planet.
		Also reads and plots the planetary orbit of Earth from previous exercises
		which is used to compare the orbits with and without escape velocity.
		"""
		fig = plt.figure()
		self.read_data_escape_vel("Escape_velocity_system.txt")
		self.read_data_ESJ_sys("Earth_Sun_sys_verlet_largerdt.txt", False)
		plt.plot(self.Earth_pos_escape_vel[0][:], self.Earth_pos_escape_vel[1][:])
		plt.hold("on")
		plt.plot(self.Sun_pos_escape_vel[0][:], self.Sun_pos_escape_vel[1][:])
		plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:])
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		plt.legend(['Planet $v = v_{escape}$', 'Sun', r'Planet $v \neq v_{escape}$'])
		plt.title('Escape velocity of a planet. \n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))		
		if self.savefile:
			fig.savefig('../Plots/Escape_velocity.pdf')
		else:		
			plt.show()


	def plotting_2D(self):
		""" Function that orbit of all planets in 2 dimensions. """
		fig = plt.figure()
		self.read_data("SolarSys_All_planets.txt")
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
		plt.title('2D Orbits of the planets. N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		if self.savefile:
			fig.savefig('../Plots/All_planets_2D_plot.pdf')
		else:
			plt.show()

	def plotting_3D(self):	
		""" 
		Function that plots all the planets in 3 dimensions.
		Creates a scatter plot to show the final position of the planets. 
		(Scatter plot size not scaled to the real planet sizes)
		Also plots the first four planets in it's own plot to better see their orbits.
		"""	
		fig1 = plt.figure()
		self.read_data("SolarSys_All_planets.txt")
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
		ax.set_title('3D Orbits of all planets. \n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		ax.scatter(self.Jupiter_pos[0][0], self.Jupiter_pos[1][0], self.Jupiter_pos[2][0], color='green', s=200)
		ax.scatter(self.Saturn_pos[0][0], self.Saturn_pos[1][0], self.Saturn_pos[2][0], color='green', s=150)
		ax.scatter(self.Neptune_pos[0][0], self.Neptune_pos[1][0], self.Neptune_pos[2][0], color='green', s=90)
		ax.scatter(self.Uranus_pos[0][0], self.Uranus_pos[1][0], self.Uranus_pos[2][0], color='green', s=100)
		ax.scatter(self.Pluto_pos[0][0], self.Pluto_pos[1][0], self.Pluto_pos[2][0], color='green', s=20)


		fig2 = plt.figure()
		self.read_data("SolarSys_All_planets_First4.txt")
		ax2 = fig2.gca(projection='3d')
		ax2.plot(self.Earth_pos[0][1:100], self.Earth_pos[1][1:100], self.Earth_pos[2][1:100], 'b-')
		plt.hold("on")
		ax2.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:], 'r-')
		ax2.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:], self.Mercury_pos[2][:], 'k-')
		ax2.plot(self.Venus_pos[0][:], self.Venus_pos[1][:], self.Venus_pos[2][:], 'm-')
		ax2.plot(self.Mars_pos[0][:], self.Mars_pos[1][:], self.Mars_pos[2][:], 'c-')
		ax2.set_xlabel('X - [AU]')
		ax2.set_ylabel('Y - [AU]')
		ax2.set_zlabel('Z - [AU]')
		ax2.legend(['Earth', 'Sun', 'Mercury', 'Venus', 'Mars'])
		ax2.set_title('3D Plot inner planets. \n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		ax2.scatter(self.Sun_pos[0][0], self.Sun_pos[1][0], self.Sun_pos[2][0], color='green', s=300)
		ax2.scatter(self.Mercury_pos[0][0], self.Mercury_pos[1][0], self.Mercury_pos[2][0], color='green', s=30)
		ax2.scatter(self.Venus_pos[0][0], self.Venus_pos[1][0], self.Venus_pos[2][0], color='green', s=65)
		ax2.scatter(self.Earth_pos[0][0], self.Earth_pos[1][0], self.Earth_pos[2][0], color='green', s=80)
		ax2.scatter(self.Mars_pos[0][0], self.Mars_pos[1][0], self.Mars_pos[2][0], color='green', s=70)

		if self.savefile:
			fig1.savefig('../Plots/All_planets_3D_plot.pdf')
			fig2.savefig('../Plots/First_planets_3D_plot.pdf')
		else:
			plt.show()
	def plot_mercury_GR(self):
		""" Function plotting the effect of general relativity correction for Mercury """
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
solve = Plotter(True)
solve.Earth_Sun_sys()
solve.plot_escape_velocity()
solve.ESJ_System()
#solve.plotting_3D()
#solve.animate()
#solve.plot_mercury_GR()