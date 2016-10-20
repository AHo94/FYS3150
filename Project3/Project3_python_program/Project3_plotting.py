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
			ax = fig.add_subplot(1,1,1)
			self.read_data_ESJ_sys(filename_method[i], False)
			plt.hold("on")
			line1, = ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
			line2, = ax.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
			plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
			plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
			ax.set_xlabel('X - [AU]')
			ax.set_ylabel('Y - [AU]')
			ax.legend()
			plt.title('Earth and Sun system, '+method[i]+'.\n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
			plt.axis('equal')

			if self.savefile:
				fig.savefig('../Plots/Earth_Sun_'+method[i]+'.pdf')
			else:		
				plt.show()

		larger_dt_method = ['Euler_method', 'Verlet_method']
		larger_dt_filename_method = ["Earth_Sun_sys_euler_largerdt.txt", "Earth_Sun_sys_verlet_largerdt.txt"]	
		for i in range(0, len(larger_dt_method)):
			fig = plt.figure()
			ax = fig.add_subplot(1,1,1)
			self.read_data_ESJ_sys(larger_dt_filename_method[i], False)
			plt.hold("on")
			line1, = ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
			line2, = ax.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
			plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
			plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
			ax.set_xlabel('X - [AU]')
			ax.set_ylabel('Y - [AU]')
			ax.legend()
			plt.title('Earth and Sun system, '+larger_dt_method[i]+'. \n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
			plt.axis('equal')

			if self.savefile:
				fig.savefig('../Plots/Earth_Sun_'+larger_dt_method[i]+'_larger_dt'+'.pdf')
			else:		
				plt.show()

	def ESJ_System(self):
		""" Function that plots the system consisting of Earth, Sun and Jupiter. """
		# 2D plot
		fig = plt.figure()
		ax1 = fig.add_subplot(1,1,1)
		self.read_data_ESJ_sys("ESJ_sys.txt", True)
		plt.hold("on")
		line1, = ax1.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
		line2, = ax1.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
		line3, = ax1.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], label='Jupiter')
		plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
		plt.plot(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], 'o', color=line3.get_color())
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax1.legend()
		plt.title('Earth, Jupiter and Sun system. N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
		plt.axis('equal')

		# Sun and Earth plot
		fig2 = plt.figure()
		ax2 = fig2.add_subplot(1,1,1)
		plt.hold("on")
		line1, = ax2.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
		line2, = ax2.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
		plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax2.legend()
		plt.title('Plot of Earth and Sun in the Earth-Sun-Jupiter system.\n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
		ax2.set_xlim(-2,2)
		ax2.set_ylim(-2,2)
		ax2.set_aspect('equal')

		# 3D Plot
		fig3 = plt.figure()
		ax = fig3.gca(projection='3d')
		plt.hold("on")
		ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], self.Earth_pos[2][:])
		ax.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:])
		ax.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], self.Jupiter_pos[2][:])
		ax.set_xlabel('X - [AU]')
		ax.set_ylabel('Y - [AU]')
		ax.set_zlabel('Z - [AU]')
		ax.legend(['Earth', 'Sun', 'Jupiter'], bbox_to_anchor=(-0.1, 0.9), loc=2)
		ax.set_title('3D Orbits of the Earth-Jupiter-Sun system. N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		ax.scatter(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], self.Jupiter_pos[2][-1], color='green', s=150)
		ax.scatter(self.Sun_pos[0][-1], self.Sun_pos[1][-1], self.Sun_pos[2][-1], color='green', s=300)
		ax.scatter(self.Earth_pos[0][-1], self.Earth_pos[1][-1], self.Earth_pos[2][-1], color='green', s=90)

		# 10MJ
		fig4 = plt.figure()
		ax3 = fig4.add_subplot(1,1,1)
		self.read_data_ESJ_sys("ESJ_sys_10MJ.txt", True)
		plt.hold("on")
		line1, = ax3.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
		line2, = ax3.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
		line3, = ax3.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], label='Jupiter')
		plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
		plt.plot(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], 'o', color=line3.get_color())
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax3.legend()
		plt.title('Plot Earth-Jupiter-Sun system, with $10M_{Jupiter}$. \n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
		plt.axis('equal')

		# 1000MJ
		fig5 = plt.figure()
		ax4 = fig5.add_subplot(1,1,1)
		self.read_data_ESJ_sys("ESJ_sys_1000MJ.txt", True)
		plt.hold("on")
		line1, = ax4.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
		line2, = ax4.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
		line3, = ax4.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], label='Jupiter')
		plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
		plt.plot(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], 'o', color=line3.get_color())
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax4.legend()
		plt.title('Plot Earth-Jupiter-Sun system, with $1000M_{Jupiter}$. \
			\n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
		plt.axis('equal')

		fig6 = plt.figure()
		ax5 = fig6.gca(projection='3d')
		plt.hold("on")
		ax5.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], self.Earth_pos[2][:])
		ax5.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:])
		ax5.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], self.Jupiter_pos[2][:])
		ax5.set_xlabel('X - [AU]')
		ax5.set_ylabel('Y - [AU]')
		ax5.set_zlabel('Z - [AU]')
		ax5.legend(['Earth', 'Sun', 'Jupiter'], bbox_to_anchor=(-0.1, 0.9), loc=2)
		ax5.set_title('3D Orbits of the Earth-Jupiter-Sun system, with $1000M_{Jupiter}$. \
				 	\n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		ax5.scatter(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], self.Jupiter_pos[2][-1], color='green', s=150)
		ax5.scatter(self.Sun_pos[0][-1], self.Sun_pos[1][-1], self.Sun_pos[2][-1], color='green', s=300)
		ax5.scatter(self.Earth_pos[0][-1], self.Earth_pos[1][-1], self.Earth_pos[2][-1], color='green', s=90)

		fig7 = plt.figure()
		ax6 = fig7.add_subplot(1,1,1)
		self.read_data_ESJ_sys("ESJ_sys_largerdt.txt", True)
		plt.hold("on")
		line1, = ax6.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label='Earth')
		line2, = ax6.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], label='Sun')
		line3, = ax6.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], label='Jupiter')
		plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos[0][-1], self.Sun_pos[1][-1], 'o', color=line2.get_color())
		plt.plot(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], 'o', color=line3.get_color())

		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax6.legend()
		plt.title('Earth, Jupiter and Sun system. N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
		plt.axis('equal')
				
		if self.savefile:
			fig.savefig('../Plots/Earth_Sun_Jupiter.pdf')
			fig2.savefig('../Plots/ESJ_EarthandSun.pdf')
			fig3.savefig('../Plots/Earth_Sun_Jupiter_3D.pdf')
			fig4.savefig('../Plots/Earth_Sun_Jupiter_10MJ.pdf')
			fig5.savefig('../Plots/Earth_Sun_Jupiter_1000MJ.pdf')
			fig6.savefig('../Plots/Earth_Sun_Jupiter_1000MJ_3D.pdf')
			fig7.savefig('../Plots/Earth_Sun_Jupiter_largerdt.pdf')
		else:
			plt.show()

	def plot_escape_velocity(self):
		""" 
		Function that plots the escape velocity of a planet.
		Also reads and plots the planetary orbit of Earth from previous exercises
		which is used to compare the orbits with and without escape velocity.
		"""
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		self.read_data_escape_vel("Escape_velocity_system.txt")
		self.read_data_ESJ_sys("Earth_Sun_sys_verlet_largerdt.txt", False)
		plt.hold("on")
		line1, = ax.plot(self.Earth_pos_escape_vel[0][:], self.Earth_pos_escape_vel[1][:], label = 'Planet $v = v_{escape}$')
		line2, = ax.plot(self.Sun_pos_escape_vel[0][:], self.Sun_pos_escape_vel[1][:], label = 'Sun')
		line3, = ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], label=r'Planet $v \neq v_{escape}$')
		plt.plot(self.Earth_pos_escape_vel[0][-1], self.Earth_pos_escape_vel[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos_escape_vel[0][-1], self.Sun_pos_escape_vel[1][-1], 'o', color=line2.get_color())
		plt.plot(self.Earth_pos[0][-1], self.Earth_pos[1][-1], 'o', color=line3.get_color())
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax.legend()
		plt.title('Escape velocity of a planet. \n N=%.f, dt=%.g, years = %.f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))	
		plt.axis('equal')	
		if self.savefile:
			fig.savefig('../Plots/Escape_velocity.pdf')
		else:		
			plt.show()


	def plotting_2D(self):
		""" Function that orbit of all planets in 2 dimensions. """
		fig = plt.figure()
		plt.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], 'r--')
		self.read_data("SolarSys_All_planets.txt")
		plt.hold("on")
		plt.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:], 'k--')
		plt.plot(self.Venus_pos[0][:], self.Venus_pos[1][:], 'm--')
		plt.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], 'b--')
		plt.plot(self.Mars_pos[0][:], self.Mars_pos[1][:], 'c--')
		plt.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], 'g-')
		plt.plot(self.Saturn_pos[0][:], self.Saturn_pos[1][:], 'b-')
		plt.plot(self.Uranus_pos[0][:], self.Uranus_pos[1][:], 'r-')
		plt.plot(self.Neptune_pos[0][:], self.Neptune_pos[1][:], 'm-')
		plt.plot(self.Pluto_pos[0][:], self.Pluto_pos[1][:], 'k-')
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		plt.legend(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
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
		ax.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:], 'r--')
		plt.hold("on")
		ax.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:], self.Mercury_pos[2][:], 'k--')
		ax.plot(self.Venus_pos[0][:], self.Venus_pos[1][:], self.Venus_pos[2][:], 'm--')
		ax.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], self.Earth_pos[2][:], 'b--')
		ax.plot(self.Mars_pos[0][:], self.Mars_pos[1][:], self.Mars_pos[2][:], 'c--')
		ax.plot(self.Jupiter_pos[0][:], self.Jupiter_pos[1][:], self.Jupiter_pos[2][:], 'g-')
		ax.plot(self.Saturn_pos[0][:], self.Saturn_pos[1][:], self.Saturn_pos[2][:], 'b-')
		ax.plot(self.Uranus_pos[0][:], self.Uranus_pos[1][:], self.Uranus_pos[2][:], 'r-')
		ax.plot(self.Neptune_pos[0][:], self.Neptune_pos[1][:], self.Neptune_pos[2][:], 'm-')
		ax.plot(self.Pluto_pos[0][:], self.Pluto_pos[1][:], self.Pluto_pos[2][:], 'k-')
		ax.set_xlabel('X - [AU]')
		ax.set_ylabel('Y - [AU]')
		ax.set_zlabel('Z - [AU]')
		ax.legend(['Sun', 'Earth', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'], bbox_to_anchor=(-0.1, 0.9), loc=2)
		ax.set_title('3D Orbits of all planets. \n N=%.f, dt=%.g, years = %.1f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		ax.scatter(self.Jupiter_pos[0][-1], self.Jupiter_pos[1][-1], self.Jupiter_pos[2][-1], color='green', s=200)
		ax.scatter(self.Saturn_pos[0][-1], self.Saturn_pos[1][-1], self.Saturn_pos[2][-1], color='green', s=150)
		ax.scatter(self.Neptune_pos[0][-1], self.Neptune_pos[1][-1], self.Neptune_pos[2][-1], color='green', s=90)
		ax.scatter(self.Uranus_pos[0][-1], self.Uranus_pos[1][-1], self.Uranus_pos[2][-1], color='green', s=100)
		ax.scatter(self.Pluto_pos[0][-1], self.Pluto_pos[1][-1], self.Pluto_pos[2][-1], color='green', s=20)
		ax.set_aspect('equal')


		fig2 = plt.figure()
		self.read_data("SolarSys_All_planets_First4.txt")
		ax2 = fig2.gca(projection='3d')
		ax2.plot(self.Sun_pos[0][:], self.Sun_pos[1][:], self.Sun_pos[2][:], 'r-')
		plt.hold("on")
		ax2.plot(self.Mercury_pos[0][:], self.Mercury_pos[1][:], self.Mercury_pos[2][:], 'k-')
		ax2.plot(self.Venus_pos[0][:], self.Venus_pos[1][:], self.Venus_pos[2][:], 'm-')
		ax2.plot(self.Earth_pos[0][:], self.Earth_pos[1][:], self.Earth_pos[2][:], 'b-')
		ax2.plot(self.Mars_pos[0][:], self.Mars_pos[1][:], self.Mars_pos[2][:], 'c-')
		ax2.set_xlabel('X - [AU]')
		ax2.set_ylabel('Y - [AU]')
		ax2.set_zlabel('Z - [AU]')
		ax2.legend(['Sun', 'Mercury', 'Earth', 'Venus', 'Mars'], bbox_to_anchor=(-0.1, 0.9), loc=2)
		ax2.set_title('3D Plot of the first four planets. \n N=%.f, dt=%.g, years = %.1f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))

		ax2.scatter(self.Sun_pos[0][-1], self.Sun_pos[1][-1], self.Sun_pos[2][-1], color='green', s=400)
		ax2.scatter(self.Mercury_pos[0][-1], self.Mercury_pos[1][-1], self.Mercury_pos[2][-1], color='green', s=30)
		ax2.scatter(self.Venus_pos[0][-1], self.Venus_pos[1][-1], self.Venus_pos[2][-1], color='green', s=65)
		ax2.scatter(self.Earth_pos[0][-1], self.Earth_pos[1][-1], self.Earth_pos[2][-1], color='green', s=80)
		ax2.scatter(self.Mars_pos[0][-1], self.Mars_pos[1][-1], self.Mars_pos[2][-1], color='green', s=70)
		ax2.set_aspect('equal')

		if self.savefile:
			fig1.savefig('../Plots/All_planets_3D_plot.pdf')
			fig2.savefig('../Plots/First_planets_3D_plot.pdf')
		else:
			plt.show()
	def plot_mercury_GR(self):
		""" Function plotting the effect of general relativity correction for Mercury """
		self.read_data_mercury_GR("Mercury_GR.txt")
		fig2 = plt.figure()
		ax2 = fig2.add_subplot(1,1,1)
		plt.hold("on")
		line1, = ax2.plot(self.Mercury_pos_GR[0][:], self.Mercury_pos_GR[1][:], label='Mercury')
		line2, = ax2.plot(self.Sun_pos_GR[0][:], self.Sun_pos_GR[1][:], label='Sun')
		plt.plot(self.Mercury_pos_GR[0][-1], self.Mercury_pos_GR[1][-1], 'o', color=line1.get_color())
		plt.plot(self.Sun_pos_GR[0][-1], self.Sun_pos_GR[1][-1], 'o', color=line2.get_color())
		plt.xlabel('X - [AU]')
		plt.ylabel('Y - [AU]')
		ax2.legend()
		plt.title('Orbital path of Mercury with GR correction. \n N=%.f, dt=%.g, years = %.1f' %(self.Nsteps, self.dt, self.Nsteps*self.dt))
		plt.axis('equal')

		if self.savefile == True:
			fig2.savefig('../Plots/Mercury_orbit_GR.pdf')
		else:
			plt.show()


solve = Plotter(True)
#solve.Earth_Sun_sys()
#solve.plot_escape_velocity()
solve.ESJ_System()
#solve.plotting_3D()
#solve.plot_mercury_GR()