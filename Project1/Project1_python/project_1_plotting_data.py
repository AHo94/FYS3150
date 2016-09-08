import matplotlib.pyplot as plt
import numpy as np
import os

file_directory = 'C:/Users/Alex/Documents/FYS3150/FYS3150_projects/build-Project_1_cpp-Desktop_Qt_5_7_0_MinGW_32bit-Debug'
filename = 'project_1_data.txt'
f = open(os.path.join(file_directory, filename), 'r')

i = 0
data = []
for line in f:
	n_value = line.split()
	if i != 1:
		data.append(n_value)
	i += 1
f.close()
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
plt.legend(['Exact solution','Numerical solution'])
plt.title('Plot of velocities, n = %g' %(n))
plt.show()