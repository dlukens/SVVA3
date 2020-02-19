import numpy as np
import math
import matplotlib.pyplot as plt
import interpolationA03
import tools


    
data = np.loadtxt('aeroload.dat',dtype='float', delimiter=',')


#####Cordinates of aero force####

Nz = 81
Nx = 41
Ca = 0.515
la = 2.691

data_z = np.zeros(Nz) #Z-coordinats of the grid
data_x = np.zeros(Nx) #X-coordinates of the grid

def theta(i,N):
    t =  math.pi*(i-1)/N
    return t


for i in range(Nz):
    data_z[i] = -Ca/4*(2 - math.cos(theta(i, Nz)) - math.cos(theta(i+1, Nz)))
    

for i in range(Nx):
    data_x[i] = la/4*(2 - math.cos(theta(i, Nx)) - math.cos(theta(i+1, Nx)))
    
plt.plot(data, data_z[0])
plt.show()


########Aero Load Integration to find resultant force, location and moment##########
    
#print(interpolationA03.S)
    


    

