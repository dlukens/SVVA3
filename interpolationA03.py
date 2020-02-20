import numpy as np
import math
import matplotlib.pyplot as plt

#import data set
data = np.loadtxt('aeroload.dat',dtype='float', delimiter=',')

#### Cordinates of aero force ####

Nz = 81
Nx = 41
Ca = 0.515
la = 2.691

grid_z = np.zeros(Nz)
grid_x = np.zeros(Nx)

def theta(i,N):
    t =  math.pi*(i-1)/N
    return t

#Grid on chord - Z-Axis
for i in range(1, Nz+1):
    grid_z[i-1] = -Ca/4*(2 - math.cos(theta(i, Nz)) - math.cos(theta(i+1, Nz)))
    
#Grid on span - X-Axis
for i in range(1, Nx+1):
    grid_x[i-1] = la/4*(2 - math.cos(theta(i, Nx)) - math.cos(theta(i+1, Nx)))

######### interpolation ##########

aeroforce_z = data.T
aeroforce_x = data


def interpolate(force, grid):
    C0 = np.zeros((len(force)-1, 4, len(force.T)-1))

    for r in range(len(force)-1):

        n = len(grid)-1
        h = np.zeros(n)
        a = np.zeros(n)
        b = np.zeros(n)
        c = np.zeros(n)
        d = force[r,:] #d variable solved
        
        
        for i in range(0,n):           
            h[i] = (grid[i+1] - grid[i])
        
        #create matrix A
        A = np.zeros((n+1, n+1))
        
        #non-zero corners of matrix A (top left, bottom right)
        A[0,0] = 1
        A[n,n] = 1
        
        #rest of non-zero values of matrix A
        for j in range(1,n):
            A[j,j] = (h[j-1] + h[j])/3
            A[j,j-1] = h[j-1]/6
            A[j,j+1] = h[j]/6
            
        #create v_matrix
        v = np.zeros(n+1)  
            
        for k in range(0,n-1):
            v[k+1] = 3*((d[k+2]-d[k+1])/h[k+1] - (d[k+1]-d[k])/h[k])
        
        #solve A*b = v
        #b are the values M(0) --> M(n)
        b = np.linalg.solve(A, v)
        
        #Solve for a and c   
        for l in range(0, n):
            c[l] = (d[l+1]-d[l])/h[l] - h[l]*(2*b[l]+b[l+1])/3        #c variable solved
            a[l] = (b[l+1] - b[l])/(3*h[l])                           #a variable solved
            
        #Remove last element to have 80 intervals
        b = b[:-1]
        d = d[:-1]
        
        #Now 3D array with all C0icients for each section and each chord
            #Index 0: Spanwise chord section - along X-axis
            #Index 1: C0icient value (a, b, c, d)
            #Index 2: Chordwise interval - along -Z-Axis
        
        C0[r, 0, :] = a
        C0[r, 1, :] = b
        C0[r, 2, :] = c
        C0[r, 3, :] = d
        
      
    return C0

def interplot(delta, grid, C0):
    #Use negative step size
    n = len(grid)-1

    S = []
    z = []
    
    #S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d
    
    for t in range(0,n-1):
        for u in np.arange(grid[t], grid[t+1], delta):
            S.append(C0[0, t]*(u - grid[t])**3 + C0[1, t]*(u - grid[t])**2 + C0[2, t]*(u - grid[t]) + C0[3, t])
            z.append(u)
            
    plt.plot(z, S)
    plt.show()


#Calling Functions

C0_z = interpolate(aeroforce_z, grid_z)
C0_x = interpolate(aeroforce_x, grid_x)


interplot(0.0001, grid_x, C0_x[5,:,:])
interplot(-0.0001, grid_z, C0_z[5,:,:])



    