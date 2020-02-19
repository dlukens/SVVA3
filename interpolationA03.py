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

data_z = np.zeros(Nz)
data_x = np.zeros(Nx)

def theta(i,N):
    t =  math.pi*(i-1)/N
    return t

for i in range(1, Nz+1):
    data_z[i-1] = -Ca/4*(2 - math.cos(theta(i, Nz)) - math.cos(theta(i+1, Nz)))
    

for i in range(1, Nx+1):
    data_x[i-1] = la/4*(2 - math.cos(theta(i, Nx)) - math.cos(theta(i+1, Nx)))

######### interpolation ##########

#new 2x81 matrix of data_z and aero force 
aeroforce_z = np.array([data_z,data.T[0,:]])


n = len(data_z)-1
a = np.zeros(n)
b = np.zeros(n)
c = np.zeros(n)
d = aeroforce_z[1,:]                                                       #d variable solved
h = np.zeros(n)


#define variable function for h, dependent on z data              
for i in range(0,n):           
    h[i] = (data_z[i+1] - data_z[i])

#create matrix A
A = np.zeros((n+1, n+1))

#non-zero corners of matrix A (top left, bottom right)
A[0,0] = 1
A[n,n] = 1

#rest of non-zero values of matrix A
for j in range(1,n):
    A[j,j] = (h[j-1] + h[j])/3                               #values along diagonal of matrix
    A[j,j-1] = h[j-1]/6                                                #values parallell to diagonal below diagonal
    A[j,j+1] = h[j]/6                                           #values parallel to diagonal above diagonal
    
#create v_matrix
v = np.zeros(n+1)  
    
for k in range(0,n-1):
    v[k+1] = 3*((d[k+2]-d[k+1])/h[k+1] - (d[k+1]-d[k])/h[k])

#solve A*b = v (where A, b and v are matrices)
#b are the values values M(0) --> M(n)
b = np.linalg.solve(A, v)                                           #b variable solved

#d = y values of data_z 
#calculate other coefficients         
     
for l in range(0, n-1):
    c[l] = (d[l+1]-d[l])/h[l] - h[l]*(2*b[l]+b[l+1])/3        #c variable solved
    a[l] = (b[l+1] - b[l])/(3*h[l])                           #a variable solved

    
#make a plot of interpolation
#step size, step size is negative as x coordinates become increasingly negative towards TE
#if step size becomes smaller than -0.01, function becomes unbounded
delta = -0.0001                                           
S = []
z = []

#S(x) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d

for t in range(0,n-1):
    for u in np.arange(data_z[t], data_z[t+1], delta):
        S.append(a[t]*(u - data_z[t])**3 + b[t]*(u - data_z[t])**2 + c[t]*(u - data_z[t]) + d[t])
        z.append(u)
    
    
print(a, b, c, d)

plt.plot(z, S)
plt.show()


    