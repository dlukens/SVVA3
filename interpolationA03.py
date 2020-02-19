#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:29:40 2020

@author: frederik
"""

import numpy as np
import scipy as sc
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd

#import data set
data = np.loadtxt('/Users/frederik/GitHub/SVVA3/aeroload.dat',dtype='float', delimiter=',')

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

for i in range(Nz):
    data_z[i] = -Ca/4*(2 - math.cos(theta(i, Nz)) - math.cos(theta(i+1, Nz)))

for i in range(Nx):
    data_x[i] = la/4*(2 - math.cos(theta(i, Nx)) - math.cos(theta(i+1, Nx)))

#### interpolation ####

#new 2x81 matrix of data_z and aero force 
aeroforce_z = np.array([data_z,data.T[0,:]])


n = len(data_z)
a = np.zeros(n)
b = np.zeros(n)
c = np.zeros(n)
d = aeroforce_z[1,:]                                                       #d variable solved
h = np.zeros(n)


#define variable function for h, dependent on z data              
for i in range(0,n-1):           
    h[i] = (data_z[i+1] - data_z[i])

#create matrix A
A = np.zeros((n+1, n+1))

#non-zero corners of matrix A (top left, bottom right)
A[0,0] = 1
A[n,n] = 1

#rest of non-zero values of matrix A
for j in range(0,n-1):
    A[j+1,j+1] = 2*(h[j] + h[j+1])                                   #values along diagonal of matrix
    A[j+1,j] = h[j]                                                  #values parallell to diagonal below diagonal
    A[j+1,j+2] = h[j+1]                                              #values parallel to diagonal above diagonal
    
#create v_matrix
v = np.zeros(n+1)  
    
for k in range(0,n-2):
    v[k+1] = (3*((d[k+2]-d[k+1]/h[k+1] - (d[k+1]-d[k])/h[k])))

#solve A*b = v (where A, b and v are matrices)
#b are the values values M(0) --> M(n)
b = np.linalg.solve(A, v)                                           #b variable solved

#d = y values of data_z 
#calculate other coefficients         
     
for l in range(0, n-1):
    a[l] = (b[l+1] - b[l])/(6*h[l])                                 #a variable solved, not sure if the 6 should be six or 3
    c[l] = (d[l+1]-d[l])/h[l] - h[l]*b[l]/3 - h[l]*b[l+1]/6         #c variable solved
    
#make a plot of interpolation
#step size, step size is negative as x coordinates become increasingly negative towards TE
#if step size becomes smaller than -0.01, function becomes unbounded
delta = -0.01                                           
S = []
z = []

#S(x) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d

for t in range(3,n-1):
    for u in np.arange(data_z[t], data_z[t+1], delta):
        S.append(a[t]*abs(u - data_z[t])**3 + b[t]*abs(u - data_z[t])**2 + c[t]*abs(u - data_z[t]) + d[t])
        z.append(u)
    

plt.plot(z, S)
plt.show()


    