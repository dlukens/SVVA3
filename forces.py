import numpy as np
import tools as t
from interpolationA03 import C0_z, C0_x, data_z, data_x, data

#Integrating once the interpolated force function

def polyintegrate(C0):
    C1 = np.zeros((len(C0[:,0,0]), 5, len(C0[0,0,:])))
    #S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d
    #to
    #S(z) = 1/4*a(z-z_i)^4 + 1/3*b(z-z_i)^3 + 1/2*c(z-z_i)^2 + dz + e
    for i in range(len(C0_z[:,0,0])):
    
        C1[i,0,:] = 1/4*C0[i,0,:]
        C1[i,1,:] = 1/3*C0[i,1,:]
        C1[i,2,:] = 1/2*C0[i,2,:]
        C1[i,3,:] = C0[i,3,:]
        C1[i,4,:] = 1
        
    return(C1)

C1_z = polyintegrate(C0_z)
C1_x = polyintegrate(C0_x)

#Integrating twice the interpolated force function

def polyintegrate2(C1):
    C2 = np.zeros((len(C1[:,0,0]), 6, len(C1[0,0,:])))
    #S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d
    #to
    #S(z) = 1/5*a(z-z_i)^5 + 1/4*b(z-z_i)^4 + 1/3*c(z-z_i)^3 + dz^2 + ez
    for i in range(len(C0_z[:,0,0])):
    
        C2[i,0,:] = 1/5*C1[i,0,:]
        C2[i,1,:] = 1/4*C1[i,1,:]
        C2[i,2,:] = 1/3*C1[i,2,:]
        C2[i,3,:] = 1/2*C1[i,3,:]
        C2[i,4,:] = C1[i,4,:]
        C2[i,5,:] = 1
        
    return(C2)

C2_z = polyintegrate(C1_z)
C2_x = polyintegrate(C1_x)