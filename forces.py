import numpy as np
from interp import C0_z, C0_x, grid_z, grid_x, data
import interp
import matplotlib.pyplot as plt
import tools as t

#Integrating once the interpolated force function
#S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d


def polyintegrate(C0):
    C1 = np.zeros((len(C0[:,0,0]), 4, len(C0[0,0,:])))
    #S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d
    #to
    #S(z) = 1/4*a(z-z_i)^4 + 1/3*b(z-z_i)^3 + 1/2*c(z-z_i)^2 + dz
    for i in range(len(C0[:,0,0])):
    
        C1[i,0,:] = 1/4*C0[i,0,:]
        C1[i,1,:] = 1/3*C0[i,1,:]
        C1[i,2,:] = 1/2*C0[i,2,:]
        C1[i,3,:] = C0[i,3,:]
        
    return(C1)


def integrate(C0, grid):
    I = np.zeros((len(C0[:,0,0]), len(grid)))
    I_sum = np.zeros(len(grid))
    
    for i in range(len(C0[:,0,0])):
        for j in range(len(C0[0,0,:])):
            
            def f(z):
                return 1/4*C0[i,0,j]*(z - grid[j])**4 + 1/3*C0[i,1,j]*(z - grid[j])**3 + 1/2*C0[i,2,j]*(z - grid[j])**2 + C0[i,3,j]*z
            
            I[i, j+1] = I[i, j] + f(grid[j+1]) - f(grid[j])
            
    I_sum = I.sum(axis=1)

    return I, I_sum



def integrate2(C1, grid):
    I = np.zeros((len(C1[:,0,0]), len(grid)))
    I_sum = np.zeros(len(grid))
    
    for i in range(len(C1[:,0,0])):
        for j in range(len(C1[0,0,:])):
            
            def f(z):
                return 1/5*C1[i,0,j]*(z - grid[j])**5 + 1/4*C1[i,1,j]*(z - grid[j])**4 + 1/3*C1[i,2,j]*(z - grid[j])**3 + 1/2*C1[i,3,j]*z**2
            
            I[i, j+1] = I[i, j] + f(grid[j+1]) - f(grid[j])
            
    I_sum = I.sum(axis=1)

    return I, I_sum


I_z, I_zsum = integrate(C0_z, grid_z)
I_x, I_xsum = integrate(C0_x, grid_x)

A_coeff = interp.interpolate1d(I_zsum, grid_x)    

#Coefficients from first integral

C1_z = polyintegrate(C0_z)
C1_x = polyintegrate(C0_x)


#Second integration
II_z, II_zsum = integrate2(C1_z, grid_z)
II_x, II_xsum = integrate2(C1_x, grid_x)


##cp location
Cp_x = -II_zsum/I_zsum
Cp_coeff = interp.interpolate1d(Cp_x, grid_x)

# interp.interplot(0.001, grid_x[:-2], Cp_coeff)


###Plotting###
def plot2d(gridx, gridz, data, title):
    #2D contour plot
    
    X, Z = np.meshgrid(gridx, gridz)
    Y = data
    
    plt.figure()
    cp = plt.contourf(X, Z, Y)
    plt.colorbar(cp)
    
    plt.title('{}'.format(title))
    plt.xlabel('X along wingspan [m]')
    plt.ylabel('Z along chord [m]')
    plt.plot(grid_x, Cp_x)
    plt.show()
    

# plot2d(grid_x, grid_z, data, 'Distributed force q')
# plot2d(grid_x, grid_z, I_x, 'First integral')
# plot2d(grid_x, grid_z, II_x, 'Second Integral')

# plot2d(grid_x, grid_z, data, 'Distributed force q')
# plot2d(grid_x, grid_z, I_z.T, 'First integral')
# plot2d(grid_x, grid_z, II_z.T, 'Second Integral')