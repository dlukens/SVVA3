# -*- coding: utf-8 -*-
import math
from sectionproperties import SCz, x1, x2, x3, xa, R1x, R1y, R2, R2x, R2y, R3, R3x, R3y, P
import numpy as np
from forces import A_coeff, Cp_coeff
from interp import grid_x, grid_z


###variables
alpha = math.radians(25)

###Funcions

def step(x, x1, exp): #Mcaulay step function
    y = x - x1
    if y <= 0: return 0 
    else: return y**exp

def A_SC_int(x):
    idx = np.searchsorted(grid_x, x)-1
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx]
    e = Cp_coeff[0, idx]
    f = Cp_coeff[1, idx]
    g = Cp_coeff[2, idx]
    h = Cp_coeff[3, idx]
    return x*(a*SCz/4*(x-grid_x[idx])**3 + b*SCz/3*(x-grid_x[idx])**2 + c*SCz/2*(x-grid_x[idx]) + d*SCz   \
                        -  a/7*(x-grid_x[idx])**6                                  \
                        - (a*f + b*e)/6*(x-grid_x[idx])**5                         \
                        - (a*g + b*f + c*e)/5*(x-grid_x[idx])**4                   \
                        - (a*h + b*g + c*f + d*e)/4*(x-grid_x[idx])**3             \
                        - (b*h + c*g + d*f)/3*(x-grid_x[idx])**2                   \
                        - (c*h + d*g)/2*(x-grid_x[idx])                            \
                        -  d*h)

                        
def A_SC_doubleint(x):
    idx = np.searchsorted(grid_x, x)-1
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx]
    e = Cp_coeff[0, idx]
    f = Cp_coeff[1, idx]
    g = Cp_coeff[2, idx]
    h = Cp_coeff[3, idx]
    return x**2*(a*SCz/(4*5)*(x-grid_x[idx])**3 + b*SCz/(3*4)*(x-grid_x[idx])**2 + c*SCz/(2*3)*(x-grid_x[idx]) + d*SCz/2 \
                        -  a/(8*7)*(x-grid_x[idx])**6 \
                        - (a*f + b*e)/(6*7)*(x-grid_x[idx])**5 \
                        - (a*g + b*f + c*e)/(5*6)*(x-grid_x[idx])**4 \
                        - (a*h + b*g + c*f + d*e)/(4*5)*(x-grid_x[idx])**3 \
                        - (b*h + c*g + d*f)/(3*4)*(x-grid_x[idx])**2 \
                        - (c*h + d*g)/(2*3)*(x-grid_x[idx]) \
                        -  d*h/2)
    
    
def A_int(x):
    idx = np.searchsorted(grid_x, x)-1
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx]    
    return x*(a/4*(x-grid_x[idx])**3      \
            + b/3*(x-grid_x[idx])**2      \
            + c/2*(x-grid_x[idx])         \
            + d)
    
def A_doubleint(x):
    idx = np.searchsorted(grid_x, x)-1
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx]  
    return x**2*(a/(4*5)*(x-grid_x[idx])**3     \
               + b/(3*4)*(x-grid_x[idx])**2     \
               + c/(2*3)*(x-grid_x[idx])        \
               + d/2)

def A_quadint(x):
    idx = np.searchsorted(grid_x, x)-1
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx] 
    return x**4*(a/(4*5*6*7)*(x-grid_x[idx])**3     \
               + b/(3*4*5*6)*(x-grid_x[idx])**2     \
               + c/(2*3*4*5)*(x-grid_x[idx])        \
               + d/(2*3*4))

    
#Reaction forces
    
xI=x2-0.5*xa

# Rxn= [[0,-SCz,0,-SCz,0,-SCz,-SCz*m.sin(alpha),0,0,0,0,0],                                                                                                   #T(la)
#        [La-x1,0,La-x2,0,La-x3,0,-m.cos(alpha)*(La-x2+0.5*xa),0,0,0,0,0],                                                                                     #My(la)
#        [0,-(La-x1),0,(La-x2),0,(La-x3),-m.sin(alpha)*(La-x2+xa*0.5),0,0,0,0,0],                                                                              #Mz(la)
#        [1,0,1,0,1,0,-m.cos(alpha),0,0,0,0,0],                                                                                                                #Sz(la)
#        [0,-1,0,-1,0,-1,-m.sin(alpha),0,0,0,0,0],                                                                                                             #Sy(la)
#        [0,0,0,0,0,0,0,0,0,x1**3,1,0],                                                                                                                        #w(x1)
#        [(x2-x1)**3,0,0,0,0,0,-m.cos(alpha)*(xa/2)**3,0,0,x2**3,1,0],                                                                                         #w(x2)
#        [(x3-x1)**3,0,(x3-x2)**3,0,0,0,0,0,0,x3**3,1,0],                                                                                                      #w(x3)
#        [-m.cos(alpha)*(xI-x1)**3/(E*Izz),-m.sin(alpha)*(xI-x1)**3/(E*Iyy) + SCz**2*(xI-x1)*m.sin(alpha)/(G*J),0,0,0,0,0,xI*m.sin(alpha)/(E*Iyy),m.sin(alpha)/(E*Iyy),-xI*cos(alpha)/(E*Izz),-m.cos(alpha)/(E*Izz),-SCz*m.sin(alpha)/(G*J)],                    #w'(xI)=0
#        [0,0,0,0,0,0,0,x1/(-E*Izz),1/(-E*Izz),0,0,SCz/(G*J)],                                                                                                                           #v(x1)+theta(x1)
#        [0,-(x2-x1)**3/(-E*Izz)-SCz**2*(x2-x1)/(G*J),0,0,0,0,-m.sin(alpha)*(xa/2)**3/(-E*Izz) - SCz**2*m.sin(alpha)*(xa/2)/(G*J),x2/(-E*Izz),1/(-E*Izz),0,0,SCz],                                               #v(x2)+theta(x2)
#        [0,-(x3-x1)**3/(-E*Izz)-SCz**2*(x3-x1)/(G*J),0,-(x3-x2)**3/(-E*Izz)-SCz**2*(x3-x2)/(G*J),0,0,-m.sin(alpha)*(x3-x2+0.5*xa)**3/(-E*Izz) - SCz**2*m.sin(alpha)*(x3-x2+xa*0.5)/(G*J),x3/(-E*Izz),1/(-E*Izz),0,0,SCz/(G*J)]]      #v(x3)+theta(x3)
    

# F=np.transpose([R1z,R1y,R2z,R2y,R3z,R3y,RI,C1,C2,C3,C4,C5])

# Bc= [[-SingleIntegralSC -SCz*P*m.sin(alpha)],               #T(la)
#       [-P*m.cos(alpha)*(La-x2-0.5*xa)],                     #My(la)
#       [-doubleIntegralA - P*m.sin(alpha)*(La-x2-0.5*xa)],    #Mz(la)
#       [-P*m.cos(alpha)],                                    #Sz(la)
#       [-SingleIntegralA - P*m.sin(alpha)],                   #Sy(la)
#       [-d1*m.sin(alpha)*E*Iyy],                             #w(x1)
#       [0],                                                  #w(x2)
#       [-d3*m.sin(alpha)*E*Iyy-P*m.cos(alpha)*(x3-x2-0.5*xa)**3],        #w(x3)
#       [0],                                                          #w(xI)
#       [d1*m.cos(alpha)+(4thIntegral/(E*Izz))-(SCz*DoubleintegralSC/(G*J))], #vertical deflection at x1
#       [4thintegral/(E*Izz) -SCz*DoubleinegrealSC/(G*J)],
#       [d3*m.cos(alpha)+(4thIntegral/(E*Izz))-(SCz*DoubleintegralSC/(G*J)) +P*m.sin(alpha)*(x3-x2-0.5*xa)**3/(E*Izz) -SCz**2 *P*m.sin(alpha)*(x3-x2-0.5xa)/(G*J)]]

#Torque X-axis
def Tx(x): return A_SC_int(x) - SCz*R1y*step(x,x1,0) - SCz*R1*math.sin(alpha)*step(x, x2-xa/2, 0) - SCz*R2y*step(x, x2, 0) + SCz*P*math.sin(alpha)*(x, x2 + xa/2, 0) - SCz*R3y*step(x, x3, 0)

#Moment in y-axis
def My(x): return R1z*step(x, x1, 1) - R1*math.cos(alpha)*(step(x, x2 - xa/2, 1)) + R2z*step(x, x2, 1) + P*math.cos(alpha)*step(x, x2 + xa/2, 1) + R3z*step(x,x3, 1)

#Moment in X-axis                                  
def Mz(x): return A_SC_doubleint(x) - R1y*step(x, x1, 1) - R1*math.sin(alpha)*(x, x2 - xa/2, 1) - R2y*step(x, x2, 1) + P*math.sin(alpha)*step(x, x2+xa/2, 1) - R3y*step(x, x3, 1)

#Shear y-axis
def Sy(x): return A_SC_int(x) - R1y*step(x, x1, 0) - R1*math.sin(alpha)*(x, x2 - xa/2, 0) - R2y*step(x, x2, 0) + P*math.sin(alpha)*step(x, x2 + xa/2, 0) - R3y*math.sin(alpha)*step(x, x3, 0)

#Shear Z-Axis
def Sz(x): return R1z*step(x, x1, 0) - R1*math.cos(alpha)*(step(x, x2 - xa/2, 0)) + R2z*step(x, x2, 0) + P*math.cos(alpha)*step(x, x2 + xa/2, 0) + R3z*step(x, x3, 0)

#Deflection in Y-axis
def v(x): return (-1/(E*Izz))*(A_quadint(x) - R1y*step(x, x1, 3) - R1*math.sin(alpha)*(x, x2 - xa/2, 3) - R2y*step(x, x2, 3) + P*math.sin(alpha)*step(x, x2+xa/2, 3) - R3y*step(x, x3, 3) + C1*step(x,0, 1) + C2)


def w(x): return (-1/(E*Iyy))*(R1z*step(x, x1, 3) - R1*math.cos(alpha)*(step(x, x2 - xa/2, 3)) + R2z*step(x, x2, 3) + P*math.cos(alpha)*step(x, x2 + xa/2, 3) + R3z*step(x,x3, 3) + C3*step(x,0, 1) + C4)

#Twist 
def theta(x): return (1/(G*J))*(A_SC_doubleint(x) - SCz*R1y*step(x,x1,1) - SCz*R1*math.sin(alpha)*step(x, x2-xa/2, 1) - SCz*R2y*step(x, x2, 1) + SCz*P*math.sin(alpha)*(x, x2 + xa/2, 1) - SCz*R3y*step(x, x3, 1) + C5)

