# -*- coding: utf-8 -*-
import math as m
from sectionproperties import h, SCz, x1, x2, x3, xa, P, G, J, E, Izz, Iyy, La, d1, d3
import numpy as np
from forces import A_coeff, Cp_coeff
import interp
from interp import grid_x
import tools as t


###variables
alpha = m.radians(25) # angle of attack
eta = SCz + h/2       # distance midpoint spar - shear center
xI=x2-0.5*xa
xII=x2+0.5*xa
###Funcions

def step(x, x1, exp): #Mcaulay step function
    y = x - x1
    if y <= 0: return 0 
    else: return y**exp
    

    
def A_SC_int(x):
    idx = np.searchsorted(grid_x, x)-1
    if x > grid_x[-1]:  #edge cases
        idx = len(grid_x)-2 
        x = grid_x[-1]
    if x < grid_x[0]: 
        idx = 0
        x = grid_x[0]
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
    if x > grid_x[-1]: #edge cases
        idx = len(grid_x)-2 
        x = grid_x[-1]
    if x < grid_x[0]: 
        idx = 0
        x = grid_x[0]
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
    if x > grid_x[-1]:  #edge cases
        idx = len(grid_x)-2 
        x = grid_x[-1]
    if x < grid_x[0]: 
        idx = 0
        x = grid_x[0]
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx]    
    return x*(a/4*(x-grid_x[idx])**3      \
            + b/3*(x-grid_x[idx])**2      \
            + c/2*(x-grid_x[idx])         \
            + d)

S = A_int(4)
    
def A_doubleint(x):
    idx = np.searchsorted(grid_x, x)-1
    if x > grid_x[-1]:  #edge cases
        idx = len(grid_x)-2 
        x = grid_x[-1]
    if x < grid_x[0]: 
        idx = 0
        x = grid_x[0]
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
    if x > grid_x[-1]:  #edge cases
        idx = len(grid_x)-2 
        x = grid_x[-1]
    if x < grid_x[0]: 
        idx = 0
        x = grid_x[0]
    a = A_coeff[0, idx]
    b = A_coeff[1, idx]
    c = A_coeff[2, idx]
    d = A_coeff[3, idx] 
    return x**4*(a/(4*5*6*7)*(x-grid_x[idx])**3     \
               + b/(3*4*5*6)*(x-grid_x[idx])**2     \
               + c/(2*3*4*5)*(x-grid_x[idx])        \
               + d/(2*3*4))

    
#Reaction forces
    
Rxn= [[0,-eta,0,-eta,0,-eta,(m.sin(alpha)*(eta+0.5*h) +m.cos(alpha)*0.5*h)*(La-xI),0,0,0,0,0],                                                                                                   #T(la)
        [La-x1,0,La-x2,0,La-x3,0,-m.cos(alpha)*(La-xI),0,0,0,0,0],                                                                                     #My(la)
        [0,(La-x1),0,(La-x2),0,(La-x3),-m.sin(alpha)*(La-xI),0,0,0,0,0],                                                                              #Mz(la)
        [1,0,1,0,1,0,-m.cos(alpha),0,0,0,0,0],                                                                                                                #Sz(la)
        [0,1,0,1,0,1,-m.sin(alpha),0,0,0,0,0],                                                                                                             #Sy(la)
        [0,0,0,0,0,0,0,0,0,x1,1,0],                                                                                                                        #w(x1)
        [(x2-x1)**3,0,0,0,0,0,-m.cos(alpha)*(x2-xI)**3,0,0,x2,1,0],                                                                                         #w(x2)
        [(x3-x1)**3,0,(x3-x2)**3,0,0,0,-m.cos(alpha)*(x3-xI)**3,0,0,x3,1,0],                                                                                                      #w(x3)
        [-m.cos(alpha)*(xI-x1)**3/(6*E*Iyy),m.sin(alpha)*((xI-x1)**3/(6*E*Izz) +eta**2*(xI-x1)**3/(G*J)),0,0,0,0,0,-xI*m.sin(alpha),-m.sin(alpha),xI*m.cos(alpha),m.cos(alpha),-eta*m.sin(alpha)],                    #w'(xI)=0
        [0,0,0,0,0,0,0,x1,1,0,0,eta],                                                                                                                           #v(x1)+theta(x1)
        [0,-(x2-x1)**3/(6*E*Izz)-eta**2*(x2-x1)/(G*J),0,0,0,0,m.sin(alpha)*(x2-xI)**3/(6*E*Izz) + eta**2*m.sin(alpha)*(x2-xI)/(G*J),x2,1,0,0,eta],                                               #v(x2)+theta(x2)
        [0,-(x3-x1)**3/(6*E*Izz)-eta**2*(x3-x1)/(G*J),0,-(x3-x2)**3/(6*E*Izz)-eta**2*(x3-x2)/(G*J),0,0,m.sin(alpha)*(x3-xI)**3/(6*E*Izz) - eta**2*m.sin(alpha)*(x3-xI)/(G*J),x3,1,0,0,eta]]      #v(x3)+theta(x3)
  

Bc= [[-A_SC_int(La) - eta*P*m.sin(alpha)],                              #T(la) 
       [P*m.cos(alpha)*(La - x2 - 0.5*xa)],                             #My(la)
       [A_SC_doubleint(La) + P*m.sin(alpha)*(La-x2 - 0.5*xa)],          #Mz(la)
       [P*m.cos(alpha)],                                                #Sz(la)
       [A_int(La) + P*m.sin(alpha)],                                    #Sy(la)
       [d1*m.sin(alpha)],                                               #w(x1)
       [0],                                                             #w(x2)
       [-d3*m.sin(alpha)*E*Iyy - P*m.cos(alpha)*(x3 - x2 - 0.5*xa)**3], #w(x3)
       [0],                                                             #w(xI)
       [d1*m.cos(alpha) - (A_quadint(x1)/(6*E*Izz)) - (eta*A_SC_doubleint(x1)/(G*J))], #vertical deflection at x1
       [-A_quadint(x2)/(6*E*Izz) - eta*A_SC_doubleint(x2)/(G*J)],
       [d3*m.cos(alpha) + (A_quadint(x3)/(E*Izz)) - (eta*A_SC_doubleint(x3)/(G*J)) + P*m.sin(alpha)*(x3 - x2 - 0.5*xa)**3/(E*Izz) - eta**2 *P*m.sin(alpha)*(x3 - x2 - 0.5*xa)/(G*J)]]


F = np.linalg.solve(Rxn, Bc)

# #F=np.transpose([R1z,R1y,R2z,R2y,R3z,R3y,RI,C1,C2,C3,C4,C5])

# #Torque X-axis
# def Tx(x): return A_SC_int(x) - eta*R1y*step(x,x1,0) - eta*R1*m.sin(alpha)*step(x, x2-xa/2, 0) - eta*R2y*step(x, x2, 0) + eta*P*m.sin(alpha)*(x, x2 + xa/2, 0) - eta*R3y*step(x, x3, 0)

# #Moment in y-axis
# def My(x): return R1z*step(x, x1, 1) - R1*m.cos(alpha)*(step(x, x2 - xa/2, 1)) + R2z*step(x, x2, 1) + P*m.cos(alpha)*step(x, x2 + xa/2, 1) + R3z*step(x,x3, 1)

# #Moment in X-axis                                  
# def Mz(x): return A_SC_doubleint(x) - R1y*step(x, x1, 1) - R1*m.sin(alpha)*(x, x2 - xa/2, 1) - R2y*step(x, x2, 1) + P*m.sin(alpha)*step(x, x2+xa/2, 1) - R3y*step(x, x3, 1)

# #Shear y-axis
# def Sy(x): return A_SC_int(x) - R1y*step(x, x1, 0) - R1*m.sin(alpha)*(x, x2 - xa/2, 0) - R2y*step(x, x2, 0) + P*m.sin(alpha)*step(x, x2 + xa/2, 0) - R3y*m.sin(alpha)*step(x, x3, 0)

# #Shear Z-Axis
# def Sz(x): return R1z*step(x, x1, 0) - R1*m.cos(alpha)*(step(x, x2 - xa/2, 0)) + R2z*step(x, x2, 0) + P*m.cos(alpha)*step(x, x2 + xa/2, 0) + R3z*step(x, x3, 0)

# #Deflection in Y-axis
# def v(x): return (-1/(E*Izz))*(A_quadint(x) - R1y*step(x, x1, 3) - R1*m.sin(alpha)*(x, x2 - xa/2, 3) - R2y*step(x, x2, 3) + P*m.sin(alpha)*step(x, x2+xa/2, 3) - R3y*step(x, x3, 3) + C1*step(x,0, 1) + C2)

# #Deflection in Z-axis
# def w(x): return (-1/(E*Iyy))*(R1z*step(x, x1, 3) - R1*m.cos(alpha)*(step(x, x2 - xa/2, 3)) + R2z*step(x, x2, 3) + P*m.cos(alpha)*step(x, x2 + xa/2, 3) + R3z*step(x,x3, 3) + C3*step(x,0, 1) + C4)

# #Twist 
# def theta(x): return (1/(G*J))*(A_SC_doubleint(x) - eta*R1y*step(x,x1,1) - eta*R1*m.sin(alpha)*step(x, x2-xa/2, 1) - eta*R2y*step(x, x2, 1) + eta*P*m.sin(alpha)*(x, x2 + xa/2, 1) - eta*R3y*step(x, x3, 1) + C5)