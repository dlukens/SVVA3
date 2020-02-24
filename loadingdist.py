# -*- coding: utf-8 -*-
import math
from sectionproperties import *
import forces
from forces import A_coeff, Cp_coeff
from interp import grid_x, grid_z


###variables
alpha = math.radians(25)

###Funcions

def step(x, x1, exp): #Mcaulay step function
    y = x - x1
    if y <= 0:
        return 0 
    else:
        return y**exp
    
def integrator(f, grid):
    I = np.zeros(len(grid))
    for i in range(len(grid)):
        I[i] = abs(f(i)-f(i-1))
    return I
        

def A_SC_int(x, A_coeff = A_coeff, SCz = SCz, Cp_coeff = Cp_coeff):
    return A_coeff[0, x]*SCz/4*x^4 + A_coeff[1, x]*SCz/3*x**3 + A_coeff[3, x]*SCz/2*x**2 +  A_coeff[4, x]*SCz*x**2 \
                        -  A_coeff[0, x]/7*x**7 \
                        - (A_coeff[0, x]*Cp_coeff[1, x] + A_coeff[1, x]*Cp_coeff[0, x])/6*x**6 \
                        - (A_coeff[0, x]*Cp_coeff[2, x] + A_coeff[1, x]*Cp_coeff[1, x] + A_coeff[2, x]*Cp_coeff[0, x])/5*x**5 \
                        - (A_coeff[0, x]*Cp_coeff[3, x] + A_coeff[1, x]*Cp_coeff[2, x] + A_coeff[2, x]*Cp_coeff[1, x] + A_coeff[3, x]*Cp_coeff[0, x])/4*x**4 \
                        - (A_coeff[1, x]*Cp_coeff[3, x] + A_coeff[2, x]*Cp_coeff[2, x] + A_coeff[3, x]*Cp_coeff[1, x])/3*x**3 \
                        - (A_coeff[2, x]*Cp_coeff[3, x] + A_coeff[3, x]*Cp_coeff[2, x])/2*x**2 \
                        -  A_coeff[3, x]*Cp_coeff[3, x]*x
         
P = integrator(A_SC_int, grid_x)               
                        
def A_SC_doubleint(x, A_coeff, SCz, Cp_coeff):
    return A_coeff[0, x]*SCz/(4*5)*x**5 + A_coeff[1, x]*SCz/(3*4)*x**4 + A_coeff[3, x]*SCz/(2*3)*x**3 +  A_coeff[4, x]*SCz/(2)*x**2 \
                        -  A_coeff[0, x]/(8*7)*x**8 \
                        - (A_coeff[0, x]*Cp_coeff[1, x] + A_coeff[1, x]*Cp_coeff[0, x])/(6*7)*x**7 \
                        - (A_coeff[0, x]*Cp_coeff[2, x] + A_coeff[1, x]*Cp_coeff[1, x] + A_coeff[2, x]*Cp_coeff[0, x])/(5*6)*x**6 \
                        - (A_coeff[0, x]*Cp_coeff[3, x] + A_coeff[1, x]*Cp_coeff[2, x] + A_coeff[2, x]*Cp_coeff[1, x] + A_coeff[3, x]*Cp_coeff[0, x])/(4*5)*x**5 \
                        - (A_coeff[1, x]*Cp_coeff[3, x] + A_coeff[2, x]*Cp_coeff[2, x] + A_coeff[3, x]*Cp_coeff[1, x])/(3*4)*x**4 \
                        - (A_coeff[2, x]*Cp_coeff[3, x] + A_coeff[3, x]*Cp_coeff[2, x])/(2*3)*x**3 \
                        - A_coeff[3, x]*Cp_coeff[3, x]/2*x**2
    
    
def A_int(x, A_coeff):
    return A_coeff[0, x]/4*x**4         \
            + A_coeff[1, x]/3*x**3      \
            + A_coeff[2, x]/2*x**2      \
            + A_coeff[3, x]*x**1
            
    
def A_doubleint(x, A_coeff):
    return A_coeff[0, x]/(4*5)*x**5     \
            + A_coeff[1, x]/(3*4)*x**4  \
            + A_coeff[2, x]/(2*3)*x**3  \
            + A_coeff[3, x]/2*x**2 

def A_quadint(x, A_coeff):
    return A_coeff[0, x]/(4*5*6*7)*x**7     \
            + A_coeff[1, x]/(3*4*5*6)*x**6  \
            + A_coeff[2, x]/(2*3*4*5)*x**5  \
            + A_coeff[3, x]/(2*3*4)*x**4

    
    
###Reaction forces
    

# Rxn= [[0,-SCz,0,-SCz,0,-SCz,-SCz*m.sin(alpha),0,0,0,0,0],                                                                                                   #T(la)
#       [La-x1,0,La-x2,0,La-x3,0,-m.cos(alpha)*(La-x2+0.5*xa),0,0,0,0,0],                                                                                     #My(la)
#       [0,-(La-x1),0,(La-x2),0,(La-x3),-m.sin(alpha)*(La-x2+xa*0.5),0,0,0,0,0],                                                                              #Mz(la)
#       [1,0,1,0,1,0,-m.cos(alpha),0,0,0,0,0],                                                                                                                #Sz(la)
#       [0,-1,0,-1,0,-1,-m.sin(alpha),0,0,0,0,0],                                                                                                             #Sy(la)
#       [0,0,0,0,0,0,0,0,0,x1**3,1,0],                                                                                                                        #w(x1)
#       [(x2-x1)**3,0,0,0,0,0,-m.cos(alpha)*(xa/2)**3,0,0,x2**3,1,0],                                                                                         #w(x2)
#       [(x3-x1)**3,0,(x3-x2)**3,0,0,0,0,0,0,x3**3,1,0],                                                                                                      #w(x3)
#       [0,0,0,0,0,0,0,0,0,0,0,1],                                                                                                                            #theta(0)
#       [0,0,0,0,0,0,0,x1/(-E*Izz),1/(-E*Izz),0,0,SCz/(G*J)],                                                                                                                           #v(x1)+theta(x1)
#       [0,-(x2-x1)**3/(-E*Izz)-SCz**2*(x2-x1)/(G*J),0,0,0,0,-m.sin(alpha)*(xa/2)**3/(-E*Izz) - SCz**2*m.sin(alpha)*(xa/2)/(G*J),x2/(-E*Izz),1/(-E*Izz),0,0,SCz],                                               #v(x2)+theta(x2)
#       [0,-(x3-x1)**3/(-E*Izz)-SCz**2*(x3-x1)/(G*J),0,-(x3-x2)**3/(-E*Izz)-SCz**2*(x3-x2)/(G*J),0,0,-m.sin(alpha)*(x3-x2+0.5*xa)**3/(-E*Izz) - SCz**2*m.sin(alpha)*(x3-x2+xa*0.5)/(G*J),x3/(-E*Izz),1/(-E*Izz),0,0,SCz/(G*J)]]      #v(x3)+theta(x3)
    


Bc= []


#Torque X-axis
def Tx(x): return A_SC_doubleint(x, A_coeff, SCz, Cp_coeff) - SCz*R1y*step(x,x1,0) - SCz*RI*math.sin(alpha)*step(x, x2+xa/2, 0) - SCz*R2y*step(x, x2, 0) + SCz*P*math.sin(alpha)*(x, x2 - xa/2, 0) - SCz*R3y*step(x, x3, 0)


#Moment in y-axis
def My(x): return R1z*step(x, x1, 1) - R1*math.cos(alpha)*(step(x, x2 + xa/2, 1)) + R2z*step(x, x2, 1) + P*math.cos(alpha)*step(x, x2 - xa/2, 1) + R3z*step(x,x3, 1)

#Shear Z-Axis
def Sz(x): return R1z*step(x, x1, 0) - R1*math.cos(alpha)*(step(x, x2 + xa/2, 0)) + R2z*step(x, x2, 0) + P*math.cos(alpha)*step(x, x2 - xa/2, 0) + R3z*step(x,x3, 0)


#Moment in X-axis                                  
def Mz(x): return doubleintegrate(Ax, dx) - R1y*step(x, x1, 1) - RI*math.sin(alpha)*(x, x2 + xa/2, 1) - R2y*step(x, x2, 1) + P*math.sin(alpha)*step(x, x2-xa/2, 1) - R3y*step(x, x3, 1)

#Shear y-axis
def Sy(x): return integrate(Ax, dx) - R1y*step(x, x1, 0) - R1*math.sin(alpha)*(x, x2 + xa/2, 0) + R2y*step(x, x2, 0) + P*math.sin(alpha)*step(x, x2-xa/2, 0) - R3y*math.sin(alpha)*step(x, x3, 0)


