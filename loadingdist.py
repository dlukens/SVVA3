# -*- coding: utf-8 -*-
import math
import tools as t


###variables
alpha = math.radians(25)

###Funcions

def step(x, x1, exp): #Mcaulay step function
    y = x - x1
    if y <= 0:
        return 0 
    else:
        return y**exp
    
def A(x):
    pass
    
    
###Reaction forces
    
#Torque X-axis
Tx = integral(Ax*(SCz - Cp_x), dx) - SCz*R1y*step(x,x1,0) - SCz*R1*math.sin(alpha)*step(x, x2+xa/2, 0) - SCz*R2y*step(x, x2, 0) + SCz*P*math.sin(alpha)*(x, x2 - xa/2, 0) + SCz*R3y*step(x, x3, 0)


#Moment in y-axis
My = R1z*step(x, x1, 1) - R1*math.cos(alpha)*(step(x, x2 + xa/2, 1)) + R2z*step(x, x2, 1) + P*math.cos(alpha)*step(x, x2 - xa/2, 1) + R3z*step(x,x3, 1)

#Shear Z-Axis
Sz = R1z*step(x, x1, 0) - R1*math.cos(alpha)*(step(x, x2 + xa/2, 0)) + R2z*step(x, x2, 0) + P*math.cos(alpha)*step(x, x2 - xa/2, 0) + R3z*step(x,x3, 0)


#Moment in X-axis                                  
Mz = doubleintegrate(Ax, dx) - R1y*step(x, x1, 1) - R1*math.sin(alpha)*(x, x2 + xa/2, 1) + R2y*step(x, x2, 1) + P*math.sin(alpha)*step(x, x2-xa/2, 1) - R3y*math.sin(alpha)*step(x, x3, 1)

#Shear y-axis
Sy = integrate(Ax, dx) - R1y*step(x, x1, 0) - R1*math.sin(alpha)*(x, x2 + xa/2, 0) + R2y*step(x, x2, 0) + P*math.sin(alpha)*step(x, x2-xa/2, 0) - R3y*math.sin(alpha)*step(x, x3, 0)


