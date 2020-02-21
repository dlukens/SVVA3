# -*- coding: utf-8 -*-
import math
import tools as t


###variables
alpha = math.radians(25)

###Funcions

def step(x, x1): #Mcaulay step function
    y = x - x1
    if y <= 0:
        return 0 
    else:
        return y
    
def A(x):
    pass
    
    
###Reaction forces


My = R1z*step(x, x1) - R1*math.cos(alpha)*(step(x, x2 + xa/2)) + 
                                    R2z*step(x, x2) + P*math.cos(alpha)*step(x, x2 - xa/2) + R3z*step(x,x3)
                                    
Mz = t.integral(A, )