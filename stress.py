# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:06:53 2020

@author: sarta
"""
import math as m
import matplotlib as plt
from sectionproperties import *
import numpy as np


#Make a numerical model of the cross section 
#Make a numerical model of the whole aileron





Npoints=100                     #Number of points on half of Aileron cross section
delta= Lsk/(2*(Npoints-1))      #difference between successive points
circfrac= 0.5*m.pi*h/Lsk        #fraction of total skin lenght that comes from arc
circpoint=circfrac*Npoints      #Number of points on the arc
deltacirc=m.pi/(2*circpoint)    #

circpoints=int(circpoint//1)

zpoints=np.zeros(2*Npoints-1)
ypoints=np.zeros(2*Npoints-1)

for i in range(circpoints):
    zpoints[i]=-0.124+0.5*h*m.cos(i*deltacirc)
    
    ypoints[i]=0.5*h*m.sin(i*deltacirc)

remcirc=0.5*m.pi-circpoints*deltacirc 
lindif=delta-remcirc*h*0.5

zpoints[circpoints]= -0.5*h-lindif*m.cos(angle_Lssk)
ypoints[circpoints]=0.5*h-lindif*m.sin(angle_Lssk)

for i in range(circpoints+1,Npoints):
    zpoints[i]=zpoints[i-1]-delta*m.cos(angle_Lssk)
    ypoints[i]=ypoints[i-1]-delta*m.sin(angle_Lssk)
    
for i in range( Npoints-1):
    zpoints[Npoints+i]=zpoints[Npoints-2-i]
    ypoints[Npoints+i]=-ypoints[Npoints-2-i]
    
    
beginnode=0.01
endnode=2.01
nodesnumber=100
deltax=(endnode-beginnode)/(nodesnumber-1)

bendingstress=np.zeros(nodesnumber,2*Npoints-1)

for i in range(nodesnumber):
    x=beginnode+deltax*i
    for j in range(2*Npoints-1):
        sigma = My(x)*z/Iyy + Mz(x)*y/Izz
        bendingstress(i,j)=sigma
        
        



