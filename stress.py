# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:06:53 2020

@author: sarta
"""
import math as m
import matplotlib.pyplot as plt
from sectionproperties import *
import numpy as np
from tools import *

#Make a numerical model of the cross section 
#Make a numerical model of the whole aileron



#################################################
#Correct stringer positions
striz=[z_st0,z_st1,z_st2,z_st3,z_st4,z_st5,z_st5,z_st4,z_st3,z_st2,z_st1,z_st0]
striy=[y_st0,y_st1,y_st2,y_st3,y_st4,y_st5,-y_st5,-y_st4,-y_st3,-y_st2,-y_st1,y_st0]

zstri=np.zeros(11) #Stringer z position in correct coordinate system
ystri=np.zeros(11) #Stringer y position in correct coordinate system

for i in range(Nst):
    zstri[i]=striz[i]+0.5*h
    ystri[i]=striy[i]



#####################################
#FUNCTIONS
def fz01(angle):
    return 0.5*h*m.cos(angle)

def fy01(angle):    
    return 0.5*h*m.sin(angle)
    

def fz12(length):
    return -length*m.cos(angle_Lssk)

def fy12(length):
    return 0.5*h-length*m.sin(angle_Lssk)


def fz23(length):
    return -(Ca-0.5*h)+length*m.cos(angle_Lssk)

def fy23(length):
    return -length*m.sin(angle_Lssk)

def fz30(angle):
    return 0.5*h*m.sin(angle)

def fy30(angle):
    return-0.5*h*m.cos(angle)

def fz34(length):
    return 0

def fy34(length):
    return -0.5*h+length
    
def fz41(length):
    return  0

def fy41(length):
    return length
#####################################################
#Discretise aileron cross section

#Arc points 
arcsec=20
dtheta=m.pi/(2*arcsec)
deltatheta=np.arange(0,m.pi/2+dtheta,dtheta)
circpoint=np.size(deltatheta)

#Points on the arc 01
z01=np.zeros(arcsec+1)    
y01=np.zeros(arcsec+1)
for i in range(circpoint):
    z01[i]=fz01(deltatheta[i])
    y01[i]=fy01(deltatheta[i])

linsec=40    
dl=Lssk/linsec
l=np.arange(0,Lssk+dl,dl)
lpoints=np.size(l)
#Points on the straight 12
z12=np.zeros(linsec+1)
y12=np.zeros(linsec+1)
for i in range(lpoints):
    z12[i]=fz12(l[i])
    y12[i]=fy12(l[i])
    
#Points on the straight 23
z23=np.zeros(linsec+1)
y23=np.zeros(linsec+1)    
for i in range(lpoints):
    z23[i]=fz23(l[i])
    y23[i]= fy23(l[i])
    

#Points on the arc 30
z30=np.zeros(arcsec+1)
y30=np.zeros(arcsec+1)

for i in range(circpoint):
    z30[i]=fz30(deltatheta[i])
    y30[i]=fy30(deltatheta[i])

#Points on the stiffener segment 34
stiffsec=5
dlstiff=0.5*h/stiffsec
lstiff=np.arange(0,0.5*h+dlstiff,dlstiff)
stiffpoints=np.size(lstiff)

z34=np.zeros(stiffsec+1)
y34=np.zeros(stiffsec+1)

for i in range(stiffpoints):
    z34[i]=fz34(lstiff[i])
    y34[i]=fy34(lstiff[i])

#Points on the stiffener segment 41    
z41=np.zeros(stiffsec+1)
y41=np.zeros(stiffsec+1)
for i in range(stiffpoints):
    z41[i]= fz41(lstiff[i])
    y41[i]= fy41(lstiff[i])

###################################################################
#Calculate the shear flow distribtuon in the aileron 
#q01z2int=np.zeros(circpoint-1)        #Integrated value for each section due to z

#for i in range(1,circpoint-1):
    #q01z2int[0]=0.5*h*Tsk*dintegration
    #q01z2int[i]=0.5*h*Tsk*dintegration+q01z2int[i-1]
    
##################################################################    
#Discretise the span of the aileron
beginnode=0.001
endnode=La-0.002
nodesnumber=100
deltax=(endnode-beginnode)/(nodesnumber)




        


plt.scatter(z01,y01,color='blue')
plt.scatter(z12,y12,color='blue')
plt.scatter(z23,y23,color='blue')
plt.scatter(z30,y30,color='blue')
plt.scatter(z34,y34,color='blue')
plt.scatter(z41,y41,color='blue')
plt.scatter(zstri,ystri,color='red')
plt.show





