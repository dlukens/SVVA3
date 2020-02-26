# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:06:53 2020

@author: sarta
"""

# imports
import math as m
import matplotlib.pyplot as plt
from sectionproperties import z, y, Nst, h, Lssk, angle_Lssk, Ca, La, SCz, Izz, Iyy, Tsk, Tsp, Ast, A1, A2
import numpy as np
from tools import integral, dintegral
#---

#Make a numerical model of the cross section 
#Make a numerical model of the whole aileron

#################################################
#Correct stringer positions
z_st0, z_st1, z_st2, z_st3, z_st4, z_st5 = z
y_st0, y_st1, y_st2, y_st3, y_st4, y_st5 = y
striz=[z_st0,z_st1,z_st2,z_st3,z_st4,z_st5,z_st5,z_st4,z_st3,z_st2,z_st1]
striy=[y_st0,y_st1,y_st2,y_st3,y_st4,y_st5,-y_st5,-y_st4,-y_st3,-y_st2,-y_st1]

zstri=np.zeros(11) #Stringer z position in correct coordinate system
ystri=np.zeros(11) #Stringer y position in correct coordinate system

for i in range(Nst):
    zstri[i]=striz[i]+0.5*h
    ystri[i]=striy[i]

#####################################
#FUNCTIONS

# functions to relate theta or s to y and z coordinates
def fz01(angle):
    return 0.5*h*m.cos(angle)

def fy01(angle):    
    return 0.5*h*m.sin(angle)
    
def fz12(s):
    return -s*m.cos(angle_Lssk)

def fy12(s):
    return 0.5*h-s*m.sin(angle_Lssk)

def fz23(s):
    return -(Ca - 0.5*h)+s*m.cos(angle_Lssk)

def fy23(s):
    return -s*m.sin(angle_Lssk)

def fz30(angle):
    return 0.5*h*m.sin(angle)

def fy30(angle):
    return -0.5*h*m.cos(angle)

def fz34(s):
    return 0

def fy34(s):
    return -0.5*h+s
    
def fz41(s):
    return  0

def fy41(s):
    return s
#---
    
# step functions for the booms
def step01(theta,y_st1) :
    if fy01(theta) < y_st1 : return 0
    else: return 1
    

    

#####################################################
#Discretise aileron cross section

#Arc points 
arcsec = 20 # sections per quarter arc
dtheta = m.pi/(2*arcsec)
deltatheta = np.arange(0, m.pi/2 + dtheta, dtheta)
circpoint = np.size(deltatheta) # points per quarter arc

#Points on the arc 01
z01 = np.zeros(arcsec + 1)    
y01 = np.zeros(arcsec + 1)

for i in range(circpoint) :
    z01[i] = fz01(deltatheta[i])
    y01[i] = fy01(deltatheta[i])

linsec = 40    
dl = Lssk/linsec
l = np.arange(0, Lssk + dl, dl)
lpoints = np.size(l)

#Points on the straight 12
z12 = np.zeros(linsec + 1)
y12 = np.zeros(linsec + 1)

for i in range(lpoints) :
    z12[i] = fz12(l[i])
    y12[i] = fy12(l[i])
    
#Points on the straight 23
z23 = np.zeros(linsec + 1)
y23 = np.zeros(linsec + 1) 
   
for i in range(lpoints) :
    z23[i] = fz23(l[i])
    y23[i] = fy23(l[i])
    
#Points on the arc 30
z30 = np.zeros(arcsec + 1)
y30 = np.zeros(arcsec + 1)

for i in range(circpoint) :
    z30[i] = fz30(deltatheta[i])
    y30[i] = fy30(deltatheta[i])

#Points on the stiffener segment 34
stiffsec = 5
dlstiff = 0.5*h/stiffsec
lstiff = np.arange(0, 0.5*h + dlstiff, dlstiff)
stiffpoints = np.size(lstiff)

z34 = np.zeros(stiffsec+1)
y34 = np.zeros(stiffsec+1)

for i in range(stiffpoints) :
    z34[i] = fz34(lstiff[i])
    y34[i] = fy34(lstiff[i])

#Points on the stiffener segment 41    
z41 = np.zeros(stiffsec + 1)
y41 = np.zeros(stiffsec + 1)
for i in range(stiffpoints) :
    z41[i] = fz41(lstiff[i])
    y41[i] = fy41(lstiff[i])
    
##################################################################    
# Discretise the span of the aileron
beginnode = 0.001
endnode = La - 0.002
nodesnumber = 100
deltax = (endnode - beginnode)/(nodesnumber)

##################################################################
#Plot the aileron
#plt.scatter(z01,y01,color='blue',label="nodes")
#plt.scatter(z12,y12,color='blue')
#plt.scatter(z23,y23,color='blue')
#plt.scatter(z30,y30,color='blue')
#plt.scatter(z34,y34,color='blue')
#plt.scatter(z41,y41,color='blue')
#plt.scatter(zstri,ystri,color='red',marker="^",label="stringer locations")
#plt.legend(), plt.show()

#############################################################
#Shear flow calculation for a specified cross-section (specified x-coordinate)
Vy, Vz = 1, 1        # [N]
T = 1                # [Nm]
eta = abs(SCz + h/2) # [m] -> distance spar-shear center

# base shear flow gradients without booms
def d_qb01(theta,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsk*fy01(theta)) - (Vz/Iyy)*(Tsk*fz01(theta))
    return dq

def d_qb12(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsk*fy12(s)) - (Vz/Iyy)*(Tsk*fz12(s))
    return dq

def d_qb23(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsk*fy23(s)) - (Vz/Iyy)*(Tsk*fz23(s))
    return dq

def d_qb30(theta,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsk*fy30(theta)) - (Vz/Iyy)*(Tsk*fz30(theta))
    return dq

def d_qb34(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsp*fy34(s)) - (Vz/Iyy)*(Tsp*fz34(s))
    return dq

def d_qb41(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsp*fy41(s)) - (Vz/Iyy)*(Tsp*fz41(s))
    return dq

# base shear flows with booms

qb01 = integral(d_qb01, 0, m.pi/2) - Ast*((Vy/Izz)*striy[1] + (Vz/Iyy)*striz[1])                      # [N/m]
qb12 = integral(d_qb12, 0, Lssk) - Ast*((Vy/Izz)*sum(striy[2:6]) + (Vz/Iyy)*sum(striz[2:6])) + qb01   # [N/m]
qb23 = integral(d_qb23, 0, Lssk) - Ast*((Vy/Izz)*sum(striy[6:10]) + (Vz/Iyy)*sum(striz[6:10])) + qb12 # [N/m]
qb30 = integral(d_qb30, 0, m.pi/2) - Ast*((Vy/Izz)*striy[-1] + (Vz/Iyy)*striz[-1]) + qb23             # [N/m]
qb34 = integral(d_qb34, 0, h/2) + qb23                                                                # [N/m]
qb41 = integral(d_qb41, 0, h/2) + qb34                                                                # [N/m]
print(qb01)
print(qb12)
print(qb23)
print(qb30)
print(qb34)
print(qb41)

# redundant shear flows (defined counterclockwise positive)

X = np.zeros((3,3))
x1, x2, x3 = 2*A1, 2*A2, 0
x4, x5, x6 = 0, 0, 0
x7, x8, x9 = 0, 0, 0

Y = np.zeros((3,1))
y1 = Vy*eta + T
y2 = 0
y3 = 0






