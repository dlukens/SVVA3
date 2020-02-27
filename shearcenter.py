# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 17:31:59 2020

@author: ellio
"""

# imports
import math as m
from sectionproperties import z, y, Nst, h, Lssk, angle_Lssk, Ca, Izz, Iyy, Tsk, Tsp, Ast, delta_str
import numpy as np
from tools import integral, integral2
#---

#################################################
#Correct stringer positions
z_st0, z_st1, z_st2, z_st3, z_st4, z_st5 = z
y_st0, y_st1, y_st2, y_st3, y_st4, y_st5 = y
striz = [z_st0,z_st1,z_st2,z_st3,z_st4,z_st5,z_st5,z_st4,z_st3,z_st2,z_st1]
striy = [y_st0,y_st1,y_st2,y_st3,y_st4,y_st5,-y_st5,-y_st4,-y_st3,-y_st2,-y_st1]

zstri = np.zeros(11) #Stringer z position in correct coordinate system
ystri = np.zeros(11) #Stringer y position in correct coordinate system

for i in range(Nst):
    zstri[i] = striz[i]+0.5*h
    ystri[i] = striy[i]

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
   
# functions for the distances along which the booms act
def dB01(zstri) :
    theta = m.acos(2*zstri[1]/h)
    dist = (m.pi/2 - theta)*h/2
    return dist

def dB12(zstri,ystri) :
    ze, ye = -(Ca - h/2), 0
    d1 = m.sqrt((zstri[2] - ze)**2 + (ystri[2] - ye)**2)
    d2 = m.sqrt((zstri[3] - ze)**2 + (ystri[3] - ye)**2)
    d3 = m.sqrt((zstri[4] - ze)**2 + (ystri[4] - ye)**2)
    d4 = m.sqrt((zstri[5] - ze)**2 + (ystri[5] - ye)**2)
    return np.array([d1, d2, d3, d4])

def dB23(zstri,ystri) :
    ze, ye = 0, -h/2
    d1 = m.sqrt((zstri[6] - ze)**2 + (ystri[6] - ye)**2)
    d2 = m.sqrt((zstri[7] - ze)**2 + (ystri[7] - ye)**2)
    d3 = m.sqrt((zstri[8] - ze)**2 + (ystri[8] - ye)**2)
    d4 = m.sqrt((zstri[9] - ze)**2 + (ystri[9] - ye)**2)
    return np.array([d1, d2, d3, d4])
    
def dB30(zstri) :
    theta = m.asin(2*zstri[-1]/h)
    dist = (m.pi/2 - theta)*h/2
    return dist
#---
    
# step functions for the booms
def stepz01(theta,zstri,oint) :
    if theta < m.acos(2*zstri[1]/h) : return 0
    elif theta > m.acos(2*zstri[1]/h) and not oint : return 1
    else : return dB01(zstri)

def stepy01(theta,ystri,oint) :
    if theta < m.asin(2*ystri[1]/h) : return 0
    elif theta > m.cos(2*ystri[1]/h) and not oint : return 1
    else : return dB01(zstri)

def stepz12(s,zstri,oint) :
    if   s < -zstri[5]/m.cos(angle_Lssk) : return np.array([0,0,0,0])
    elif s > -zstri[5]/m.cos(angle_Lssk) and not oint : return np.array([1,1,1,1])
    else : return dB12(zstri,ystri)
    
def stepy12(s,ystri,oint) :
    if   s < (h/2 -ystri[5])/m.sin(angle_Lssk) : return np.array([0,0,0,0])
    elif s > (h/2 -ystri[5])/m.sin(angle_Lssk) and not oint : return np.array([1,1,1,1])
    else : return dB12(zstri,ystri)
    
def stepz23(s,zstri,oint) :
    if   s < (zstri[9] + Ca - h/2)/m.cos(angle_Lssk) : return np.array([0,0,0,0])
    elif s > (zstri[9] + Ca - h/2)/m.cos(angle_Lssk) and not oint : return np.array([1,1,1,1])
    else : return dB23(zstri,ystri)
    
def stepy23(s,ystri,oint) :
    if   s < -ystri[9]/m.sin(angle_Lssk) : return np.array([0,0,0,0])
    elif s > -ystri[9]/m.sin(angle_Lssk) and not oint : return np.array([1,1,1,1])
    else : return dB23(zstri,ystri)

def stepz30(theta,zstri,oint) :
    if theta < m.asin(2*zstri[-1]/h) : return 0
    elif theta > m.asin(2*zstri[-1]/h) and not oint : return 1
    else : return dB30(zstri)
    
def stepy30(theta,ystri,oint) :
    if theta < m.acos(-2*ystri[-1]/h) : return 0    
    elif theta > m.acos(-2*ystri[-1]/h) and not oint : return 1
    else : return dB30(zstri)
#---

#############################################################
#Shear flow calculation for a specified cross-section (specified x-coordinate)
Vy, Vz = 1, 0               # [N]
d = 0.5*h*m.cos(angle_Lssk) # [m]
Delta = delta_str
A1 = 0.5*m.pi*(h/2)**2      # area of cell 1
A2 = (Ca - h/2)*h/2         # area of cell 2

# base shear flow gradients without booms
def d_qb01(theta,Vy=Vy,Vz=Vz) :
    dq = -(Vy*h/(2*Izz))*(Tsk*fy01(theta)) - (Vz*h/(2*Iyy))*(Tsk*fz01(theta))
    return dq

def d_qb12(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsk*fy12(s)) - (Vz/Iyy)*(Tsk*fz12(s))
    return dq

def d_qb23(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsk*fy23(s)) - (Vz/Iyy)*(Tsk*fz23(s))
    return dq

def d_qb30(theta,Vy=Vy,Vz=Vz) :
    dq = -(Vy*h/(2*Izz))*(Tsk*fy30(theta)) - (Vz*h/(2*Iyy))*(Tsk*fz30(theta))
    return dq

def d_qb34(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsp*fy34(s)) - (Vz/Iyy)*(Tsp*fz34(s))
    return dq

def d_qb41(s,Vy=Vy,Vz=Vz) :
    dq = -(Vy/Izz)*(Tsp*fy41(s)) - (Vz/Iyy)*(Tsp*fz41(s))
    return dq
#---
    
# shear flows due to the booms
def qB01(theta,oint,Vy=Vy,Vz=Vz) :
    q = -Ast*((Vy/Izz)*stepy01(theta,ystri,oint)*ystri[1] + (Vz/Iyy)*stepz01(theta,zstri,oint)*zstri[1])
    return q

def qB12(s,oint,Vy=Vy,Vz=Vz) :
    yarr = np.array(ystri[2:6]).reshape((4,1))
    zarr = np.array(zstri[2:6]).reshape((4,1))
    q = -Ast*((Vy/Izz)*np.dot(stepy12(s,ystri,oint),yarr)[0] + (Vz/Iyy)*np.dot(stepz12(s,zstri,oint),zarr)[0])
    return q

def qB23(s,oint,Vy=Vy,Vz=Vz) :
    yarr = np.array(ystri[6:10]).reshape((4,1))
    zarr = np.array(zstri[6:10]).reshape((4,1))
    q = -Ast*((Vy/Izz)*np.dot(stepy23(s,ystri,oint),yarr)[0] + (Vz/Iyy)*np.dot(stepz23(s,zstri,oint),zarr)[0])
    return q

def qB30(theta,oint,Vy=Vy,Vz=Vz) :
    q = -Ast*((Vy/Izz)*stepy30(theta,ystri,oint)*ystri[-1] + (Vz/Iyy)*stepz30(theta,zstri,oint)*zstri[-1])
    return q

# base shear flows with booms

qb01 = integral(d_qb01, 0, m.pi/2) + qB01(m.pi/2, False)           # [N/m]
qb41 = integral(d_qb41, 0, h/2)                                    # [N/m]
qb12 = integral(d_qb12, 0, Lssk) + qB12(Lssk, False) + qb01 + qb41 # [N/m]
qb23 = integral(d_qb23, 0, Lssk) + qB23(Lssk, False) + qb12        # [N/m]
qb30 = integral(d_qb30, 0, m.pi/2) + qB30(m.pi/2, False) + qb23    # [N/m]
qb34 = integral(d_qb34, 0, h/2) + qb23                             # [N/m]


# redundant shear flows (defined clockwise positive like the base flows)

X = np.zeros((2,2))
x1, x2 = ((h*m.pi)/(2*Tsk) + h/Tsp), -h/Tsp
x3, x4 = -h/Tsp, ((2*Lssk)/Tsk + h/Tsp)

X[0,:] = x1, x2
X[1,:] = x3, x4

Y = np.zeros((2,1))

int_r11 = -(1/Tsk)*((h/2)*integral2(d_qb01, 0, m.pi/2) + qB01(m.pi/2, True) + (h/2)*integral2(d_qb30, 0, m.pi/2) + qB30(m.pi/2, True) + qb23*h*m.pi/4)
int_r12 = (1/Tsp)*(integral2(d_qb34, 0, h/2) + qb23*h/2 + integral2(d_qb41, 0, h/2))
y1 = int_r11 + int_r12

int_r21 =  -(1/Tsk)*(integral2(d_qb12, 0, Lssk) + qB12(Lssk, True) + (qb01 + qb41)*Lssk + integral2(d_qb23, 0, Lssk) + qB23(Lssk, True) + qb12*Lssk)
int_r22 = -(1/Tsp)*(integral2(d_qb34, 0, h/2) + qb23*h/2 + integral2(d_qb41, 0, h/2))
y2 = int_r21 + int_r22

Y[:,0] = y1, y2

Q = np.linalg.solve(X,Y)

q0I, q0II = Q[:,0]

Mi = (h/2)*((h/2)*integral2(d_qb01, 0, m.pi/2) + qB01(m.pi/2, True) + (h/2)*integral2(d_qb30, 0, m.pi/2) + qB30(m.pi/2, True) + qb23*h*m.pi/4)
Mi += d*(integral2(d_qb12, 0, Lssk) + qB12(Lssk, True) + (qb01 + qb41)*Lssk + integral2(d_qb23, 0, Lssk) + qB23(Lssk, True) + qb12*Lssk)
Mi += 2*(A1*q0I + A2*q0II)

eta = -Mi/Vy
SCz = -eta - h/2
#print(SCz)
