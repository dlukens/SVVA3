# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:06:53 2020

@author: sarta
"""

# imports
import math as m
import matplotlib.pyplot as plt
from loadingdist import Sz, Sy, T, Mz, My
from sectionproperties import Cz, z, y, Nst, h, Lssk, angle_Lssk, Ca, La, eta, Izz, Iyy, Tsk, Tsp, Ast, A1, A2, G, d
import numpy as np
from tools import integral, integral2
#---

#############################################
# Function to find shear flow, shear and bending stress distributions
#==========================================
def stress_dist(Vz,Vy,Mz,My,T,x) :
    #FUNCTIONS
    #======================================
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
    
    # step functions for the booms
    def stepz01(theta,zstri) :
        if theta < m.acos(2*zstri[1]/h) : return 0
        else : return 1
    
    def stepy01(theta,ystri) :
        if theta < m.asin(2*ystri[1]/h) : return 0
        else : return 1
    
    def stepz12(s,zstri) :
        if   s > -zstri[5]/m.cos(angle_Lssk) : return np.array([1,1,1,1])
        elif s > -zstri[4]/m.cos(angle_Lssk) : return np.array([1,1,1,0])
        elif s > -zstri[3]/m.cos(angle_Lssk) : return np.array([1,1,0,0])
        elif s > -zstri[2]/m.cos(angle_Lssk) : return np.array([1,0,0,0])
        else : return np.array([0,0,0,0])
        
    def stepy12(s,ystri) :
        if   s > (h/2 -ystri[5])/m.sin(angle_Lssk) : return np.array([1,1,1,1])
        elif s > (h/2 -ystri[4])/m.sin(angle_Lssk) : return np.array([1,1,1,0])
        elif s > (h/2 -ystri[3])/m.sin(angle_Lssk) : return np.array([1,1,0,0])
        elif s > (h/2 -ystri[2])/m.sin(angle_Lssk) : return np.array([1,0,0,0])
        else : return np.array([0,0,0,0])
        
    def stepz23(s,zstri) :
        if   s > (zstri[9] + Ca - h/2)/m.cos(angle_Lssk) : return np.array([1,1,1,1])
        elif s > (zstri[8] + Ca - h/2)/m.cos(angle_Lssk) : return np.array([1,1,1,0])
        elif s > (zstri[7] + Ca - h/2)/m.cos(angle_Lssk) : return np.array([1,1,0,0])
        elif s > (zstri[6] + Ca - h/2)/m.cos(angle_Lssk) : return np.array([1,0,0,0])
        else : return np.array([0,0,0,0])
        
    def stepy23(s,ystri) :
        if   s > -ystri[9]/m.sin(angle_Lssk) : return np.array([1,1,1,1])
        elif s > -ystri[8]/m.sin(angle_Lssk) : return np.array([1,1,1,0])
        elif s > -ystri[7]/m.sin(angle_Lssk) : return np.array([1,1,0,0])
        elif s > -ystri[6]/m.sin(angle_Lssk) : return np.array([1,0,0,0])
        else : return np.array([0,0,0,0])
    
    def stepz30(theta,zstri) :
        if theta < m.asin(2*zstri[-1]/h) : return 0
        else : return 1
        
    def stepy30(theta,ystri) :
        if theta < m.acos(-2*ystri[-1]/h) : return 0
        else : return 1
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
        
    # shear flows due to the booms
    def qB01(theta,Vy=Vy,Vz=Vz) :
        q = -Ast*((Vy/Izz)*stepy01(theta,ystri)*ystri[1] + (Vz/Iyy)*stepz01(theta,zstri)*zstri[1])
        return q
    
    def qB12(s,Vy=Vy,Vz=Vz) :
        yarr = np.array(ystri[2:6]).reshape((4,1))
        zarr = np.array(zstri[2:6]).reshape((4,1))
        q = -Ast*((Vy/Izz)*np.dot(stepy12(s,ystri),yarr)[0] + (Vz/Iyy)*np.dot(stepz12(s,zstri),zarr)[0])
        return q
    
    def qB23(s,Vy=Vy,Vz=Vz) :
        yarr = np.array(ystri[6:10]).reshape((4,1))
        zarr = np.array(zstri[6:10]).reshape((4,1))
        q = -Ast*((Vy/Izz)*np.dot(stepy23(s,ystri),yarr)[0] + (Vz/Iyy)*np.dot(stepz23(s,zstri),zarr)[0])
        return q
    
    def qB30(theta,Vy=Vy,Vz=Vz) :
        q = -Ast*((Vy/Izz)*stepy30(theta,ystri)*ystri[-1] + (Vz/Iyy)*stepz30(theta,zstri)*zstri[-1])
        return q
    #---
        
    # shear forces due to the booms
    def SB01() :
        S = -Ast*((Vy/Izz)*dB01(zstri)*ystri[1] + (Vz/Iyy)*dB01(zstri)*zstri[1])
        return S
    
    def SB12() :
        yarr = np.array(ystri[2:6]).reshape((4,1))
        zarr = np.array(zstri[2:6]).reshape((4,1))
        S = -Ast*((Vy/Izz)*np.dot(dB12(zstri,ystri),yarr)[0] + (Vz/Iyy)*np.dot(dB12(zstri,ystri),zarr)[0]) 
        return S
    
    def SB23() :
        yarr = np.array(ystri[6:10]).reshape((4,1))
        zarr = np.array(zstri[6:10]).reshape((4,1))
        S = -Ast*((Vy/Izz)*np.dot(dB23(zstri,ystri),yarr)[0] + (Vz/Iyy)*np.dot(dB23(zstri,ystri),zarr)[0])
        return S
    
    def SB30() :
        S = -Ast*((Vy/Izz)*dB30(zstri)*ystri[-1] + (Vz/Iyy)*dB30(zstri)*zstri[-1])
        return S
        
    ###################################################
        
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
    #Plot the aileron
    #plt.scatter(z01,y01,color='blue',label="nodes")
    #plt.scatter(z12,y12,color='blue')
    #plt.scatter(z23,y23,color='blue')
    #plt.scatter(z30,y30,color='blue')
    #plt.scatter(z34,y34,color='blue')
    #plt.scatter(z41,y41,color='blue')
    #plt.scatter(zstri,ystri,color='red',marker="^",label="stringer locations")
    #plt.legend(), plt.show()    
    #####################################################
    
    #Shear flow calculation for a specified cross-section (specified x-coordinate)
    
    # base shear flows with booms
    qb01 = integral(d_qb01, 0, m.pi/2) + qB01(m.pi/2)          # [N/m]
    qb41 = integral(d_qb41, 0, h/2)                            # [N/m]
    qb12 = integral(d_qb12, 0, Lssk) + qB12(Lssk) + qb01 +qb41 # [N/m]
    qb23 = integral(d_qb23, 0, Lssk) + qB23(Lssk) + qb12       # [N/m]
    #qb30 = integral(d_qb30, 0, m.pi/2) + qB30(m.pi/2) + qb23   # [N/m]
    #qb34 = integral(d_qb34, 0, h/2) + qb23                     # [N/m]

    # redundant shear flows (defined clockwise positive like the base flows)
    X = np.zeros((3,3))
    x1, x2, x3 = 2*A1, 2*A2, 0
    x4, x5, x6 = (1/(2*A1))*((h*m.pi)/(2*Tsk) + h/Tsp), -(1/(2*A1))*h/Tsp, -G
    x7, x8, x9 = -(1/(2*A2))*h/Tsp, (1/(2*A2))*((2*Lssk)/Tsk + h/Tsp), -G
    
    X[0,:] = x1, x2, x3
    X[1,:] = x4, x5, x6
    X[2,:] = x7, x8, x9
    
    Y = np.zeros((3,1))
    
    int_r1 = -(h/2)*((h/2)*integral2(d_qb01, 0, m.pi/2) + SB01() + (h/2)*integral2(d_qb30, 0, m.pi/2) + SB30() + qb23*h*m.pi/4)
    int_r1 -= d*(integral2(d_qb12, 0, Lssk) + SB12() + (qb01 + qb41)*Lssk + integral2(d_qb23, 0, Lssk) + SB23()+ qb12*Lssk)
    y1 = -Vy*eta - T + int_r1 
    
    int_r2 =  -(1/(2*Tsk*A1))*((h/2)*integral2(d_qb01, 0, m.pi/2) + SB01() + (h/2)*integral2(d_qb30, 0, m.pi/2) + SB30() + qb23*h*m.pi/4)
    int_r2 += (1/(2*Tsp*A1))*(integral2(d_qb34, 0, h/2) + qb23*h/2 + integral2(d_qb41, 0, h/2))
    y2 = int_r2
     
    int_r3 = -(1/(2*Tsk*A2))*(integral2(d_qb12, 0, Lssk) + SB12() + (qb01 + qb41)*Lssk + integral2(d_qb23, 0, Lssk) + SB23() + qb12*Lssk)
    int_r3 -= (1/(2*Tsp*A2))*(integral2(d_qb34, 0, h/2) + qb23*h*m.pi/4 + integral2(d_qb41, 0, h/2))
    y3 = int_r3
    
    Y[:,0] = y1, y2, y3
    
    Q = np.linalg.solve(X,Y)
    
    q0I, q0II, dtheta_dx = Q[:,0]
    
    # shear flow distribution functions
    def get_q01(theta) :
        q01 = integral(d_qb01, 0, theta) + qB01(theta) + q0I
        return q01
    
    def get_q12(s) :
        q12 = integral(d_qb12, 0, s) + qB12(s) + qb01 +qb41 + q0II
        return q12
    
    def get_q23(s) :
        q23 = integral(d_qb23, 0, s) + qB23(s) + qb12 + q0II
        return q23
    
    def get_q30(theta) :
        q30 = integral(d_qb30, 0, theta) + qB30(theta) + qb23 + q0I
        return q30
    
    def get_q34(s) :
        q34 = integral(d_qb34, 0, h/2) + qb23 - q0I + q0II
        return q34
    
    def get_q41(s) :
        q41 = integral(d_qb41, 0, h/2) - q0I + q0II
        return q41
    
    # shear stress distribution functions
    def get_tau01(theta) :
        tau01 = get_q01(theta)/Tsk
        return tau01
    
    def get_tau12(s) :
        tau12 = get_q12(s)/Tsk
        return tau12
    
    def get_tau23(s) :
        tau23 = get_q23(s)/Tsk
        return tau23
    
    def get_tau30(theta) :
        tau30 = get_q30(theta)/Tsk
        return tau30
    
    def get_tau34(s) :
        tau34 = get_q34(s)/Tsp
        return tau34
    
    def get_tau41(s) :
        tau41 = get_q41(s)/Tsp
        return tau41
    
    # functions for the bending stresses
    def sigx01(theta) :
        sigx = (1/Izz)*Mz*fy01(theta) + (1/Iyy)*My*(fz01(theta) - (Cz + h/2))
        return sigx
    
    def sigx12(s) :
        sigx = (1/Izz)*Mz*fy12(s) + (1/Iyy)*My*(fz12(s) - (Cz + h/2))
        return sigx
    
    def sigx23(s) :
        sigx = (1/Izz)*Mz*fy23(s) + (1/Iyy)*My*(fz23(s) - (Cz + h/2))
        return sigx
    
    def sigx30(theta) :
        sigx = (1/Izz)*Mz*fy30(theta) + (1/Iyy)*My*(fz30(theta) - (Cz + h/2))
        return sigx
    
    def sigx34(s) :
        sigx = (1/Izz)*Mz*fy34(s) + (1/Iyy)*My*(fz34(s) - (Cz + h/2))
        return sigx
    
    def sigx41(s) :
        sigx = (1/Izz)*Mz*fy41(s) + (1/Iyy)*My*(fz41(s) - (Cz + h/2))
        return sigx
    
    # function for the Von Mises stresses
    def vm01(theta) :
        vm = m.sqrt((sigx01(theta))**2 + 3*(get_tau01(theta))**2)
        return vm
    
    def vm12(s) :
        vm = m.sqrt((sigx12(s))**2 + 3*(get_tau12(s))**2)
        return vm
    
    def vm23(s) :
        vm = m.sqrt((sigx23(s))**2 + 3*(get_tau23(s))**2)
        return vm
    
    def vm30(theta) :
        vm = m.sqrt((sigx30(theta))**2 + 3*(get_tau30(theta))**2)
        return vm
    
    def vm34(s) :
        vm = m.sqrt((sigx34(s))**2 + 3*(get_tau34(s))**2)
        return vm
    
    def vm41(s) :
        vm = m.sqrt((sigx41(s))**2 + 3*(get_tau41(s))**2)
        return vm

    # max von mises stress in a cross-section
    def maxvm(x) :
        vmlst = []
        for th in np.arange(0,m.pi/2 + 0.01,0.01) :
            vmlst.append(vm01(th))
            vmlst.append(vm30(th))
        vmlst = [max(vmlst)]
        for s in np.arange(0,Lssk + 0.01,0.01) :
            vmlst.append(vm12(s))
            vmlst.append(vm23(s))
        vmlst = [max(vmlst)]
        for s in np.arange(0,h/2 + 0.01,0.01) :
            vmlst.append(vm34(s))
            vmlst.append(vm41(s))
        return max(vmlst)
    
    return maxvm(x)
        
X = np.arange(0.01,La - 0.01, 0.01) 
vmlst = []

for x in X :
    Vz, Vy, M_z, M_y, Tor = Sz(x), Sy(x), Mz(x), My(x), T(x)
    vmlst.append(stress_dist(Vz,Vy,M_z,M_y,Tor,x))
    
VMmax = max(vmlst)
print(VMmax)
        
    
