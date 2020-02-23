#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 23:10:17 2020

@author: frederik
"""

import math as m
import matplotlib.pyplot as plt

#### section properties ####

#definition of aileron geometry variables
Ca = 0.515      #aileron chord length [m]
La = 2.691      #aileron span [m]
x1 = 0.174       #x-coordinate hinge 1 [m]
x2 = 1.051       #x-coordinate hinge 2 [m]
x3 = 2.512       #x-coordinate hinge 3 [m]
xA = 0.3       #spacing between actuator I & II [m]
h = 0.248      #aileron profile thickness (height) [m]
Tsp = 0.0022   #spar thickness [m]
Tst = 0.0012   #stiffener thickness [m]
Tsk = 0.0011   #skin thickness [m]
Hst = 0.015    #stiffener height [m]
Wst = 0.030      #stiffener width [m]
Nst = 11      #number of stiffeners
delta_y1 = 1.034 #vertical displacement of hinge 1 [cm]
delta_y3 = 2.066 #vertical displacement of hinge 3 [cm]
theta = 25      #maximum upward deflection [deg]
P = 20.6        #maximum load in actuator 2 [kN]

#calculate z_coordinates stiffeners
Lssk = m.sqrt(0.5*0.5*h*h + (Ca - 0.5*h)*(Ca - 0.5*h)) #length of straight skin (one half)
Lsk = 2*Lssk + 0.5*m.pi*h #length of skin [m]

delta_str = Lsk/Nst                 #distance between stringers
angle_Lssk = m.asin((0.5*h)/Lssk)   #angle of straight skin with symmetry line

z_st0 = 0

rad_st1 = 2*delta_str/h

#upper half
z_st1 = z_st0 - (0.5*h - m.cos(rad_st1)*0.5*h)  #z_coordinate first stringer (on bent part)
z_st2 = z_st0 - Ca + m.cos(angle_Lssk)*(Lssk - abs(0.25*m.pi*h - 2*delta_str))
z_st3 = z_st0 - Ca + m.cos(angle_Lssk)*(Lssk - delta_str - abs(0.25*m.pi*h - 2*delta_str))
z_st4 = z_st0 - Ca + m.cos(angle_Lssk)*(Lssk - 2*delta_str - abs(0.25*m.pi*h - 2*delta_str))
z_st5 = z_st0 - Ca + m.cos(angle_Lssk)*(Lssk - 3*delta_str - abs(0.25*m.pi*h - 2*delta_str))

#calculate y_location stiffeners (upper half)
y_st0 = 0
y_st1 = m.sin(rad_st1)*0.5*h                    #y_coordinate first stringer (on bent part)
y_st2 = y_st0 + m.sin(angle_Lssk)*(Lssk - abs(0.25*m.pi*h - 2*delta_str))
y_st3 = y_st0 + m.sin(angle_Lssk)*(Lssk - delta_str - abs(0.25*m.pi*h - 2*delta_str))
y_st4 = y_st0 + m.sin(angle_Lssk)*(Lssk - 2*delta_str - abs(0.25*m.pi*h - 2*delta_str))
y_st5 = y_st0 + m.sin(angle_Lssk)*(Lssk - 3*delta_str - abs(0.25*m.pi*h - 2*delta_str))

Ast = Hst * Tst + Wst * Tst                     #stiffener area [mm^2]

Cy = 0 # y_coordinate of centroid of cross section [m]

#centroid z_coordinate (y_coordinate is on symmetry axis)

z_bar_Lssk = m.cos(angle_Lssk)*0.5*Lssk - Ca                   #centroid Lssk
z_bar_arc = h/m.pi - 0.5*h                                   #centroid semi-arc
z_bar_spar = -0.5*h                                               #centroid spar
A_tot = Lsk*Tsk + 11*Ast + Tsp*h                                         #total area [m^2]

Cz = (2*Ast*(z_st1+z_st2+z_st3+z_st4+z_st5) + Tsp*h*z_bar_spar  + 2*z_bar_Lssk*Lssk*Tsk + z_bar_arc*0.5*m.pi*h*Tsk)/A_tot #aileron centroid [m]

print(Cz)

z = [z_st0, z_st1, z_st2, z_st3, z_st4,z_st5]
y = [y_st0, y_st1, y_st2, y_st3, y_st4, y_st5]

# Moments of inertia
Izz_ss = (1/12)*Tsk*Lssk*Lssk*Lssk*(m.sin(angle_Lssk))**2 #MOI one straight skin
Ad_zss = Lssk*Tsk*(0.5*Lssk*m.sin(angle_Lssk))**2 #steiner term straight skin

Izz_sp = (1/12)*Tsp*h*h*h #MOI spar, no steiner term

Izz_arc = (m.pi/2)*(h/2)**3*Tsk #MOI arc, no steiner term

Izz = 0
for i in y :
    steiner = 2*Ast*i*i
    Izz += steiner

Izz += 2*(Izz_ss + Ad_zss) + Izz_sp + Izz_arc # [m^4]
print(Izz)

Iyy_ss = (1/12)*Tsk*Lssk*Lssk*Lssk*(m.sin(0.5*m.pi - angle_Lssk))**2
Ad_yss = Lssk*Tsk*(0.5*Lssk*m.cos(angle_Lssk) - Cz - Ca)**2

Iyy_sp = Tsp*h*(Cz - 0.5*h)**2 # only steiner term

Iyy_arc = Izz_arc
Ad_yarc = 0.5*m.pi*h*Tsk*(h/m.pi - 0.5*h - Cz)**2

Iyy = 0
for i in z[1:] :
    steiner = 2*Ast*(i - Cz)**2
    Iyy += steiner
    
Iyy += Ast*Cz*Cz + 2*(Iyy_ss + Ad_yss) + Iyy_sp + Iyy_arc + Ad_yarc # [m^4]
print(Iyy)

# Shear center

Sy = 1 # [N] unit shear load

# counterclockwise positive
q12 = 0.0
q23 =  Sy/Izz * (Ast*(y_st2+y_st3+y_st4+y_st5) + h*Tsk*Lssk/4) 
q31 = -Sy/Izz * (Ast*(y_st2+y_st3+y_st4+y_st5) + h*Tsk*Lssk/4)
q21 = 0
qs0_1 = Sy*h*h/(6*Izz) * (1/(m.pi/(2*Tsk) + 1/Tsp))
qs0_2 = Sy/(3*Izz) * ((h*h*h/4 + h*Lssk*Lssk)/(h/Tsp + 2*Lssk/Tsk))

# final shear flows
q12 += qs0_1
q23 += qs0_2
q31 += qs0_2
q21 += qs0_1 + qs0_2

# shear center
d = h*Lssk/(4*Lssk*Lssk + h*h) * m.sqrt(h*h + 4*Lssk*Lssk) # moment arm for q23 and q31

Mi = q12*m.pi*h*h/4 + q23*Lssk*d + q31*Lssk*d
SCz = -Mi/Sy # shear center [m]
print(SCz)
