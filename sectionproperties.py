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
Tsp = 2.2      #spar thickness [mm]
Tst = 1.2      #stiffener thickness [mm]
Tsk = 1.1      #skin thickness [mm]
Hst = 150      #stiffener height [mm]
Wst = 300      #stiffener width [mm]
Nst = 11      #number of stiffeners
delta_y1 = 1.034 #vertical displacement of hinge 1 [cm]
delta_y3 = 2.066 #vertical displacement of hinge 3 [cm]
theta = 25      #maximum upward deflection [deg]
P = 20.6        #maximum load in actuator 2 [kN]

#calculate z_coordinates stiffeners
Lssk = m.sqrt((0.5*h)**2 + (Ca-0.5*h)**2) #length of straight skin (one half)
Lsk = 2*Lssk + 0.5*m.pi*h #length of skin [m]

delta_str = Lsk/Nst                 #distance between stringers
angle_Lssk = m.asin((0.5*h)/Lssk)   #angle of straight skin with symmetry line

z_st0 = 0

rad_st1 = delta_str/(h*0.25*m.pi)

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


#centroid z_coordinate (y_coordinate is on symmetry axis)

z_bar_Lssk = (m.cos(angle_Lssk)*0.5*Lssk*1000 - Ca)*Tsk*Lssk*1000   #centroid Lssk
z_bar_arc = (h*1000)/m.pi - 0.5*h*1000
z_bar_spar = -Tsp*(h*1000)*0.5*(h*1000)
A_tot = Lsk*Tsk*1000 + 11*Ast * Tsp*h*1000

z_bar = (2*Ast*1000*(z_st1+z_st2+z_st3+z_st4+z_st5) + z_bar_spar  + 2*z_bar_Lssk + z_bar_arc)/A_tot

print(z_bar)

z = [z_st0, z_st1, z_st2, z_st3, z_st4,z_st5]
y = [y_st0, y_st1, y_st2, y_st3, y_st4, y_st5]

# plt.plot(z, y)
# plt.show()