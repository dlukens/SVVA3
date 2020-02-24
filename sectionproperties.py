import math as m
import numpy as np

#### section properties ####

#definition of aileron geometry variables
Ca = 0.515      #aileron chord length [m]
La = 2.691      #aileron span [m]
x1 = 0.174       #x-coordinate hinge 1 [m]
x2 = 1.051       #x-coordinate hinge 2 [m]
x3 = 2.512       #x-coordinate hinge 3 [m]
xa = 0.3       #spacing between actuator I & II [m]
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
E = 73.1 *10**9   #Young's modulus [Pa]

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

# Centroid z_coordinate (y_coordinate is on symmetry axis)

z_bar_Lssk = m.cos(angle_Lssk)*0.5*Lssk - Ca                   #centroid Lssk
z_bar_arc = h/m.pi - 0.5*h                                   #centroid semi-arc
z_bar_spar = -0.5*h                                               #centroid spar
A_tot = Lsk*Tsk + 11*Ast + Tsp*h                                         #total area [m^2]

Cz = (2*Ast*(z_st1+z_st2+z_st3+z_st4+z_st5) + Tsp*h*z_bar_spar  + 2*z_bar_Lssk*Lssk*Tsk + z_bar_arc*0.5*m.pi*h*Tsk)/A_tot #aileron centroid [m]


z = [z_st0, z_st1, z_st2, z_st3, z_st4,z_st5]
y = [y_st0, y_st1, y_st2, y_st3, y_st4, y_st5]

# Moments of Inertia
Izz_ss = ((1/12)*Tsk*Lssk**3)*(m.sin(angle_Lssk))**2 #MOI one straight skin #checked
Ad_zss = Lssk*Tsk*(0.5*Lssk*m.sin(angle_Lssk))**2 #steiner term straight skin #checked

Izz_sp = (1/12)*Tsp*h**3 #MOI spar, no steiner term #checked

Izz_arc = Tsk*(m.pi/16)*(h)**3 #MOI arc, no steiner term 

Izz = 0
for i in y :
    steiner = 2*Ast*i**2
    Izz += steiner

Izz += 2*(Izz_ss + Ad_zss) + Izz_sp + Izz_arc # [m^4]

Iyy_ss = (1/12)*Tsk*Lssk**3*(m.cos(angle_Lssk))**2
Ad_yss = Lssk*Tsk*(Ca+Cz-0.5*Lssk*m.cos(angle_Lssk))**2

Iyy_sp = Tsp*h*(-Cz - 0.5*h)**2 # only steiner term


Iyy_arc=Izz_arc
Ad_yarc = 0.5*m.pi*h*Tsk*(Cz**2 + 0.5*h*(0.5*h - 2 * -Cz + 4 * -Cz / m.pi - 4 * 0.5*h / m.pi))

Iyy = 0
for i in z[1:] :
    steiner = 2*Ast*(Cz-i)**2
    Iyy += steiner
    
Iyy += Ast*Cz**2 + 2*(Iyy_ss + Ad_yss) + Iyy_sp + Iyy_arc +Ad_yarc # [m^4]

# Shear Center

Sy = 1 # [N] unit shear load

# clockwise positive
q12 = 0
q23 = -Sy/Izz * (Ast*(y_st2+y_st3+y_st4+y_st5) + h*Tsk*Lssk/4) 
q31 =  Sy/Izz * (Ast*(y_st2+y_st3+y_st4+y_st5) + h*Tsk*Lssk/4)
q142 = 0
qs0_1 = Sy*h**3/(3*Izz*(m.pi*h/(2*Tsk) + h/Tsp))
qs0_2 = Sy*h*(Lssk**2 - h**2)/(12*Izz*(2*Lssk/Tsk + h/Tsp))

# final shear flows
q12 += qs0_1
q23 += qs0_2
q31 += qs0_2
q142 += qs0_1 + qs0_2

# taking moment about midpoint of the spar and clockwise positive
d = h*m.cos(angle_Lssk)/2 # moment arm for q23 and q31
Mi = h/2*(m.pi*h*qs0_1/2 - Sy*h**3*Tsk/(4*Izz)) + d*(2*qs0_2*Lssk - Sy*h*Lssk**2*Tsk/(12*Izz))

Ksi = -Mi/Sy # distance midpoint spar - shear center (positive)
SCz = -h/2 - Ksi # shear center [m]

# Torsional Stiffness

T = 1 # [Nm] unit torsional load
G = 28e9 # Aluminium 2024-T3 shear modulus [Pa]

A1 = 0.5*m.pi*(h/2)**2 # area of cell 1
A2 = (Ca - h/2)*h/2 # area of cell 2

# setting up the 3x3 system
b = np.matrix([[T],
               [0],
               [0]])

X = np.zeros((3,3))
x1, x2, x3 = 2*A1, 2*A2, 0
x4, x5, x6 = 1/(2*A1)*(m.pi*h/(2*Tsk) + h/Tsp), -1/(2*A1)*h/Tsp, -G
x7, x8, x9 = -1*h/(2*Tsp*A2), 1*(h/Tsp + 2*Lssk/Tsk)/(2*A2), -G


X[0,:] = x1, x2, x3
X[1,:] = x4, x5, x6
X[2,:] = x7, x8, x9

Q = np.linalg.solve(X, b)
dtheta_dz = Q[2,0] # rate of twist [rad/m]

J = T/(G*dtheta_dz) # torsional stiffness [m^4]