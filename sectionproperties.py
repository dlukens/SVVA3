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
d1 = 0.01034 #vertical displacement of hinge 1 [m]
d2 = 0 #[m]
d3 = 0.02066 #vertical displacement of hinge 3 [m]
theta = 25      #maximum upward deflection [deg]
P = -20600        #maximum load in actuator 2 [N]
E = 73.1 *10**9 #[Pa]

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

Ast = Hst * Tst + Wst * Tst                     #stiffener area [m^2]

Cy = 0 # y_coordinate of centroid of cross section [m]

# Centroid z_coordinate (y_coordinate is on symmetry axis)

z_bar_Lssk = m.cos(angle_Lssk)*0.5*Lssk - Ca                   #centroid Lssk
z_bar_arc = h/m.pi - 0.5*h                                   #centroid semi-arc
z_bar_spar = -0.5*h                                               #centroid spar
A_tot = Lsk*Tsk + 11*Ast + Tsp*h                                         #total area [m^2]

Cz = (2*Ast*(z_st1+z_st2+z_st3+z_st4+z_st5) + Tsp*h*z_bar_spar  + 2*z_bar_Lssk*Lssk*Tsk + z_bar_arc*0.5*m.pi*h*Tsk)/A_tot #aileron centroid [m]

print(Cz)

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
print(Izz)

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
print(Iyy)

# Shear Center

Sy = 1 # [N] unit shear load

# base shear flows
q1 = -(Sy*Tsk*h**2)/(4*Izz) - (Sy*Ast*y_st1)/Izz  # [N/m]
q2 = -(Sy*Tsp*h**2)/(8*Izz) # [N/m]
q3 = -(Sy*h*Lssk*Tsk)/(4*Izz) - ((Sy*Ast)/Izz)*(y_st2+y_st3+y_st4+y_st5) + q1 + q2 # [N/m]
q4 = (Sy*h*Lssk*Tsk)/(4*Izz) + ((Sy*Ast)/Izz)*(y_st2+y_st3+y_st4+y_st5) + q3 # [N/m] 
q5 = (Sy*Tsp*h**2)/(8*Izz) + q4 # [N/m]
q6 = (Sy*Tsk*h**2)/(4*Izz) + (Sy*Ast*y_st1)/Izz + q4

# distances along which the boom shears act
e1 = Lssk - 2*delta_str + m.pi*h/4 # length of straight skin between boom 2 and the TE
e2 = 5*delta_str - m.pi*h/4 # length of straight skin between boom 5 and the start of the arc
E1 = np.array([e1, e1 - delta_str, e1 - 2*delta_str, e1 - 3*delta_str])
E2 = np.array([e2 - 3*delta_str, e2 - 2*delta_str, e2 - delta_str, e2])

# base shear forces
S1 = (Sy*Tsk*(h**3)*(1 - m.pi/2))/(8*Izz) - (Sy*Ast*y_st1*(h*m.pi/4 - delta_str))/Izz # [N]
S2 = -(Sy*Tsp*h**3)/(16*Izz) # [N]
S3 = -(Sy*Tsk*h*Lssk**2)/(6*Izz) - ((Sy*Ast)/(Izz))*(y_st2*E1[0]+y_st3*E1[1]+y_st4*E1[2]+y_st5*E1[3]) + (q1 + q2)*Lssk # [N]
S4 = (Sy*Tsk*h*Lssk**2)/(12*Izz) + ((Sy*Ast)/(Izz))*(y_st2*E2[0]+y_st3*E2[1]+y_st4*E2[2]+y_st5*E2[3]) + q3*Lssk # [N]
S5 = -(Sy*Tsp*h**3)/(12*Izz) + q4*h/2 # [N]
S6 = -(Sy*Tsk*h**2)/(4*Izz) + (Sy*Ast*y_st1*delta_str)/(Izz) + q4*h*m.pi/4

# redundant shear flows: clockwise positive
A = np.zeros((2,2))
b = np.zeros((2,1))
a1, a2, b1 = ((m.pi*h)/(2*Tsk) + h/Tsp), -h/Tsp, -(S1 + S6)/Tsk + (S2 + S5)/Tsp
a3, a4, b2 = -h/Tsp, (h/Tsp + (2*Lssk)/Tsk), -(S3 + S4)/Tsk - (S2 + S5)/Tsp
A[0,:] = a1, a2
A[1,:] = a3, a4
b[0,:] = b1
b[1,:] = b2

# results of redundant shear flows
Y = np.linalg.solve(A,b)
q0_1, q0_2 = Y[0,0], Y[1,0]

# taking moment about midpoint of the spar and clockwise positive
A1 = 0.5*m.pi*(h/2)**2 # area of cell 1
A2 = (Ca - h/2)*h/2 # area of cell 2
d = h*m.cos(angle_Lssk)/2 # moment arm for q23 and q31
Mi = (h/2)*(S1 + S6) + d*(S3 + S4) + 2*(A1*q0_1 + A2*q0_2)

Ksi = -Mi/Sy # distance midpoint spar - shear center (positive)
SCz = -h/2 - Ksi # shear center [m]
print(SCz)

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
x4, x5, x6 = (1/(2*A1))*(m.pi*h/(2*Tsk) + h/Tsp), (-1/(2*A1))*h/Tsp, -G
x7, x8, x9 = (-1/(2*A2))*h/Tsp, (1/(2*A2))*(h/Tsp + 2*Lssk/Tsk), -G

X[0,:] = x1, x2, x3
X[1,:] = x4, x5, x6
X[2,:] = x7, x8, x9

Q = np.linalg.solve(X, b)
dtheta_dx = Q[2,0] # rate of twist [rad/m]

J = T/(G*dtheta_dx) # torsional stiffness [m^4]
print(J)




