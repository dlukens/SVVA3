import math as m
import numpy as np
from tools import integral, integral2

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
P = 20600        #maximum load in actuator 2 [N]
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
#print(Cz)

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
#print(Izz)

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
#print(Iyy)

# Shear Center
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
X1, X2 = ((h*m.pi)/(2*Tsk) + h/Tsp), -h/Tsp
X3, X4 = -h/Tsp, ((2*Lssk)/Tsk + h/Tsp)

X[0,:] = X1, X2
X[1,:] = X3, X4

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
X1, X2, X3 = 2*A1, 2*A2, 0
X4, X5, X6 = (1/(2*A1))*(m.pi*h/(2*Tsk) + h/Tsp), (-1/(2*A1))*h/Tsp, -G
X7, X8, X9 = (-1/(2*A2))*h/Tsp, (1/(2*A2))*(h/Tsp + 2*Lssk/Tsk), -G

X[0,:] = X1, X2, X3
X[1,:] = X4, X5, X6
X[2,:] = X7, X8, X9

Q = np.linalg.solve(X, b)
dtheta_dx = Q[2,0] # rate of twist [rad/m]

J = T/(G*dtheta_dx) # torsional stiffness [m^4]
#print(J)




