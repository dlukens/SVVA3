# -*- coding: utf-8 -*-
import math as m
from sectionproperties import h, SCz, x1, x2, x3, xa, P, G, J, E, Izz, Iyy, La, d1, d3
import numpy as np
from forces import A_coeff, Cp_coeff
import interp
from interp import grid_x
import tools as t
import matplotlib.pyplot as plt

###variables
alpha = m.radians(25) # angle of attack
eta = np.abs (SCz + h/2)       # distance midpoint spar - shear center
xI=x2-0.5*xa
xII=x2+0.5*xa
###Funcions

def step(x, x1, exp): #Mcaulay step function
    y = x - x1
    if y < 0: return 0 
    else: return y**exp
    

    
def A_SC_int(x):
    idx = np.searchsorted(grid_x, x)-1
    I = np.zeros(idx)
    
    for i in range(idx-1):
        a = A_coeff[0, i]
        b = A_coeff[1, i]
        c = A_coeff[2, i]
        d = A_coeff[3, i]
        e = Cp_coeff[0, i]
        f = Cp_coeff[1, i]
        g = Cp_coeff[2, i]
        h = Cp_coeff[3, i]
            
        def func(z):
            return z*(a*SCz/4*(z-grid_x[i])**3 + b*SCz/3*(z-grid_x[i])**2 + c*SCz/2*(z-grid_x[i]) + d*SCz   \
                        -  a/7*(z-grid_x[i])**6                                  \
                        - (a*f + b*e)/6*(z-grid_x[i])**5                         \
                        - (a*g + b*f + c*e)/5*(z-grid_x[i])**4                   \
                        - (a*h + b*g + c*f + d*e)/4*(z-grid_x[i])**3             \
                        - (b*h + c*g + d*f)/3*(z-grid_x[i])**2                   \
                        - (c*h + d*g)/2*(z-grid_x[i])                            \
                        -  d*h)
        
        if i != idx-2:
            I[i+1] = I[i] + func(grid_x[i + 1]) - func(grid_x[i])
        else:
            I[i+1] = I[i] + func(x) - func(grid_x[i])
    return I[-1]                     

                        
def A_SC_doubleint(x):
    idx = np.searchsorted(grid_x, x)-1
    I = np.zeros(idx)
    
    for i in range(idx-1):
        a = A_coeff[0, i]
        b = A_coeff[1, i]
        c = A_coeff[2, i]
        d = A_coeff[3, i]
        e = Cp_coeff[0, i]
        f = Cp_coeff[1, i]
        g = Cp_coeff[2, i]
        h = Cp_coeff[3, i]
            
        def func(z):
            return z**2*(a*SCz/(4*5)*(z-grid_x[i])**3 + b*SCz/(3*4)*(z-grid_x[i])**2 + c*SCz/(2*3)*(z-grid_x[i]) + d*SCz/2 \
                        -  a/(8*7)*(z-grid_x[i])**6 \
                        - (a*f + b*e)/(6*7)*(z-grid_x[i])**5 \
                        - (a*g + b*f + c*e)/(5*6)*(z-grid_x[i])**4 \
                        - (a*h + b*g + c*f + d*e)/(4*5)*(z-grid_x[i])**3 \
                        - (b*h + c*g + d*f)/(3*4)*(z-grid_x[i])**2 \
                        - (c*h + d*g)/(2*3)*(z-grid_x[i]) \
                        -  d*h/2)
        
        if i != idx-2:
            I[i+1] = I[i] + func(grid_x[i + 1]) - func(grid_x[i])
        else:
            I[i+1] = I[i] + func(x) - func(grid_x[i])
    return I[-1]
    
    
def A_int(x):
    idx = np.searchsorted(grid_x, x)-1
    I = np.zeros(idx)
    
    for i in range(idx-1):
        a = A_coeff[0, i]
        b = A_coeff[1, i]
        c = A_coeff[2, i]
        d = A_coeff[3, i]
            
        def f(z):
            return 1/4*a*(z - grid_x[i])**4 + 1/3*b*(z - grid_x[i])**3 + 1/2*c*(z - grid_x[i])**2 + d*z
        
        if i != idx-2:
            I[i+1] = I[i] + f(grid_x[i + 1]) - f(grid_x[i])
        else:
            I[i+1] = I[i] + f(x) - f(grid_x[i])
    return I[-1]

    
def A_doubleint(x):
    idx = np.searchsorted(grid_x, x)-1
    I = np.zeros(idx)
    
    for i in range(idx-1):
        a = A_coeff[0, i]
        b = A_coeff[1, i]
        c = A_coeff[2, i]
        d = A_coeff[3, i]
            
        def f(z):
            return z**2*(a/(4*5)*(z-grid_x[i])**3     \
               + b/(3*4)*(z-grid_x[i])**2     \
               + c/(2*3)*(z-grid_x[i])        \
               + d/2)
        
        if i != idx-2:
            I[i+1] = I[i] + f(grid_x[i + 1]) - f(grid_x[i])
        else:
            I[i+1] = I[i] + f(x) - f(grid_x[i])
    return I[-1]

def A_quadint(x):
    idx = np.searchsorted(grid_x, x)-1
    I = np.zeros(idx)
    
    for i in range(idx-1):
        a = A_coeff[0, i]
        b = A_coeff[1, i]
        c = A_coeff[2, i]
        d = A_coeff[3, i]
            
        def f(z):
            return z**4*(a/(4*5*6*7)*(z-grid_x[i])**3     \
               + b/(3*4*5*6)*(z-grid_x[i])**2     \
               + c/(2*3*4*5)*(z-grid_x[i])        \
               + d/(2*3*4))
        
        if i != idx-2:
            I[i+1] = I[i] + f(grid_x[i + 1]) - f(grid_x[i])
        else:
            I[i+1] = I[i] + f(x) - f(grid_x[i])
    return I[-1]


    
#Reaction forces
    
Rxn=np.zeros((12,12))

Rxn[0]= [0,-eta*step(La,x1,0),0,-eta*step(La, x2, 0),0,-eta*step(La, x3, 0),-( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(La, xI, 0),0,0,0,0,0] #T(la) =0

Rxn[1]= [step(La, x1, 1), 0, step(La, x2, 1), 0, step(La,x3, 1), 0, -m.cos(alpha)*(step(La, xI, 1)), 0, 0, 0, 0, 0]  #My(la)=0

Rxn[2]= [0, step(La, x1, 1), 0, step(La, x2, 1), 0 , step(La, x3, 1) , -m.sin(alpha)*step(La, xI, 1), 0, 0 , 0, 0, 0] #Mz(la)

Rxn[3]= [0, step(La, x1, 0), 0, step(La, x2, 0), 0 , step(La, x3, 0) , -m.sin(alpha)*step(La, xI, 0), 0, 0 , 0, 0, 0] #Sy(la)

Rxn[4]= [step(La, x1, 0), 0, step(La, x2, 0), 0, step(La,x3, 0), 0, -m.cos(alpha)*(step(La, xI, 0)), 0, 0, 0, 0, 0]   #Sz(la)

Rxn[5]=  [-step(x1, x1, 3)/(6*E*Iyy), 0 , -step(x1, x2, 3)/(6*E*Iyy) , 0 , -step(x1, x3, 3)/(6*E*Iyy), 0 , m.cos(alpha)*step(x1, xI, 3)/(6*E*Iyy) , 0 , 0, step(x1,0,1) , 1, 0] #w(x1)

Rxn[6]=  [-step(x2, x1, 3)/(6*E*Iyy), 0 , -step(x2, x2, 3)/(6*E*Iyy) , 0 , -step(x2, x3, 3)/(6*E*Iyy), 0 , m.cos(alpha)*step(x2, xI, 3)/(6*E*Iyy) , 0 , 0, step(x2,0,1) , 1, 0] #w(x2)

Rxn[7]=  [-step(x3, x1, 3)/(6*E*Iyy), 0 , -step(x3, x2, 3)/(6*E*Iyy) , 0 , -step(x3, x3, 3)/(6*E*Iyy), 0 , m.cos(alpha)*step(x3, xI, 3)/(6*E*Iyy) , 0 , 0, step(x3,0,1) , 1, 0] #w(x3)

Rxn[8]= [-step(xI, x1, 3)/(6*E*Iyy)*m.cos(alpha), -m.sin(alpha)*(step(xI,x1,3)/(6*E*Izz) - eta**2*step(xI,x1,1)/(G*J)) , -step(xI, x2, 3)/(6*E*Iyy)*m.cos(alpha) , 0 , -step(xI, x3, 3)/(6*E*Iyy)*m.cos(alpha), 0 , m.cos(alpha)*step(xI, xI, 3)/(6*E*Iyy)*m.cos(alpha) , step(xI,0,1)*m.sin(alpha) , m.sin(alpha), step(xI,0,1)*m.cos(alpha) , m.cos(alpha), -(eta+h/2)*m.sin(alpha)]  #w(xI)=0    

Rxn[9]= [0,step(x1,x1,3)/(6*E*Izz) + eta**2*step(x1,x1,1)/(G*J), 0, step(x1,x2,3)/(6*E*Izz) + eta**2*step(x1,x2,1)/(G*J), 0 ,step(x1,x3,3)/(6*E*Izz) + eta**2*step(x1,x3,1)/(G*J), -m.sin(alpha)*step(x1,xI,3)/(6*E*Izz) + eta*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x1,xI,1)/(G*J), step(x1,0,1),1,0,0,-eta] #d(x1)

Rxn[10]= [0,step(x2,x1,3)/(6*E*Izz) + eta**2*step(x2,x1,1)/(G*J), 0, step(x2,x2,3)/(6*E*Izz) + eta**2*step(x2,x2,1)/(G*J), 0 ,step(x2,x3,3)/(6*E*Izz) + eta**2*step(x2,x3,1)/(G*J), -m.sin(alpha)*step(x2,xI,3)/(6*E*Izz) + eta*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x2,xI,1)/(G*J), step(x2,0,1),1,0,0,-eta] #d(x2)

Rxn[11]= [0,step(x3,x1,3)/(6*E*Izz) + eta**2*step(x3,x1,1)/(G*J), 0, step(x3,x2,3)/(6*E*Izz) + eta**2*step(x3,x2,1)/(G*J), 0 ,step(x3,x3,3)/(6*E*Izz) + eta**2*step(x3,x3,1)/(G*J), -m.sin(alpha)*step(x3,xI,3)/(6*E*Izz) + eta*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x3,xI,1)/(G*J), step(x3,0,1),1,0,0,-eta] #d(x3)


#Bc1= [[-A_SC_int(La) - P*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)],  #T(la) #g
#       [P*m.cos(alpha)*(La - xII)],                             #My(la) #g
#       [A_doubleint(La) + P*m.sin(alpha)*(La-xII)],          #Mz(la) #g
#        [P*m.sin(alpha) +A_int(La) ],           #g                         #Sy(la) 
#       [P*m.cos(alpha)],                                                #Sz(la) #g
#       [d1*m.sin(alpha) ],                                               #w(x1)
#       [0],                                                             #w(x2)
#       [d3*m.sin(alpha) -P*m.cos(alpha)*(x3-xII)**3/(6*E*Iyy)],   #w(x3)
#       [0-A_quadint(xI)*m.sin(alpha)/(E*Izz) + eta*A_SC_doubleint(xI)*m.sin(alpha)/(G*J) ],   #w'(xI)
#       [d1*m.cos(alpha) - (A_quadint(x1)/(6*E*Izz)) - (eta*A_SC_doubleint(x1)/(G*J))], #vertical deflection at x1
#       [-A_quadint(x2)/(E*Izz) - eta*A_SC_doubleint(x2)/(G*J)],
#       [d3*m.cos(alpha) - (A_quadint(x3)/(6*E*Izz)) - eta*A_SC_doubleint(x3)/(G*J) - P*m.sin(alpha)*(x3 - xII)**3/(6*E*Izz) - P*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*(x3-xII)*eta/(G*J)]]


Bc= np.zeros(12)

Bc=[[-A_SC_int(La) - P*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(La, xII, 0)], #T(la)=0
    [P*m.cos(alpha)*step(La, xII, 1)], #My(la)=0
    [A_doubleint(La) +P*m.sin(alpha)*step(La, xII, 1)], #Mz(la)
    [A_int(La) + P*m.sin(alpha)*step(La, xII, 0)], #Sy(la)=0
    [ P*m.cos(alpha)*step(La, xII, 0) ], #Sz(la)
    [-d1*m.sin(alpha) -P*m.cos(alpha)*step(x1, xII, 3)/(6*E*Iyy) ], #w(x1)
    [0-P*m.cos(alpha)*step(x2, xII, 3)/(6*E*Iyy)], #w(x2)
    [-d3*m.sin(alpha)-P*m.cos(alpha)*step(x3, xII, 3)/(6*E*Iyy) ],  #w(x3))
    [0-A_quadint(xI)/(E*Izz)*m.sin(alpha) + m.sin(alpha)*eta*A_SC_doubleint(xI)/(G*J)],
    [d1*m.cos(alpha) - (A_quadint(x1)/(E*Izz)) + eta*A_SC_doubleint(x1)/(G*J) - P*m.sin(alpha)*step(x1,xII,3)/(6*E*Izz) + P*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x1,xII,1)*eta/(G*J)],
    [0 - (A_quadint(x2)/(E*Izz)) + eta*A_SC_doubleint(x2)/(G*J) - P*m.sin(alpha)*step(x2,xII,3)/(6*E*Izz) + P*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x2,xII,1)*eta/(G*J)],
    [d3*m.cos(alpha) - (A_quadint(x3)/(E*Izz)) + eta*A_SC_doubleint(x3)/(G*J) - P*m.sin(alpha)*step(x3,xII,3)/(6*E*Izz) + P*( m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x3,xII,1)*eta/(G*J)]]


F= np.linalg.solve(Rxn, Bc)


R1z=float(F[0])
R1y=float(F[1])
R2z=float(F[2])
R2y=float(F[3])
R3z=float(F[4])
R3y=float(F[5])
RI=float(F[6])
C1=float(F[7])
C2=float(F[8])
C3=float(F[9])
C4=float(F[10])
C5=float(F[11])

 #Torque X-axis
def T(xx): return A_SC_int(xx) - eta*R1y*step(xx,x1,0) +RI*(m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(xx,xI,0) - eta*R2y*step(xx,x2,0) +P*(m.sin(alpha)*(eta+h/2)-m.cos(alpha)*h*0.5)*step(xx,xII,0) -eta*R3y*step(xx,x3,0)

 #Moment in y-axis
def My(x): return R1z*step(x, x1, 1) - RI*m.cos(alpha)*(step(x, xI, 1)) + R2z*step(x, x2, 1) - P*m.cos(alpha)*step(x, xII, 1) + R3z*step(x,x3, 1)

 #Moment in X-axis                                  
def Mz(x): return A_doubleint(x) - R1y*step(x, x1, 1) + RI*m.sin(alpha)*step(x, xI, 1) - R2y*step(x, x2, 1) + P*m.sin(alpha)*step(x, xII, 1) - R3y*step(x, x3, 1)

 #Shear y-axis
def Sy(x): return A_int(x) - R1y*step(x, x1, 0) + RI*m.sin(alpha)*step(x, xI, 0) - R2y*step(x, x2, 0) + P*m.sin(alpha)*step(x, xII, 0) - R3y*step(x, x3, 0)

 #Shear Z-Axis
def Sz(x): return R1z*step(x, x1, 0) - RI*m.cos(alpha)*(step(x, xI, 0)) + R2z*step(x, x2, 0) - P*m.cos(alpha)*step(x, xII, 0) + R3z*step(x,x3, 0)

 #Deflection in Y-axis
def v(x): return (1/(6*E*Izz))*(6*A_quadint(x) - R1y*step(x, x1, 3) + RI*m.sin(alpha)*step(x, xI, 3) - R2y*step(x, x2, 3) + P*m.sin(alpha)*step(x, xII, 3) - R3y*step(x, x3, 3)) +C1*step(x,0,1) + C2


 #Deflection in Z-axis
def w(x): return (1/(6*E*Iyy))*(-R1z*step(x, x1, 3) + RI*m.cos(alpha)*step(x, xI, 3) - R2z*step(x, x2, 3) + P*m.cos(alpha)*step(x, xII, 3) - R3z*step(x,x3, 3)) + C3*step(x,0,1) + C4

 #Twist around X-axis
def theta(x): return (1/(G*J))*(A_SC_doubleint(x) - eta*R1y*step(x,x1,1) +RI*(m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x,xI,1) - eta*R2y*step(x,x2,1) +P*(m.sin(alpha)*(eta+0.5*h)-m.cos(alpha)*0.5*h)*step(x,xII,1) -eta*R3y*step(x,x3,1)) + C5              


beginnode = 0.01
endnode = 2.68
nodesnumber = 1000
deltax = (endnode - beginnode)/(nodesnumber)

lst=np.arange(beginnode,endnode+deltax,deltax)
thetalist=np.zeros(len(lst))
vdeflectionlist=np.zeros(len(lst))
wdeflectionlist=np.zeros(len(lst))
Torquelist=np.zeros(len(lst))
Mylist = np.zeros(len(lst))
Mzlist = np.zeros(len(lst))
Sylist = np.zeros(len(lst))
Szlist = np.zeros(len(lst))



for i in range(len(lst)):
    thetalist[i]=theta(lst[i])
    vdeflectionlist[i]=v(lst[i]) - eta*theta(lst[i])
    wdeflectionlist[i]=w(lst[i])
    Torquelist[i]=T(lst[i])
    Mylist[i] = My(lst[i])
    Mzlist[i] = Mz(lst[i])
    Sylist[i] = Sy(lst[i])
    Szlist[i] = Sz(lst[i])
plt.subplot(331)    
plt.plot(lst,thetalist, color='red')
plt.xlabel('x[m]')
plt.ylabel('[degrees]')
plt.title('Twist')

plt.subplot(332)
plt.plot(lst,vdeflectionlist, color='green')
plt.xlabel('x[m]')
plt.ylabel('v[m]')
plt.title('Deflection in y')

#plt.subplot(333)
#plt.plot(lst,wdeflectionlist, color='blue')
#plt.xlabel('x[m]')
#plt.ylabel('w[m]')
#plt.title('Deflection in z')
#plt.show()

plt.subplot(334)    
plt.plot(lst,Torquelist, color='red')
plt.xlabel('x[m]')
plt.ylabel('T[Nm]')
plt.title('Torque')  

#plt.subplot(335)
#plt.plot(lst,Mylist, color='green')
#plt.xlabel('x[m]')
#plt.ylabel('My[m^4]')
#plt.title('Bending Moment-Y')

#plt.subplot(336)
#plt.plot(lst,Mzlist, color='blue')
#plt.xlabel('x[m]')
#plt.ylabel('Mz[m^4]')
#plt.title('Bending Moment-Z')
#plt.show()

plt.subplot(337)
plt.plot(lst,Sylist, color='red')
plt.xlabel('x[m]')
plt.ylabel('Sy[N]')
plt.title('Shear Y')
#
#plt.subplot(338)
#plt.plot(lst,Szlist, color='blue')
#plt.xlabel('x[m]')
#plt.ylabel('Sz[N]')
#plt.title('Shear Z')
#plt.show()