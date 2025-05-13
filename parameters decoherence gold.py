# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

vel=2*1E7
leng=0.0015
tf=leng/vel
h=6.63*1E-34
z=0.85*1E-6
Pi=3.14
e=1.6E-19
kB=1.38E-23
T=300
pho=1E-7
dx=4E-7
C1=1
C2=1
eps0=8.85E-12
c=3E8
tau=4*h**2*z**3/(Pi*e**2*kB*T*pho*dx**2*C1*C2)
#tau=4*h**2*z**3*(4*Pi*eps0)/(Pi*e**2*kB*T*pho*dx**2*C1*C2)
Rz=tf/tau
print(Rz)


P=e**2*pho*vel**2/(16*Pi*z**3)/e/(4*Pi*eps0)
dE=P*tf
print(dE)


me = 9.11E-31
kB = 1.38E-23
T = 300
hbar = 1.055E-34
qe = 1.60E-19
E = 1.67E3*qe
v0 = np.sqrt(2*E/me)
Pi=3.14 
tflight = 1.5E-3/v0
rho = 1.9E-7
z0 = 0.85E-6
C1C2 = 1
dx = 400E-9
decoh_rate = Pi*qe**2*kB*T*rho*tflight*dx**2/(2*(2*Pi*hbar)**2*z0**3)
print(decoh_rate)

print((2/3)*1/(16*(30E-9)**(1.5))/(tflight/(z0**3)))
