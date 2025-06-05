import math

#resGC=1.5*10**-2 #5 #2.5*10**6 resistivity for dark GaAs
#resGC=1.5*10**-6
resGC=18*8
print("rho=",resGC,"Ohm*m")
#height=1.7*10**-6  #for 12 micron tall beam, maybe use 1.5 micron instead
#height=1.0*10**-6
height=2.0*10**-6 # change this from 1 to 2 microns
lambda_l=532.0*10**-9
lplate=0.00004 # plate length

 
h=6.087*10**-34
e=1.602*10**-19
#En=380.0*e
En=2.5*10**3*e
kB=1.38*10**-23
Temp=300.0
Pi=3.14
m=9.1*10.0**-31
v=math.sqrt(2*En/m)
#flightT=1.34*10**-12  #for approx 1 cm long strip at E=1.6 keV
flightT=lplate/v
dx=200*10**-9 # <----------- transverse coherence length

dphi=(4*10**-6)/0.24

dpx=m*v*dphi

L_coh=h/(2.0*Pi*dpx)

print("dx=",dx)
print("L_coh=",L_coh)

#dx=400.0*10**-9
C1=(dx/height)**2
#C1=1.0
#w0=1.0*10**-6
w0=50*10**-9
length=0.002  #gap from second slit to gaas edge + half gaas surface
ldB=h/(m*v)
w=w0*math.sqrt(1+(length*ldB/(Pi*w0*w0))**2)
#C2=(dx/w)**2   #assuming that dx is smaller than w
C2=1.0
tauZ=(height**3)/((dx**2)*(C1*C2))*(4*h**2/(Pi*e**2*kB*Temp*resGC))


L_nf=(lambda_l/2)**2/ldB #talbot length=2L_nf 

print("L_nf=",L_nf,"m")
 
dec_events=flightT/tauZ
print("C1=",C1,"C2=",C2)
print("tau_Zurek=",tauZ,"s")
print("flightT=",flightT,"s") # length of 10 microns
print("v=",v,"m/s")
print("dec_amount=",dec_events)
	
P=(e**2*resGC*v**2)/(16*Pi*height**3)
#tauR=e/P
E_loss=P*flightT/e
print("E_loss=",E_loss)
