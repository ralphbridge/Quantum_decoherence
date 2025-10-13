clear all
clc

e=1.602e-19;
m=9.1e-31;
E=700*e; % Energy in Joules
v=sqrt(2*E/m);
theta=1e-3;
vT=v*theta;
rhoAu=2.44e-8; % Resistivity of Gold at room temperature
% z=1e-6; % Height of the electron with respect to the surface
z=1e-9;
L=2e-3; % Length of the surface
% flightT=L/vT; % Flight time
flightT=3e-15;

P=(e^2*rhoAu*v^2)/(16*pi*z^3);
E_loss=P*flightT/e