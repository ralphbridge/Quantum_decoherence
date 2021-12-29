tic

%ws=[5.5e-6,5.6e-6,5.7e-6,5.8e-6,5.9e-6,6e-6,7e-6,10e-6,20e-6,21e-6,22e-6,23e-6,24e-6,25e-6];
%ws=[6.1e-6,6.2e-6,6.3e-6,6.4e-6,6.5e-6,6.6e-6,6.7e-6,6.8e-6,6.9e-6,7.2e-6,7.4e-6,7.6e-6,7.8e-6,8e-6,8.2e-6,8.4e-6,8.6e-6,8.8e-6,9e-6];
%ws=[5e-6,20e-6];
%ws=[6.7e-6,22e-6,7.9e-6,8e-6,25e-6];

% for index=1:length(ws)
% clearvars -except ws index
% w1=ws(index);

% This code computes the KD diffraction pattern
% following the double slit (full) code
% with a point source.
% Everything in SI units.
clear all
clc
format long

global m;
m=9.10938356e-31;
q=1.602e-19;
global hbar;
hbar=1.0545718e-34;
c=299792458;
eps0=8.85e-12;

E=2.5e3; % 380 electron energy in eV
E_loss=545; % (68 for 22 um and 545 for 100 um) Energy loss in eV

E=E-E_loss;

D=125e-6; % Laser beam waist
lamL=532e-9; % laser wavelength
k=2*pi/lamL;
wL=k*c;
I=18e14; % 18e14 up to 10^16 W/m^2

w1=400e-6; % (6.7-22 microns) FWHM of incoherent gaussian slit width
sigma1=w1/(2*sqrt(2*log(2)));
w2=1e-6; % FWHM of collimating gaussian slit width
sigma2=w2/(2*sqrt(2*log(2)));
w3=5e-6; % FWHM of detection gaussian slit width
sigma3=w3/(2*sqrt(2*log(2)));
l0=24e-2; % distance from incoherent source slit to collimating slit
l1=1e-2; % 5.5 cm distance from second slit to grating
l2=24e-2; % distance from grating to screen

N0=401; % Number of points for the incoherent source slit x-coordinate
if N0==1
    x0=0;
    dx0=w1/N0;
else
    x0=linspace(-w1,w1,N0);
    dx0=w1/(N0-1);
    x0=x0';
end

x00=linspace(-w1,w1,100*N0);

% figure
% subplot(1,3,1)
% stem(x0*1e6,conj(exp(-(x0).^2/(2*sigma1^2))).*exp(-(x0).^2/(2*sigma1^2)))
% hold on
% plot(x00*1e6,conj(exp(-(x00).^2/(2*sigma1^2))).*exp(-(x00).^2/(2*sigma1^2)))
% grid on
% xlabel('Incoherent source position $x_0\ \mu m$','fontsize',15,'interpreter','latex')
% ylabel('$P_{inc}=\left|\Psi_{inc}(x_0)\right|^2$','fontsize',15,'interpreter','latex')

N01=500; % 500 Number of points for the collimating second slit x-coordinate
x01=linspace(-w2,w2,N01);
dx01=w2/(N01-1);
x01=x01';

xmax1=1e-6;
N1=1001; % 10001 Number of points for the grating (laser) x-coordinate
x1=linspace(-xmax1,xmax1,N1);
dx1=xmax1/(N1-1);
x1=x1';

xmax2=600e-6; % 100 microns to see only first and second diffraction orders
N2=6001; % 6001 Number of points for the detection screen x-coordinate
x2=linspace(-xmax2,xmax2,N2);
dx2=xmax2/(N2-1);
x2=x2';

V0=q^2*I/(2*m*eps0*c*wL^2); % Ponderomotive potential amplitude

Vpond=V0*(cos(k*x1)).^2; % ponderomotive potential

P_detect=zeros(length(x2),length(x0)); % initializing the detection probability matrix

rho=zeros(length(x1),length(x1)); % initializing density matrix
rho_d=zeros(length(x2),length(x2)); % initializing density matrix

% This section computes the propagation of the source wavefunction from the source to the one on the screen

for iter=1:length(x0)
    v=sqrt(2*(E*q)/m); % electron velocity ~ 1.1e7 m/s
    t=D/v; % total interaction time
    
    Psi_inc=zeros(size(x0)); % initializign wavefunction at incoherent source slit
    Psi_slit=zeros(size(x01)); % initializign wavefunction at collimating slit
    Psi=zeros(size(x1)); % initializing wavefunction before laser interaction
    Psi_detect=zeros(size(x2)); % initializing wavefunction at detection screen
    
    rhok=zeros(length(x1),length(x1)); % initializing density matrix for k-th incoherent source
        
    % This section computes the source wavefunction (wf) using a incoherent slit
    
    Psi_inc(iter)=exp(-(x0(iter))^2/(2*sigma1^2));
    
%     if iter==1||iter==length(x0)
%         figure
%         subplot(2,1,1)
%         stem(x0*1e6,conj(Psi_inc).*Psi_inc)
%         hold on
%         plot(x00*1e6,conj(exp(-(x00).^2/(2*sigma1^2))).*exp(-(x00).^2/(2*sigma1^2)))
%         xlabel('Incoherent source position $x_0\ \mu m$','fontsize',15,'interpreter','latex')
%         ylabel('$P_{inc}=\left|\Psi_{inc}(x_0)\right|^2$','fontsize',15,'interpreter','latex')
%     end
    
    % This section computes the wf propagation from the first (incoherent source) gaussian slit of FWHM w1 to the second (collimating) gaussian slit of FWHM w2

    Psi_slit=exp(-(x01).^2/(2*sigma2^2)).*propagate(Psi_inc,x0,dx0,x01,l0,v);
        
    % This section computes the wf propagation from the second (collimating) slit to the laser

    Psi=propagate(Psi_slit,x01,dx01,x1,l1,v);
    
    for i=1:length(x1) % Right before the laser beam
        for j=1:length(x1)
            rhok(i,j)=conj(Psi(i))*Psi(j);
            rho(i,j)=rho(i,j)+conj(Psi_inc(iter))*Psi_inc(iter)*rhok(i,j);
        end
    end

    % This section computes the wf after laser interaction
    
%     v=0.9*v;
%     t=D/v;
    Psi_grate=Psi.*exp(-1i*Vpond*t/hbar);

    % This section computes the wf propagation from the laser to the screen

    Psi_detect=propagate(Psi_grate,x1,dx1,x2,l2,v);
    
    P_detect(:,iter)=conj(Psi_detect).*Psi_detect; % This line computes the detection probability for each delta
    
%     if iter==1||iter==length(x0)
%         subplot(2,1,2)
%         plot(x2*1e6,P_detect(:,iter))
%         xlabel('Screen position $x_2\ \mu m$','fontsize',15,'interpreter','latex')
%         ylabel('$P_{final}=\left|\Psi_{final}(x_2)\right|^2$','fontsize',15,'interpreter','latex')
%     end
end

% This section sums over the incoherent source

P_detincoh=sum(P_detect,2)*dx0;

Dslit=zeros(size(x2)); % detection slit
for i=1:length(x2)
    Dslit(i)=exp(-(x2(i))^2/(2*sigma3^2));
end

P_final=conv(Dslit,P_detincoh,'same'); % this line convolutes the detection slit with the probability pattern
norm=trapz(x2,P_final);
P_final=P_final/norm;

subplot(1,2,1)
hold on
plot(x2*1e6,P_final)
xlabel('Screen position $x_2\ \mu m$','fontsize',15,'interpreter','latex')
ylabel('$P_{final}=\left|\Psi_{final}(x_2)\right|^2$','fontsize',15,'interpreter','latex')
% axis([min(x2*1e6) max(x2*1e6) min(P_final) 1.2*max(P_final)])
grid on

adiagonal=zeros(size(x1));

for i=1:length(x1)
    adiagonal(i)=conj(rho(length(x1)+1-i,i))*rho(length(x1)+1-i,i);
end

% adiagonal=adiagonal/max(adiagonal);
% adiagonal=adiagonal;

subplot(1,2,2)
hold on
plot(x1*1e9,adiagonal)
axis([min(x1)*1e9 max(x1)*1e9 min(adiagonal) 1.2*max(adiagonal)])
xlabel('Grating (laser) x-coordinate $x_1\ nm$','interpreter','latex','fontsize',16)
ylabel('Density matrix anti-diagonal profile (arb. units)','interpreter','latex','fontsize',16)
% text(200,1.03,'$w_1=$','interpreter','latex','fontsize',20)
% text(500,1.03,sprintf('%g',w1*1e6),'fontsize',15)
% text(700,1.03,'$\mu m$','interpreter','latex','fontsize',20)
grid on

% figure
% imagesc(x1*1e6,x1*1e6,abs(rho))
% axis square
% xlabel('Grating (laser) x-coordinate $x_1\ \mu m$','fontsize',15,'interpreter','latex')
% ylabel('Grating (laser) x-coordinate $x_1\ \mu m$','fontsize',15,'interpreter','latex')

for iter=1:length(x0)
    rhok_d=zeros(length(x2),length(x2)); % initializing density matrix for k-th incoherent source
    for i=1:length(x2)
        for j=1:length(x2)
            rhok_d(i,j)=P_final(j);
            rho_d(i,j)=rho_d(i,j)+conj(Psi_inc(iter))*Psi_inc(iter)*rhok_d(i,j);
        end
    end
end

diagonal=zeros(size(x2));

for i=1:length(x2)
    diagonal(i)=conj(rho_d(length(x2)+1-i,i))*rho_d(length(x2)+1-i,i);
end

% end

timeElapsed=string(toc) + ' seconds' % execution time in seconds

function Psi_out=propagate(Psi_in,x_in,dx,x_out,l,v)
    global hbar;
    global m;
    lamdB=2*pi*hbar/(m*v); % deBroglie wavelength
	Psi_out=zeros(size(x_out));
	for i=1:length(x_out)
		for j=1:length(x_in)
			Kij=exp(1j*2*pi*sqrt((x_out(i)-x_in(j))^2+l^2)/lamdB);
			Psi_out(i)=Psi_out(i)+Kij*Psi_in(j)*dx;
        end
    end
end