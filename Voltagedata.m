N=10; % Number of voltage measurements per sweeping step
bw=2; % Bin width

Vdata=readmatrix("data/voltage_042125.xlsx");

t=Vdata(:,1);
V=Vdata(:,2);

% plot(t,V,'o')
% title('Raw data')
% xlabel('Relative time t s','FontSize',15)
% ylabel('RFA Voltage V','FontSize',15)
% grid on

% figure
% plot(t(941:1220),V(941:1220),'o')
% % axis([min(t(941:1260))-5 max(t(941:1260)) -950.5 -1005.3])
% grid on

Vave=zeros(length(t)/N,1);
tave=zeros(length(t)/N,1);
sigma=zeros(length(t)/N,1);

for i=1:length(Vave)
    Vave(i)=mean(V((i-1)*N+1:i*N));
    tave(i)=mean(t((i-1)*N+1:i*N));
    tmp=0;
    for j=1:N
        tmp=tmp+(Vave(i)-V((i-1)*N+j))^2;
    end
    sigma(i)=sqrt(tmp/(N-1));
end

measurement=1:length(t)/N;

% figure
% errorbar(measurement,Vave,sigma,'both','o')
% title('Averaged data')
% xlabel('Measurement number n','FontSize',15)
% ylabel('RFA Voltage V','FontSize',15)
% grid on

% figure
% errorbar(measurement(960/N:1220/N),Vave(960/N:1220/N),sigma(960/N:1220/N),'both','o')
% % axis([min(measurement(960/N:1220/N))-1 max(measurement(960/N:1220/N))+1 -1005.4 -1005.3])
% grid on

Cdata=readmatrix("data/Trigger_042125.txt");

bin=Cdata(:,1)/bw;
cts=Cdata(:,2);

tfix=zeros(size(cts));
Vfix=zeros(size(cts));

tfix(1)=tave(1);
Vfix(1)=Vave(1);
j=0;
for i=2:length(tave)
    if (tave(i)-tave(i-1))<1.2*bw
        tfix(i+j)=tave(i);
        Vfix(i+j)=Vave(i);
    else
        for j=0:(round(tave(i)-tave(i-1))/2)-1
            tfix(i+j)=tave(i-1)+(j+1)*bw;
            Vfix(i+j)=min(Vave);
        end
    end
end

Vfix(find(tfix==0),:)=[];
tfix(find(tfix==0),:)=[];
measurementfix=1:length(tfix);

Vfix=Vfix+max(Vfix);

figure
plot(measurementfix,Vfix,'o')
title('Averaged data')
xlabel('Measurement number n','FontSize',15)
ylabel('RFA Voltage V','FontSize',15)
grid on

figure
plot(bin,cts,'o')
grid on
