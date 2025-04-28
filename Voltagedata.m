N=10; % Number of voltage measurements per sweeping step
bw=2; % Bin width
vrange=20-0;
vstep=0.05;
steps=vrange/vstep;

Vdata=readmatrix("data/voltage_042425.xlsx");

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

Cdata=readmatrix("data/Trigger_042425.txt");

bin=Cdata(:,1)/bw;
cts=Cdata(:,2);

m=round(length(t)/N/steps);

for i=2:length(tave)-m
    if abs(tave(i)-tave(i-1))>2*bw
        extras=round((tave(i)-tave(i-1)))/bw;
        tave(i-1)=[];
        Vave(i-1)=[];
    end
end

measurement=1:length(tave);

tfix=zeros(size(cts));
Vfix=zeros(size(cts));

% figure
% subplot(1,2,1)
% plot(measurement,Vave,'o')
% grid on

j=0;
for i=1:m
    Vfix((i-1)*steps+(i-1)*extras+1:i*steps+(i-1)*extras)=Vave((i-1)*steps+1:i*steps);
    tfix((i-1)*steps+(i-1)*extras+1:i*steps+(i-1)*extras)=tave((i-1)*steps+1:i*steps);
    if i<m
        Vfix(i*steps+(i-1)*extras+1:i*steps+i*extras)=Vave(i*steps+1);
        tfix(i*steps+(i-1)*extras+1:i*steps+i*extras)=tave(i*steps+1);
    end
end

cts(find(tfix==0),:)=[];
bin(find(tfix==0),:)=[];
Vfix(find(tfix==0),:)=[];
tfix(find(tfix==0),:)=[];
measurementfix=1:length(tfix);

bin=bin+1;

Vfix=abs(Vfix);
Vfix=Vfix-min(Vfix);

figure
yyaxis left;
plot(measurementfix,Vfix,'o')
title('Averaged data')
xlabel('Measurement number n','FontSize',15)
ylabel('RFA Voltage V','FontSize',15)
grid on

yyaxis right;
plot(bin,cts,'o')
ylabel('Counts','FontSize',15)
grid on

xline(400)
xline(425)
xline(825)
xline(850)
xline(1250)
xline(1275)
xline(1675)

figure
yyaxis left;
plot(measurementfix(1:400),Vfix(1:400),'o')
title('Averaged data')
xlabel('Measurement number n','FontSize',15)
ylabel('RFA Voltage V','FontSize',15)
grid on

yyaxis right;
plot(bin(1:400),cts(1:400),'o')
ylabel('Counts','FontSize',15)
grid on

hold on
yyaxis right;
plot(bin(1:400),cts(426:825),'bo')
ylabel('Counts','FontSize',15)
grid on

yyaxis right;
plot(bin(1:400),cts(851:1250),'mo')
ylabel('Counts','FontSize',15)
grid on

yyaxis right;
plot(bin(1:400),cts(1276:1675),'ko')
ylabel('Counts','FontSize',15)
grid on
