N=20;

data=readmatrix("voltage_041625.xlsx");

t=data(:,1);
V=data(:,2);

plot(t,V,'o')
title('Raw data')
xlabel('Relative time t s','FontSize',15)
ylabel('RFA Voltage V','FontSize',15)
grid on

figure
plot(t(941:1220),V(941:1220),'o')
axis([min(t(941:1260))-5 max(t(941:1260)) -1005.4 -1005.3])
grid on

Vave=zeros(1,length(t)/N);
sigma=zeros(1,length(t)/N);

for i=1:length(Vave)
    Vave(i)=mean(V((i-1)*N+1:i*N));
    tmp=0;
    for j=1:N
        tmp=tmp+(Vave(i)-V((i-1)*N+j))^2;
    end
    sigma(i)=sqrt(tmp/(N-1));
end

measurement=1:length(t)/N;

figure
errorbar(measurement,Vave,sigma,'both','o')
title('Averaged data')
xlabel('Measurement number n','FontSize',15)
ylabel('RFA Voltage V','FontSize',15)
grid on

figure
errorbar(measurement(960/N:1220/N),Vave(960/N:1220/N),sigma(960/N:1220/N),'both','o')
axis([min(measurement(960/N:1220/N))-1 max(measurement(960/N:1220/N))+1 -1005.4 -1005.3])
grid on