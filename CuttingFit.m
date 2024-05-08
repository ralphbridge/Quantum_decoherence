clear all
clc

data=readmatrix('data/2024-05-03_fullbeam_2ndslit.xlsx');

xdata=data(:,1)*1e-3;
ydata=data(:,2);

if data(1,2)<data(length(ydata),2) % 1 for increasing data (erf), -1 for decreasing data (erfc)
    slope=1;
else
    slope=-1;
end

if slope==1
    ydata=ydata-(max(ydata)+min(ydata))/2; % Adjusting the vertical data to matlabs erf(x)
    lb=[max(ydata)/4,xdata(abs(ydata)==min(abs(ydata)))/2,0.01]; % lower bounds for the fit (they might change depending on the units)
    ub=[max(ydata),xdata(abs(ydata)==min(abs(ydata)))*2,100]; % upper bounds for the fit
    fun=@(x,xdata)x(1)*erf((xdata-x(2))/x(3));
    x0=[max(ydata)/2,max(xdata)+min(xdata),0.0001]; % Initial guess
elseif slope==-1
    ydata=ydata-min(ydata);
    lb=[max(ydata)/3,min(xdata(abs(ydata)==min(abs(ydata)))/2),0.01];
    ub=[max(ydata),min(xdata(abs(ydata)==min(abs(ydata)))*2),100];
    fun=@(x,xdata)x(1)*erfc((xdata-x(2))/x(3));
    x0=[max(ydata)/2,(max(xdata)+min(xdata))/2,0.0001];
end

x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);

FWHM=2*sqrt(2*log(2))*x(3)/sqrt(2);

xteo=linspace(min(xdata),max(xdata),1000);

figure
subplot(2,1,1)
plot(xdata,ydata,'ko',xteo,fun(x,xteo),'b-')%,linspace(min(xdata),max(xdata)),max(ydata)*erf((linspace(min(xdata),max(xdata))-xdata(abs(ydata)==min(abs(ydata))))/7),'r-')
legend('Experimental data','Fitted $A$ erf$\left(\frac{x-x_0}{\sigma\sqrt 2}\right)$','interpreter','latex','FontSize',15)
xlabel('$x\ mm$','interpreter','latex','FontSize',15)
ylabel('Counts')
grid on
title('Data and Fitted Curve')
subplot(2,1,2)
plot(linspace(-5*x(3),5*x(3),2000),max(ydata)*exp(-(linspace(-5*x(3),5*x(3),2000)).^2/(x(3)^2)))
xlabel('$x\ mm$','interpreter','latex','FontSize',15)
ylabel('Arb. units')
text(x(3),max(ydata)/2,sprintf('FWHM=%f mm',FWHM),"FontSize",15)
grid on
