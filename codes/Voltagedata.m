clear all
clc

nmeas=6; % Number of measurements per position (9 point method)
d=50; % Wire diameter in microns
npoints=400; % Maker sure to measure 400 points in a full sweep

ctsmatrix=zeros(npoints,nmeas);
Vmatrix=zeros(npoints,nmeas);
gap=zeros(nmeas,1);
Defl=zeros(nmeas,2);

syms x % Defining position in the channel's horizontal micrometer
f=4.7388*x-41.453; % This is the linear fit obtained experimentally for the gap size as function of position (x in mm)

syms y % Defining channel's tilt angle theta
g=y/675; % This is the angle in radians for L=6.75 meters of laser travel for the reflection (y in cm)

angle=zeros(nmeas,1);

for idx=1:nmeas
    % Vdata=readmatrix(['data/Post/NoDrift/9pointtest/',int2str(idx),'.xlsx']);
    % Cdata=readmatrix(['data/Post/NoDrift/9pointtest/',int2str(idx),'.txt']);
    Vdata=readmatrix(['data/Post/mesh/measHerman/zoomin/',int2str(idx),'.xlsx']);
    Cdata=readmatrix(['data/Post/mesh/measHerman/zoomin/',int2str(idx),'.txt']);

    N=Vdata(15,4); % Number of voltage measurements per sweeping step
    bw=Vdata(11,4); % Bin width
    vrange=Vdata(12,4);
    vstep=Vdata(13,4);
    steps=vrange/vstep;
    pos=Vdata(4,4);

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

    measurement=1:length(tave);

    bin=Cdata(:,1)/bw;
    cts=Cdata(:,2);

    m=round(length(t)/N/steps);

    extras=(length(cts)-steps*m)/m;

    cfinal=zeros(steps,1);
    Vfinal=zeros(steps,1);
    sigmafinal=zeros(steps,1);

    for i=1:steps
        Vfinal(i)=Vave(i);
        sigmafinal(i)=sigma(i);
    end

    cfinal=Cdata(:,2);

    % figure
    % yyaxis left;
    % plot(1:length(Vfinal),Vfinal,'o')
    % title('Averaged data')
    % xlabel('Measurement number n','FontSize',15)
    % ylabel('RFA Voltage V','FontSize',15)
    % title(['Gap position ',num2str(pos),' mm'])
    % grid on
    % 
    % yyaxis right;
    % plot(1:length(Vfinal),cfinal,'o')
    % ylabel('Counts','FontSize',15)
    % grid on

    Vfinal=Vfinal';
    Vfinal=abs(fliplr(Vfinal));
    Vmatrix(:,idx)=Vfinal;

    cfinal=cfinal';
    cfinal=fliplr(cfinal);
    ctsmatrix(:,idx)=cfinal;

    gap(idx)=subs(f,x,pos);
    Defl(idx,1)=Vdata(8,4); % V_E
    Defl(idx,2)=Vdata(9,4); % V_P
    angle(idx)=Vdata(16,4)-Vdata(17,4);
    angle(idx)=subs(g,y,angle(idx))/2;
end
%% Fitting section
% datatemp=ctsmatrix;
data=ctsmatrix;

%gap=[9.24,10,9.5,9.3,12.5,13,18,20.5];

cutoff=zeros(size(gap));
DeltaV_cutoff=zeros(size(gap));
order=[5,6,4,1,7,8,9,3,2]; % Make sure that the order of measurements is the usual one
epsilon=0.00005;
for i=1:length(gap)
    V=Vmatrix(:,i);
    vec1=[];
    y=[];
    V1=[];
    x=[];
    for j=1:length(V)
        if V(j)<(max(V)-min(V))/8+min(V)
            V1=horzcat(V1,V(j));
            vec1=horzcat(vec1,data(j,i));
        end
    end
    P1=polyfit(V1,vec1,1);
    cond=0;
    for j=length(V):-1:3
        % if data(j,i)<(P1(1)*V(j)+P1(2)) && data(j,i)>0.1*max(data(:,i))
        %     if data(j,i)<(P1(1)*V(j)+P1(2))
        %     yp2=(data(j-2,i)-data(j-1,i))/(2*(max(ctsmatrix(:,i))-min(ctsmatrix(:,i))));
        %     yp1=(data(j-1,i)-data(j,i))/(2*(max(ctsmatrix(:,i))-min(ctsmatrix(:,i))));
        %     if abs(yp2-yp1)/(abs(yp2)+abs(yp1))<1.5
        %         V2=horzcat(V(j),V2);
        %         vec2=horzcat(data(j,i),vec2);
        %         cond=1;
        %     elseif cond==1
        %         if length(V2)<6
        %             cond=0;
        %             V2=[];
        %             vec2=[];
        %         else
        %             break;
        %         end
        %     end
        % end 
        if data(j,i)>20*mean(data(length(data(:,i))-20:end,i)) && data(j,i)<0.8*max(data(:,i))
            if isempty(y)
                y=horzcat(data(j,i),y);
                x=horzcat(V(j),x);
            else
                y=horzcat(data(j,i),y);
                y1=data(j-1,i);
                % x=V(j,i); % Create Voltage matrix instead of same voltage for all of the scans
                x=horzcat(V(j),x);
                x1=V(j-1);
                P=polyfit(x,y,1);
                P2=polyfit([x1 x(end)],[y1 y(end)],1);
                % plot(V(:),data(:,i),'o')
                % hold on
                % plot(x,P(1)*x+P(2),'green','linewidth',1.5)
                % plot(x,P2(1)*x+P2(2),'red','linewidth',1.5)
                % xlim([972 974.5])
                % hold off
                if abs(1-P2(1)/P(1))>1+epsilon || abs(1-P2(1)/P(1))<1-epsilon
                    y(end)=[];
                    x(end)=[];
                end
                y=horzcat(y1,y);
                x=horzcat(x1,x);
                if y1>0.75*max(data(:,i))
                    break
                end
            end
        end
    end
    V2=x;
    vec2=y;
    P2=P;

    figure
    p=plot(V,data(:,i),'-o','color','black','LineWidth',1.5);
    p.MarkerFaceColor=[0 0 0];
    p.MarkerSize=3;
    hold on
    plot(V,P1(1)*V+P1(2),'color','blue','LineWidth',2)
    plot(V,P2(1)*V+P2(2),'color','red','LineWidth',2)
    % txt=['\leftarrow Cutoff: ' num2str(cutoff(i)) ' V'];
    % t=text(cutoff(i),(-P2(1)*P1(2)+P1(1)*P2(2))/(P1(1)-P2(1))/2,txt);
    % t.FontSize = 14;
    grid on
    xlabel('$V_{RFA}$ (V)','Interpreter','latex','FontSize',20)
    ylabel('$\# e^-$ (Counts/s)','Interpreter','latex','FontSize',20)
    title(['gap size of approximately ',num2str(gap(i)),' microns'])
    %title(['Gap position ',num2str(gap(i)),' mm for DeflX=',num2str(Defl(1,i)),'V and DeflY=',int2str(Defl(2,i)),' V'])
    axis([min(V) max(V) 0 1.2*max(data(:,i))])

    % if i<10
    %     subplot(3,3,order(i))
    %     % subplot(1,2,i)
    %     p=plot(V,data(:,i),'-o','color','black','LineWidth',1.5);
    %     p.MarkerFaceColor=[0 0 0];
    %     p.MarkerSize=3;
    %     hold on
    %     plot(V,P1(1)*V+P1(2),'color','blue','LineWidth',2)
    %     plot(V,P2(1)*V+P2(2),'color','red','LineWidth',2)
    %     txt=['\leftarrow Cutoff: ' num2str(cutoff(i)) ' V'];
    %     t=text(cutoff(i),(-P2(1)*P1(2)+P1(1)*P2(2))/(P1(1)-P2(1))/2,txt);
    %     t.FontSize = 14;
    %     grid on
    %     xlabel('$V_{RFA}$ (V)','Interpreter','latex','FontSize',20)
    %     ylabel('$\# e^-$ (Counts/s)','Interpreter','latex','FontSize',20)
    %     title(['V_E=',num2str(Defl(i,1)),'V and V_P=',num2str(Defl(i,2))])
    %     %title(['Gap position ',num2str(gap(i)),' mm for DeflX=',num2str(Defl(1,i)),'V and DeflY=',int2str(Defl(2,i)),' V'])
    %     axis([min(V) max(V) 0 1.2*max(data(:,i))])
    % else
    %     figure
    %     p=plot(V,data(:,i),'-o','color','black','LineWidth',1.5);
    %     p.MarkerFaceColor=[0 0 0];
    %     p.MarkerSize=3;
    %     hold on
    %     plot(V,P1(1)*V+P1(2),'color','blue','LineWidth',2)
    %     plot(V,P2(1)*V+P2(2),'color','red','LineWidth',2)
    %     txt=['\leftarrow Cutoff: ' num2str(cutoff(i)) ' V'];
    %     t=text(cutoff(i),(-P2(1)*P1(2)+P1(1)*P2(2))/(P1(1)-P2(1))/2,txt);
    %     t.FontSize = 14;
    %     grid on
    %     xlabel('$V_{RFA}$ (V)','Interpreter','latex','FontSize',20)
    %     ylabel('$\# e^-$ (Counts/s)','Interpreter','latex','FontSize',20)
    %     title(['V_E=',num2str(Defl(i,1)),'V and V_P=',num2str(Defl(i,2))])
    %     %title(['Gap position ',num2str(gap(i)),' mm for DeflX=',num2str(Defl(1,i)),'V and DeflY=',int2str(Defl(2,i)),' V'])
    %     axis([min(V) max(V) 0 1.2*max(data(:,i))])
    % end

    %%%%%% Wayne's idea to understand bump
    % bumpC=zeros(length(V),1);

    % for j=1:length(V)
    %     if data(j,i)-(P1(1)*V(j)+P1(2))>=0
    %         bumpC(j)=data(j,i)-(P1(1)*V(j)+P1(2));
    %     end
    % end

    % subplot(2,1,2)
    % plot(V,bumpC,'.')
    % grid on
    %%%%%%

    cutoff(i)=-P2(2)/P2(1);
    DeltaV_cutoff(i)=(P2(1)*P1(2)-P1(1)*P2(2))/((P1(1)-P2(1))*P2(1));
end

%% Mean and standard deviation (for multiple measurements only)
[G, gapmean] = findgroups(gap);      % group indices + unique values
cutoffmean  = splitapply(@mean, cutoff, G);
cutoffstd  = splitapply(@std, cutoff, G);

%% Plotting Voltage cutoff as a function of gap size, their mean and standard deviation
figure
plot(gap,cutoff,'or','linewidth',1.5)
grid on
xlabel('Gap size $\mu$m','Interpreter','latex','FontSize',20)
ylabel('Cutoff voltage V','FontSize',20)
title('Gap size vs. Cutoff voltage','FontSize',20)
xticks=0:1:ceil(max(gap)+2);
set(gca,'XTick',xticks);
axis([min(xticks) max(xticks) 1000*Vdata(1,4) 1000*Vdata(1,4)+5])
hold on
errorbar(gapmean,cutoffmean,cutoffstd,'xb','linewidth',1.5)

%% Plotting Voltage cutoff as a function of angle
figure
plot(angle*1e3,cutoff,'o','linewidth',1.5)
grid on
xlabel('Tilt angle mrad','FontSize',20)
ylabel('Cutoff voltage V','FontSize',20)
title(['Tilt angle vs. Cutoff voltage for gap size ~',num2str(gap(i)),' microns'],'FontSize',20)
% xticks=0:1:ceil(max(gap)+2);
% set(gca,'XTick',xticks);
axis([floor(min(angle*1e3)-1) ceil(max(angle*1e3)+1) 1000*Vdata(1,4) 1000*Vdata(1,4)+5])

%%
VE=unique(Defl(:,1));
VP=unique(Defl(:,2));

C = [cutoff(4), cutoff(9), cutoff(8);           % Counts matrix (3x3) (replace with your actual values)
     cutoff(3), cutoff(1), cutoff(2);
     cutoff(5), cutoff(6), cutoff(7)];
C(2,1)=704;
C(3,3)=703.6;

%C(3,2)=974.5;

% Create meshgrid
[VE_grid, VP_grid] = meshgrid(VE, VP);

[VE_fine, VP_fine] = meshgrid(linspace(min(VE), max(VE), 100), ...
                              linspace(min(VP), max(VP), 100));
C_fine = interp2(VE_grid, VP_grid, C', VE_fine, VP_fine, 'spline');

[maxC_fine, idx_fine] = max(C_fine(:));
[row_fine, col_fine] = ind2sub(size(C_fine), idx_fine);

% Coordinates of the maximum
V1_max = VE_fine(1, col_fine);
V2_max = VP_fine(row_fine, 1);

figure;
surf(VE_fine, VP_fine, C_fine, 'EdgeColor', 'none');
xlabel('V_E V','FontSize',20);
ylabel('V_P V','FontSize',20);
zlabel('V_{cutoff} V','FontSize',20);
formattedText = sprintf('Maximum = %.3f at VE = %.3f V, VP = %.3f V\n', maxC_fine, V1_max, V2_max);
text(2, 0.8, formattedText, 'FontSize', 20, 'Color', 'blue');
title(formattedText)

colormap jet;
colorbar;
shading interp;

%% Proper fit
[xData, yData, zData] = prepareSurfaceData( VE_fine, VP_fine, C_fine );

% Set up fittype and options.
ft = fittype('a*x.^2 + b*y.^2 + c*x + d*y + e', ...
                         'dependent', 'z', 'independent', {'x', 'y'}, ...
                         'coefficients', {'a', 'b','c','d','e'});

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 1', 'C_fine vs. VE_fine, VP_fine', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'VE_fine', 'Interpreter', 'none' );
ylabel( 'VP_fine', 'Interpreter', 'none' );
zlabel( 'C_fine', 'Interpreter', 'none' );
grid on
view( -14.7, -5.8 );

%%
figure
p=plot(gap,cutoff,'o');
p.MarkerFaceColor=[0 0 0];
p.MarkerSize=5;
hold on
% for i=1:length(gap)
%     txt=['\DeltaV: ' num2str(DeltaV_cutoff(i)) ' V'];
%     t=text(gap(i)-1,cutoff(i)-0.05,txt);
%     t.FontSize = 14;
% end
xlabel('Gap height $z\ \mu m$','Interpreter','latex','FontSize',30)
ylabel('$V_{cutoff}\ V$','Interpreter','latex','FontSize',30)
% title(['Gap height ',num2str(gap(1)),'\mu m'],'FontSize',30)
% xTickLocations=[0 5 10 15 20 25 50];
yTickLocations=702:0.5:706;
% yticks=([1000 1010 1020]);
% set(gca,'XTick', xTickLocations);
set(gca,'YTick', yTickLocations);
axis([0 60 min(yTickLocations) max(yTickLocations)])
hold on
errorbar(gap(1),mean(cutoff),std(cutoff),'rx')
grid on

%% Plotting as a function of time

%time=linspace(0,118,nmeas);
t0=5*60+42;
time=[t0,5*60+55,6*60+09,6*60+33,6*60+46,7*60+00,7*60+14,7*60+28,7*60+41]-t0;
figure
plot(time,cutoff,'o','LineWidth',2)
xlabel('time min','FontSize',20)
ylabel('RFA voltage V','FontSize',20)
grid on
hold on
% P=polyfit([time(1) time(2)],[cutoff(1) cutoff(2)],1);
% plot(linspace(0,max(time),100),P(1)*linspace(0,max(time),100)+P(2))

%% Long mesurement of e-gun power supply

metering=[11*60+08,67.720,-3.337,-2.758,33.909;...
            11*60+15,68.092,-3.124,-2.886,33.715;...
            11*60+20,69.893,-1.754,-1.974,34.767;...
            11*60+25,69.901,-1.755,-1.975,34.771;...
            12*60+07,69.932,-1.758,-1.979,34.775;...
            12*60+12,69.910,-1.824,-2.045,34.759;...
            12*60+35,69.919,-1.825,-2.049,34.761;...
            13*60+18,69.931,-1.827,-2.048,34.767;...
            13*60+40,69.933,-1.829,-2.048,34.769;...
            14*60+00,69.937,-1.826,-2.052,34.768;...
            14*60+09,69.941,-1.819,-2.041,34.766;...
            14*60+22,69.942,-1.819,-2.040,34.766;...
            15*60+16,69.948,-1.820,-2.041,34.769;...
            15*60+31,69.949,-1.820,-2.041,34.767;...
            15*60+44,69.950,-1.821,-2.042,34.768;...
            15*60+59,69.951,-1.820,-2.041,34.769;...
            16*60+16,69.952,-1.820,-2.041,34.770;...
            16*60+30,69.952,-1.821,-2.042,34.768;...
            16*60+55,69.951,-1.821,-2.042,34.769;...
            18*60+19,69.952,-1.821,-2.042,34.770];

figure
yyaxis left;
plot(metering(:,1),metering(:,2),'o','linewidth',2)
ylabel('Energy voltage V','fontsize',20)
yyaxis right;
plot(metering(:,1),metering(:,5),'o','linewidth',2)
ylabel('Focus voltage V','fontsize',20)
xlabel('Time t min','fontsize',20)
title('Source voltages as function of time','fontsize',20)
grid on

figure
yyaxis left;
plot(metering(:,1),metering(:,3),'o','linewidth',2)
ylabel('Deflection X V','fontsize',20)
yyaxis right;
plot(metering(:,1),metering(:,4),'o','linewidth',2)
ylabel('Deflection Y V','fontsize',20)
xlabel('Time t min','fontsize',20)
title('Deflection voltages as function of time','fontsize',20)
grid on