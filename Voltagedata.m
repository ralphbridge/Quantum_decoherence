ctsmatrix=zeros(16/0.05,8);
Vmatrix=zeros(16/0.05,8);
gap=[];
Defl=zeros(2,8);

for idx=1:8
    Vdata=readmatrix(['data/voltage_050525_',int2str(idx),'.xlsx']);
    Cdata=readmatrix(['data/counts_050525_',int2str(idx),'.txt']);

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

    for i=1:m
        cfinal=cfinal+cts((i-1)*steps+(i-1)*extras+1:i*steps+(i-1)*extras);
    end

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

    gap=horzcat(gap,pos);
    Defl(1,idx)=Vdata(8,4);
    Defl(2,idx)=Vdata(9,4);
end

%% Fitting section
% datatemp=ctsmatrix;
data=ctsmatrix;

%gap=[9.24,10,9.5,9.3,12.5,13,18,20.5];

cutoff=zeros(size(gap));
DeltaV_cutoff=zeros(size(gap));
for i=1:length(gap)
    V=Vmatrix(:,i);
    vec1=[];
    vec2=[];
    V1=[];
    V2=[];
    for j=1:length(V)
        if V(j)<(max(V)-min(V))/4+min(V)
            V1=horzcat(V1,V(j));
            vec1=horzcat(vec1,data(j,i));
        end
    end
    P1=polyfit(V1,vec1,1);
    if i==length(gap)
        1;
    end
    for j=length(V):-1:1
        if data(j,i)>0.1*(P1(1)*V(j)+P1(2)) && data(j,i)>0.05*max(data(:,i))
            if data(j,i)<0.8*(P1(1)*V(j)+P1(2))
                V2=horzcat(V(j),V2);
                vec2=horzcat(data(j,i),vec2);
            else
                break
            end
        end
    end
    P2=polyfit(V2,vec2,1);

    figure
    subplot(2,1,1)
    p=plot(V,data(:,i),'-o','color','black','LineWidth',1.5);
    p.MarkerFaceColor=[0 0 0];
    p.MarkerSize=3;
    hold on
    plot(V,P1(1)*V+P1(2),'color','blue','LineWidth',2)
    plot(V,P2(1)*V+P2(2),'color','red','LineWidth',2)
    txt=['\leftarrow Cutoff: ' num2str(cutoff(i)) ' V'];
    t=text(cutoff(i),(-P2(1)*P1(2)+P1(1)*P2(2))/(P1(1)-P2(1))/2,txt);
    t.FontSize = 14;
    grid on
    xlabel('$V_{RFA}$ (V)','Interpreter','latex','FontSize',16)
    ylabel('$\# e^-$ (Counts/s)','Interpreter','latex','FontSize',16)
    title(['Gap height ',num2str(round((gap(i)-8.57)*50/11.2)),' \mu m'])
    %title(['Gap position ',num2str(gap(i)),' mm for DeflX=',num2str(Defl(1,i)),'V and DeflY=',int2str(Defl(2,i)),' V'])
    axis([min(V) max(V) 0 1.2*max(data(:,i))])

    %%%%%% Wayne's idea to understand bump
    bumpC=zeros(length(V),1);

    for j=1:length(V)
        if data(j,i)-(P1(1)*V(j)+P1(2))>=0
            bumpC(j)=data(j,i)-(P1(1)*V(j)+P1(2));
        end
    end

    subplot(2,1,2)
    plot(V,bumpC,'.')
    grid on
    %%%%%%

    cutoff(i)=-P2(2)/P2(1);
    DeltaV_cutoff(i)=(P2(1)*P1(2)-P1(1)*P2(2))/((P1(1)-P2(1))*P2(1));
end

% figure
% p=plot((gap-8.57)*50/11.2,cutoff,'o');
% p.MarkerFaceColor=[0 0 0];
% p.MarkerSize=5;
% hold on
% for i=1:length(gap)
%     txt=['\DeltaV: ' num2str(DeltaV_cutoff(i)) ' V'];
%     t=text(gap(i)-1,cutoff(i)-0.05,txt);
%     t.FontSize = 14;
% end
% xlabel('Gap height $z\ \mu m$','Interpreter','latex','FontSize',30)
% ylabel('$V_{cutoff}\ V$','Interpreter','latex','FontSize',30)
% xTickLocations=[0 5 10 15 20 25];
% yTickLocations=[1000 1010 1020];
% yticks=([1000 1010 1020]);
% plot(gap,1017.7-(0.3*(1.5)^3)./gap)
% set(gca,'XTick', xTickLocations);
% set(gca,'YTick', yTickLocations);
% axis([0 60 952 955])
% grid on
% yline(953.75)
% yline(953.25)
