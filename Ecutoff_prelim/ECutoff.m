clear all
clc

data1=readmatrix("NewAuChannel_cutoff_1.txt");
data1b=readmatrix("NewAuChannel_cutoff_1b.txt");
data2=readmatrix("NewAuChannel_cutoff_2.txt");
data2b=readmatrix("NewAuChannel_cutoff_2b.txt");
data3=readmatrix("NewAuChannel_cutoff_3.txt");
data3b=readmatrix("NewAuChannel_cutoff_3b.txt");
data4=readmatrix("NewAuChannel_cutoff_4.txt");
data5=readmatrix("NewAuChannel_cutoff_5.txt");
data6=readmatrix("NewAuChannel_cutoff_6.txt");
data7=readmatrix("NewAuChannel_cutoff_7.txt");
data8=readmatrix("NewAuChannel_cutoff_8.txt");

datatemp=horzcat(data1,data1b,data2,data2b,data3,data3b,data4,data5,data6,data7,data8);

data=zeros(size(datatemp,1),size(datatemp,2)/2);

for j=1:size(data,2)
    for i=1:size(data,1)
        data(i,j)=datatemp(i,2*j);
    end
end

gap=[20.77674024,20.77674024,17.59337861,17.59337861,12.28777589,12.28777589,5.921052632,4.859932088,3.798811545,2.737691002,1.591680815];
crate=[733.6,755.8,504.2,703.4,547,487.8,198.8,113.2,164,106,26.6];

V=linspace(2025*0.495,2075*0.495,length(data));
cutoff=zeros(size(gap));
DeltaV_cutoff=zeros(size(gap));
for i=1:length(gap)
    vec1=[];
    vec2=[];
    V2=[];
    if i==1
        V1=[];
    end
    for j=1:length(V)
        if V(j)>1003 && V(j)<1012
            if i==1
                V1=horzcat(V1,V(j));
            end
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
    % figure
    % plot(V,data(:,i),'o')
    % hold on
    % plot(V,P1(1)*V+P1(2))
    % plot(V,P2(1)*V+P2(2))
    % txt=['\leftarrow Cutoff: ' num2str(cutoff(i)) ' V'];
    % t=text(cutoff(i),(-P2(1)*P1(2)+P1(1)*P2(2))/(P1(1)-P2(1))/2,txt);
    % t.FontSize = 14;
    % grid on
    % xlabel('RFA Voltage $V$','Interpreter','latex','FontSize',14)
    % ylabel('Counts/s','Interpreter','latex','FontSize',14)
    % title(['Gap height ',num2str(gap(i)),' um for Average count rate ',num2str(crate(i)),' C/s'])
    % axis([min(V) max(V) 0 1.2*max(data(:,i))])
    cutoff(i)=-P2(2)/P2(1);
    DeltaV_cutoff(i)=(P2(1)*P1(2)-P1(1)*P2(2))/((P1(1)-P2(1))*P2(1));
end

figure
plot(gap,cutoff,'o')
for i=1:length(gap)
    txt=['\DeltaV: ' num2str(DeltaV_cutoff(i)) ' V'];
    t=text(gap(i)-1,cutoff(i)-0.05,txt);
    t.FontSize = 14;
end
xlabel('Gap size $\mu m$','Interpreter','latex','FontSize',14)
ylabel('Cutoff voltage $V$','Interpreter','latex','FontSize',14)
grid on