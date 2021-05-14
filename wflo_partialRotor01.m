function [fobj,f,g]=wflo_partialRotor01(x0,ploton)
%
%%% WFLO problem 1
%%% Code from the paper;
%%% Kunakote, T., Sabangban, N., Kumar, S., Tejani, G. G., Panagant, N., Pholdee, N., S.Bureerat & Yildiz, A. R. (2021), 
%%% Comparative Performance of Twelve Metaheuristics for Wind Farm Layout Optimisation, Archives of Computational Methods in Engineering, 1-14.

if (nargin == 1)
   ploton=1; 
end


x0=round(x0);
if max(x0)>0
ns=find(x0==1);
x=ns;

WF.v0=12;%m/s, inflow wind speed
WF.Wdirection=10:10:360;%degree, wind direction angle w.r.t. +x (East) 
WF.RoseDiagram=0;% to be defined later
WF.h=60;%m, WT hub height
WF.h0=0.3;% surface roughness
WF.Dr=40;%m, HAWT rotor diameter
WF.CT=0.88;%thrust coefficient
WF.efficiency=0.4;% turbine efficientcy
WF.rho_air=1.225;%kg/m^3, air density
WF.k=0.5/log(WF.h/WF.h0);% wake expansion rate
WF.xmax=1000;

% Based on Katic, I., J?rgen H?jstrup, and Niels Otto Jensen. "A simple 
% model for cluster efficiency." European wind energy association conference
% and exhibition. A. Raguzzi, 1987.
% the problem is how to place 39 wind turbines one those grid to get the
% maximum wind power
[xwf,ywf]=meshgrid(linspace(0,2000,11));% divide the 2000x2000 m^2 area to be 100 cells.
[xc0,yc0]=meshgrid(linspace(100,2000-100,10));% centres of the 100 cells.
% figure(1),clf,hold on
% surf(xwf,ywf,0*xwf)
% plot3(xc0,yc0,0*xc0,'d')
% xlabel('x'),ylabel('y')
WF.xc0=reshape(xc0,100,1);
WF.yc0=reshape(yc0,100,1);
% WF.nWT=39;

% figure(1),clf,hold on
% for i=1:100
%     plot(xc0(i),yc0(i),'d')
%     text(xc0(i),yc0(i),num2str(i))
% end
% pause


% x=rand(100,1);% 100 design variables
% [xs,nsort]=sort(-x);% decoding to be the places of 39 wind turbines 
% x=[1 5 10 15 20]';     %  wind turbine position No.
WF.nWT=length(x); 
% WT's positions
% pre-processing for wind turbines' data
nWT=WF.nWT;
xcO=WF.xc0(x,1);
ycO=WF.yc0(x,1);

% plot(xcO,ycO,'sr')

for k=1:length(WF.Wdirection)
    theta=WF.Wdirection(k);
    WT.xc=xcO*cosd(theta)+ycO*sind(theta);
    WT.yc=-xcO*sind(theta)+ycO*cosd(theta);
    % check wind turbines' positions and compute modified inflow speed
    WShadow=0*ones(nWT,nWT);
    for i=1:nWT
        for j=(i+1):nWT
            [WShadow(i,j),vk(i,j),Ak(i,j)]=inWakeCheck(WT.xc(i),WT.yc(i),...
                WT.xc(j),WT.yc(j),WF);
            [WShadow(j,i),vk(j,i),Ak(j,i)]=inWakeCheck(WT.xc(j),WT.yc(j),...
                WT.xc(i),WT.yc(i),WF);
        end
    end
    % compute total power all the turbines can harvest
    for i2=1:nWT
        if sum(WShadow(:,i2))>0
            Iwakes=find(WShadow(:,i2)==1);
            vij=vk(Iwakes,i2);
            VD=sum(vij.^2);
            v(i2)=WF.v0*(1-sqrt(VD));
        else
            v(i2)=WF.v0;
        end
    end
    P_windFarm(k)=WF.efficiency*0.5*WF.rho_air*(pi*WF.Dr^2/4)*sum(v.^3);
end
% P_max=WF.efficiency*0.5*WF.rho_air*(pi*WF.Dr^2/4)*nWT*WF.v0^3;

Cost=nWT*(2/3+1/3*exp(-0.00174*nWT^2));
fobj=Cost*1000/mean(P_windFarm);
f=fobj;
g=0;
else
    fobj=1e10;
    f=1e10;
    g=1e10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
if ploton~=1
    figure,clf, hold on
    
    for i=1:nWT
        wakeWindMap1(WT,WF,i)
        text(WT.xc(i)-10,WT.yc(i)-10,['\bf' num2str(i)],'fontsize',12)
    end
    axis equal
    title(['Power ' num2str(mean(P_windFarm)/1e3) ' kW, Cost = ' num2str(Cost)])
%     xlabel(['Possible maximum power is ' num2str(P_max/1e3) ' kW'])
hold off
end
%%%%%%%%%%%sub-functions%%%%%%%%%%%%
function [Wij,vkij,Aij]=inWakeCheck(xci,yci,xcj,ycj,WF)
% to check whether wind turbine j is in the wake region of wind turbine i
dxij=xci-xcj;dyij=yci-ycj;dij=abs(dyij);
Dr=WF.Dr;% wind turbine diameter
Dw=Dr+2*WF.k*abs(dxij);
A0=pi*Dr^2/4;
Cr=1-sqrt(1-WF.CT);
if dxij<0&&(dij-Dr/2)<(Dw/2)% wind turbine i is shadowed by wind turbine j
    Wij=1;
    rk=Dw/2;rj=Dr/2;
    Aij=real(rk^2*acos((dij^2+rk^2-rj^2)/2/dij/rk)+rj^2*acos((dij^2+rj^2-rk^2)/2/dij/rj)-...
        0.5*sqrt((-dij+rk+rj)*(dij-rk+rj)*(dij+rk-rj)*(dij+rk+rj)));
    vkij=Cr*(Dr^2/Dw^2)*(Aij/A0);%velocity deficit weighted with partial area
else
    Wij=0;
    vkij=0;
    Aij=pi*Dr^2/4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [m1,m2]=crossLines(p1,p2,s1,s2)
% m=inv([s1 -s2])*(p2-p1);
% m1=m(1);m2=m(2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XY=WTwakeRegion(iwt,WT,WF)
xc=WT(iwt).xc;yc=WT(iwt).yc;
Dr=WF.Dr;k=WF.k;xmax=WF.xmax;
Dw=Dr+2*WF.k*xmax;
theta=WF.Wdirection*pi/180;% wind direction
x1=0;x2=xmax;x3=xmax;x4=0;
y1=-Dr/2;y2=-Dw/2;y3=Dw/2;y4=Dr/2;
A=[cos(theta)  -sin(theta)
   sin(theta) cos(theta)];
XY=A*[x1 x2 x3 x4;y1 y2 y3 y4];
XY(1,:)=xc+XY(1,:);XY(2,:)=yc+XY(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vx=windAfterWake(x,WF)
% x = distance from upstream turbines
% vx = velocity
vx=WF.v0*(1-(1-sqrt(1-WF.CT))/(1+x/(log(WF.h/WF.h0)*WF.Dr))^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function vx=findWakeShadow(WT,i,j,WF,xmax)
% xmax = maximum distance or distance from WT to the edge of the fam
% xc,yc = WT coordinate
% compare turbines i and j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wakeWindMap1(WT,WF,nt)
xc=WT.xc(nt);yc=WT.yc(nt);
Dr=WF.Dr;k=WF.k;xmax=WF.xmax;
Dw=Dr+2*WF.k*xmax;
% theta=WF.Wdirection*pi/180;% wind direction
x1=0;x2=xmax;
y1=Dr/2;y2=Dw/2;
nint=50;sc=linspace(0,1,nint);
x0=linspace(x1,x2,nint+1);
y0=linspace(y1,y2,nint+1);


plot(xc,yc,'o','markerfacecolor','b')
% for i=1:nint    
%     x1=x0(i);x2=x0(i+1);x3=x0(i+1);x4=x0(i);
%     y1=-y0(i);y2=-y0(i+1);y3=y0(i+1);y4=y0(i);
%     vxi=windAfterWake(x1,WF);
%     
%     A=[cos(theta) -sin(theta)
%        sin(theta) cos(theta)];
%     XY=A*[x1 x2 x3 x4;y1 y2 y3 y4];
%     XY(1,:)=xc+XY(1,:);XY(2,:)=yc+XY(2,:);
%     
%     surface([XY(1,1) XY(1,2);XY(1,4) XY(1,3)],...
%         [XY(2,1) XY(2,2);XY(2,4) XY(2,3)],...
%         [0 0;0 0],[vxi vxi;vxi vxi],'edgecolor','none',...
%         'facealpha',0.5)
% end
% colorbar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function drawWake(WT,WF,xmax,Colorstyle)
% xc=WT.xc;yc=WT.yc;
% Dr=WF.Dr;k=WF.k;
% Dw=Dr+2*WF.k*xmax;
% theta=WF.Wdirection*pi/180;% wind direction
% x1=0;x2=xmax;x3=xmax;x4=0;
% y1=-Dr/2;y2=-Dw/2;y3=Dw/2;y4=Dr/2;
% A=[cos(theta)  -sin(theta)
%    sin(theta) cos(theta)];
% XY=A*[x1 x2 x3 x4;y1 y2 y3 y4];
% XY(1,:)=xc+XY(1,:);XY(2,:)=yc+XY(2,:);
% hold on
% plot(XY(1,[1 4]),XY(2,[1 4]),'k','linewidth',2)
% plot(xc,yc,'o','markerfacecolor',Colorstyle)
% fill([XY(1,:) XY(1,1)],[XY(2,:) XY(2,1)],Colorstyle,'facealpha',0.25)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 

