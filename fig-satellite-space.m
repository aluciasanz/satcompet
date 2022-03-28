% Adriana Sanz 14/0672017
% this script is to visualice parameter space and bifurcation diagrams from
% satcompet.c
clear
clc
%% Import data
name='./Data/e_gdx_'; % commensalism
eh=importdata(strcat(name,'h.dat')); % H (d_x,g)*
ex=importdata(strcat(name,'x.dat')); % X (d_x,g)*
ey=importdata(strcat(name,'y.dat')); % Y (d_x,g)*
es=importdata(strcat(name,'s.dat')); % S (d_x,g)*

name='./Data/c_gdx_'; % coexistence
ch=importdata(strcat(name,'h.dat')); % H (d_x,g)*
cx=importdata(strcat(name,'x.dat')); % X (d_x,g)*
cy=importdata(strcat(name,'y.dat')); % Y (d_x,g)*
cs=importdata(strcat(name,'s.dat')); % S (d_x,g)*

name='./Data/b_gdx_'; % bistability
bh=importdata(strcat(name,'h.dat')); % H (d_x,g)*
bx=importdata(strcat(name,'x.dat')); % X (d_x,g)*
by=importdata(strcat(name,'y.dat')); % Y (d_x,g)*
bs=importdata(strcat(name,'s.dat')); % S (d_x,g)*
%% Parameters
pars.g=1;	
pars.d=0.05;	
pars.px=0.1;	
pars.py=0.1;	
pars.dx=0.25;	
pars.dy=0.25;	
pars.ps=0.8;
pars.ds=0.3;
pars.psy=0.04;

%% Create the canvas
size=length(es);
scale=linspace(0,0.5,size);
lim=0.0001;
Eh=(eh>lim).*7;
Ex=(ex>lim).*13;
Ey=(ey>lim).*20;
Es=(es>lim).*5;

Bh=(bh>lim).*7;
Bx=(bx>lim).*13;
By=(by>lim).*20;
Bs=(bs>lim).*5;

Ch=(ch>lim).*7;
Cx=(cx>lim).*13;
Cy=(cy>lim).*20;
Cs=(cs>lim).*5;
%% Figure 1 - Comensalism
pars.pY=0.076666;
figure(1);
set(gcf, 'Position',  [200, 200, 1200, 450])% set position, width and height of plot
subplot(1,2,1);
contourf(scale,scale,(Eh+Ex+Ey+Es));ylabel('d_{x}');xlabel('g');title('y ic');caxis([0 1].*40);colormap('jet')

hold on;
% vlim dx=0.25
plot(scale,0.25.*ones(1,size),'--r','linewidth',2)
%hlim g/d=d+dx/px & g/d=d+dy/py
plot(linspace(pars.d^2/pars.px,0.15,10),linspace(pars.d^2/pars.px,0.15,10).*pars.px./pars.d -pars.d, '--k','linewidth',2)
plot(ones(1,10).*0.15,linspace(0.25,0.5,10), '--k','linewidth',2)
%hlim g=
g_hlim=((pars.d+pars.ds-pars.psy*(pars.d+pars.dy)/pars.py)*pars.py/pars.ps+pars.d)*(pars.d+pars.dy)/pars.py;
plot(ones(1,size).*g_hlim,scale,'--b','linewidth',2)
hold off;

xlabel('g')
ylabel('d_x')
title('(a)R_c=R_y')
set(gca, 'fontname','helvetica','fontsize', 16,'linewidth',1.5)

subplot(1,2,2);
hold on;
plot(scale, eh(:,end),'Color',[0.3,0.3,0.8],'linewidth',2);
plot(scale, ex(:,end),'Color',[0.7,1.0,0.3],'linewidth',2);
plot(scale, ey(:,end),'Color',[1.0,0.8,0.0],'linewidth',2);
plot(scale, es(:,end),'Color',[1.0,0.4,0.0],'linewidth',2);
hold off;
xlim([0.1,0.5]);
ylim([0,3.5]);
box on;
xlabel('d_x')
ylabel('Count')
title('(b) bifurcation map')
legend('H^*','X^*','Y^*','S^*','Location','Best');
set(gca, 'fontname','helvetica','fontsize', 16,'linewidth',1.5)
saveas(gcf,'./Figures/commensalism.pdf','pdf');
%% Figure 2 - Mutualism (Bistability)
figure(2);
pars.pY=0.2;

set(gcf, 'Position',  [200, 200, 1200, 450])% set position, width and height of plot

subplot(1,2,1);
contourf(scale,linspace(0,0.5,481),(Bh+Bx+By+Bs));ylabel('d_{x}');xlabel('g');title('y ic');caxis([0,1]*40);colormap('jet')
hold on;

% vlim dx=0.25
plot(scale,0.25.*ones(1,size),'--r','linewidth',2)
%hlim g/d=d+dx/px & g/d=d+dy/py
plot(linspace(pars.d^2/pars.px,0.15,10),linspace(pars.d^2/pars.px,0.15,10).*pars.px./pars.d -pars.d, '--k','linewidth',2)
plot(ones(1,10).*0.15,linspace(0.25,0.5,10), '--k','linewidth',2)
%hlim g=
g_hlim=((pars.d+pars.ds-pars.psy*(pars.d+pars.dy)/pars.py)*pars.py/pars.ps+pars.d)*(pars.d+pars.dy)/pars.py;
plot(ones(1,size).*g_hlim,scale,'--b','linewidth',2)
hold off;
xlabel('g')
ylabel('d_x')
title('a R_c > R_y')
set(gca, 'fontname','helvetica','fontsize', 16,'linewidth',1.5)

subplot(1,2,2);
hold on;
plot(linspace(0,0.5,481), bh(:,end),'Color',[0.3,0.3,0.8],'linewidth',2);
plot(linspace(0,0.5,481), bx(:,end),'Color',[0.7,1.0,0.3],'linewidth',2);
plot(linspace(0,0.5,481), by(:,end),'Color',[1.0,0.8,0.0],'linewidth',2);
plot(linspace(0,0.5,481), bs(:,end),'Color',[1.0,0.4,0.0],'linewidth',2);
hold off;
xlim([0.1,0.5])
ylim([0,3.5])
box on;
xlim([0.1,0.5])
xlabel('d_x')
ylabel('Count')
title('b bifurcation map')
legend('H^*','X^*','Y^*','S^*','Location','Best');
set(gca, 'fontname','helvetica','fontsize', 16,'linewidth',1.5)
saveas(gcf,'./Figures/mutualism.pdf','pdf');
%% Figure 3 - Parasitism (coexistence)
figure(3);
pars.pY=0.016;

set(gcf, 'Position',  [200, 200, 1200, 450])% set position, width and height of plot

subplot(1,2,1);
contourf(scale,scale,(Ch+Cx+Cy+Cs));ylabel('d_{x}');xlabel('g');title('y ic');caxis([0 1].*40);colormap('jet')

hold on;
% vlim dx=0.25
plot(scale,0.25.*ones(1,size),'--r','linewidth',2)
%hlim g/d=d+dx/px & g/d=d+dy/py
plot(linspace(pars.d^2/pars.px,0.15,10),linspace(pars.d^2/pars.px,0.15,10).*pars.px./pars.d -pars.d, '--k','linewidth',2)
plot(ones(1,10).*0.15,linspace(0.25,0.5,10), '--k','linewidth',2)
%hlim g=
g_hlim=((pars.d+pars.ds-pars.psy*(pars.d+pars.dy)/pars.py)*pars.py/pars.ps+pars.d)*(pars.d+pars.dy)/pars.py;
plot(ones(1,size).*g_hlim,scale,'--b','linewidth',2)
hold off;

xlabel('g')
ylabel('d_x')
title('(a) R_c < R_y')
set(gca, 'fontname','helvetica','fontsize', 16,'linewidth',1.5)

subplot(1,2,2);
hold on;
plot(scale, ch(:,end),'Color',[0.3,0.3,0.8],'linewidth',2);
plot(scale, cx(:,end),'Color',[0.7,1.0,0.3],'linewidth',2);
plot(scale, cy(:,end),'Color',[1.0,0.8,0.0],'linewidth',2);
plot(scale, cs(:,end),'Color',[1.0,0.4,0.0],'linewidth',2);
hold off;
box on;
xlim([0.1,0.5])
ylim([0,5.5])
xlabel('d_x')
ylabel('Count')
title('(b) bifurcation map')
legend('H^*','X^*','Y^*','S^*','Location','Best');
set(gca, 'fontname','helvetica','fontsize', 16,'linewidth',1.5);
saveas(gcf,'./Figures/parasitism.pdf','pdf');