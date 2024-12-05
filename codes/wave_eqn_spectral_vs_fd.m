%FFT based method for wave equation with periodic BCs
clear all;close all

%define a spatial grid
xmin = 0;
xmax = 50;
N = 2^7;
%Nx is the actual number of points, Nx is the doubled grid for the BCs
x = linspace(xmin,xmax,N+1);
x=x(1:end-1);
dx=x(2)-x(1);

% define wave speed
c=1; c2=c*c;

% For FD
xfd = linspace(xmin,xmax,N+1);
dxfd=xfd(2)-xfd(1); dxfd2=dxfd*dxfd;
e=ones(N+1,1);
Dxx = spdiags([e -2*e e], -1:1, N+1, N+1);
Dxx=(1/dxfd2)*Dxx;
% for periodic BCs
Dxx(1,end)=Dxx(1,2);
Dxx(N+1,1)=Dxx(N+1,N);

%make initial condition
%u0 = [-sech((xn+25)) sech((x-25))];
u0= sech(x-0.5*xmax);
un=u0;
up = un;
uf = zeros(size(un));
unfd=sech(xfd-0.5*xmax)';
upfd=unfd;
uffd=zeros(N+1,1);

%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2=ks.*ks; ks3=ks2.*ks;

%time step and number of steps
dt = 0.01; dt2=dt*dt;
numsteps=15;
numouts=250;

figure(1)
clf
plot(x,u0,'k-','linewidth',2)
grid on
xlabel('x','fontweight','bold','fontsize',12);
ylabel('u','fontweight','bold','fontsize',12);
axis([0 xmax -1 1])
drawnow
t=0;
% second order time stepping
for ii=1:numouts
 for jj=1:numsteps
    t=t+dt;
% spectral
    uf = 2*un-up+dt2*c2*ifft(-ks2.*fft(un)); 
    % strictly speaking we don't need to transform back
% FD
    uffd = 2*unfd-upfd+dt2*c2*Dxx*unfd;
    % update
    up=un; un=uf;
    upfd=unfd; unfd=uffd;
 end
    
figure(1)
clf
 set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
plot(x,un,'k-',xfd,unfd,'r--','linewidth',2)
grid on
xlabel('x','fontweight','bold','fontsize',12);
ylabel('u','fontweight','bold','fontsize',12);
title(['time = ' num2str(t,2)]);
axis([0 xmax -0.05 1.05])
legend('spectral','FD')
drawnow
end
