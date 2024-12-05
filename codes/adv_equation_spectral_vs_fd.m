%FFT and FD based methods for advection equation with periodic BCs
clear all;close all

%define a spatial grid
xmin = 0;
xmax = 50;
N = 2^8;
%Nx is the actual number of points, Nx is the doubled grid for the BCs
x = linspace(xmin,xmax,N+1);
x=x(1:end-1);
dx=x(2)-x(1);

% define wave speed
c=1; cfl=dx/c;

% For FD
xfd = linspace(xmin,xmax,N+1);
dxfd=xfd(2)-xfd(1); dxfd2=dxfd*dxfd;
e=ones(N+1,1);
Dx = spdiags([-e e], -1:0, N+1, N+1);
Dx=(1/dxfd)*Dx;
% for periodic BCs
Dx(1,end)=-Dx(1,1);

%make initial condition
%u0 = [-sech((xn+25)) sech((x-25))];
u0= sech(x-0.5*xmax);
un=u0;
unfd=sech(xfd-0.5*xmax)';

%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2=ks.*ks; ks3=ks2.*ks;

%time step and number of steps
dt = 0.1*cfl; dt2=dt*dt;
numsteps=15;
numouts=50;

figure(1)
clf
plot(x,u0,'k-','linewidth',2)
grid on
xlabel('x','fontweight','bold','fontsize',12);
ylabel('u','fontweight','bold','fontsize',12);
axis([0 xmax -0.05 1.05])
drawnow
t=0;
% first order time stepping
kappa_num=1e-2; % this is artificial viscosity to stabilize the method
diff_fact=1./(1+dt*kappa_num*ks2); % this is the implict diffusion factor
unf=fft(un); % evolve in spectral space
for ii=1:numouts
 for jj=1:numsteps
    t=t+dt;
% spectral
    unf = (unf-dt*c*sqrt(-1)*ks.*unf).*diff_fact;
% FD
    unfd = unfd-dt*c*Dx*unfd;
 end
 un=real(ifft(unf)); % transform back to physical space only for graphics
figure(1)
clf
 set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
plot(x,un,'k-',xfd,unfd,'r--',x,u0,'b-','linewidth',2)
grid on
xlabel('x','fontweight','bold','fontsize',12);
ylabel('u','fontweight','bold','fontsize',12);
title(['time = ' num2str(t,2)]);
axis([0 xmax -0.05 1.05])
legend('spectral','FD','ICs','Location','northwest')
drawnow
end
