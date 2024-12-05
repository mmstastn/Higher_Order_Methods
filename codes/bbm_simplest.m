%FFT based method for BBM equation.
% Uses explicit scheme

clear all; close all
%Physical parameters
H=100;
h1=20; %upper layer depth
h2=H-h1; %lower layer depth
g=9.81; drho=0.005; gp=g*drho;
ctwolayer=sqrt(gp*h1*h2/H);
betatwolayer=ctwolayer*h1*h2/6;
alphatwolayer=1.5*ctwolayer*(h1-h2)/(h1*h2);

%define a spatial grid
xmin = -0.5e4;
xmax = 0.5e4;
N = 2^10;
x = linspace(xmin,xmax,N+1); x=x(1:end-1);
dx = x(2)-x(1);

%make initial condition
B0 = -0.1*H*sech((x+0.5*xmax)/(0.035*xmax)).^2;
B1 = B0;
%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2=ks.*ks; ks3=ks2.*ks;

%time step and number of steps
t=0; cfl=dx/ctwolayer;
dt = 0.0025*cfl;
numstps=2000;
numouts=200;

figure(1)
clf
% this sets the thick line plotting parameters useful for hard copies and journals
 set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');

bbmfact=1./(1+(betatwolayer/ctwolayer).*ks2);
%start with Euler time stepping
for ii=1:numouts
   for jj=1:numstps 
    t=t+dt;
    % explicit method for BBM
    B1lin=-ctwolayer*sqrt(-1)*ks.*fft(B1);
    B1nl=-0.5*alphatwolayer*sqrt(-1)*ks.*fft(B1.^2);
    B1 = B1+dt*real(ifft(bbmfact.*(B1lin+B1nl)));
   end 
   figure(1)
   plot(x,B1,'b-',x,B0)
   grid on
   xlabel('x');
   ylabel('B');
   title(['time = ' num2str(t,2)]);
   axis([xmin xmax -20 1])
   drawnow
end
