%FFT based method for Burgers equation.
% Uses implicit scheme
clear all;

tic;

%define physical paramaters for the model
nu=5e-2; % both stay stable for this value
%nu=2e-1; % explicit goes unstable


%define a spatial grid
xmin = -10;
xmax = 10;
N = 2^12;
x = linspace(xmin,xmax,N);
dx = x(2)-x(1);

%make initial condition
u1 = 10*sech(x/2).*cos(5*pi*x/xmax).*cos(0.5*pi*x/xmax);
u2 = u1;%.*sin(2*pi*x/2);
u0=u1;


%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2=ks.*ks; ks3=ks2.*ks; ks4=ks2.^2;

%time step and number of steps
t=0;
dt = 1e-5; % try 1e-4 to see that explicit fails
impop=1./(1+dt*nu*ks2);
%impop=1./(1+dt*nu*ks4);

expop=-dt*ks2*nu;
numstps=500;
numouts=200;
figure(1)
clf
% this sets the thick line plotting parameters useful for hard copies and journals
 set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
subplot(2,1,1)
plot(x,u0,'k-')
subplot(2,1,2)
u0f=fft(u0);spec0=log10(u0f.*conj(u0f));
plot(fftshift(ks),fftshift(spec0),'k.-')

%start with Euler time stepping
for ii=1:numouts
   for jj=1:numstps 
    t=t+dt;
    u1 = real(ifft(impop.*(fft(u1) - dt*i*0.5*ks.*fft(u1.*u1))));
    
   end 
   figure(1)
   clf
   subplot(2,1,1)
   plot(x,u0,'k-',x,u1,'b-',x,u2,'r--')
   grid on
   xlabel('x');
   ylabel('u');
   title(['time = ' num2str(t,2)]);
   axis([xmin xmax -max(max(u0)) 1.4*max(max(u0))])
   legend('initial','implicit','explicit')
   subplot(2,1,2)
   u1f=fft(u1);spec1=log10(u1f.*conj(u1f));u2f=fft(u2);spec2=log10(u2f.*conj(u2f));
   plot(fftshift(ks),fftshift(spec0),'k.',fftshift(ks),fftshift(spec1),'b.',fftshift(ks),fftshift(spec2),'r.')
   xlabel('k');
   ylabel('log10 of spectrum');
   drawnow 
end
