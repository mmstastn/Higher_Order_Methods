%% This script shows some sample filters for the inviscid Burgers
clear all,close all
%define a spatial grid
xmin = -1;
xmax = 1;
Lx=xmax-xmin;
N = 2^10; % try 2^8
x = linspace(xmin,xmax,N);
dx = x(2)-x(1);

%make initial condition
% packet
%u0 = 5*sech(x/(0.1*xmax)).*cos(4*pi*x/xmax);
% soliton-like
%u0 = 5*sech(x/(0.1*xmax)).^8;
% boxcar like
u0=0.5*(1+tanh((x+0.25*Lx)/(0.05*Lx)));
u0=u0-0.5*(1+tanh((x-0.25*Lx)/(0.05*Lx)));
u1n=u0;u2n=u0;u3n=u0;
%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;

% define the filter
% Here I use hypervisocsity so the order should be even
% The parameter can be somewhat automated but here I tune by hand
 filtord=8; eps=1e-6;
 hypervisc=(1-eps)/(eps*max(abs(ks)).^filtord);
 myfiltu1=1./(1+hypervisc*ks.^filtord);
valend=1e-12;
filtalpha=-log(valend);
knyq=max(abs(ks));
kcut=0.35*knyq;
filtbeta=8;
dummy=ones(size(ks));
myfiltu2=dummy.*(abs(ks)<kcut)+exp(-filtalpha*(((abs(ks)-kcut)/(knyq-kcut)).^filtbeta)).*(abs(ks)>=kcut);
valend=1e-12;
filtalpha=-log(valend);
knyq=max(abs(ks));
kcut=0.35*knyq;
filtbeta=4;
dummy=ones(size(ks));
myfiltu3=dummy.*(abs(ks)<kcut)+exp(-filtalpha*(((abs(ks)-kcut)/(knyq-kcut)).^filtbeta)).*(abs(ks)>=kcut);

dt=1e-4;t=0;
numstps=200;
numouts=30;
% one step of backward Euler
u1p = real(ifft((fft(u1n) + 0.5*dt*i*ks.*fft(u1n.*u1n))));
u2p = real(ifft((fft(u2n) + 0.5*dt*i*ks.*fft(u2n.*u2n))));
u3p = real(ifft((fft(u3n) + 0.5*dt*i*ks.*fft(u3n.*u3n))));

%leapfrog time stepping
for ii=1:numouts
   for jj=1:numstps 
    t=t+dt;
    u1f = real(ifft(myfiltu1.*(fft(u1p) - dt*i*ks.*fft(u1n.*u1n))));
    u2f = real(ifft(myfiltu2.*(fft(u2p) - dt*i*ks.*fft(u2n.*u2n))));
    u3f = real(ifft(myfiltu3.*(fft(u3p) - dt*i*ks.*fft(u3n.*u3n))));
    % rotate
    u1p=u1n;u1n=u1f;
    u2p=u2n;u2n=u2f;
    u3p=u3n;u3n=u3f;   
   end 
   figure(1)
   clf
   betterplots
   subplot(2,1,1)
   plot(x,u1n,'k-',x,u2n,'b-',x,u3n,'r-',x,u0,'g-')
   grid on
   xlabel('x');
   ylabel('u');
   title(['time = ' num2str(t,2)]);
   axis([xmin xmax -max(max(u0)) 1.4*max(max(u0))])
   subplot(2,1,2)
   plot(x,u1n,'k-',x,u2n,'b-',x,u3n,'r-')
   grid on
   xlabel('x');
   ylabel('u (detail)');
   %title(['time = ' num2str(t,2)]);
   axis([0.7 0.9 -0.2 1.2])

   %legend('initial','implicit','explicit')
   % subplot(2,1,2)
   % u1f=fft(u1);spec1=log10(u1f.*conj(u1f));u2f=fft(u2);spec2=log10(u2f.*conj(u2f));
   % plot(fftshift(ks),fftshift(spec0),'k.',fftshift(ks),fftshift(spec1),'b.',fftshift(ks),fftshift(spec2),'r.')
   % xlabel('k');
   % ylabel('log10 of spectrum');
   drawnow 
end