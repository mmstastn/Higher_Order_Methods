%% This script shows some sample filters
%define a spatial grid
xmin = -1;
xmax = 1;
Lx=xmax-xmin;
N = 2^9; % try 2^8
x = linspace(xmin,xmax,N);
dx = x(2)-x(1);

%make initial condition
% boxcar
u=0.5*(1+tanh((x+0.25*Lx)/(0.05*Lx)));
u=u-0.5*(1+tanh((x-0.25*Lx)/(0.05*Lx)))
% packet
%u = 5*sech(x/(0.1*xmax)).*cos(4*pi*x/xmax);
% soliton-like
%u = 5*sech(x/(0.1*xmax)).^8;
% boxcar like
u=0.5*(1+tanh((x+0.25*Lx)/0.01*Lx));
u=u-0.5*(1+tanh((x-0.25*Lx)/0.01*Lx));

%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;

% define the filter
% Here I use hypervisocsity so the order should be even
% The parameter can be somewhat automated but here I tune by hand
 valend=1e-2;
 filtord=8; eps=1e-6;
 hypervisc=(1-eps)/(eps*max(abs(ks)).^filtord);
 myfiltu1=1./(1+hypervisc*ks.^filtord)
% filtord=4; eps=1e-6;
% hypervisc=(1-eps)/(eps*max(abs(ks)).^filtord);
% myfiltT=1./(1+hypervisc*ks.^filtord)
% filtord=8; eps=1e-7;
% hypervisc=(1-eps)/(eps*max(abs(ks)).^filtord);
% myfiltS=1./(1+hypervisc*ks.^filtord)
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

% plots
figure(1)
clf
betterplots
subplot(3,1,1)
plot(x,u,'b.'),xlabel('x'),ylabel('u')
subplot(3,1,2)
plot(ks,myfiltu1,'k.',ks,myfiltu2,'b.',ks,myfiltu3,'r.')
xlabel('k'),ylabel('filter')
legend('hyperviscosity','exponential 1','exponential 2')
ufilt1=real(ifft(myfiltu1.*fft(u)));
ufilt2=real(ifft(myfiltu2.*fft(u)));
ufilt3=real(ifft(myfiltu3.*fft(u)));
subplot(3,1,3)
plot(x,ufilt1,'k',x,ufilt2,'b:',x,ufilt3,'r--')
xlabel('x'),ylabel('u fitered (detail)')
axis([0.15*Lx 0.35*Lx -0.1 1.1])

