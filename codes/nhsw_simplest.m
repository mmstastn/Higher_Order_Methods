% 1D SW equations in periodic deomain in physical formulation
% not in conservative form, simplest case
clear all,close all
g=0.01*9.81; %reduced gravity
H=10; % effective layer depth
H2o6=H^2/6; % the NH factor
clin=sqrt(g*H); % linear wave speed
f=1e-4; % rotation
N=8*1024; % num points
L=8*2.048e3; % domain half-width, domain is [-L,L]

% grid and Fourier
x=((0.5:N-0.5)/N)*2*L-L; dx=x(2)-x(1);
dk=pi/L;
ksvec=zeros(size(x));
ksvec(1)=0; ksvec(N/2+1)=0;
for ii=2:(N/2)
   ksvec(ii)=ii-1;
   ksvec(N/2+ii)=-N/2 + ii -1;
end
k=ksvec'*dk; k2=k.*k; ik=sqrt(-1)*k;
nhfact=1./(1+k2*H*H/6);

%Build our `radial' filter
filtorder=8;
kmax=max(k(:));
cutoff=0.66;
kcrit=kmax*cutoff;
alpha1=1;
alpha2=1; 
kmag=sqrt(k.*k);
myfilter=exp(-alpha1*(kmag/(alpha2*kcrit)).^filtorder);

% time stepping info
cfl=dx/clin;
numsteps=200;
numouts=100;
dt=0.5*cfl;

% ICs
um1=zeros(N,1);
vm1=zeros(N,1);
zm1=H+0.1*H*sech(x/(0.1*L)).^2';
t=0;
% One step of Euler
t=t+dt;
uf1=fft(um1);zf1=fft(zm1);
uhef1=fft(um1.*zm1);uuf1=0.5*fft(um1.^2);
vf1=fft(vm1);

un1=um1+dt*(real(ifft((-i*k.*uuf1-i*g*k.*zf1+f*vf1))));
vn1=vm1+dt*(um1.*real(ifft(-i*k.*vf1))-f*um1);
zn1=zm1+dt*(real(ifft(-i*k.*uhef1)));
% Main loop uses leapfrog for now
% notation is up1 is u at n+1, un1 is u now, um1 is u at n-1
for ii=1:numouts
 for jj=1:numsteps
     t=t+dt;
     uf1=fft(un1);zf1=fft(zn1);
     uhef1=fft(un1.*(H+(zn1-H)));uuf1=0.5*fft(un1.^2);
     vf1=fft(vn1);
     phterm=real(ifft(i*g*k.*zf1));
     up1=um1+2*dt*(real(ifft(nhfact.*(-i*k.*uuf1-i*g*k.*zf1+f*vf1))));
     vp1=vm1+2*dt*(un1.*real(ifft(-i*k.*vf1))-f*un1);
     zp1=zm1+2*dt*(real(ifft(-i*k.*uhef1)));
     dummy=fft((up1-um1)/(2*dt));
     pnhterm=real(ifft(-H2o6*k2.*dummy));
     um1=un1;vm1=vn1;zm1=zn1;
     un1=real(ifft(myfilter.*fft(up1)));
     vn1=real(ifft(myfilter.*fft(vp1)));
     zn1=real(ifft(myfilter.*fft(zp1)));
 end
 figure(1)
 clf
 betterplots
 subplot(3,1,1)
 plot(x,(zn1-H)/H,'k'),grid on
 axis([-L L -0.15 0.15])
 ylabel('eta/H')
 title(['t = ' num2str(t,4) ' s']);
 subplot(3,1,2)
 plot(x,un1/clin,'b',x,vn1/clin,'r'),grid on
 axis([-L L -0.15 0.15])
 ylabel('(u,v)/c')
 subplot(3,1,3)
 plot(x,phterm,'b',x,pnhterm,'r'),grid on
 xlabel('x')
 ylabel('ph_x and pnh_x')
 axis([-L L -0.01 0.01])
 drawnow
end