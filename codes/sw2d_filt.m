%FFT based method for 2D SW constant depth
% with different filters

clear all;

%define physical paramaters for the model
g=9.81; H=10;

%define a spatial grid
xmin = -100;
xmax = 100;
Nx = 2^8;
x = linspace(xmin,xmax,Nx);
dx = x(2)-x(1);
ymin = -100;
ymax = 100;
Ny = 2^8;
y = linspace(ymin,ymax,Ny);
dy = y(2)-y(1);
[xx,yy]=meshgrid(x,y);

%make initial condition
% Sech for traditional shock/rarefaction calculation
r=sqrt(xx.^2+yy.^2); myamp=0.8*H;
eta = myamp*sech(r/(0.1*max(r(:)))).^8;
u=zeros(size(eta));v=zeros(size(eta));
eta2 = myamp*sech(r/(0.1*max(r(:)))).^8;
u2=zeros(size(eta));v2=zeros(size(eta));

u0=u; eta0=eta;v0=v;

%make wave numbers
nyquist_freqx = 2*pi/(xmax-xmin);
ks=[0:Nx/2-1 0 -Nx/2+1:-1]*nyquist_freqx;
nyquist_freqy = 2*pi/(ymax-ymin);
ls=[0:Ny/2-1 0 -Ny/2+1:-1]*nyquist_freqy;
[kk,ll]=meshgrid(ks,ls);
kk2=kk.*kk; ll2=ll.*ll;

% define the filter functions
f_bump = @(a,b,k,x) 1.*(abs(x)<=k) + exp(a*(1-1./(1-((abs(x)-k)./(1-k)).^b))).*(abs(x)>k);
f_spec = @(a,b,k,x) 1.*(abs(x)<=k) + exp(-a*(abs(x)-k).^b./(1-k).^b).*(abs(x)>k);

valend=1e-12;
filtalpha=-log(valend);
knyq=max(abs(ks));lnyq=max(abs(ls));
kcut=0.35*knyq;lcut=0.35*lnyq;
%kcut=0.25*knyq;
filtbeta=2;
dummy=ones(size(ks));
myfiltax = f_bump(0.5, filtbeta, kcut/knyq, ks/knyq);
myfiltay = f_bump(0.5, filtbeta, lcut/lnyq, ls/lnyq);
myfilta=repmat(myfiltay',[1 Nx]).*repmat(myfiltax,[Ny 1]);
myfilta(Ny/2+1,:)=0;
myfilta(:,Nx/2+1)=0;
myfiltbx = f_spec(filtalpha, filtbeta, kcut/knyq, ks/knyq);
myfiltby = f_spec(filtalpha, filtbeta, lcut/lnyq, ls/lnyq);
myfiltb=repmat(myfiltby',[1 Nx]).*repmat(myfiltbx,[Ny 1]);
myfiltb(Ny/2+1,:)=0;
myfiltb(:,Nx/2+1)=0;


%time step and number of steps
t=0;
dt = 1e-3; twodt=2*dt;

numstps=400;
numouts=50;
%numouts=48;


%start with one step of Euler time stepping
up=u0;up2=u0;ep=eta0;ep2=eta0;vp=v0;vp2=v0;
upf=fft2(up);upf2=fft2(up2);
upx=ifft2(i*kk.*upf);upx2=ifft2(i*kk.*upf);
upy=ifft2(i*ll.*upf);upy2=ifft2(i*ll.*upf2);
vpf=fft2(vp);vpf2=fft2(vp2);
vpx=ifft2(i*kk.*vpf);vpx2=ifft2(i*kk.*vpf);
vpy=ifft2(i*ll.*vpf);vpy2=ifft2(i*ll.*vpf2);
epf=fft2(ep);epf2=fft2(ep2);
epx=ifft2(i*kk.*epf);epx2=ifft2(i*kk.*epf2);
epy=ifft2(i*ll.*epf);epy2=ifft2(i*ll.*epf2);
eprhs=ifft2(i*kk.*fft2((H+ep).*up)+i*ll.*fft2((H+ep).*vp));
eprhs2=ifft2(i*kk.*fft2((H+ep2).*up2)+i*ll.*fft2((H+ep2).*vp2));


un=up+dt*(-up.*upx-vp.*upy-g*epx);
un2=up2+dt*(-up2.*upx2-vp2.*upy2-g*epx2);
vn=vp+dt*(-up.*vpx-vp.*vpy-g*epy);
vn2=vp2+dt*(-up2.*vpx2-vp2.*vpy2-g*epy2);
en=ep+dt*(-eprhs);
en2=ep2+dt*(-eprhs2);


for ii=1:numouts
   for jj=1:numstps 
    t=t+dt;
    % the main loop uses leapfrog
    unf=fft2(un);unf2=fft2(un2);
    unx=ifft2(i*kk.*unf);unx2=ifft2(i*kk.*unf2);
    uny=ifft2(i*ll.*unf);uny2=ifft2(i*ll.*unf2);
    vnf=fft2(vn);vnf2=fft2(vn2);
    vnx=ifft2(i*kk.*vnf);vnx2=ifft2(i*kk.*vnf2);
    vny=ifft2(i*ll.*vnf);vny2=ifft2(i*ll.*vnf2);
    enf=fft2(en);enf2=fft2(en2);
    enx=ifft2(i*kk.*enf);enx2=ifft2(i*kk.*enf2);
    eny=ifft2(i*ll.*enf);eny2=ifft2(i*ll.*enf2);
    enrhs=ifft2(i*kk.*fft2((H+en).*un)+i*ll.*fft2((H+en).*vn));
    enrhs2=ifft2(i*kk.*fft2((H+en2).*un2)+i*ll.*fft2((H+en2).*vn2));
     
    uf=up+twodt*(-un.*unx-vn.*uny-g*enx);
    uf2=up2+twodt*(-un2.*unx2-vn2.*uny2-g*enx2);
    vf=vp+twodt*(-un.*vnx-vn.*vny-g*eny);
    vf2=vp2+twodt*(-un2.*vnx2-vn2.*vny2-g*eny2);
    ef=ep+twodt*(-enrhs);
    ef2=ep2+twodt*(-enrhs2);
    
    % now rotate and filter as you go
    up=un; up2=un2; ep=en; ep2=en2; vp=vn; vp2=vn2;    
    un=real(ifft2(myfilta.*fft2(uf)));
    un2=real(ifft2(myfiltb.*fft2(uf2)));
    vn=real(ifft2(myfilta.*fft2(vf)));
    vn2=real(ifft2(myfiltb.*fft2(vf2)));
    en=real(ifft2(myfilta.*fft2(ef)));
    en2=real(ifft2(myfiltb.*fft2(ef2)));
   end 
   
   figure(1)
   clf
    set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
      colormap darkjet
   subplot(2,2,1)
   pcolor(xx,yy,un),shading flat
   title(num2str(t,4))
   subplot(2,2,2)
   pcolor(xx,yy,un2),shading flat
   subplot(2,2,3)
   pcolor(xx,yy,en),shading flat
   subplot(2,2,4)
   pcolor(xx,yy,en2),shading flat
   drawnow   

end
