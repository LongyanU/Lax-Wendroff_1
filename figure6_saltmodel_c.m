%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%            High-order time discretizations in seismic modeling
%%%%%%%%%                作者：梁文全  Lax-Wendroff methods   %%%%%%%%%%%%%%%%%

% 时间已过 173.654854 秒。
% 时间已过 170.721254 秒。
clear
clc %%%%%%%
close all
nt=2500;    % number of time steps
eps=.6;     % stability
isnap=30;    % snapshot sampling

load('vv.mat')
c1=flipud(c);
c=c1;

[nz,nx]=size(c1);
v=c1;


nz=nz+50;

vv=zeros(nz,nx);
for ii=1:nz-50
    for jj=1:nx
        vv(ii+50,jj)=v(ii,jj);
    end
end

for ii=1:50  %%top
    for jj=1:nx
        vv(ii,jj)=vv(51,jj);
    end
end



clear v
v=vv;


p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.002; % calculate time step from stability criterion
tau=dt;
% r2=v.*v.*dt*dt/dx/dx;
v2=v.*v;

f0=38
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff(diff(src))/dx^2;					% time derivative to obtain gaussian



seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
% zs=nz/2;
% xs=nx/2;

zs=51;
xs=400;







p=zeros([nz nx]);
coeff1= [ 1.83806, -0.355029, 0.100949, -0.0291963, 0.00730575, -0.00135064, 0.000135486];
temp1= -3.12175;

% coeff2=[3.17198, -0.983605, 0.303994, -0.0846446, 0.0191256, -0.00304942, 0.000254116];
% temp2= -4.8481;

coeff2=[3.25708, -1.04773, 0.343324, -0.103701, 0.0260243, -0.00471958, 0.00045775];
temp2=-4.94149;

% temp2= -4
tic
for it=1:nt-2,
    
    p(zs,xs)= p(zs,xs)+src(it)*dt^2;
    d2px11=(circshift(p,[0 -1])+circshift(p,[0 1]));
    d2px12=(circshift(p,[0 -2])+circshift(p,[0 2]));
    d2px13=(circshift(p,[0 -3])+circshift(p,[0 3]));
    d2px14=(circshift(p,[0 -4])+circshift(p,[0 4]));
    d2px15=(circshift(p,[0 -5])+circshift(p,[0 5]));
    d2px16=(circshift(p,[0 -6])+circshift(p,[0 6]));
    d2px17=(circshift(p,[0 -7])+circshift(p,[0 7]));
    
    d2pz11=(circshift(p,[-1 0])+circshift(p,[1 0]));
    d2pz12=(circshift(p,[-2 0])+circshift(p,[2 0]));
    d2pz13=(circshift(p,[-3 0])+circshift(p,[3 0]));
    d2pz14=(circshift(p,[-4 0])+circshift(p,[4 0]));
    d2pz15=(circshift(p,[-5 0])+circshift(p,[5 0]));
    d2pz16=(circshift(p,[-6 0])+circshift(p,[6 0]));
    d2pz17=(circshift(p,[-7 0])+circshift(p,[7 0]));
    
    d2px=coeff1(1)*d2px11+coeff1(2)*d2px12+coeff1(3)*d2px13+coeff1(4)*d2px14+coeff1(5)*d2px15+coeff1(6)*d2px16...
        +coeff1(7)*d2px17;
    d2pz=coeff1(1)*d2pz11+coeff1(2)*d2pz12+coeff1(3)*d2pz13+coeff1(4)*d2pz14+coeff1(5)*d2pz15+coeff1(6)*d2pz16...
        +coeff1(7)*d2pz17;
    
    
    
    d2pxz=v2.*(d2px+d2pz+2*temp1*p)/h^2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lax-Wendroff methods
    
    d2px2=coeff2(1)*d2px11+coeff2(2)*d2px12+coeff2(3)*d2px13+coeff2(4)*d2px14+coeff2(5)*d2px15+coeff2(6)*d2px16...
        +coeff2(7)*d2px17;
    d2pz2=coeff2(1)*d2pz11+coeff2(2)*d2pz12+coeff2(3)*d2pz13+coeff2(4)*d2pz14+coeff2(5)*d2pz15+coeff2(6)*d2pz16...
        +coeff2(7)*d2pz17;
    
    d2pxz2=v2.*(d2px2+d2pz2+2*temp2*p)/h^2;
    
    d2px4=(circshift(d2pxz2,[0 -1])+circshift(d2pxz2,[0 1]));
    d2pz4=(circshift(d2pxz2,[-1 0])+circshift(d2pxz2,[1 0]));
    
    d2pxz3=v2.*(d2px4+d2pz4-4*d2pxz2)*dt^2/12/h^2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lax-Wendroff methods
    pnew=2*p-pold+(d2pxz*dt^2+d2pxz3*dt^2);
    
    pold=p;
    p=pnew;
    
%     [pold,p]=spongeABC(pold,p,nx,nz,50,50,0.005);
     [pold,p]=spongeABC(pold,p,nx,nz,50,50,0.007);
    seis_record(it,:)=p(zs,:);
    
    if rem(it,isnap)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    
    
end
toc
save('figure6c.mat')
