clear;
clc

M=7;

close all

kh=linspace(0.001,0.62*pi,M+1);


AA=zeros(M+1,M+1);
b=zeros(M+1,1);

% %二阶导数的系数
for ii=1:M+1 %kx=1
    AA(ii,1)=( 2*cos(kh(ii))-2);
    for kk=1:M
        AA(ii,kk+1)=2*cos(kk*kh(ii))*( 2*cos(kh(ii))-2);
    end
    
    b(ii)=(kh(ii))^4;
end


c=AA\b;%求系数
length(c)
digits(6)
vpa(c)'

figure;
kh=linspace((pi)/(100),(pi),100);


for  jj=1:4
    xita=15*(jj-1)*pi/180;
    
    a=c(1)*( 2*cos(kh*sin(xita))-2);
    for m=2:length(c)
        a=a+2*c(m)*cos((m-1)*kh*cos(xita)).*( 2*cos(kh*sin(xita))-2);
    end
    
    vv=a-kh.^4*sin(xita)^2*cos(xita)^2;
    if jj==1
        hold on; plot(vv,'b','LineWidth',2);
    elseif jj==2
        hold on;plot(vv,'m','LineWidth',2);
    elseif jj==3
        hold on;plot(vv,'r','LineWidth',2);
    elseif jj==4
        hold on;plot(vv,'k','LineWidth',2);
    end
    
    
end
grid on
legend('\theta=0','\theta=\pi/12','\theta=\pi/6','\theta=\pi/4')
xlabel('percentage of kh')
ylabel('{\itE}_3')
grid on
box on




figure;
for  jj=1:4
    xita=15*(jj-1)*pi/180;
    
    a=c(1);
    for m=2:length(c)
        a=a+2*c(m)*cos((m-1)*kh*sin(xita));
    end
    
    vv=a.*( 2*cos(kh*cos(xita))-2)-kh.^4*sin(xita)^2*cos(xita)^2;
    if jj==1
        hold on; plot(vv,'b-.','LineWidth',2);
    elseif jj==2
        hold on;plot(vv,'m-.','LineWidth',2);
    elseif jj==3
        hold on;plot(vv,'r-.','LineWidth',2);
    elseif jj==4
        hold on;plot(vv,'k-.','LineWidth',2);
    end
    
    
end
grid on
legend('\theta=0','\theta=\pi/12','\theta=\pi/6','\theta=\pi/4')
xlabel('percentage of kh')
ylabel('{\itE}_4')
grid on
box on



figure;
for  jj=1:40
    xita=15*(jj-1)*pi/180;
    
    a=c(1);
    for m=2:length(c)
        a=a+2*c(m)*cos((m-1)*kh*sin(xita));
    end
    
    
      b=c(1);
    for m=2:length(c)
        b=b+2*c(m)*cos((m-1)*kh*cos(xita));
    end
    
    vv=a.*( 2*cos(kh*cos(xita))-2)+ b.*( 2*cos(kh*sin(xita))-2)-2* kh.^4*sin(xita)^2*cos(xita)^2;
    vv=-vv./2;
    
    
    if jj==1
        hold on; plot(vv,'b-.','LineWidth',2);
    elseif jj==2
        hold on;plot(vv,'m-.','LineWidth',2);
    elseif jj==3
        hold on;plot(vv,'r-.','LineWidth',2);
    elseif jj==4
        hold on;plot(vv,'k-.','LineWidth',2);
    end
    
    
end
grid on
legend('\theta=0','\theta=\pi/12','\theta=\pi/6','\theta=\pi/4')
xlabel('percentage of kh')
% ylabel('(E3+E4)/2')
ylabel('({\itE}_3+{\itE}_4)/2')
grid on
box on
axis([0 100 -7 1])
% axis([0 100 -0.007 0.001])