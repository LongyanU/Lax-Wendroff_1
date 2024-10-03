clear;
clc

M=7;

close all

kh=linspace(0.001,0.62*pi,M+1);


AA=zeros(M+1,M+1);
b=zeros(M+1,1);

% %二阶导数的系数
for ii=1:M+1 %kx=1
    AA(ii,1)=1;
    for kk=1:M
        AA(ii,kk+1)=2*cos(kk*kh(ii));
    end
    
    b(ii)=-(kh(ii))^2;
end


c=AA\b;%求系数
length(c)
digits(6)
vpa(c)'

figure;
kh=linspace((pi)/(100),(pi),100);
a=c(1);
for m=2:length(c)
    a=a+2*c(m)*cos((m-1)*kh);
end



vv=a.*a-kh.^4;
hold on; plot(vv,'r','LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kh=linspace(0.001,0.62*pi,M+1);


AA=zeros(M+1,M+1);
b=zeros(M+1,1);

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
% c=[ -4.8481, 3.17198, -0.983605, 0.303994, -0.0846446, 0.0191256, -0.00304942, 0.000254116];

kh=linspace((pi)/(100),(pi),100);


a=c(1)*( 2*cos(kh)-2);

for m=2:length(c)
    a=a+2*c(m)*cos((m-1)*kh).*( 2*cos(kh)-2);
end

vv=a-kh.^4;

hold on; plot(vv,'k','LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xlabel('percentage of kh')
ylabel('{\itE}_1')
grid on
box on
legend('Tra FD scheme for Lax-Wendroff({\itM}_1={\itM}_2=7)','Proposed FD scheme for Lax-Wendroff({\itM}_1=1,{\itM}_2=7)')
axis([0 100 -0.007 0.0019])