
clear;
clc
close all;

load('figure3a.mat')
seis_recordaa=seis_record;

plotimage(p1(51:end,50:end-50))
xlabel('x/dx')
ylabel('z/dz')
title('')

load('figure3b.mat')
seis_recordbb=seis_record;
plotimage(p1(51:end,50:end-50))
xlabel('x/dx')
ylabel('z/dz')
title('')



load('figure3c.mat')
seis_recordcc=seis_record;

NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

plotimage(p1(51:end,50:end-50))
xlabel('x/dx')
ylabel('z/dz')
title('')



plotimage((1:1000)*2,seis_recordaa(:,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

plotimage((1:1000)*2,seis_recordbb(:,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

plotimage((1:1000)*2,seis_recordcc(:,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')


figure;
plot((1:1000)*2,seis_recordaa(:,300+50),'b')


hold on;plot((1:1000)*2,I *seis_recordaa(:,300+50),'b-.')

hold on;plot((1:1000)*2,seis_recordbb(:,300+50),'k')
hold on;plot((1:1000)*2,seis_recordcc(:,300+50),'r')

legend('Tra FD sheme(Without using Lax-Wendroff method)','TDE(Time Dispersion Elimination)','Tra FD scheme for Lax-Wendroff({\itM}_1={\itM}_2=7)','Proposed FD scheme for Lax-Wendroff({\itM}_1=7,{\itM}_2=1)')
% 

grid on
xlabel('Time(ms)')
ylabel('Amp')