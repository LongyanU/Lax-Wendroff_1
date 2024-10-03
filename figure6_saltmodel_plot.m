
clear;
clc
close all;

load('figure6a.mat')
seis_recordaa=seis_record;

load('figure6b.mat')
seis_recordbb=seis_record;

load('figure6c.mat')
seis_recordcc=seis_record;
% 
% load('figure6d.mat')
% seis_recorddd=seis_record;
% 
% load('figure6e.mat')
% seis_recordee=seis_record;
% 
% load('figure6f.mat')
% seis_recordff=seis_record;





% NT=nt-2;
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix
% u = I * record';   % <- Transforming the 1xN solution vector v

% % % 
plotimage((1:2500)*2,seis_recordaa(:,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

plotimage((1:2500)*2,seis_recordbb(:,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

plotimage((1:2500)*2,seis_recordcc(:,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')


figure;
plot((1:2500)*2,seis_recordaa(:,300+50),'b')


hold on;plot((1:2500)*2,I *seis_recordaa(:,300+50),'b-.')
% 
hold on;plot((1:2500)*2,seis_recordbb(:,300+50),'k')
hold on;plot((1:2500)*2,seis_recordcc(:,300+50),'r')

legend('Tra FD sheme(Without using Lax-Wendroff method)','TDE(Time Dispersion Elimination)','Tra FD scheme for Lax-Wendroff({\itM}_1={\itM}_2=7)','Proposed FD scheme for Lax-Wendroff({\itM}_1=7,{\itM}_2=1)')
% 
grid on
xlabel('Time(ms)')
ylabel('Amp')
