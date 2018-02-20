%generating parameters distribution for control case
clear all;close all;clc;
addpath('C:\Users\Alexander\Documents\MATLAB\CuMcmcControl\mex');
addpath('C:\Users\Alexander\Documents\MATLAB\CuMcmcControl\mcmcstat');
load chain12
load results12
load modeloutput12
data.repetitions = 4;
data.experiments = 600;

 %  8 15 - overall mortality
 % 79 80 - exit
 % 53 22  -perc fed mosq
data.obs =[ 8 15 79 80 53 22];

%theta0 = [3*1e-5 1e-4 9*1e-5 1.5*1e-4 0.9995 3.1e-04 4];
%theta0 = [1e-6 2*1e-4 1e-5 8*1e-4 0.9995 3.1*1e-4];
theta0 = chain12(end,:);
res = hut_exp_ss(theta0,data)
sigma2 = 1;%std for experimental data
%%%% the MCMC part
model.ssfun=@hut_exp_ss;
model.sigma2=sigma2;
model.N = 1;

% the parameter structure: name, init.val, min, max
params = {
{'\sigma^g_{acc}(1)',theta0(1),0,0.1}
{'\sigma^g_{acc}(2)',theta0(2),0,0.1}
{'\sigma^a_{acc}(1)',        theta0(3),0,0.1}      
{'\sigma^a_{acc}(2)',        theta0(4),0,0.1} 
{'p_{net}', theta0(5), 0,1}
{'p_{hut}', theta0(6), 0,1}
};
% MCMC options
options.nsimu = 30000; % number of samples
%v = [1e-2 1e-2 1e-2 1e-2 1e-3 1e-3];
%options.qcov = diag(v.^2); % (initial) proposal covariance
%options.method = 'am'; % method: DRAM
%options.adaptint = 100; % adaptation interval

global modeloutput modeloutpc
% calling MCMCRUN
[results12_2,chain12_2,s2chain12_2,sschain12_2, hchain12_2] = mcmcrun(model,data,params,options,results12);

save('results12_2','results12_2');
save('chain12_2','chain12_2');
save('s2chain12_2','s2chain12_2');
save('sschain12_2','sschain12_2');
save('hchain12_2','hchain12_2');
modeloutput12_2 = modeloutput;
save('modeloutput12_2','modeloutput12_2');

chain = [chain12;chain12_2];
figure;
mcmcplot(chain,[],results12_2.names);
figure;
mcmcplot(chain,[],results12_2.names,'pairs');


modeloutput = [modeloutput12;modeloutput12_2];
figure;
N=length(modeloutput(:,3));
subplot(2,2,1);
hold on
 plot(1:N, modeloutput(:,3), 'r-') ; 
 plot(1:N, data.obs(3), 'k-') ;
 hold off
 title('Exit Gambiae');
 axis tight
 subplot(2,2,2); 
hold on
plot(1:N, modeloutput(:,4), 'r-') ;
 plot(1:N, data.obs(4), 'k-') ;
 
 hold off
  title('Exit Arabiensis');
  axis tight
   subplot(2,2,3);
hold on

 plot(1:N, modeloutput(:,5), 'r-') ;
 plot(1:N, data.obs(5), 'k-') ;
 hold off
  title('Fed Gambiae');
  axis tight
     subplot(2,2,4);
hold on
 
 plot(1:N, modeloutput(:,6), 'r-') ;
 plot(1:N, data.obs(6), 'k-') ;
 hold off
  title('Fed Arabiensis');
  axis tight