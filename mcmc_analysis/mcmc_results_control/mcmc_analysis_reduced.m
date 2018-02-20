%generating parameters distribution for control case
clear all;close all;clc;
addpath('C:\Users\Alexander\Documents\MATLAB\CuMcmcControl\mex');
addpath('C:\Users\Alexander\Documents\MATLAB\CuMcmcControl\mcmcstat');

data.repetitions = 4;
data.experiments = 600;

 %  10 - natural overall mortality
 % 80 - exit
 % 53  -perc fed mosq
data.obs =[10 80 53];

theta0 = [1e-7 1e-4 0.99981 6.1*1e-4];
%theta0 = [1e-6 2*1e-4 1e-5 8*1e-4 0.9995 3.1*1e-4];

%res = hut_exp_ss_reduced(theta0,data)
sigma2 = 10;%std for experimental data
%%%% the MCMC part
model.ssfun=@hut_exp_ss_reduced;
model.sigma2=sigma2;
model.N = 1;

% the parameter structure: name, init.val, min, max
params = {
{'\sigma^g_{acc}(1)',theta0(1),0,1e-4}
{'\sigma^g_{acc}(2)',theta0(2),0,1e-2}
{'p_{net}', theta0(3), 0.9,1}
{'p_{hut}', theta0(4), 0,1}
};
% MCMC options
options.nsimu = 30000; % number of samples
v = [1e-5 1e-4 1e-4 1e-5];
options.qcov = diag(v.^2); % (initial) proposal covariance
options.method = 'am'; % method: DRAM
options.adaptint = 100; % adaptation interval

global modeloutput modeloutpc
% calling MCMCRUN
[resultsReduced1,chainReduced1,s2chainReduced1,sschainReduced1, hchainReduced1] = mcmcrun(model,data,params,options);

save('resultsReduced1','resultsReduced1');
save('chainReduced1','chainReduced1');
save('s2chainReduced1','s2chainReduced1');
save('sschainReduced1','sschainReduced1');
save('hchainReduced1','hchainReduced1');
modeloutputReduced1 = modeloutput;
save('modeloutputReduced1','modeloutputReduced1');


figure;
mcmcplot(chainReduced1,[],resultsReduced1.names);
figure;
mcmcplot(chainReduced1,[],resultsReduced1.names,'pairs');


figure;

N=length(modeloutput(:,3));

subplot(2,2,1);
hold on
plot(1:N, modeloutput(:,1), 'r-') ; 
plot(1:N, data.obs(1), 'k-') ;
hold off
title('Natural mortality');
axis tight

subplot(2,2,2); 
hold on
plot(1:N, modeloutput(:,2), 'r-') ;
plot(1:N, data.obs(2), 'k-') ;
hold off
title('Exit rate in control case');
axis tight

subplot(2,2,3);
hold on
plot(1:N, modeloutput(:,3), 'r-') ;
plot(1:N, data.obs(3), 'k-') ;
hold off
title('Fed rate');
axis tight
