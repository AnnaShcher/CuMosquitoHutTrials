%generating parameters distribution for control case
clear all;
close all;clc;

%addpath('C:\StorageDriveData\Dropbox\MATLAB\CuMosquitoHutTrialsControl\CuMosquitoHutTrialsControl\mcmcstat');

addpath('C:\Users\china\Downloads\CuMosquitoHutTrialsControl\CuMosquitoHutTrialsControl\CuMosquitoHutTrialsControl\mcmcstat');


data.repetitions = 4;
data.experiments = 600;

 %  10 - natural overall mortality
 % 80 - exit
 % 53  -perc fed mosq
 data.obs =[10 80 53];
%data.obs =[10 70.7 68.6];

%theta0 = [0.99978 6.5*1e-4];
%theta0 = [0.9965 4*1e-2 6];
  
theta0 = [0.99955 4.3*1e-4];
theta0 = [0.9996 4.3*1e-4]; 
%theta0 = [0.999699 4.2*1e-4];
%theta0 = [0.999699 4.2*1e-4];%
%theta0 = [0.9962 2*1e-2];
%theta0 = [1e-6 2*1e-4 1e-5 8*1e-4 0.9995 3.1*1e-4];

tic
res = hut_exp_ss_reduced(theta0,data)
toc
sigma2 = 0.9;%std for experimental data
%%%% the MCMC part
model.ssfun=@hut_exp_ss_reduced;
model.sigma2=sigma2;
model.N = 1;

% the parameter structure: name, init.val, min, max
params = {
{'p_{net}', theta0(1), 0.9,1}
{'p_{hut}', theta0(2), 0,1}
};
% MCMC options
options.nsimu = 20000; % number of samples
v = [1e-5 1e-5];
options.qcov = diag(v.^2); % (initial) proposal covariance
options.method = 'am'; % method: DRAM
options.adaptint = 100; % adaptation interval

global modeloutput modeloutpc
% calling MCMCRUN
[results_2param_resamp,chain_2param_resamp,s2chain_2param_resamp,sschain_2param_resamp, hchain_2param_resamp] = mcmcrun(model,data,params,options);


save('results_2param_resamp','results_2param_resamp');
save('chain_2param_resamp','chain_2param_resamp');
save('s2chain_2param_resamp','s2chain_2param_resamp');
save('sschain_3param_resamp','sschain_2param_resamp');
save('hchain_2param_resamp','hchain_2param_resamp');
modeloutput_2param_resamp = modeloutput;
save('modeloutput_2param_resamp','modeloutput_2param_resamp');


figure;
mcmcplot(chain_2param_resamp,[],results_2param_resamp.names);
figure;
mcmcplot(chain_2param_resamp,[],results_2param_resamp.names,'pairs');


figure;

N=length(modeloutput(:,3));

subplot(2,2,1);
hold on
plot(1:N, modeloutput(:,1), 'r-') ; 
plot(1:N, repmat(data.obs(1),1,N), 'k-') ;
hold off
title('Natural mortality');
axis tight

subplot(2,2,2); 
hold on
plot(1:N, modeloutput(:,2), 'r-') ;
plot(1:N, repmat(data.obs(2),1,N), 'k-') ;
hold off
title('Exit rate in control case');
axis tight 

subplot(2,2,3);
hold on
plot(1:N, modeloutput(:,3), 'r-') ;
plot(1:N, repmat(data.obs(3),1,N), 'k-') ;
hold off
title('Fed rate');
axis tight
