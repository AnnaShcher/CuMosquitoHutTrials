function modoutp = modelfun_reduced(theta1, theta2, repetitions, nexperiments)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%initializing experimental conditions

global modeloutpc

exp_params.pnet = theta1;
exp_params.phut = theta2;
exp_params.xlim = [-3 3];
exp_params.ylim = [-3 3];
exp_params.eps = 0.65;
exp_params.tmax = 5;

%both An. Gambiae and An. Arabiensis
exp_params.mu = 0.1;
exp_params.sig_acc = [1e-3 1.2500e-5];
[ in, dead, trap, fed, ~ ] = ...
MATLAB_interface(repetitions, nexperiments, exp_params);
in_Gamb = sum(in);
pdead = mean(sum(dead)./in_Gamb);
ptrap = mean(sum(trap)./in_Gamb);
pfed = mean(sum(fed)./in_Gamb);



modoutp = 100*[pdead ptrap pfed]; 
modeloutpc = modoutp;
end

