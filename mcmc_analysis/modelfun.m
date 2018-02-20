function modoutp = modelfun(theta1, theta2, theta3, theta4, theta5, theta6, repetitions, nexperiments)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%initializing experimental conditions

global modeloutpc

exp_params.pnet = theta5;
exp_params.phut = theta6;
exp_params.xlim = [-3 3];
exp_params.ylim = [-3 3];
exp_params.eps = 0.2;
exp_params.tmax = 4;

%An. Gambiae
exp_params.mu = 0.1;
exp_params.sig_acc = [theta1 theta2];
[ in, dead, trap, fed, ~ ] = ...
CuMosquitoRepellentMEX(repetitions, nexperiments, exp_params);
in_Gamb = sum(in);
pdead_Gamb = mean(sum(dead)./in_Gamb);
ptrap_Gamb = mean(sum(trap)./in_Gamb);
pfed_Gamb = mean(sum(fed)./in_Gamb);

%An. Arabiensis
exp_params.mu = 0.1;
exp_params.sig_acc = [theta3 theta4];
[ in, dead, trap, fed, ~ ] = ...
CuMosquitoRepellentMEX(repetitions, nexperiments, exp_params);
in_Arab = sum(in);
pdead_Arab = mean(sum(dead)./in_Arab);
ptrap_Arab  = mean(sum(trap)./in_Arab);
pfed_Arab = mean(sum(fed)./in_Arab);

modoutp = ceil(100*[pdead_Gamb pdead_Arab ptrap_Gamb ptrap_Arab pfed_Gamb pfed_Arab]); 
modeloutpc = modoutp;
end

