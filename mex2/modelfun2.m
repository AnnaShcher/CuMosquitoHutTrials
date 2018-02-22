function modoutp = modelfun2(theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, repetitions, nexperiments)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%initializing experimental conditions

global modeloutpc

exp_params.mu = 0.1;%natural death rate
%exp_params.sig_acc = [1e-10 1e-2/80];
%exp_params.pnet = 0.9998;%IconMAxx
exp_params.pnet = theta7;%Alphacypermethrin
%exp_params.pnet = 0.99975;
%exp_params.pnet = 0.0;
%exp_params.phut = 6.6177*1e-4;%Iconmaxx
exp_params.phut = theta8;%Alphacypermethrin
%exp_params.phut = 4.8*1e-4;
%exp_params.phut = 1;
exp_params.xlim = [-1.5 1.5];
exp_params.ylim = [-1.5 1.5];
exp_params.eps = 0.65;
exp_params.Cmax = Inf;
exp_params.d50 = theta3;
exp_params.r = theta4;
exp_params.s = 0.01;
exp_params.alpha_p = theta5;
exp_params.tmax = theta6;%fixed time, from the control case parameters
exp_params.alpha_d = 0;%not used

%An. Gambiae
exp_params.sig_acc = [1e-4 1.2500e-5];
%exp_params.tmax = theta4;%fixed time, from the control case parameters
[ in, dead, trap, fed, unfd_gamb, Ctot] = ...
CuMosquitoHutTrialsITNmex(repetitions, nexperiments, exp_params);
in_Gamb = sum(in);
pdead_Gamb =(mean(sum(dead)./in_Gamb) - 0.1)/(1 - 0.1);
ptrap_Gamb = mean(sum(trap)./in_Gamb);
pfed_Gamb = mean(sum(fed)./in_Gamb);
punfd_gamb = mean(sum(unfd_gamb)./in_Gamb);
Ctot_Gamb  = mean(mean(Ctot));

%An. Arabiensis
exp_params.sig_acc = [theta1 theta2];
%exp_params.tmax = theta5;
[ in, dead, trap, fed, unfd_arab, Ctot ] = ...
CuMosquitoHutTrialsITNmex(repetitions, nexperiments, exp_params);
in_Arab = sum(in);
pdead_Arab = (mean(sum(dead)./in_Arab) - 0.1)/(1 - 0.1);
ptrap_Arab  = mean(sum(trap)./in_Arab);
pfed_Arab = mean(sum(fed)./in_Arab);
punfd_arab = mean(sum(unfd_arab)./in_Arab);
Ctot_Arab  = mean(mean(Ctot));

modoutp = 100*[pdead_Gamb pdead_Arab ptrap_Gamb ptrap_Arab...
                    pfed_Gamb pfed_Arab punfd_gamb punfd_arab Ctot_Gamb/100 Ctot_Arab/100]; 
modeloutpc = modoutp;
end

