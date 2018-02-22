function modoutp = modelfun_logistic_death_rate(theta, repetitions, nexperiments)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%initializing experimental conditions

global modeloutpc

exp_params.mu = 0.1;%natural death rate
exp_params.sig_acc = [theta(3) 1.2500e-5];

exp_params.pnet = 0.9998899;%Reference 5 from Kitau paper, funestus
% exp_params.pnet = 0.99978;%IconMAxx
% exp_params.phut = 6.1*1e-4;%Icon Maxx
exp_params.phut = 8.5*1e-4;%Reference 5 from Kitau paper, funestus
%exp_params.phut = 4.3*1e-4;%Alphacypermethrin
%exp_params.phut = 1;
exp_params.xlim = [-1.5 1.5];
exp_params.ylim = [-1.5 1.5];
exp_params.eps = 0.65;
exp_params.d50 = 0.755;
exp_params.r = theta(1);
exp_params.s = 0.01;
exp_params.alpha_p = theta(2);
exp_params.alpha_d = 1;%not used
exp_params.s_NetCont = 1;
exp_params.tmax = 5;

exp_params.d50_NetCont = theta(4);
exp_params.tmax = theta(5);
%An. Gambiae
[ in, dead, trap, fed, unfd_gamb, Ctot_Gamb, NetCont_Gamb] = ...
CuMosquitoHutTrialsITNmex(repetitions, nexperiments, exp_params);
in_Gamb = sum(in);
pdead_Gamb =(mean(sum(dead)./in_Gamb) - 0.1)/(1 - 0.1);
ptrap_Gamb = mean(sum(trap)./in_Gamb);
pfed_Gamb = mean(sum(fed)./in_Gamb);
punfd_gamb = mean(sum(unfd_gamb)./in_Gamb);
mCtot_Gamb  = squeeze(mean(mean(Ctot_Gamb)));
mNetCont_Gamb  = squeeze(mean(mean(NetCont_Gamb)));
%An. Arabiensis
exp_params.d50_NetCont = theta(6);
exp_params.tmax = theta(7);
%An. arab
[ in, dead, trap, fed, unfd_arab, Ctot_Arab, NetCont_Arab] = ...
CuMosquitoHutTrialsITNmex(repetitions, nexperiments, exp_params);
in_Arab = sum(in);
pdead_Arab = (mean(sum(dead)./in_Arab) - 0.1)/(1 - 0.1);
ptrap_Arab  = mean(sum(trap)./in_Arab);
pfed_Arab = mean(sum(fed)./in_Arab);
punfd_arab = mean(sum(unfd_arab)./in_Arab);
mCtot_Arab  = squeeze(mean(mean(Ctot_Arab)));
mNetCont_Arab  = squeeze(mean(mean(NetCont_Arab)));
% 
% modoutp = 100*[pdead_Gamb pdead_Arab ptrap_Gamb ptrap_Arab...
%                     pfed_Gamb pfed_Arab punfd_gamb punfd_arab];
                
modoutp = 100*[pdead_Gamb pdead_Arab ptrap_Gamb ptrap_Arab...
                    pfed_Gamb pfed_Arab]; 
modeloutpc = modoutp;

% data_NCgamb = reshape(NetCont_Gamb,800*4,18000);
%  data_NCarab = reshape(NetCont_Arab,800*4,18000);
%   data_Ctotarab = reshape(Ctot_Arab,800*4,18000);
%    data_Ctotgamb = reshape(Ctot_Gamb,800*4,18000);
%    figure;subplot(3,1,1);boxplot(data_NCgamb(:,1:500:end));
%    ylabel('Net Contacts');
%    title('Number of net contacts (An. gambiae)');
%    subplot(3,1,2);boxplot(data_NCarab(:,1:500:end));
%    hold on
%    ylabel('Net Contacts');
%    xlabel('time (seconds)');
%    title('Number of net contacts (An. Arabiensis)');
%    subplot(3,1,3);
%       plot(1:500:length(mNetCont_Gamb),mNetCont_Gamb(1:500:end),'b');hold on
%    plot(1:500:length(mNetCont_Gamb),mNetCont_Arab(1:500:end),'r');
%    legend('Av. num. net cont (Gamb)','Av. num. net cont (Arab)')
%       ylabel('Net Contacts');
%    xlabel('time (seconds)');
%    
%       figure;subplot(3,1,1);boxplot(data_Ctotgamb(:,1:500:end));
%    hold on
%    ylabel('C_{tot}');
%    title('Total concentration (An. Gambiae)');
%    subplot(3,1,2);boxplot(data_Ctotarab(:,1:500:end));
%    ylabel('C_{tot}');
%    xlabel('time (seconds)');
%    title('Total concentration (An. Arabiensis)');
%    subplot(3,1,3);
%       plot(1:500:length(mCtot_Gamb),mCtot_Gamb(1:500:end),'b');hold on
%    plot(1:500:length(mCtot_Gamb),mCtot_Arab(1:500:end),'r');
%    legend('Av. total conc. (Gamb)','Av. total conc. (Arab)')
%       ylabel('C_{tot}');
%    xlabel('time (seconds)');
%    
%    figure;
%    hold on
%    plot(mNetCont_Gamb,mCtot_Gamb);
%    plot(mNetCont_Arab,mCtot_Arab,'r');
%    xlabel('Average net contacts');
%    ylabel('Ctot');
%    legend('Av. total conc. (Gamb)','Av. total conc. (Arab)')
%    title('Consumed concentration conditioned on number of contacts');
end

