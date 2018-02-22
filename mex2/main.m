clear all;close all;clc;


repetitions = 6;
nexperiments = 600;
%theta0 = [6*1e-4 1e-4 0.9983 0.00031];
%theta0 = [1e-5 1e-4 9*1e-5 1.5*1e-4 0.99945 1-0.9997];
%theta0 = [0.0, 1*1e+7 , 0.2, 0.08, 1e-3, 1e-2];
% d50 s alpha_p tmax_A
%theta0 = [0.73 0.08 1.2*1e-10 5.5];
%theta0 = [0.73 0.08 6*1e-13 3.0];
%theta0 = [0.69 0.02 2.5*1e-8 3.0 1.8];%cubic power
%theta0 = [0.78 0.85 1.5*1e-9 3.0 1.2];%IconMaxx
%theta0 = [0.78 0.4 1.6*1e-9 3.2 0.7];%Olyset
%theta0 = [0.78 0.4 1.6*1e-9 3.2 0.7];%KOTAB
%theta0 = [0.756 0.87 3.2*1e-7 2.2 0.2];%IconMAxx
theta0 = [0.755 0.87 4.8*1e-8 4.5 1.2];%IconMAxx
theta0 = [0.755 0.945 1014 77.8763 5.0 2.8];%IconMAxx, logistic death rate
 theta0 = [0.755 0.03 105 4.8763 2700 10.0 400.0 3 4.0 3.4];%Olyset, logistic death rate
  theta0 = [0.1 1.9*1e-7 3*1e-4 1.3*1e-5 2.1 2.1*1e-4 4.8];%Olyset, logistic death rate
% theta0 = [0.755 1e-1 2.1*1e-9 1.2 0.3 0.89];%Alphacypermetrin LN funestus
% theta0 = [0.755 1e-3 1.5*1e-5 1.2 0.5 0.89];%Alphacypermetrin funestus
% theta0 = [0.755 1e-3 4*1e-6 3.8 0.7 0.89];%Olyset LN funestus
% theta0 = [0.755 1e-2 4*1e-8 3.8 0.7 0.89];%Olyset LN funestus
%theta0 = [0.84 0.008 4.5*1e-10 2.8 1.5];
%theta0 = [1.0 0.01 4*1e-10 3.5 2.7];
%theta0 = [1.0 0.25 3.1*1e-11 0.021];
tic

%modoutp =modelfun_logistic_death_rate(theta0, repetitions, nexperiments)
data = [72 15 93 84 16 22]
modoutp = modelfun_logistic_death_rate(theta0, repetitions, nexperiments)
% toc
% modelfun1 = @(theta)modelfun_logistic_death_rate(theta, repetitions, nexperiments)
% ssfun = @(x)sum((data - modelfun1([theta0(1) x])).^2)
% % ssfun(theta0(2:end))4
% PSoptions = optimoptions(@patternsearch,'Display','iter');
% LB = zeros(1,10);
% UB = [1 Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];
% [Xps,Fps] = patternsearch(ssfun,theta0(2:end),[],[],[],[],LB,UB,PSoptions)

% display_ITN;

% neg_Olyset = [66 11 8 15 85 73];
%  pos_Olyset = [78 18 29 30 97 91];
% neg_IconMaxx = [59 37 2 2 62 77];
% pos_IconMaxx = [89 51 34 21 89 92];
% 
% initg_Olyset = [70 19 25 17 93 96];
% initg_IconMaxx = [74 47 6 5 80 89];
% 
% mean_Olyset = [72 15 16 22 93 84];
% mean_IconMaxx = [74 44 9 7 79 87];
% 
% figure;
% subplot(2,1,1);
% hold on
% errorbar(1:6,mean_Olyset,-neg_Olyset+mean_Olyset,pos_Olyset-mean_Olyset);
% plot(1:6, initg_Olyset,'*');
% title('Olyset, initial guess versus Kitau paper data');
% subplot(2,1,2);
% hold on
% errorbar(1:6,mean_IconMaxx,-neg_IconMaxx+mean_IconMaxx,pos_IconMaxx-mean_IconMaxx);
% plot(1:6, initg_IconMaxx,'*');
% title('IconMaxx, initial guess versus Kitau paper data');
