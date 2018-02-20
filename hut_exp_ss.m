function res = hut_exp_ss(theta,data)
%ss - least square sum
% data - vector [8,2]
%we compare:
%1) overall mortality (corrected for control)
%2) blood feeding in control arm
%3) exophily in control arm

global modeloutpc

modoutp= modelfun(theta(1),theta(2),theta(3),theta(4),theta(5),theta(6),data.repetitions,data.experiments);
% fed_diff_obs = data.obs(end-1) - data.obs(end);
% fed_diff_mod = modoutp(end-1) - modoutp(end);
% ex_diff_obs = data.ydata(end-2) - data.ydata(end-3);
% ex_diff_mod = modoutp(end-2) - modoutp(end-3);
res=sum((data.obs(3:end)-modoutp(3:end)).^2);
end