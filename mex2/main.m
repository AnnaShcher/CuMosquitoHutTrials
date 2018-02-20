clear all;close all;clc;


repetitions = 4;
nexperiments = 600;

%theta0 = [1e-6 2*1e-4 1e-5 8*1e-4 0.9995 3.1*1e-4];
%theta0 = [1e-6 2*1e-4 1e-5 8*1e-4 0 1];
%theta0 = [4.675824121377307e-05   2.004469446332999e-04 4.195019954476745e-05 9.180249122576662e-04 0.999362425306867 3.301410328230170e-04];
%theta0 = [4.675824121377307e-05   2.004469446332999e-04 4.195019954476745e-05 9.180249122576662e-04 0.9995 4.6e-04];
%theta0 = [ 0.0001    0.0002    0.0006    0.0011    0.9991    0.0003];
%theta0 = [1e-7 1e-4 0.99981 6.1*1e-4];
%theta0 = [1e-10 1e-2/80 1e-10 1e-2/80 0.9998 6.5*1e-4];
theta0= linspace(1e-4,1e-4,200) ;
theta1=linspace(1.2500e-5,1.2500e-5,200);
theta2=linspace(1e-4,1e-4,200);
theta3=linspace(1.2500e-5,1.2500e-5,200);


theta4=0.998:0.00001:0.9999999;
%theta4=linspace(0.99899,0.99899,200);
%theta5=linspace( 0.004701, 0.004701,200);
theta5= 0.000001:0.00005:0.01;
theta0 = [1e-4 1.25e-5 1e-4 1.25e-5 0.99978 6.1e-4];
modoutp = modelfun(theta0(1), theta0(2), theta0(3), theta0(4), ...
    theta0(5), theta0(6), repetitions, nexperiments)
%theta0 = [1e-4 1.2500e-5 1e-4 1.2500e-5 0.999899 1.0*1e-3];%initial guess
%for ref 5 in Kitau.
%experimental parameters: [sigma_acc_1 sigma_acc2 sigma_acc_1 sigma_acc2
%p_net p_hut] - initial guess of the mcmc
% for x=1:200
% 
% modoutp(x).model = modelfun(theta0(x), theta1(x), theta2(x), theta3(x), theta4(x), theta5(x), repetitions, nexperiments);
% 
% 
% %data = [10 10 80 79 53 22]%data from Kitau paper
% end