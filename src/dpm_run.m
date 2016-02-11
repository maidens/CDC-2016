clear all 
close all
clc

r = 0.5;      % growth rate
K = 1000;     % carrying capacity

% create state grid
grd.Nx{1}    = 50;
grd.Xn{1}.lo = 0;
grd.Xn{1}.hi = 1;
grd.Nx{2}    = 50;
grd.Xn{2}.lo = -1;
grd.Xn{2}.hi = 1;

% create input grid 
grd.Nu{1}    = 1000;
grd.Un{1}.lo = 0;
grd.Un{1}.hi = 1;

% set initial state
grd.X0{1}    = 1;
grd.X0{2}    = 0; 

% set final state constraints
grd.XN{1}.hi = 1;
grd.XN{1}.lo = 0;
grd.XN{2}.hi = 1;
grd.XN{2}.lo = -1;


% define problem
prb.Ts = 1;
prb.N  = 25;

% set options
options = dpm(); 
% options.UseLine  = 1; 
options.Waitbar = 'on'; 
options.SaveMap  = 1; 
options.MyInf  = 10000; 
options.FixedGrid = 0;
options.BoundaryMethod = 'none'; 

par = [r, K]; 

% load('results.mat')
tic 
[res, dyn] = dpm(@logistic_nondimensionalized, par, grd, prb, options);
time = toc


%% 
u_opt = res.U{1}; 
figure 
plot(u_opt, '-o') 
xlabel('t')
ylabel('u_t')

x1_opt = res.X{1}; 
x2_opt = res.X{2}; 
figure
hold on 
plot(K*x1_opt, '-or')
xlabel('t')
ylabel('x_t')
% plot(x2_opt)
hold off