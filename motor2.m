clear all
close all
clc

%% Defaults
FS = 18;
LW = 1.5;
set(0, 'defaultAxesFontSize', FS)
set(0, 'defaultLineLineWidth', LW)
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')

%% Data
D = load('datos.mat');

f = 50;
i1 = 0.5;
i2 = 0;
R1 = 1;
L1 = 0.011;
R2 = 2;
L2 = 0.013;
Rfe = 1603.46;
Lu = 0.515;
n1 = 3000;
n = 2860;
s0 = 0.98;
Vp = 380;

tc = Lu*i1/Vp;

tau = D.t(1:end-1)/tc;
V = D.v1;

figure, plot(tau, V)

%% Initial conditions and parameters
Lambda1 = L1/Lu;
Lambda2 = L2/Lu;
alpha = 1 + Lambda1;
beta = 1 + Lambda2;
nu = i1*R1/Vp;
rho = R2/R1;

xi1 = -1;
xi2 = 1;

F = ode45(@(t,y) motorODE(t, y, alpha, beta, nu, rho, s0, V, tau),...
    tau, [xi1 xi2]);

figure, plot(F.x*tc, F.y*i1)



