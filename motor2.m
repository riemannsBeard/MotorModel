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

i1 = 0.5;
f = 50;
R1 = 1;
L1 = 0.011;
R2 = 2;
L2 = 0.013;
Rfe = 1603.46;
Lu = 0.515;
n1 = 3000;
n = 2860;
Vp = 380;

tc = Lu*i1/Vp;

tt = D.t/tc;
V = D.v1;

N = 5;
V = repmat(V, [N, 1]);
tt = linspace(tt(1), N*tt(end), length(V))';

%% Initial conditions and parameters
Lambda1 = L1/Lu;
Lambda2 = L2/Lu;
alpha = 1 + Lambda1;
beta = 1 + Lambda2;
nu = i1*R1/Vp;
rho = R2/R1;
p = 3;
ns = 60*f/p;
Tc = (2*pi*ns/60)/(Vp*i1);
Tr = 5e3/Tc;
ji = 3e6;

xi1 = -1;
xi2 = 1;
s0 = 1;

%% Solution

% Numerical solution using an RK4/5 method
F = ode23s(@(t,y) motorODE(t, y, alpha, beta, nu, rho, Tr, ji, V, tt),...
    tt, [xi1 xi2 s0]);

tau = F.x';
xi1 = F.y(1,:)';
xi2 = F.y(2,:)';
s = F.y(3,:)';

% Torque
T = 3*rho./s.*xi2.^2;

% Rotor velocity
nr = ns*(1 - s);

%% Plots
if exist('figs', 'dir') == 0
    mkdir('figs');
end

figure,
plot(tt*tc, V*Vp)
xlabel('$t$ (s)')
ylabel('$V_1$ (V)')
saveas(gcf, './figs/V1_vs_t', 'jpg')

figure,
plot(tau*tc, [xi1 xi2]*i1)
legend('$i_1$', '$i_2$')
xlabel('$t$ (s)')
ylabel('$I$ (A)')
saveas(gcf, './figs/I_vs_t', 'jpg')

figure,
plot(tau*tc, T*Tc, '-', tau*tc, T*0 + Tr*Tc, '--')
legend('$T_i$', '$T_r$')
xlabel('$t$ (s)')
ylabel('$T$ (Nm)')
saveas(gcf, './figs/T_vs_t', 'jpg')

figure,
subplot(211), plot(tau*tc, s)
% xlabel('$t$ (s)')
ylabel('$s$')
% saveas(gcf, './figs/s_vs_t', 'jpg')

subplot(212), plot(tau*tc, nr)
xlabel('$t$ (s)')
ylabel('$n_r$ (rpm)')
saveas(gcf, './figs/s_nr_vs_t', 'jpg')
