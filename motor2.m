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

i1 = 10*0.5;
f = 50;
R1 = 2;
L1 = 0.011;
R2 = 2;
L2 = 0.013;
Rfe = 1603.46;
Lu = 0.515;
n1 = 3000;
n = 2860;
Vp = 380;
p = 1;
ns = 60*f/p;

% Internal mechanical torque (N m)
Mc = 3*R2*i1^2/(2*pi*ns/60);
% Moment of inertia (kg m^2)
J = 0.02;

% Nondimensional resistive torque (Tr/Mmi)
Mr = 0.125;
% Caracteristic (starting) time
tc = J*2*pi/60*ns/(Mc*(1 - Mr));

tt = D.t/tc;
V = D.v1;

N = 15;
V = repmat(V, [N, 1]);
tt = linspace(tt(1), N*tt(end), length(V))';

% The prescribed voltage signal can be replaced by
% V = Vp*cos(2*pi*f*tt*tc);
% It might be useful to compute the permanent state solution analytically.

%% Initial conditions and dimensionless parameters
alpha = (L1 + Lu)*i1/Vp/tc;
beta = Lu*i1/Vp/tc;
gamma = (L2 + Lu)*i1/Vp/tc;

nu = i1*R1/Vp;
rho = R2/R1;

xi1 = -1;
xi2 = 1;
s0 = 1;

%% Solution

% Numerical solution using an RK4/5 method
% F = ode23s(@(t,y) motorODE(t, y, alpha, beta, nu, rho, Mr, ji, V, tt),...
%     tt, [xi1 xi2 s0]);

F = ode23s(@(t,y) motorODE2(t, y, alpha, beta, gamma, nu, rho, Mr, V, tt),...
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
plot(tt*tc, V*Vp, tt*tc, Vp*sin(2*pi*f*tt*tc))
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
plot(tau*tc, T*Mc, '-', tau*tc, T*0 + Mr*Mc, '--')
legend('$T_i$', '$T_r$')
xlabel('$t$ (s)')
ylabel('$T$ (Nm)')
saveas(gcf, './figs/T_vs_t', 'jpg')

figure,
subplot(211), plot(tau*tc, s, tau*tc, s*0, '--')
ylabel('$s$')
subplot(212), plot(tau*tc, nr, tau*tc, nr*0 + ns, '--')
xlabel('$t$ (s)')
ylabel('$n_r$ (rpm)')
saveas(gcf, './figs/s_nr_vs_t', 'jpg')
