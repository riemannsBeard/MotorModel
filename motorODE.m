function dydt = motorODE(t, y, alpha, beta, nu, rho, s0, Vt, tau)

dydt = zeros(2,1);

dydt(1) = -(beta*s0*interp1(tau, Vt, t) - beta*nu*s0*y(1) - nu*rho*y(2))/(s0*(1 - alpha*beta));
dydt(2) = (-s0*interp1(tau, Vt, t) + nu*s0*y(1) + alpha*rho*nu*y(2))/(s0*(1 - alpha*beta));

end
