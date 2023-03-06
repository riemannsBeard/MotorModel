function dydt = motorODE(t, y, alpha, beta, nu, rho, Tr, ji, Vt, tau)

dydt = zeros(3,1);

dydt(1) = -(beta*y(3)*interp1(tau, Vt, t) - beta*nu*y(3)*y(1) - nu*rho*y(2))/(y(3)*(1 - alpha*beta));
dydt(2) = (-y(3)*interp1(tau, Vt, t) + nu*y(3)*y(1) + alpha*rho*nu*y(2))/(y(3)*(1 - alpha*beta));
dydt(3) = (Tr - 3*rho*y(2)^2/y(3))/ji;

end
