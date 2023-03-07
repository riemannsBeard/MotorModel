function dydt = motorODE2(t, y, alpha, beta, gamma, nu, rho, Mr, Vt, tau)

dydt = zeros(3,1);

dydt(1) = -(gamma*y(3)*interp1(tau, Vt, t) - gamma*nu*y(3)*y(1) - beta*nu*rho*y(2))/...
    (y(3)*(beta^2 - alpha*gamma));
dydt(2) = (-beta*y(3)*interp1(tau, Vt, t) + beta*nu*y(3)*y(1) + alpha*nu*rho*y(2))/...
    (y(3)*(beta^2 - alpha*gamma));
dydt(3) = -y(2)^2/y(3) + Mr;

end
