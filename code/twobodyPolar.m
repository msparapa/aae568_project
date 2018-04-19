function X = twobodyPolar(X0, tspan, dt, Chaser, Cov, alpha, alpha_t)
%...Pass the initial conditions and time interval to ode113,
%...which calculates the position and velocity at discrete
%...times t, returning the solution in the column vector X.
%...ode113 uses the m-function 'prop_UT_odes' to evaluate the
%...acceleration at each integration time step.
tol = 1e-10;
options = odeset('RelTol',tol,'AbsTol',tol);
odes = @(t, X) prop_estim_odes(t, X, Chaser, Cov, alpha, alpha_t);
[~,X] = ode113(odes, tspan(1):dt:tspan(2), X0, options);