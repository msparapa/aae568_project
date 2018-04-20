function F = indirect_fsolver(X, Chaser, Target, Nav, t0, fixed_x0, Opt_sol)

X0 = [fixed_x0; X(1:4)];
tf_rel = X(5);
slack = X(6);

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
odes = @(tau, X) indirect_odes(tau, X, tf_rel, Chaser);
[tau,Y] = ode113(odes, [0 1], X0, options);

Yf = Y(end,:)';

PSI = indirect_bcs(X0, Yf, tf_rel, Chaser, Target, Nav, t0, slack);

if Opt_sol
    F = PSI(5:10);
else
    F = PSI([5:8,10]);
end
return
