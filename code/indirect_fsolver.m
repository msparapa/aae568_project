function F = indirect_fsolver(X,Chaser,Target,Nav,t0,x0)

X0 = [x0; X(1:4)];
tf_rel = X(5);

options = odeset('RelTol',1e-9,'AbsTol',1e-9);
odes = @(tau, X, tf_rel) indirect_odes(tau, X, tf_rel, Chaser);
[tau,Y] = ode45(@(tau, X) indirect_odes(tau, X, tf_rel, Chaser),[0 1],X0,options);

Yf = Y(end,:)';

PSI = indirect_bcs(X0, Yf, tf_rel, Chaser, Target, Nav, t0);

F = PSI(5:9);
return
