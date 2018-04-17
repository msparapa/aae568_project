function Actual = prop_actual(Actual,Chaser,Cov,alpha,alpha_t,t0,tf);


X0 = [Actual.r0; Actual.theta0; Actual.rdot0; Actual.thetadot0];

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,X] = ode113(@(t,X)prop_actual_odes(t,X,Chaser,Cov,alpha,alpha_t),[t0 tf],X0,options);

Actual.X = X(end,:)';
return