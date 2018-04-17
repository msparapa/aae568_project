function Nav = prop_EKF(Nav,Chaser,Cov,alpha,alpha_t,t0,tf);


X0 = [Nav.r0; Nav.theta0; Nav.rdot0; Nav.thetadot0];
for i = 1:4
    for j = 1:4
        X0(4*i+j,1) = Nav.P0(i,j);
    end
end

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,X] = ode45(@(t,X)prop_EKF_odes(t,X,Chaser,Cov,alpha,alpha_t),[t0 tf],X0,options);

Nav.r = X(end,1);
Nav.theta = X(end,2);
Nav.rdot = X(end,3);
Nav.thetadot = X(end,4);
for i = 1:4
    for j = 1:4
        Nav.P(i,j) = X(end,4*i+j);
    end
end

return