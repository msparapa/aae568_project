function Nav = prop_EKF(Nav,Chaser,Cov,alpha,alpha_t,t0,tf)
% PROP_EKF Propagate the state and covariance in the EKF
%
%   Nav = prop_EKF(Nav, Chaser, Cov, alpha, alpha_t, t0, tf)
%
%   Inputs:
%
%       - Nav: Structure containing information about the estimated Chaser 
%       s/c state
%       - Chaser: Structure containing informationa bout the Chaser s/c
%       - Cov:  Structure containing informationa bout the Covariance
%       - alpha: Optimal control history
%       - alpha_t: time associated with alpha, relative to mission start
%       - t0: initial time on segment, relative to mission start
%       - tf: final time on segment, relative to mission start
%
%   Outputs:
%
%       - Nav: Updated Nav structure

X0 = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
for i = 1:4
    for j = 1:4
        X0(4*i+j,1) = Nav.P0(i,j);
    end
end

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
odes = @(t,X)prop_EKF_odes(t,X,Chaser,Cov,alpha,alpha_t);
[t,X] = ode113(odes, [t0 tf], X0, options);

Nav.X_history{end+1} = X;
Nav.t_history{end+1} = t;

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