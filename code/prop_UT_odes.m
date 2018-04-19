function dX_dt = prop_UT_odes(t, X, Chaser, Cov, alpha, alpha_t)
% PROP_UT_ODES Propagate dynamics for the UT
%
%   dX_dt = prop_UT_odes(t, X, Chaser, Cov, alpha, alpha_t)
%
%   Inputs:
%   
%       - t: time relative to mission start
%       - X: State vector
%       - Chaser: Structure containing information about the Chaser s/c
%       - Cov: Structure containing information about the Covariance
%           (Unused)
%       - alpha: Optimal control
%       - alpha_t: time associated with alpha; relative to mission start
%
%   Outputs:
%
%       - dX_dt: derivative of X w.r.t. t

r = X(1);
theta = X(2);
rdot = X(3);
thetadot = X(4);

m = 1 - Chaser.mdot*t;
T = Chaser.T;

alpha_i = interp1(alpha_t, alpha, t);

dX_dt = zeros(size(X));
dX_dt(1) = rdot;
dX_dt(2) = thetadot;
dX_dt(3) = r*thetadot - 1/r^2 + T/m*cos(alpha_i - theta);
dX_dt(4) = -2*rdot*thetadot/r + T/m/r*sin(alpha_i - theta);

end
