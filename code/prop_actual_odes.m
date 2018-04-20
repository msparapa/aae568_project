function dX_dt = prop_actual_odes(t, X, Chaser, Cov, alpha, alpha_t, w)
% PROP_ACTUAL_ODES Compute the state derivative of the spacecraft by
% leveraging the true, or actual, dynamics
%
%   dX_dt = prop_actual_odes(t, X, Chaser, Cov, alpha, alpha_t, w)
%
%   Inputs:
%
%       - t: propagation time relative to mission start (i.e., t = 0)
%       - X: State vector, X = [r, theta, rdot, thetadot]
%       - Chaser: Structure containing information about the Chaser s/c
%       - Cov: Structure containing covariance information
%       - alpha: optimal control history
%       - alpha_t: time history associated with alpha, relative to mission
%       start
%       - w: noise that affects the "true" dynamics, also corresponds to
%       the time values stored in alpha_t
%
%   Outputs:
%
%       - dX_dt: derivative of X w.r.t. t
%
%   Author: Collin York
%   Version: April 18, 2018

r = X(1);
theta = X(2);
rdot = X(3);
thetadot = X(4);

m = 1 - Chaser.mdot*t;
T = Chaser.T;

alpha_i = interp1(alpha_t, alpha, t);
w_i = zeros(4,1);
for i = 1:4
    w_i(i) = interp1(alpha_t, w(i,:), t);
end

cos_gamma = cos(alpha_i-theta);
sin_gamma = sin(alpha_i-theta);

dX_dt = zeros(size(X));
dX_dt(1) = rdot;
dX_dt(2) = thetadot;
dX_dt(3) = r*thetadot - 1/r^2 + T/m*cos_gamma;
dX_dt(4) = -2*rdot*thetadot/r + T/m/r*sin_gamma;

dX_dt = dX_dt + Cov.R*w_i;




