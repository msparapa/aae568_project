function dX_dtau = indirect_odes(tau, X, tf_rel, Chaser)
% INDIRECT_ODES Compute the state derivatives for the indirect optimization
% method
%
%   dX_dtau = indirect_odes(tau, X, tf, Chaser) 
%
%   Inputs:
%
%       - tau: nondimensional propagation time, relative to the beginning
%       of the segment: tau = [0, ..., 1]
%       - X: state vector, 
%           X = [r, theta, rdot, thetadot, lambda1, ..., lambda4]
%       - tf_rel: time of flight for the segment
%       - Chaser: structure that stores information about the Chaser s/c
%
%   Outputs:
%
%       - dX_dtau: derivative of the state vector w.r.t. tau
%
% Collin York
% AAE 508
% Final Project
% Minimum Time for Circular Orbit Phase Change ODEs

r = X(1);
theta = X(2);
rdot = X(3);
thetadot = X(4);
lambda1 = X(5);
lambda2 = X(6);
lambda3 = X(7);
lambda4 = X(8);
m = Chaser.m0 - Chaser.mdot*tau*tf_rel;

cos_gamma = -lambda3*r/sqrt(lambda3^2 * r^2 + lambda4^2);
sin_gamma = -lambda4/sqrt(lambda3^2 * r^2 + lambda4^2);

dX_dt = zeros(size(X));

dX_dt(1) = X(3);
dX_dt(2) = X(4);
dX_dt(3) = r*thetadot^2 - 1/r^2 + Chaser.T/m*cos_gamma;
dX_dt(4) = -2*rdot*thetadot/r + Chaser.T/(m*r)*sin_gamma;
dX_dt(5) = -lambda3*(thetadot^2 + 2/r^3) - lambda4*(2*rdot*thetadot/r^2 - Chaser.T/(m*r^2)*sin_gamma);
dX_dt(6) = -2*Chaser.T/m * lambda3*sin_gamma;
% dX_dt(6) = -lambda3*Chaser.T/m*sin_gamma - lambda4*Chaser.T/m/r*(-cos_gamma);
dX_dt(7) = -lambda1 + lambda4*2*thetadot/r;
dX_dt(8) = -lambda2 - 2*lambda3*r*thetadot + lambda4*2*rdot/r;
% dX_dt(8) = -lambda2 - lambda3*r - lambda4*(-2*rdot/r);

dX_dtau = tf_rel*dX_dt;
return