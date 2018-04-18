function dX_dt = indirect_odes_tf(t,X,Chaser)
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
m = Chaser.m0-Chaser.mdot*t;
T = Chaser.T;

cos_gamma = -lambda3*r/sqrt(lambda3^2 * r^2 + lambda4^2);
sin_gamma = -lambda4/sqrt(lambda3^2 * r^2 + lambda4^2);

dX_dt(1) = X(3);
dX_dt(2) = X(4);
dX_dt(3) = r*thetadot - 1/r^2 + T/m*cos_gamma;
dX_dt(4) = -2*rdot*thetadot/r + T/m/r*sin_gamma;
dX_dt(5) = -lambda3*(thetadot + 2/r^3) - lambda4*(2*rdot*thetadot/r^2 - T/m/r^2*sin_gamma);
dX_dt(6) = -lambda3*T/m*sin_gamma - lambda4*T/m/r*(-cos_gamma);
dX_dt(7) = -lambda1 - lambda4*(-2*thetadot/r);
dX_dt(8) = -lambda2 - lambda3*r - lambda4*(-2*rdot/r);

dX_dt = dX_dt';
return