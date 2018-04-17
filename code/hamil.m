% Collin York
% Indiret Opt Hamiltonian Function

function H = hamil(X,tf)

global T mdot

r = X(1);
theta = X(2);
rdot = X(3);
thetadot = X(4);
lambda1 = X(5);
lambda2 = X(6);
lambda3 = X(7);
lambda4 = X(8);
m = 1-mdot*tf;
cos_gamma = -lambda3*r/sqrt(lambda3^2 * r^2 + lambda4^2);
sin_gamma = -lambda4/sqrt(lambda3^2 * r^2 + lambda4^2);

H = lambda1*rdot + lambda2*thetadot + lambda3*(r*thetadot - 1/r^2 + T/m*cos_gamma)...
    + lambda4*(-2*rdot*thetadot/r + T/m/r*sin_gamma);