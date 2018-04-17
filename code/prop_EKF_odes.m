function dX_dt = prop_EKF_odes(t,X,Chaser,Cov,alpha,alpha_t)

r = X(1);
theta = X(2);
rdot = X(3);
thetadot = X(4);

for i = 1:4
    for j = 1:4
        P(i,j) = X(4*i+j);
    end
end

m = Chaser.m0-Chaser.mdot*t;
T = Chaser.T;

alpha_i = interp1(alpha_t,alpha,t);

cos_gamma = cos(alpha_i-theta);
sin_gamma = sin(alpha_i-theta);

dX_dt(1) = rdot;
dX_dt(2) = thetadot;
dX_dt(3) = r*thetadot - 1/r^2 + T/m*cos_gamma;
dX_dt(4) = -2*rdot*thetadot/r + T/m/r*sin_gamma;

A = [0, 0, 1, 0;...
    0, 0, 0, 1;...
    thetadot+2/r^3, T/m*sin_gamma, 0, r;...
    2*rdot*thetadot/r^2 - T/m/r^2*sin_gamma, -T/m/r*cos_gamma, -2*thetadot/r, -2*rdot/r];
C = [0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 1, 0;...
    0, 0, 0, 1];

Pdot = A*P + P*A' + C*Cov.R*C';

for i = 1:4
    for j = 1:4
        dX_dt(4*i+j) = Pdot(i,j);
    end
end

dX_dt = dX_dt';



