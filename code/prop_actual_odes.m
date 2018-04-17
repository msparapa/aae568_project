function dX_dt = prop_actual_odes(t,X,Chaser,Cov,alpha,alpha_t,w)

r = X(1);
theta = X(2);
rdot = X(3);
thetadot = X(4);

m = Chaser.m0-Chaser.mdot*t;
T = Chaser.T;

alpha_i = interp1(alpha_t,alpha,t);
for i = 1:4
    w_i(i) = interp1(alpha_t,w(i,:),t);
end

cos_gamma = cos(alpha_i-theta);
sin_gamma = sin(alpha_i-theta);

dX_dt(1) = rdot;
dX_dt(2) = thetadot;
dX_dt(3) = r*thetadot - 1/r^2 + T/m*cos_gamma;
dX_dt(4) = -2*rdot*thetadot/r + T/m/r*sin_gamma;

C = [0, 0, 0, 0;...
                  0, 0, 0, 0;...
                  0, 0, 1, 0;...
                  0, 0, 0, 1];
dX_dt = dX_dt' + C*w_i'*sqrt(Cov.R(1,1));




