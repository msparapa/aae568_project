function dX_dt = prop_EKF_odes(t, X, Chaser, Cov, alpha, alpha_t)
% PROP_EKF_ODES Propagate dynamics and covaraince for the EKF
%
%   dX_dt = prop_EKF_odes(t, X, Chaser, Cov, alpha, alpha_t)
%
%   Inputs:
%   
%       - t: time relative to mission start
%       - X: State vector
%       - Chaser: Structure containing information about the Chaser s/c
%       - Cov: Structure containing information about the Covariance
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

alpha_i = interp1(alpha_t,alpha,t);

cos_gamma = cos(alpha_i-theta);
sin_gamma = sin(alpha_i-theta);

dX_dt = zeros(size(X));
dX_dt(1) = rdot;
dX_dt(2) = thetadot;
dX_dt(3) = r*thetadot - 1/r^2 + T/m*cos_gamma;
dX_dt(4) = -2*rdot*thetadot/r + T/m/r*sin_gamma;

if(length(X) == 4 + 4*4)
    P = reshape(X(5:end), 4, 4);
    
    A = [0, 0, 1, 0;...
        0, 0, 0, 1;...
        thetadot+2/r^3, T/m*sin_gamma, 0, r;...
        2*rdot*thetadot/r^2 - T/m/r^2*sin_gamma, -T/m/r*cos_gamma, -2*thetadot/r, -2*rdot/r];
    C = [0, 0, 0, 0;...
        0, 0, 0, 0;...
        0, 0, 1, 0;...
        0, 0, 0, 1];

    Pdot = A*P + P*A' + C*Cov.R*C';

    dX_dt(5:end) = reshape(Pdot, 16, 1);
end

end



