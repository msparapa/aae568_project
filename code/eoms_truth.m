function xdot = eoms_truth(t, x, p)
% EOMS compute the state derivatives for the spacecraft
%
%   xdot = eoms(t, x, p) takes the time (t), state (x), and a
%   set of parameters (p) as inputs and outputs the time
%   derivatives of the state variables
%
%   x = [r; theta; rdot; thetadot]
%
%   p must include the fields
%
%       mu      -   Mass parameter for earth, km^3/s^2
%       T       -   Thrust magnitude, kg*km/s^2
%       u_disc  -   Discretized control history (1D, contains the angle,
%                   alpha)
%       t_disc  -   Discretized time values associated with u_disc
%       sigX    -   2x1 vector of noise standard deviations to scale a 
%                   gaussian noise addition to the accelerations
%
%   Author: Andrew Cox
%   Version: April 16, 2018

% Set size, save easy derivatives
xdot = zeros(size(x));
xdot(1:2) = x(3:4);

mu = p.mu;
T = p.T;

% Determine mass based on current time
m = p.m0 - T/(p.Isp * p.g0) * t;

% Determine state based on the discretized state and the current time
u = interpState(p.u_disc, p.t_disc, t);

% r_ddot
xdot(3) = x(1)*x(4)^2 + mu/(x(1)^2) + (T/m)*(cos(u)*cos(x(2)) + sin(u)*sin(x(2)));

% theta_ddot
xdot(4) = -2*x(3)*x(4)/x(1) + (T/(x(1)*m))*(sin(u)*cos(x(2)) - cos(u)*sin(x(2)));

% Add random noise to the accelerations
xdot(3:4) = xdot(3:4) + p.sigX.*randn(4,1);
end