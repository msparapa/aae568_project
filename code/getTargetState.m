function [X] = getTargetState(Target, t)
% GETTARGETSTATE Compute the state of the target vehicle given its initial
% condition and the current time
%
%   X = getTargetState(Target, t) takes the Target structure, which
%   includes the initial state, and the current time, t, and computes the
%   output state X = [r, theta, rdot, thetadot] of the Target vehicle
%
%   Author: Andrew Cox
%   Version: April 18, 2018

r = Target.r0;      % circular orbit, r never changes

% circular orbit: theta changes linearly with time
theta = Target.theta0 + Target.thetadot0*t;

rdot = 0;
thetadot = Target.thetadot0;

X = [r; theta; rdot; thetadot];

end