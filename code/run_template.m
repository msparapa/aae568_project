clear; clc; close all;

%% Define Parameters

g0    = 9.80665/1000;      % gravitational acceleration, km/s^2
Isp   = 1500;              % specific-impulse, s
T     = 1e-3/1000;         % thrust, kg km/s^2
mu    = 3.98600e+05;       % Earth gravitational parameter, km^3/s^2
Re    = 6378;              % earth radius, km
m     = 1000;              % spacecraft mass, kg
h     = 250;               % Initial orbit altitude, km
a     = Re + h;            % Initial orbit radius, km (need to fix for eccentric orbits)

targ_a = 10000;             % Target orbit radius
targ_w = sqrt(mu/targ_a^3); % Target orbit mean motion

%% Begin observation/propagation Loop
complete = false;
while(~complete)
    %% Using current state estimate, compute optimal trajectory to reach target
    
    %% Propagate for one "step" until next observation
    %
    %   * propagate state estimate
    %   * propagate state truth
    %   * propagate covariance
    
    %% Make Observation
    %   
    %   * Update state estimate
    %   * Update covariance
    
    %% Evaluate status
    %
    %   Has spacecraft reached the target?
end