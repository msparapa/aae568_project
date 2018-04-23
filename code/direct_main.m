% Collin York
% Initialize Indirect Optimization
clear all; clc; close all;
addpath('optim');

% Options for simulation

% Which optimization method to use
%
%   - indirect
%   - direct (TODO)
sim_opt.optim = 'direct';

% Which estimation method to use
%
%   - ekf   Extended Kalman Filter
%   - ut    Unscented Transform + Unscented Kalman Filter
sim_opt.estim = 'ekf';

% Tolerance to check if state has reached final value
sim_opt.stateTol = 1e-6;

%% Define Dimensional Initial Conditions
% Dimensional initial conditions
T = 0.2034;             % Thrust, kN; TODO - DEFINE THIS
Isp = 3000;             % Specific impulse, sec
g0 = 9.80665/1000;      % Mean Earth gravity, km/s^2
M0 = 500;               % Initial S/C Mass, kg


r0 = 6378 + 780;        % Initial orbital radius of chaser, km
theta0 = 0;             % Initial longigutde of chaser, rad
rdot0 = 0;              % Initial radial velocity of chaser, km/s

charL = r0;             % characteristic length, km; Initial Orbit Radius
charM = M0;             % characteristic mass, kg; S/C mass
muEarth = 3.986e5;      % Earth gravity parameter, km^3/s^2
charT = sqrt(charL^3/muEarth); % characteristic time, sec; sqrt(r0^3/mu)

% Compute thetadot for a circular orbit; rad/s
thetadot0 = sqrt(muEarth/r0^3);

range_err = (10/1000)/charL;            % 10 m -> nondim
rr_err = (1/1000/1000)*charT/charL;     % 1 mm/s -> nondim
            
% Chaser Nondimensional Parameters
% These will be fed in by main code
Chaser.T = T/charM/charL*charT^2;        % Non-dim thrust
Chaser.mdot = T/(Isp*g0)/charM*charT;    % Non-dim mdot
Chaser.m0 = M0/charM;                    % Non-dim mass; starts at 1; not state variable
Chaser.ts_opt = 300/charT;               % Non-dim time-of-flight for one "leg" between observations 

% Filtered Nav Initial State for Optimizer at t0
Nav.r = 1;
Nav.theta = 0;
Nav.rdot = 0;
Nav.thetadot = 1;
Nav.P = 0*1e-12*eye(4); % Initial Covariance to test EKF
Nav.X_history = {};
Nav.t_history = {};

% Actual Initial State for Optimizer at t0
Actual.X = [1; 0; 0; 1] + sqrt(Nav.P)*randn(4,1);    %Current state [r, theta, rdot, thetadot]
Actual.X_history = {};      % Each cell holds the state history for one segment
Actual.t_history = {};      % Each cell holds the time associated with the state history
Actual.alpha_history = {};
Actual.alpha_t_history = {};

% Target Actual State at t = 0 (note: not same as updated t0)
Target.r0 = 1.01;
Target.theta0 = 0.21;
Target.rdot0 = 0;
Target.thetadot0 = sqrt(1/Target.r0^3);

% Data from Montenbruk & Gill: 
%   - At LEO (e.g. Iridium at r = 63780 + 780), most significant accel pert
%   is J_2,0 w/ magnitude 1e-5 km/s^2
%   - At GEO (r = 42000 km), most signficant accel pert are J_2,0 and Lunar 
%   gravity w/ magnitude 1e-8 km/s^2
A = [zeros(2,4); zeros(2,2), eye(2)];
Cov.R = A*1e-5*(charT^2/charL);   % Acceleration Process Noise (xdot = f(x,u,t) + C*w)

switch(sim_opt.estim)
    case 'ekf'
        % Noise Covariance
        Cov.Z = [range_err, 0; 0, rr_err];
%         Cov.Z = eye(2)*1e-12; % Measurement noise (y = Hx + Gz)
    case 'ut'
        % Some arbitrary covariance matrix
        P0 = rand(4);
        P0 = P0*P0.' * 1e-9;        % Use small values to avoid larger errors that crash the Chaser into Earth
        Cov.P0 = P0;
        Cov.alpha = 1;
        Cov.beta = 2.0;
        Cov.dt = 0.05;    % Propagate step-size (nondimensional time)
end

t_now = 0;          % t_now is the current reoptimization time (not always 0)
tf_rel_guess = pi;  % pi is good guess when t_now = 0, will update as tf-t_now

% helpful IC, replaced by lambda_f when re-optimizing
lambda0_guess = [22; -7; 20; -5];


% Begin Mission Loop
gameover = false;
count = 0;

%% Using current state est., compute optimal traj. to reach target
%   
%   Optimization method should spit out the FULL optimal control
%   from the current state to the final rendezvous
%
%   We won't necessarily use the entire optimal solution, just a
%   section of it that corresponds to our flight time between
%   observations
switch(lower(sim_opt.optim))
    case 'direct'
        if(count > 0)
            lambda0_guess = lambda_seg;         % Update
            t_now = t_seg;                      % Update starting time
            tf_rel_guess = tf - t_now ;         % Update TOF guess
            if tf_rel_guess <= 0
                tf_rel_guess = 1.1*Chaser.ts_opt;
            end
        end
        [alpha, alpha_t, tf, sol] = direct_fcn(Chaser, Target, Nav, t_now, 0, tf_rel_guess);
end

r = sol.y(1,:);
theta = sol.y(2,:);
rdot = sol.y(3,:);
thetadot = sol.y(4,:);
tf = sol.parameters(1);

figure();
subplot(2,1,1);
plot(cos(theta).*r, sin(theta).*r);
xlabel('x');
ylabel('y');
title('Position State Space');

subplot(2,1,2);
stairs(((1:length(sol.control)) - 1)/(length(sol.control)-1)*tf, sol.control);
ylabel('$\alpha$','interpreter','latex');
xlabel('Time');
title('Control History');