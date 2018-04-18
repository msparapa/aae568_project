% Collin York
% Initialize Indirect Optimization
clear all; clc; close all;

% Options to Plot Indirect Optimization
plot_opt.indirect = 1;      % 1=plot, 0=no plot
plot_opt.i = 0;             % last plot figure number

% Options for simulation
sim_opt.optim = 'indirect';     % which optimization method to use
sim_opt.estim = 'ekf';          % which estimate method to use
sim_opt.stateTol = 1e-6;        % tolerance to check if state has reached final value

%% Define Dimensional Initial Conditions
charL = 6378 + 200;     % characteristic length, km; TODO - DEFINE THIS
charT = 1;              % characteristic time, sec; TODO - DEFINE THIS
charM = 500;            % characteristic mass, kg; TODO - DEFINE THIS
muEarth = 3.986e5;      % Earth gravity parameter, km^3/s^2

% Dimensional initial conditions
T = 1;                  % Thrust, Newtons; TODO - DEFINE THIS
Isp = 1500;             % Specific impulse, sec
g0 = 9.80665/1000;      % Mean Earth gravity, km/s^2

r0 = 7000;              % Initial orbital radius of chaser, km
theta0 = 0;             % Initial longigutde of chaser, rad
rdot0 = 0;              % Initial radial velocity of chaser, km/s

% Compute thetadot for a circular orbit; rad/s
thetadot0 = sqrt(muEarth/r0^3);

% Chaser Nondimensional Parameters
% These will be fed in by main code
Chaser.T = 0.05; % Non-dim thrust (Value tends to work well)
Chaser.mdot = 0.001; % Non-dim mdot (Value tends to work well)
Chaser.m0 = 1; % Non-dim mass; starts at 1; not state variable
% **** m0 would need to be updated
Chaser.ts_opt = 1; % Re-optimization non-dim sample time 

% Actual Initial State for Optimizer at t0
Actual.r0 = 1;
Actual.theta0 = 0;
Actual.rdot0 = 0;
Actual.thetadot0 = 1;

% Filtered Nav Initial State for Optimizer at t0
Nav.r0 = 1;
Nav.theta0 = 0;
Nav.rdot0 = 0;
Nav.thetadot0 = 1;
Nav.P0 = zeros(4); % Initial Covariance to test EKF

% Target Actual State at t = 0 (note: not same as updated t0)
Target.r0 = 1.01;
Target.theta0 = 0.21;
Target.rdot0 = 0;
Target.thetadot0 = sqrt(1/Target.r0^3);

% Noise Covariance
Cov.R = eye(4)*1e-4; % Acceleration Process Noise (xdot = f(x,u,t) + C*w)
Cov.Z = eye(4)*1e-4; % Measurement noise (y = x + z)

t0 = 0; % t0 is the current reoptimization time (not always 0)
tf_rel_guess = pi; % pi is good guess when t0 = 0, will update as tf-t0

% helpful IC, replaced by lambda_f when re-optimizing
lambda10 = 22; 
lambda20 = -7;
lambda30 = 20;
lambda40 = -5;
lambda0_guess = [lambda10; lambda20; lambda30; lambda40];

% Begin Mission Loop
gameover = false;
while(~gameover)
    %% Using current state est., compute optimal traj. to reach target
    switch(lower(sim_opt.optim))
        case 'indirect'
            % Indirect Optimization - Collin
            % Outputs: alpha(t), sample t's for alpha, 
            %   absolute tf (tf_rel + t0), lambdas at next re-optimization
            %   time
            [alpha, alpha_t, tf, lambda_f] = indirect_fcn(Chaser, Target,...
                Nav, t0, plot_opt, lambda0_guess, tf_rel_guess);
    end
    
    
    %% Propagate for one "step" until next observation
    
    % Propagate the true state with interpolated process noise
    Actual = prop_actual(Actual, Chaser, Cov, alpha, alpha_t, t0, tf);
    
    switch(lower(sim_opt.estim))
        case 'ekf'
            % EKF function to propagate State covariance in P in continuous
            % time with acceleration process covariance
            Nav = prop_EKF(Nav, Chaser, Cov, alpha, alpha_t, t0, tf);
    end
    
    %% Make Observation
    %   
    %   * Update state estimate
    %   * Update covariance
    
    switch(lower(sim_opt.estim))
        case 'ekf'
            H = eye(4);         % TODO - what is this?
            G = eye(4);         % TODO - what is this?
            L_k = (Nav.P * H.') / (H*Nav.P*H.' + G*Cov.Z*G.');
            y = H*Actual.X + Cov.Z*randn(4,1);  % measurement + noise
            % apriori nav estimate
            nav_pre = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
            % nav estimate post observation
            nav_post = nav_pre + L_k*(y - H*nav_pre);   % update
            
            % Update Nav state structure
            Nav.r = nav_post(1);
            Nav.theta = nav_post(2);
            Nav.rdot = nav_post(3);
            Nav.thetadot = nav_post(4);
            Nav.P = (eye(4) - L_k*H)*Nav.P;     % covariance update
    end
    
    %% Evaluate status
    %   Has spacecraft reached the target?
    
    X_targ = getTargetState(Target, t);
    X_chaser = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
    state_err = norm(X_targ - X_chaser);
    gameover = state_err < sim_opt.stateTol;
    
    if(~gameover)
        % Update Chaser state
        Chaser.m = Chaser.m - Chaser.mdot*t_ellapsed;
    end
    
end