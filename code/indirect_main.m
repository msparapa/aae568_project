% Collin York
% Initialize Indirect Optimization

clear all; clc; close all;

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

% Options to Plot Indirect Optimization
plot_opt.indirect = 1; % 1=plot, 0=no plot
plot_opt.i = 0; % last plot figure number

t0 = 0; % t0 is the current reoptimization time (not always 0)
tf_rel_guess = pi; % pi is good guess when t0 = 0, will update as tf-t0

% helpful IC, replaced by lambda_f when re-optimizing
lambda10 = 22; 
lambda20 = -7;
lambda30 = 20;
lambda40 = -5;
lambda0_guess = [lambda10; lambda20; lambda30; lambda40];

% Indirect Optimization Function
% outputs: alpha(t), sample t's for alpha, absolute tf (tf_rel+t0), lambdas
% at next reopt time
[alpha, alpha_t, tf, lambda_f] = indirect_fcn(Chaser, Target, Nav, t0,...
    plot_opt, lambda0_guess, tf_rel_guess);

% EKF Function to propagate State Covariance P in continuous time with
% Acceleration process covariance
Nav = prop_EKF(Nav,Chaser,Cov,alpha,alpha_t,t0,tf);

% Propagate actual state with interpolated process noise
Actual = prop_actual(Actual,Chaser,Cov,alpha,alpha_t,t0,tf);

% Hacked together code for update at measurement time
H = eye(4);
G = eye(4);
L_k = Nav.P*H'*(H*Nav.P*H'+G*Cov.Z*G')^-1;
y = H*Actual.X + Cov.Z*randn(4,1);
Xhat_minus = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
Xhat_plus = Xhat_minus + L_k*(y - H*Xhat_minus);
Nav.r_plus = Xhat_plus(1);
Nav.theta_plus = Xhat_plus(2);
Nav.rdot_plus = Xhat_plus(3);
Nav.thetadot_plus = Xhat_plus(4);
Nav.P_plus = (eye(4)-L_k*H)*Nav.P;