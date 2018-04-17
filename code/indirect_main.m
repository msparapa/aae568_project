% Collin York
% Initialize Indirect Optimization

clear all; clc; close all;

Chaser.T = 0.05;
Chaser.mdot = 0.001;
Chaser.m0 = 1;
Chaser.ts_opt = 1;

Actual.r0 = 1;
Actual.theta0 = 0;
Actual.rdot0 = 0;
Actual.thetadot0 = 1;

Nav.r0 = 1;
Nav.theta0 = 0;
Nav.rdot0 = 0;
Nav.thetadot0 = 1;
Nav.P0 = zeros(4);

Target.r0 = 1.01;
Target.theta0 = 0.21;
Target.rdot0 = 0;
Target.thetadot0 = sqrt(1/Target.r0^3);

Cov.R = 1e-5;

plot_opt.indirect = 1;
plot_opt.i = 0;

t0 = 0;
tf_rel_guess = pi;
lambda10 = 22;
lambda20 = -7;
lambda30 = 20;
lambda40 = -5;
lambda0_guess = [lambda10; lambda20; lambda30; lambda40];

[alpha, alpha_t, tf, lambda_f] = indirect_fcn(Chaser, Target, Nav, t0, plot_opt, lambda0_guess, tf_rel_guess);

Nav = prop_EKF(Nav,Chaser,Cov,alpha,alpha_t,t0,tf);
