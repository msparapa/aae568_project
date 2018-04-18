function [alpha, alpha_t, tf, t_seg, lambda_seg] = indirect_fcn(Chaser, Target,...
    Nav, t0, plot_opt, lambda0_guess, tf_rel_guess)
% INDIRECT_FCN Solve the indirect optimization problem
%
%   [...] = indirect_fcn(Chaser, Target, Nav, t0, plot_opt, lambda0_guess,
%                       tf_rel_guess)
%
%       Solves the indirect optimization problem of guiding the Chaser
%       spacecraft to rendezvous with the Target spacecraft.
%
%       Although the entire optimal trajectory is computed, the problem is
%       structured to only fly part of the trajectory before obtaining a
%       new observation, updating the navigation estimate, and recomputing
%       the optimal trajectory.
%
%   Inputs:
%
%       - Chaser: Structure containing information about the Chaser s/c
%       - Target: Structure containing information about the Target s/c
%       - Nav: Structure containing information about the current Nav state
%       - t0: Initial epoch associated with the current segment
%       - plot_opt: Options for plotting
%       - lambda0_guess: Initial conditions for the costates
%       - tf_rel_guess: Guess for propagation time to reach final state
%
%   Outputs: [alpha, alpha_t, tf, lambda_f]
%
%       - alpha: Control history for the entire optimal trajectory
%       - alpha_t: times associated with alpha
%       - tf: total time on the optimal trajectory
%       - t_seg: time at the end of "the segment"
%       - lambda_seg: costate values at the end of "the segment", i.e., 
%       part way through the optimal trajectory.
%
%   Author: Collin York
%   Version: April 18, 2018

yinit = [Nav.r0; Nav.theta0; Nav.rdot0; Nav.thetadot0; lambda0_guess];
Nt = 500;
tau = linspace(0, 1, Nt)'; % non-dimensional time vector

solinit = bvpinit(tau, yinit, tf_rel_guess);

bvp_opts = bvpset('Stats','on');

odes = @(tau, X, tf) indirect_odes(tau, X, tf, Chaser);
bcs = @(Y0, Yf, tf) indirect_bcs(Y0, Yf, tf, Chaser, Target, Nav, t0);

sol = bvp4c(odes, bcs, solinit, bvp_opts);
tf = sol.parameters + t0;
X0 = sol.y(:,1);

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
odes_tf = @(t, X) indirect_odes_tf(t, X, Chaser);
[alpha_t,X] = ode113(odes_tf, tau*tf, X0, options);

[r, theta, x, y, alpha, gamma, alpha_hor, i_quiv] = unpack_X(X);

i_ts_opt = find(alpha_t < Chaser.ts_opt, 1, 'last');
t_seg = alpha_t(i_ts_opt);
lambda_seg = X(i_ts_opt, 5:8).';

if plot_opt.indirect
    plot_opt.i = plot_opt.i + 1;
    figure(plot_opt.i);
    plot(x, y);
    axis equal;
    hold on;
    plot(Target.r0*cos(Target.thetadot0*(linspace(t0,tf,100))+Target.theta0),...
        Target.r0*sin(Target.thetadot0*(linspace(t0,tf,100))+Target.theta0));
    hold off;

    plot_opt.i = plot_opt.i + 1;
    figure(plot_opt.i);
    plot(r);
    hold on;
    plot(theta);
    hold off;

    plot_opt.i = plot_opt.i + 1;
    figure(plot_opt.i);
    plot(alpha/pi*180);

    close all;
    plot_opt.i = plot_opt.i + 1;
    figure(plot_opt.i);                                 
    hold on;
    plot(x,y,'b');
    plot(Target.r0*cos(Target.thetadot0*(linspace(t0,tf,100))+Target.theta0),...
        Target.r0*sin(Target.thetadot0*(linspace(t0,tf,100))+Target.theta0),'k--');
    quiver(x(i_quiv),y(i_quiv),cos(alpha(i_quiv)),sin(alpha(i_quiv)),0.5,'r');
    scatter(x(1),y(1),'r');
    scatter(x(end),y(end),'g');
    title('C. York - Min-Time Trajectory - Case 3');
    xlabel('x (DU)');
    ylabel('y (DU)');
    legend('Solution','Ref.','Thrust Direction','t_0','t_f')
    axis equal;
    hold off;

    plot_opt.i = plot_opt.i + 1;
    figure(plot_opt.i);
    plot(alpha_t, alpha_hor*180/pi, 'k');
    title('C. York - Control Time History - Case 3');
    xlabel('time (TU)');
    ylabel('$$\bar{\alpha}$$ (deg)','interpreter','latex','fontsize',12);
    return
end