% Collin York
% AAE 508
% Final Project
% Minimum Time for Circular Orbit Phase Change
%
% State and Cosstate Equations: phase_odes.m
% BCs: phase_bc.m
%
function [alpha, alpha_t, tf, lambda_f] = indirect_fcn(Chaser,Target,Nav,...
    t0,plot_opt,lambda0_guess,tf_rel_guess)


yinit = [Nav.r0; Nav.theta0; Nav.rdot0; Nav.thetadot0; lambda0_guess];
Nt = 500;
tau = linspace(0,1,Nt)'; % non-dimensional time vector
solinit = bvpinit(tau,yinit,tf_rel_guess);
bvp_opts = bvpset('Stats','on');
sol = bvp4c(@(tau,X,tf)indirect_odes(tau,X,tf,Chaser),...
    @(Y0,Yf,tf)indirect_bcs(Y0,Yf,tf,Chaser,Target,Nav,t0), solinit, bvp_opts);
tf = sol.parameters+t0;
X0 = sol.y(:,1);

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[alpha_t,X] = ode113(@(t,X)indirect_odes_tf(t,X,Chaser),[t0 tf],X0,options);

[r, theta, x, y, alpha, gamma, alpha_hor, i_quiv] = unpack_X(X);

for i = 1:length(alpha_t)
    if alpha_t(i) < Chaser.ts_opt
        i_ts_opt = i;
    else
        break
    end
end
lambda_f = X(i_ts_opt,5:8);

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