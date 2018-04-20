function [alpha, alpha_t, tf, sol] = direct_fcn(Chaser,...
    Target, Nav, t0, plot_opt, tf_rel_guess)

yinit = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
slack_guess = sqrt(tf_rel_guess);
ICsolver0 = [tf_rel_guess; slack_guess];

options = optimset('Algorithm','sqp','display','off');
options.Nodes = 50;

T = linspace(0, 1, options.Nodes);
solinit0 = bvpinit(T, yinit', ICsolver0);
solinit0.consts.Chaser = Chaser;
solinit0.consts.Target = Target;
solinit0.consts.Nav = Nav;
solinit0.consts.t0 = t0;
solinit0.control(1,:) = 0*linspace(0, 0, options.Nodes);

options.isdirect = 1;

tic;
sol = bvpmc(@direct_eoms, [], @direct_bcs, solinit0, options);
time0 = toc;

alpha = sol.control;
alpha_t = linspace(0,1,options.Nodes);
tf = sol.parameters(1);

