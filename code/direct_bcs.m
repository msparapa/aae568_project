function PSI = direct_bcs(t0, X0, u0, tf, Xf, uf, quads0, quadsf, params, consts) %Y0, Yf, tf_rel, Chaser, Target, Nav, t0, slack)

t0 = consts.t0;
Nav = consts.Nav;
Target = consts.Target;

tf_rel = params(1);
slack = params(2);

% dYf_dtf = direct_odes(0,Yf,tf_rel,Chaser)/tf_rel;

PSI = [X0(1) - Nav.r;
    X0(2) - Nav.theta;
    X0(3) - Nav.rdot;
    X0(4) - Nav.thetadot;
    Xf(1) - Target.r0;
    Xf(2) - Target.thetadot0*(tf_rel+t0) - Target.theta0;
    Xf(3) - Target.rdot0;
    Xf(4) - Target.thetadot0;
%     Yf(6) - Yf(5:8)'*dYf_dtf(1:4) - 1;
    tf_rel - slack^2];
return