function PSI = indirect_bcs(Y0,Yf,tf,Chaser,Target,Nav,t0)
%
% Collin York
% AAE 508
% Final Project
% Minimum Time for Circular Orbit Phase Change BCs


dYf_dtf = indirect_odes(0,Yf,tf,Chaser)/tf;

PSI = [Y0(1) - Nav.r0;
    Y0(2) - Nav.theta0;
    Y0(3) - Nav.rdot0;
    Y0(4) - Nav.thetadot0;
    Yf(1) - Target.r0;
    Yf(2) - Target.thetadot0*(tf+t0) - Target.theta0;
    Yf(3) - Target.rdot0;
    Yf(4) - Target.thetadot0;
    Yf(6) - Yf(5:8)'*dYf_dtf(1:4) - 1];
return