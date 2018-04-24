function PSI = indirect_bcs(Y0, Yf, tf_rel, Chaser, Target, Nav, t0, slack)
% INDIRECT_BCS Compute the boundary conditions for the indirect
% optimization method
%
%   psi = indirect_bcs(Y0, Yf, tf, Chaser, Target, Nav, t0)
%
%   Inputs:
%
%       - Y0: Initial state vector
%       - Yf: Final state vector
%       - tf_rel: final time relative to beginning of segment (i.e., tof)
%       - Chaser: Structure that stores data about the Chaser s/c
%       - Target: Structure that stores data about the Target s/c
%       - Nav: Structure that stores data about the current estimated state
%       of the Chaser s/c
%       - t0: Epoch of the beginning of the segment
%
%   Outputs:
%
%       - psi: column vector of boundary conditions; each entry should 
%       evaluate to zero when satisfied
%
% Collin York
% AAE 508
% Final Project
% Minimum Time for Circular Orbit Phase Change BCs


dYf_dtf = indirect_odes(0,Yf,tf_rel,Chaser)/tf_rel;

PSI = [Y0(1) - Nav.r;
    Y0(2) - Nav.theta;
    Y0(3) - Nav.rdot;
    Y0(4) - Nav.thetadot;
    Yf(1) - Target.r0;
    Yf(2) - Target.thetadot0*(tf_rel+t0) - Target.theta0;
    Yf(3) - Target.rdot0;
    Yf(4) - Target.thetadot0;
    Yf(6)*Target.thetadot0 - Yf(5:8)'*dYf_dtf(1:4) - 1;
    tf_rel - slack^2];
return