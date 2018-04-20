function [f,Pdot] = twobodyPolarEKF(f0,P,R,tspan,dt)
global alpha mu m T J2 Re
%...Pass the initial conditions and time interval to ode45,
%...which calculates the position and velocity at discrete
%...times t, returning the solution in the column vector f.
%...ode45 uses the m-function 'accel_polar' to evaluate the
%...acceleration at each integration time step.
tol = 1e-10;
options = odeset('RelTol',tol,'AbsTol',tol);
[~,p] = ode45('prop_EKF_odes_test', [tspan(1):dt:tspan(2)], f0,options);
f = p(end,:)'; 

A = [0, 0, 1, 0;...
    0, 0, 0, 1;...
    f0(4)+2/f0(1)^3, T/m*sin(alpha), 0, f0(1);...
    2*f0(2)*f0(4)/f0(1)^2 - T/m/f0(1)^2*sin(alpha), -T/m/f0(1)*cos(alpha), -2*f0(4)/f0(1), -2*f0(2)/f0(1)];
C = [0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 1, 0;...
    0, 0, 0, 1];

Pdot = A*P + P*A' + C*R*C';

