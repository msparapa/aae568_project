clear all; close all; clc; 

% Variances of the process and measurement noise
q = 0.001*0.001;
r = 0.01*0.01;
% Sampling interval
T = 0.02;
% Initial values of the state and covariance
Xplus = [0; 0; 0; 130/3.6; 0.01];
Pplus = eye(5);

% Nonlinear state equation (f)
f = @(x)[(x(4)/x(5))*sin(x(5)*T+x(3))-(x(4)/x(5))*sin(x(3))+x(1);
    -(x(4)/x(5))*cos(x(5)*T+x(3))+(x(4)/x(5))*cos(x(3))+x(2);
    x(5)*T+x(3);
    x(4);
    x(5)];
% Calculate jacobian matrix
x = sym('x', [5, 1]);
Ajaco = jacobian(f(x));

% Measurement equation (h)
h = @(x)[x(4);x(5)];
% Calculate Jacobian matrix
Hjaco1 = jacobian(h(x));
Hjaco = [zeros(2,3) Hjaco1];

for i = 1:500
% f is the nonlinear state equations
% Xplus is the state
% Pplus is the covariance
% h is the nonlinear measurement equation
% y is the measurement
% Q is the covariance of the process noise
% R is the covariance of the measurement noise
% Ajaco is the jacobian matrix of the state equations
% Hjaco is the jacobian matrix of the measurement
% equations
[Xplus, Pplus] = ekf(f,h,Xplus,Pplus,y,Q,R,Ajaco,Hjaco);
end

% Estimate Xplus and Pplus
function [Xplus, Pplus] = ekf(f,h,Xplus,Pplus,y,Q,R,Ajaco,Hjaco)
x = sym('x', [5, 1]);
% Propagation equations
% Xminus and Pminus are the propagated state and
covariance
Xminus = f(Xplus);
% A is the Jacobian matrix at Xplus
A = subs(Ajaco, x, Xplus);
Pminus = A*Pplus*A'+Q;

% hx is the predicted measurement
hx = h(Xminus);
% H is the Jacobian matrix at Xminus
H = subs(Hjaco, x, Xminus);
% Compute Kalman Gain
L = Pminus*H'/(H*Pminus*H'+R);
% Measurement Update equations
% Xplus and Pplus are the updated state and
covariance
Xplus = Xminus + L*(y - hx);
Pplus = Pminus - L*H*Pminus;
end