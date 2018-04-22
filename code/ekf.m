function [Xplus, Pplus] = ekf(h,Xplus,Pplus,y,Q,R)
% x = sym('x', [4, 1]);
% Propagation equations
% Xminus and Pminus are the propagated state and
% covariance
Xminus = twobodyPolarEKF(Xplus,[0 60],60);
% A is the Jacobian matrix at Xplus
% A = subs(Ajaco, x, Xplus);
[A, ~] = jacob(Xplus);
Pminus = A*Pplus*A'+Q;

% hx is the predicted measurement
hx = h(Xminus)';
% H is the Jacobian matrix at Xminus
% H = subs(Hjaco, x, Xminus);
[~, H] = jacob(Xminus);
% Compute Kalman Gain
L = Pminus*H'/(H*Pminus*H'+R);
% Measurement Update equations
% Xplus and Pplus are the updated state and
% covariance
if y == 0
    Xplus = Xminus;
    Pplus = Pminus;
else
    Xplus = Xminus + L*(y - hx);
    Pplus = Pminus - L*H*Pminus;
end
end