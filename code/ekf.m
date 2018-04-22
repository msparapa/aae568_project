function [Xplus, Pplus] = ekf(h,Xplus,Pplus,y,Q,R,Ajaco,Hjaco)
% x = sym('x', [4, 1]);
% Propagation equations
% Xminus and Pminus are the propagated state and
% covariance
Xminus = twobodyPolarEKF(Xplus,[0 60],1);
% A is the Jacobian matrix at Xplus
% A = subs(Ajaco, x, Xplus);
Pminus = Ajaco*Pplus*Ajaco'+Q;

% hx is the predicted measurement
hx = h(Xminus)';
% H is the Jacobian matrix at Xminus
% H = subs(Hjaco, x, Xminus);
% Compute Kalman Gain
L = Pminus*Hjaco'/(Hjaco*Pminus*Hjaco'+R);
% Measurement Update equations
% Xplus and Pplus are the updated state and
% covariance
if y == 0
    Xplus = Xminus;
else
    Xplus = Xminus + L*(y - hx);
end
Pplus = Pminus - L*Hjaco*Pminus;
end