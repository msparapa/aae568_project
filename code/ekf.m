function [Xplus, Pplus] = ekf(h,Xplus,Pplus,y,R,Q)
% x = sym('x', [4, 1]);
% Start by reshaping Pplus to send through integrator
Xplus_reshaped = [Xplus;reshape(Pplus, 16, 1)];
% Propagation equations
% Xminus and Pminus are the propagated state and
% covariance
% Xminus = twobodyPolarEKF(Xplus_reshaped,[0 60],60);
Xminus = twobodyPolarEKF(Xplus,[0 60],60);
% Split state and covariance
% Pminus = reshape(Xminus(5:end), 4, 4);
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
    Xplus = Xminus(1:4,1);
    Pplus = Pminus;
else
    Xplus = Xminus(1:4,1) + L*(y - hx);
    Pplus = Pminus - L*H*Pminus;
end
end