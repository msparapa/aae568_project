function [r, theta, x, y, alpha, gamma, alpha_hor, i_quiv] = unpack_X(X)
% Collin York
% AAE 508
% Final Project
% Minimum Time for Circular Orbit Phase Change
%
% State Histories Unpacking
r = X(:,1);
theta = X(:,2);

x = r.*cos(theta);
y = r.*sin(theta);

lambda3 = X(:,7);
lambda4 = X(:,8);

cos_gamma = -lambda3.*r./sqrt(lambda3.^2 .* r.^2 + lambda4.^2);
sin_gamma = -lambda4./sqrt(lambda3.^2 .* r.^2 + lambda4.^2);
gamma = unwrap(atan2(sin_gamma,cos_gamma));
alpha = unwrap(gamma + theta);

% Control History (psi rel. to local horizon)
alpha_hor = theta + pi/2 - alpha;
if alpha_hor(1) <= -pi
    alpha_hor = alpha_hor + 2*pi;
end

i_quiv = round(linspace(1, length(r), 15));

return