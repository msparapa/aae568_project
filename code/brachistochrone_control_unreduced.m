function [theta] = brachistochrone_control_unreduced(t, X, params, consts)
g = consts(1);
x = X(1,:);
y = X(2,:);
v = X(3,:);
lx = X(4,:);
ly = X(5,:);
lv = X(6,:);

theta1 = 2*atan((lx.*v - sqrt(g^2*lv.^2 - 2*g*lv.*ly.*v + lx.^2.*v.^2 + ly.^2.*v.^2))./(g*lv - ly.*v));
theta2 = 2*atan((lx.*v + sqrt(g^2*lv.^2 - 2*g*lv.*ly.*v + lx.^2.*v.^2 + ly.^2.*v.^2))./(g*lv - ly.*v));

theta = [theta1; theta2];

H1 = -g*lv.*sin(theta1) + lx.*v.*cos(theta1) + ly.*v.*sin(theta1) + 1;
H2 = -g*lv.*sin(theta2) + lx.*v.*cos(theta2) + ly.*v.*sin(theta2) + 1;

Hset = [H1; H2];

[~, ind] = min(Hset);
thetause = zeros(length(v),1);

for ii = 1:length(v)
    thetause(ii) = theta(ind(ii),ii);
end
theta = thetause';
end