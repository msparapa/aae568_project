function [dX] = brachistochrone_eoms_unreduced(t, X, params, consts)

x = X(1,:);
y = X(2,:);
v = X(3,:);
lx = X(4,:);
ly = X(5,:);
lv = X(6,:);
tf = params(1);
g = consts(1);

theta = brachistochrone_control_unreduced(t, X, params, consts);

N = length(x);
dX = zeros(6,N);

dX(1,:) = v.*cos(theta);
dX(2,:) = v.*sin(theta);
dX(3,:) = -g*sin(theta);
dX(4,:) = 0;
dX(5,:) = 0;
dX(6,:) = -lx.*cos(theta) - ly.*sin(theta);
dX = dX*tf;

end