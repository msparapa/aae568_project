function [BC] = brachistochrone_bcs_unreduced(t0, X0, u0, tf, Xf, uf, quads0, quadsf, params, consts)

x0 = X0(1);
y0 = X0(2);
v0 = X0(3);
lx0 = X0(4);
ly0 = X0(5);
lv0 = X0(6);
xf = Xf(1);
yf = Xf(2);
vf = Xf(3);
lxf = Xf(4);
lyf = Xf(5);
lvf = Xf(6);
tf = params(1);
g = consts(1);

% theta0 = brachistochrone_control_unreduced(t0, X0, params, consts);
% thetaf = brachistochrone_control_unreduced(tf, Xf, params, consts);

theta0 = u0;
thetaf = uf;

BC(1) = x0;
BC(2) = y0;
BC(3) = v0;
BC(4) = xf - 1;
BC(5) = yf + 1;

end

