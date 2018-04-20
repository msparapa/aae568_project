function dX_dtau = direct_odes(t, X, u, params, consts)

r = X(1,:);
theta = X(2,:);
rdot = X(3,:);
thetadot = X(4,:);

tf_rel = params(1);
slack = params(2);

Chaser = consts.Chaser;

m = Chaser.m0 - Chaser.mdot*t*tf_rel;

cos_gamma = cos(u);
sin_gamma = sin(u);
dX_dt = zeros(size(X));

dX_dt(1,:) = X(3);
dX_dt(2,:) = X(4);
dX_dt(3,:) = r.*thetadot.^2 - 1./r.^2 + Chaser.T./m.*cos_gamma;
dX_dt(4,:) = -2*rdot.*thetadot./r + Chaser.T./m./r.*sin_gamma;

dX_dtau = tf_rel*dX_dt;
return