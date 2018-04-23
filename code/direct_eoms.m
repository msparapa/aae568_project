function dX_dtau = direct_eoms(t, X, u, params, consts)

r = X(1,:);
theta = X(2,:);
rdot = X(3,:);
thetadot = X(4,:);

tf_rel = params(1);

Chaser = consts.Chaser;

m = Chaser.m0 - Chaser.mdot*t*tf_rel;

gamma = u - theta;

cos_gamma = cos(gamma);
sin_gamma = sin(gamma);
dX_dt = zeros(size(X));

dX_dt(1,:) = rdot;
dX_dt(2,:) = thetadot;
dX_dt(3,:) = r.*thetadot.^2 - 1./r.^2 + Chaser.T./m.*cos_gamma;
dX_dt(4,:) = -2*rdot.*thetadot./r + Chaser.T./m./r.*sin_gamma;

dX_dtau = tf_rel*dX_dt;
return