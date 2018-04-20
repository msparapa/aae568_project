function [J] = direct_cost(t, X, params, consts)

tf_rel = params(1);
Chaser = consts.Chaser;
mfinal = Chaser.m0 - Chaser.mdot*t(end)*tf_rel;

J = mfinal;
end

