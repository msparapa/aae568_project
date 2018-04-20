function dfdt = prop_EKF_odes_test(~,f)
global alpha mu m T J2 Re
r         = f( 1);
theta     = f( 2);
rDot      = f( 3);
thetaDot  = f( 4);

rDotDot = r*thetaDot^2 - mu/r^2*(1-3/2*J2*(Re/r)^2) +...
    T/m*(cos(alpha)*cos(theta)+sin(alpha)*sin(theta));
thetaDotDot = -2*rDot*thetaDot/r +...
    T/m/r*(sin(alpha)*cos(alpha)-cos(alpha)*sin(theta));

dfdt = [rDot; thetaDot; rDotDot; thetaDotDot];
