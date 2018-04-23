function dfdt = prop_UT_odes(~,f)
global alpha mu m T R Q
r         = f( 1);
theta     = f( 2);
rDot      = f( 3);
thetaDot  = f( 4);

rDotDot = r*thetaDot^2 - mu/r^2 +...
    T/m*(cos(alpha)*cos(theta)+sin(alpha)*sin(theta));
thetaDotDot = -2*rDot*thetaDot/r +...
    T/m/r*(sin(alpha)*cos(alpha)-cos(alpha)*sin(theta));

dfdt = [rDot; thetaDot; rDotDot; thetaDotDot];

if(length(f) == 4 + 4*4)
    Pplus = reshape(f(5:end), 4, 4);
    
    A = [0, 0, 1, 0;
        0, 0, 0, 1;
        thetaDot^2 + 2*mu/r^3, T/m*(-sin(theta)*cos(alpha)+sin(alpha)*cos(theta)), 0, 2*r*thetaDot;
        2*rDot*thetaDot/r^2-T/m/r^2*sin(alpha-theta), -T*cos(theta-alpha)/m/r, -(2*thetaDot)/r, -2*rDot/r];
    C = [0, 0, 0, 0;...
        0, 0, 0, 0;...
        0, 0, 1, 0;...
        0, 0, 0, 1];
    
    Pminus = A*Pplus*A'+Q + C*R*C';

    dfdt(5:20) = reshape(Pminus, 16, 1);
end