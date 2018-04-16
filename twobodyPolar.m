function f = twobodyPolar(f0,tspan,dt)
%...Pass the initial conditions and time interval to ode45,
%...which calculates the position and velocity at discrete
%...times t, returning the solution in the column vector f.
%...ode45 uses the m-function 'accel_polar' to evaluate the
%...acceleration at each integration time step.
tol = 1e-10;
options = odeset('RelTol',tol,'AbsTol',tol);
[~,p] = ode45('accel_polar', [tspan(1):dt:tspan(2)], f0,options);
f = p(:,:); 
