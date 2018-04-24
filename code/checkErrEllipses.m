function [bInPos, bInVel] = checkErrEllipses(Nav, Actual, Target, t, plot_opt)
% CHECKERRELLIPSES Construct the error ellipses and check for convergence
%
%   [bInPos, bInVel] = checkErrEllipses(Nav, Actual, Target, t, plot_opt)
%
%   Inputs:
%       
%       - Nav: Structure containing navigation information
%       - Actual: Structure containing true spacecraft state
%       - Target: Structure containing target orbit information
%       - t: Current time, relative to mission start
%       - plot_opt: Structure containing plot options
%
%   Outputs:
%
%       - bInPos: Whether or not the target position states (r, theta) are
%       within the error ellipse constructed from the navigation covariance
%       in Nav.P
%       - bInVel: Whether or not the target velocity states (rdot, thetadot) are
%       within the error ellipse constructed from the navigation covariance
%       in Nav.P
%
%   Author: Andrew Cox
%   Version: April 20, 2018

%% First, r-theta ellipse
[R, val] = eig(Nav.P(1:2,1:2));

% Semi-major and semi-minor axis sizes
a = sqrt(val(1,1));
b = sqrt(val(2,2));
x = linspace(-a, a, 5e2).';

% Compute ellipse in local coordinates
Ey_local = sqrt(b^2 - b^2*x.^2/a^2);
E_local = [x, Ey_local; flipud(x), flipud(-Ey_local)];

% Normalize eigenvectors
R(:,1) = R(:,1)/norm(R(:,1));
R(:,2) = R(:,2)/norm(R(:,2));

% Rotate ellipse by eigenvectors to put in "global" space
E_global = E_local*R;
E_global(:,1) = E_global(:,1) + Nav.r;
E_global(:,2) = E_global(:,2) + Nav.theta;

% Plot Ellipse
if(plot_opt.nav)
    figure(); hold on;
    plot(E_global(:,1), E_global(:,2), 'k', 'linewidth', 2);
    plot(Nav.r, Nav.theta, 'kd', 'markerfacecolor', 'k');
    plot(Actual.X(1), Actual.X(2), 'rd', 'markerfacecolor', 'r');
    plot(Target.r0, Target.theta0 + Target.thetadot0*t, 'bd', ...
        'markerfacecolor', 'b');

    hold off; grid on;
    legend('1 \sigma Error', 'Nav State', 'Actual State', 'Target State',...
        'location', 'best');
    xlabel('r, nondim');
    ylabel('\theta, nondim');
end

%% Second, rdot-thetadot ellipse
[R, val] = eig(Nav.P(3:4, 3:4));

% Semi-major and semi-minor axis sizes
a = sqrt(val(1,1));
b = sqrt(val(2,2));
x = linspace(-a, a, 5e2).';

% Compute ellipse in local coordinates
Ey_local = sqrt(b^2 - b^2*x.^2/a^2);
E_local = [x, Ey_local; flipud(x), flipud(-Ey_local)];

% Normalize eigenvectors
R(:,1) = R(:,1)/norm(R(:,1));
R(:,2) = R(:,2)/norm(R(:,2));

% Rotate ellipse by eigenvectors to put in "global" space
E_global2 = E_local*R;
E_global2(:,1) = E_global2(:,1) + Nav.rdot;
E_global2(:,2) = E_global2(:,2) + Nav.thetadot;

% Plot Ellipse
if(plot_opt.nav)
    figure(); hold on;
    plot(E_global2(:,1), E_global2(:,2), 'k', 'linewidth', 2);
    plot(Nav.rdot, Nav.thetadot, 'kd', 'markerfacecolor', 'k');
    plot(Actual.X(3), Actual.X(4), 'rd', 'markerfacecolor', 'r');
    plot(Target.rdot0, Target.thetadot0, 'bd', 'markerfacecolor', 'b');

    hold off; grid on;
    legend('1 \sigma Error', 'Nav State', 'Actual State', 'Target State',...
        'location', 'best');
    xlabel('$\dot{r}$, nondim', 'interpreter', 'latex');
    ylabel('$\dot{\theta}$, nondim', 'interpreter', 'latex');
end

%% Check to see if the target is within the error ellipses
bInPos = pointLiesInside([Target.r0, Target.theta0 + t*Target.thetadot0],...
    E_global);
bInVel = pointLiesInside([Target.rdot0, Target.thetadot0], E_global2);

if(bInPos)
    fprintf('Target (r, theta) is INSIDE error ellipse\n');
else
    fprintf('Target (r, theta) is OUTSIDE error ellipse\n');
end

if(bInVel)
    fprintf('Target (rdot, thetadot) is INSIDE error ellipse\n');
else
    fprintf('Target (rdot, thetadot) is OUTSIDE error ellipse\n');
end
end