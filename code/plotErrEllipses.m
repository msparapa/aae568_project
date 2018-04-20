function [] = plotErrEllipses(Nav, Actual, Target, t, plot_opt)

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

% Plot Ellipse
plot_opt.i = plot_opt.i + 1;
figure(plot_opt.i); hold on;
plot(Nav.r + E_global(:,1), Nav.theta + E_global(:,2), 'k', 'linewidth', 2);
plot(Nav.r, Nav.theta, 'kd', 'markerfacecolor', 'k');
plot(Actual.X(1), Actual.X(2), 'rd', 'markerfacecolor', 'r');
plot(Target.r0, Target.theta0 + Target.thetadot0*t, 'bd', ...
    'markerfacecolor', 'b');
% plot(Nav.r + [0, a*R(1,1)], Nav.theta + [0, a*R(2,1)]);
% plot(Nav.r + [0, b*R(1,2)], Nav.theta + [0, b*R(2,2)]);
hold off; grid on; axis equal;
legend('1 \sigma Error', 'Nav State', 'Actual State', 'Target State');
xlabel('r, nondim');
ylabel('\theta, nondim');

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
E_global = E_local*R;

% Plot Ellipse
plot_opt.i = plot_opt.i + 1;
figure(plot_opt.i); hold on;
plot(Nav.rdot + E_global(:,1), Nav.thetadot + E_global(:,2), 'k', ...
    'linewidth', 2);
plot(Nav.rdot, Nav.thetadot, 'kd', 'markerfacecolor', 'k');
plot(Actual.X(3), Actual.X(4), 'rd', 'markerfacecolor', 'r');
plot(Target.rdot0, Target.thetadot0, 'bd', 'markerfacecolor', 'b');
% plot(Nav.r + [0, a*R(1,1)], Nav.theta + [0, a*R(2,1)]);
% plot(Nav.r + [0, b*R(1,2)], Nav.theta + [0, b*R(2,2)]);
hold off; grid on; axis equal;
legend('1 \sigma Error', 'Nav State', 'Actual State', 'Target State');
xlabel('$\dot{r}$, nondim', 'interpreter', 'latex');
ylabel('$\dot{\theta}$, nondim', 'interpreter', 'latex');
%%
% keyboard;
end