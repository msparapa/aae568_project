clc; close all; clear variables;
addpath('..\optim');
g = 10;

options = optimset('Algorithm','sqp','display','off');
options.Nodes = 20;

xguess = linspace(0, 1, options.Nodes);
yguess = linspace(0, -1, options.Nodes);
vguess = linspace(0, 2.5, options.Nodes);
lxguess = -0.1*linspace(1, 1, options.Nodes);
lyguess = 0.1*linspace(1, 1, options.Nodes);
lvguess = -0.1*linspace(1, 0, options.Nodes);
uguess = 0*linspace(0, 0, options.Nodes);

tfguess = 0.5;
T = linspace(0, 1, options.Nodes);
solinit0 = bvpinit(T, [0,0,0,0,0,0], tfguess);
solinit0.y(1,:) = xguess;
solinit0.y(2,:) = yguess;
solinit0.y(3,:) = vguess;
solinit0.y(4,:) = lxguess;
solinit0.y(5,:) = lyguess;
solinit0.y(6,:) = lvguess;
solinit0.consts = [g];
solinit0.control(1,:) = uguess;

options.isdirect = 1;

tic;
sol = bvpmc(@brachistochrone_eoms_unreduced, [], @brachistochrone_bcs_unreduced, solinit0, options);
time0 = toc;
fprintf('M time : \t %.4f\n', time0);

figure();
subplot(2,1,1);
plot(sol.y(1,:), sol.y(2,:));
subplot(2,1,2);
stairs(sol.control);