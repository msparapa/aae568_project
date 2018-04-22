function [A, C] = jacob(x)
global mu alpha T m
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     x4^2 + 2*mu/x1^3, T/m*(-sin(x2)*cos(alpha)+sin(alpha)*cos(x2)), 0, 2*x1*x4;
     2*x3*x4/x1^2-T/m/x1^2*sin(alpha-x2), -T*cos(x2-alpha)/m/x1, -(2*x4)/x1, -2*x3/x1];

C = [1 0 0 0;0 0 1 0];