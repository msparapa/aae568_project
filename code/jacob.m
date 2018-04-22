function [A, C] = jacob(x)
global mu J2 Re
x1 = x(1);
x3 = x(3);
x4 = x(4);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     x4^2 - 2*mu*(x1^2-3*J2*Re^2)/x1^5, 0, 0, 2*x1*x4;
     2*x3*x4/x1^2, 0, -(2*x4)/x1, -2*x3/x1];

C = [1 0 0 0;0 0 1 0];