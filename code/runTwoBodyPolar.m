clear all; close all; clc;
global alpha g0 Isp T mu Re m J2
alpha = 1*pi/180;          % thrust angle
J2    = 1.0826269e-3; 
g0    = 9.82;              % gravitational acceleration, m/s^2
Isp   = 1000;              % specific-impulse, s
T     = 1e-3;              % thrust, N
mu    = 3.98600e+05;       % standard gravitational parameter, km^3/s^2
Re    = 6378;              % earth radius, km
m     = 100;               % spacecraft mass, kg
h     = 250;               % orbit altitude, km
a     = Re + h;            % orbit radius, km (need to fix for eccentric orbits)
P     = 2*pi*sqrt(a^3/mu); % orbit period, s
w     = 2*pi/P;            % orbit angular velocity, rad/s
dt    = 60;                % propagation step-size
revs  = 4;                 % revs to propagate

P     = floor(P/60)*60;    % even out propagation duration (for storing information)

mdot  = -T/Isp/g0;         % mass flow rate
x0    = [a;0;0;w;mdot];
f     = twobodyPolar(x0,[0 revs*P],dt);

% Plot the true orbit
figure(); polar(f(:,2),f(:,1),'b-'); grid on; hold on;

%% Some abitrary covariance matrix 
% C = [ 1.4516e-03  -5.0018e-09  -1.4306e-05  -4.0265e-10  -2.2077e-14
%      -5.0018e-09   1.0026e-09   7.1661e-09  -7.3595e-14  -5.0066e-17
%      -1.4306e-05   7.1661e-09   8.4438e-07   4.3113e-12   1.4098e-16
%      -4.0265e-10  -7.3595e-14   4.3113e-12   1.9591e-16   8.5570e-21
%      -2.2077e-14  -5.0066e-17   1.4098e-16   8.5570e-21   6.0000e-24];

% Arbitrary sigmas 
rSig    = 1e-04;
tSig    = 1e-08;
rDotSig = 1e-08;
tDotSig = 1e-08;
mDotSig = 1e-24;
C = [rSig 0 0 0 0;0 tSig 0 0 0;0 0 rDotSig 0 0;0 0 0 tDotSig 0;0 0 0 0 mDotSig];

%% Setup the Unscented Transform Options
options = struct; options.alpha = 1; options.beta = 2.0;  

% Propagate all sigma points until time of measurement
[statesOut,Wm,Wc,useSquareRoot] = unscentedTransformPolar( x0, C, options, 9, revs*P);

% Preallocation
n = revs*P/dt+1;
intMeans = zeros(n,5);
intCovars = zeros(5,5,n);
for h = 1:n
    [intMean, intCovar] = unscentedTransformNoFcn(Wm, Wc, statesOut(:,:,h), useSquareRoot);
    intMeans(h,:) = intMean;
    intCovars(:,:,h) = intCovar;  
end

polar(statesOut(2,:,1),statesOut(1,:,1),'k.'); hold on;
polar(statesOut(2,:,end),statesOut(1,:,end),'r.');

%% Store Propagated Covariances
n = length(intCovars);
storeCov = zeros(n,4);
for i = 1:n
    storeCov(i,:) = sqrt(diag(intCovars(1:4,1:4,i)));
end

%% Update Propagated Particles using a UKF 
truObs                    = [f(end,1);f(end,3)];
z                         = [3;3e-3];                       % Measurement Noise
h                         = @(j) updatePolarMeasurement(j);                                                                              
x_initial                 = statesOut(:,:,end);             % Sigma Points                   
w                         = 0.03;                           % Process Noise Standard Deviation                              
obs                       = [truObs(1,1);truObs(2,1)];      % [2x1] Measurement                     
num_iterations            = 1;
[x_update, x_updated_vec] = ukf(h,x_initial,w,z,obs,num_iterations);

%% Calculate new covariance and propagate forward
postUpdateCov = cov(x_update');

%% Propagate after update
% Propagate all sigma points 
[statesOutNew,Wm,Wc,useSquareRoot] = unscentedTransformPolar( x_updated_vec, postUpdateCov, options, 9, revs*P);

% Preallocation
n = revs*P/dt+1;
intMeansNew = zeros(n,5);
intCovarsNew = zeros(5,5,n);
for h = 1:n
    [intMean, intCovar] = unscentedTransformNoFcn(Wm, Wc, statesOutNew(:,:,h), useSquareRoot);
    intMeansNew(h,:) = intMean;
    intCovarsNew(:,:,h) = intCovar;     
end

%% Store Propagated Covariances
n = length(intCovarsNew);
storeCovNew = zeros(n,4);
for i = 1:n
    storeCovNew(i,:) = sqrt(diag(intCovarsNew(1:4,1:4,i)));
end

t = linspace(dt,revs*P,revs*P/dt);
t2 = linspace(revs*P,revs*P*2,revs*P/dt);
figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],storeCov(:,1)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeCov(:,1)/1e3,'b-'); hold on;
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t/3600],storeCov(:,3)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeCov(:,3)/1e3,'b-'); hold on;
ylabel('$\dot{r}$, m/sec','Interpreter','latex')
subplot(4,1,3)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],mod(storeCov(:,2)*180/pi,360),'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],mod(storeCovNew(:,2)*180/pi,360),'b-'); hold on;
ylabel('$\theta$, deg','Interpreter','latex')
subplot(4,1,4)
plot([0 t/3600],storeCov(:,4)*180/pi,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeCovNew(:,4)*180/pi,'b-'); hold on;
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('1\sigma Uncertainty');
set(gca,'gridlinestyle','--')