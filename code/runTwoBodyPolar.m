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
N     = 100;                % number of measurements 

P     = floor(P/60)*60;    % even out propagation duration (for storing information)

mdot  = -T/Isp/g0;         % mass flow rate
q     = 0.00;              % std dev 
x0    = [a;0;0;w;mdot];    % initial state 
s0    = [x0(1:4)+q*...     % initial state with noise
         randn(4,1);x0(5)];    
f     = twobodyPolar(x0,[0 revs*P],dt);
fTot  = twobodyPolar(x0,[0 2*revs*P+dt],dt);
meas  = twobodyPolar(f(end,:)',[0 N*dt],dt);
meas  = [meas(1:N,1),meas(1:N,3)];

% Plot the true orbit
figure(); polar(f(:,2),f(:,1),'b-'); grid on; hold on;

%% Some abitrary covariance matrix 
C = [ 1.4516e-03  -5.0018e-09  -1.4306e-05  -4.0265e-10  -2.2077e-14
     -5.0018e-09   1.0026e-09   7.1661e-09  -7.3595e-14  -5.0066e-17
     -1.4306e-05   7.1661e-09   8.4438e-07   4.3113e-12   1.4098e-16
     -4.0265e-10  -7.3595e-14   4.3113e-12   1.9591e-16   8.5570e-21
     -2.2077e-14  -5.0066e-17   1.4098e-16   8.5570e-21   6.0000e-24];

% Scale the Covariance up 
sf = 1000;
sfMatrix = [
        sf   0     0     0    0;
         0   sf    0     0    0;
         0    0   sf     0    0;
         0    0    0    sf    0;
         0    0    0     0   sf];
C = sfMatrix*C*transpose(sfMatrix);

% Draw a random sample from the scaled covariance and set this to your
% initial seed
nSamples  = 1;
s0 = repmat(x0,1,nSamples) + chol(C,'lower')*randn(5,nSamples);
 
% Arbitrary sigmas 
% rSig    = 1e-10;
% tSig    = 1e-10;
% rDotSig = 1e-12;
% tDotSig = 1e-12;
% mDotSig = 1e-20;
% C = [rSig 0 0 0 0;0 tSig 0 0 0;0 0 rDotSig 0 0;0 0 0 tDotSig 0;0 0 0 0 mDotSig];

%% Setup the Unscented Transform Options
options = struct; options.alpha = 1; options.beta = 2.0;  

% Propagate all sigma points until time of measurement
[utMeans,utCovars,utSigmaPoints,Wm,Wc] = prop_UT( s0, C, options, 9, revs*P, dt);

%% Store Propagated Sigmas
n = length(utCovars);
storeCov = zeros(n,4);
for i = 1:n
    storeCov(i,:) = sqrt(diag(utCovars(1:4,1:4,i)));
end

%% Update Propagated Particles using a UKF
    ukfCovar                  = utCovars(:,:,end);              % a-priori covariance
    ukfSigmaPoints            = utSigmaPoints(:,:,end);         % sigma points
    estMeans                  = zeros(N,5);                     % preallocate storage for estimated means
    estCovars                 = zeros(5,5,N);                   % preallocate storage for estiamte covariances
for k = 1:N
    truObs                    = meas(k,:)';                     % extract observation at Nth step 
    z                         = [5;5e-3];                       % measurement noise (a noisier measurement is less trustworthy)
    h                         = @(j) updatePolarMeasurement(j);
    w                         = 0;                              % process noise standard deviation
    obs                       = [truObs(1,1);truObs(2,1)];      % single [r,rhoDot] measurement
    num_iterations            = 1;
    [x_update,postUpdateCov]  = ukf(h,ukfSigmaPoints,ukfCovar,w,z,obs,num_iterations,Wm,Wc);
%     [x_update,postUpdateCov]  = EnKF(h,ukfSigmaPoints,w,z,obs,num_iterations);
    [meanOut,covarOut,sigmaPointsOut,~,~] = prop_UT( mean(x_update,2), postUpdateCov, options, 9, dt,dt/dt);  % Propagate measurement update to next time step
    ukfMean                   = meanOut(end,:)';
    ukfCovar                  = covarOut(:,:,end);
    ukfSigmaPoints            = sigmaPointsOut(:,:,end);
    estMeans(k,:)             = ukfMean;
    estCovars(:,:,k)          = ukfCovar;
end

%% Propagate after update
% Propagate all sigma points 
[utMeansNew,utCovarsNew,utSigmaPointsNew,Wm,Wc] = prop_UT( ukfMean, ukfCovar, options, 9, revs*P-dt*N,dt);

utMeansNew  = [estMeans;utMeansNew];

%% Store Propagated Sigmas
n = length(estCovars); l = length(utCovarsNew);
storeCovNew = zeros(n+l,4);
for i = 1:n
    storeCovNew(i,:) = sqrt(diag(estCovars(1:4,1:4,i)));
end
ind = 1;
for i = n+1:n+l
    storeCovNew(i,:) = sqrt(diag(utCovarsNew(1:4,1:4,ind)));
    ind = ind + 1;
end

t  = linspace(dt,revs*P,revs*P/dt);
t2 = linspace(revs*P,revs*P*2,revs*P/dt);
t3 = linspace(dt,revs*P*2+dt,(revs*P*2+dt)/dt);
figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],storeCov(:,1)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeCovNew(:,1)/1e3,'b-'); hold on;
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t/3600],storeCov(:,3)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeCovNew(:,3)/1e3,'b-'); hold on;
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

figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],utMeans(:,1)/1e3,'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],utMeansNew(:,1)/1e3,'b-'); hold on;
plot([0 t3/3600],fTot(:,1)/1e3, 'r--');
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t/3600],utMeans(:,3)/1e3,'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],utMeansNew(:,3)/1e3,'b-'); hold on;
plot([0 t3/3600],fTot(:,3)/1e3, 'r--');
ylabel('$\dot{r}$, m/sec','Interpreter','latex')
subplot(4,1,3)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],mod(utMeans(:,2)*180/pi,360),'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],mod(utMeansNew(:,2)*180/pi,360),'b-'); hold on;
plot([0 t3/3600],mod(fTot(:,2)*180/pi,360), 'r--');
ylabel('$\theta$, deg','Interpreter','latex')
subplot(4,1,4)
plot([0 t/3600],utMeans(:,4)*180/pi,'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],utMeansNew(:,4)*180/pi,'b-'); hold on;
plot([0 t3/3600],fTot(:,4)*180/pi, 'r--');
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('Trajectory');
legend('Traj Prior to Update','Traj After Update','True Traj');
set(gca,'gridlinestyle','--')

% Residuals 
PredTraj = [utMeans; utMeansNew]; % Predicted Trajectory
Diff = PredTraj - fTot;

figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t3/3600],Diff(:,1)*1e3, 'r-'); grid on
ylabel('$\Delta r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t3/3600],Diff(:,3)*1e3, 'r-'); grid on
ylabel('$\Delta\dot{r}$, m/sec','Interpreter','latex')
subplot(4,1,3)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t3/3600],Diff(:,2)*180/pi, 'r-'); grid on
ylabel('$\Delta\theta$, deg','Interpreter','latex')
subplot(4,1,4)
plot([0 t3/3600],Diff(:,4)*180/pi, 'r-'); grid on
ylabel('$\dot{\Delta\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('Trajectory Error');
set(gca,'gridlinestyle','--')
