clear all; close all; clc;
global alpha g0 Isp T mu Re m J2
alpha = 0;                 % thrust angle
J2    = 1.0826269e-3; 
g0    = 9.82;              % gravitational acceleration, m/s^2
Isp   = 1000;              % specific-impulse, s
T     = 0.00;              % thrust, N
mu    = 3.98600e+05;       % standard gravitational parameter, km^3/s^2
Re    = 6378;              % earth radius, km
m     = 100;               % spacecraft mass, kg
h     = 250;               % orbit altitude, km
a     = Re + h;            % orbit radius, km (need to fix for eccentric orbits)
P     = 2*pi*sqrt(a^3/mu); % orbit period, s
w     = 2*pi/P;            % orbit angular velocity, rad/s
dt    = 60;                % propagation step-size
revs  = 4;                 % revs to propagate
N     = 100;               % number of measurements 
q     = 0.01*0.01;         % variance of the process noise
r     = 0.1*0.1;           % variance of the measurement noise
Q     = eye(4)*q;
R     = eye(2)*r;
h     = @(j) updatePolarMeasurement(j);

P     = floor(P/60)*60;    % even out propagation duration (for storing information)

mdot  = -T/Isp/g0;         % mass flow rate
q     = 0.10;              % std dev 
x0    = [a;0;0;w];         % initial state 
% s0    = x0+q*randn(4,1);   % initial state with noise     
f     = twobodyPolar(x0,[0 revs*P],dt);
fTot  = twobodyPolar(x0,[0 2*revs*P+dt],dt);
true  = twobodyPolar(f(end,:)',[0 N*dt],dt);
meas  = [true(1:N,1),true(1:N,3)];

% Plot the true orbit
polar(f(:,2),f(:,1),'b-'); grid on; hold on;

%% Some abitrary covariance matrix 
C = [ 1.4516e-03  -5.0018e-09  -1.4306e-05  -4.0265e-10;
     -5.0018e-09   1.0026e-09   7.1661e-09  -7.3595e-14;
     -1.4306e-05   7.1661e-09   8.4438e-07   4.3113e-12;
     -4.0265e-10  -7.3595e-14   4.3113e-12   1.9591e-16];

% Scale the Covariance up 
sf = 50;
sfMatrix = [
        sf   0     0     0;
         0   sf    0     0;
         0    0   sf     0;
         0    0    0    sf];
C = sfMatrix*C*transpose(sfMatrix);

% Draw a random sample from the scaled covariance and set this to your
% initial seed
nSamples  = 1;
s0 = repmat(x0,1,nSamples) + chol(C,'lower')*randn(4,nSamples);
 
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
Xplus = s0; Pplus = C; y = 0;
n = length(utMeans); 
ekfMeans  = zeros(n,4);
ekfCovars = zeros(4,4,n); 
for k = 1:n
    [A, C]           = jacob(Xplus);
    [Xplus, Pplus]   = ekf(h,Xplus,Pplus,y,Q,R,A,C);
    ekfMeans(k,:)    = Xplus;
    ekfCovars(:,:,k) = Pplus;
end

polar(ekfMeans(:,2),ekfMeans(:,1),'r-'); grid on;

%% Store Propagated Sigmas
n = length(utCovars);
storeUtCov = zeros(n,4);
storeEkfCov = zeros(n,4);
for i = 1:n
    storeUtCov(i,:)  = sqrt(diag(utCovars(:,:,i)));
    storeEkfCov(i,:) = sqrt(diag(ekfCovars(:,:,i)));
end

%% Update Propagated Particles using a UKF
ukfCovar       = utCovars(:,:,end);           % a-priori covariance
ukfSigmaPoints = utSigmaPoints(:,:,end);      % sigma points
estUkfMeans    = zeros(N,4);                  % preallocate storage for estimated means
estUkfCovars   = zeros(4,4,N);                % preallocate storage for estimated covariances
estEkfMeans    = zeros(N,4);                  % preallocate storage for estimated means
estEkfCovars   = zeros(4,4,N);                % preallocate storage for estimated covariances
storeObs       = zeros(N,2);                  % preallocate storage for deviated measurements
Xplus          = ekfMeans(end,:)';
Pplus          = ekfCovars(:,:,end);
for k = 1:N
    truObs                    = meas(k,:)';                     % extract observation at Nth step 
    z                         = [5;5e-3];                       % measurement noise (a noisier measurement is less trustworthy)
    w                         = 0;                              % process noise standard deviation
    obs                       = [truObs(1,1);truObs(2,1)];      % single [r,rhoDot] measurement
    num_iterations            = 1;
    [Xout,Pout,devMeas]       = ukf(h,ukfSigmaPoints,ukfCovar,w,z,obs,num_iterations,Wm,Wc);
%     [x_update,Post_P,devMeas]  = EnKF(h,ukfSigmaPoints,w,z,obs,num_iterations);
    [meanOut,covarOut,sigmaPointsOut,Wm,Wc] = prop_UT( mean(Xout,2), Pout, options, 9, dt,dt/dt);  % Propagate measurement update to next time step
    ukfMean                   = meanOut(end,:)';
    ukfCovar                  = covarOut(:,:,end);
    ukfSigmaPoints            = sigmaPointsOut(:,:,end);
    estUkfMeans(k,:)          = ukfMean;
    estUkfCovars(:,:,k)       = ukfCovar;
    storeObs(k,:)             = devMeas;    
    [A, C]                    = jacob(Xplus);
    [Xplus, Pplus]            = ekf(h,Xplus,Pplus,devMeas,Q,0,A,C);
    estEkfMeans(k,:)          = Xplus;
    estEkfCovars(:,:,k)       = Pplus;
end

% Plots for Estimation Accuracy
numMeas = linspace(1,N,N);
figure(); 
subplot(4,1,1); plot(numMeas,estUkfMeans(:,1),'b--'); hold on; plot(numMeas,true(1:N,1),'r-'); hold on; plot(numMeas,storeObs(:,1),'gx'); grid on; 
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2); plot(numMeas,estUkfMeans(:,2),'b--'); hold on; plot(numMeas,true(1:N,2),'r-'); grid on;
ylabel('$\theta$, m/sec','Interpreter','latex')
subplot(4,1,3); plot(numMeas,estUkfMeans(:,3),'b--'); hold on; plot(numMeas,true(1:N,3),'r-'); hold on; plot(numMeas,storeObs(:,2),'gx'); grid on;
ylabel('$\dot{r}$, m/s','Interpreter','latex')
subplot(4,1,4); plot(numMeas,estUkfMeans(:,4),'b--'); hold on; plot(numMeas,true(1:N,4),'r-'); grid on;
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('Number of Measurements','Interpreter','latex')
subplot(4,1,1); title('UKF Estimation Accuracy','Interpreter','latex')

figure(); 
subplot(4,1,1); plot(numMeas,estEkfMeans(:,1),'b--'); hold on; plot(numMeas,true(1:N,1),'r-'); hold on; plot(numMeas,storeObs(:,1),'gx'); grid on; 
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2); plot(numMeas,estEkfMeans(:,2),'b--'); hold on; plot(numMeas,true(1:N,2),'r-'); grid on;
ylabel('$\theta$, m/sec','Interpreter','latex')
subplot(4,1,3); plot(numMeas,estEkfMeans(:,3),'b--'); hold on; plot(numMeas,true(1:N,3),'r-'); hold on; plot(numMeas,storeObs(:,2),'gx'); grid on;
ylabel('$\dot{r}$, m/s','Interpreter','latex')
subplot(4,1,4); plot(numMeas,estEkfMeans(:,4),'b--'); hold on; plot(numMeas,true(1:N,4),'r-'); grid on;
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('Number of Measurements','Interpreter','latex')
subplot(4,1,1); title('EKF Estimation Accuracy','Interpreter','latex')

% Post-Fit Residuals
rResidUkf    = storeObs(:,1) - estUkfMeans(:,1);
rDotResidUkf = storeObs(:,2) - estUkfMeans(:,3);
rResidEkf    = storeObs(:,1) - estEkfMeans(:,1);
rDotResidEkf = storeObs(:,2) - estEkfMeans(:,3);

figure(); 
subplot(2,1,1); plot(numMeas,rResidUkf,'bx'); grid on;
ylabel('$r$, m','Interpreter','latex')
subplot(2,1,2); plot(numMeas,rDotResidUkf,'bx'); grid on;
ylabel('$\dot{r}$, m/s','Interpreter','latex')
xlabel('Number of Measurements','Interpreter','latex')
subplot(2,1,1); title('Ukf Observation Residuals','Interpreter','latex')

figure(); 
subplot(2,1,1); plot(numMeas,rResidEkf,'bx'); grid on;
ylabel('$r$, m','Interpreter','latex')
subplot(2,1,2); plot(numMeas,rDotResidEkf,'bx'); grid on;
ylabel('$\dot{r}$, m/s','Interpreter','latex')
xlabel('Number of Measurements','Interpreter','latex')
subplot(2,1,1); title('Ekf Observation Residuals','Interpreter','latex')

% Mean Square Error Plots
storeUkfMSE = zeros(N,2); 
storeEkfMSE = zeros(N,2); 
rResidSqUkf    = rResidUkf.^2;
rDotResidSqUkf = rDotResidUkf.^2;
rResidSqEkf    = rResidEkf.^2;
rDotResidSqEkf = rDotResidEkf.^2;
for i = 1:N
    rMSEUkf = sum(rResidSqUkf(1:i))/i;
    rDotMSEUkf = sum(rDotResidSqUkf(1:i))/i;
    storeUkfMSE(i,:) = [rMSEUkf;rDotMSEUkf];
    rMSEEkf = sum(rResidSqEkf(1:i))/i;
    rDotMSEEkf = sum(rDotResidSqEkf(1:i))/i;
    storeEkfMSE(i,:) = [rMSEEkf;rDotMSEEkf];
end

figure(); 
subplot(2,1,1); plot(numMeas,storeUkfMSE(:,1),'k-'); grid on; hold on
subplot(2,1,1); plot(numMeas,storeEkfMSE(:,1),'b-'); grid on; hold on
ylabel('$r$, m','Interpreter','latex')
subplot(2,1,2); plot(numMeas,storeUkfMSE(:,2),'k-'); grid on; hold on
subplot(2,1,2); plot(numMeas,storeEkfMSE(:,2),'b-'); grid on;
ylabel('$\dot{r}$, m/s','Interpreter','latex')
xlabel('Number of Measurements','Interpreter','latex')
subplot(2,1,1); title('Mean Squared Error','Interpreter','latex')

%% Propagate after update
% Propagate all sigma points 
[utMeansNew,utCovarsNew,utSigmaPointsNew,Wm,Wc] = prop_UT( ukfMean, ukfCovar, options, 9, revs*P-dt*N,dt);
n = length(utMeansNew); 
ekfMeansNew  = zeros(n,4);
ekfCovarsNew = zeros(4,4,n); 
for k = 1:n
    [A, C]           = jacob(Xplus);
    [Xplus, Pplus]   = ekf(h,Xplus,Pplus,y,Q,R,A,C);
    ekfMeansNew(k,:)    = Xplus;
    ekfCovarsNew(:,:,k) = Pplus;
end

% Append trajectory during estimation to propagated trajectory 
utMeansNew   = [estUkfMeans;utMeansNew];
ekfMeansNew  = [estEkfMeans;ekfMeansNew];

%% Store Propagated Sigmas
l = length(utCovarsNew);
storeUtCovNew = zeros(N+l,4);
storeEkfCovNew = zeros(N+1,4);
for i = 1:N
    storeUtCovNew(i,:) = sqrt(diag(estUkfCovars(:,:,i)));
    storeEkfCovNew(i,:) = sqrt(diag(estEkfCovars(:,:,i)));
end
ind = 1;
for i = N+1:N+l
    storeUtCovNew(i,:) = sqrt(diag(utCovarsNew(:,:,ind)));
    storeEkfCovNew(i,:) = sqrt(diag(ekfCovarsNew(:,:,ind)));    
    ind = ind + 1;
end

t  = linspace(dt,revs*P,revs*P/dt);
t2 = linspace(revs*P,revs*P*2,revs*P/dt);
t3 = linspace(dt,revs*P*2+dt,(revs*P*2+dt)/dt);
figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],storeUtCov(:,1)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeUtCovNew(:,1)/1e3,'b-'); hold on;
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t/3600],storeUtCov(:,3)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeUtCovNew(:,3)/1e3,'b-'); hold on;
ylabel('$\dot{r}$, m/sec','Interpreter','latex')
subplot(4,1,3)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],mod(storeUtCov(:,2)*180/pi,360),'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],mod(storeUtCovNew(:,2)*180/pi,360),'b-'); hold on;
ylabel('$\theta$, deg','Interpreter','latex')
subplot(4,1,4)
plot([0 t/3600],storeUtCov(:,4)*180/pi,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeUtCovNew(:,4)*180/pi,'b-'); hold on;
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('UT/UKF $1\sigma$ Uncertainty','Interpreter','latex');
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
title('UT/UKF Trajectory','Interpreter','latex');
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
ylabel('$\Delta\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('UT/UKF Trajectory Error','Interpreter','latex');
set(gca,'gridlinestyle','--')

figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],storeEkfCov(:,1)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeEkfCovNew(:,1)/1e3,'b-'); hold on;
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t/3600],storeEkfCov(:,3)/1e3,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeEkfCovNew(:,3)/1e3,'b-'); hold on;
ylabel('$\dot{r}$, m/sec','Interpreter','latex')
subplot(4,1,3)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],mod(storeEkfCov(:,2)*180/pi,360),'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],mod(storeEkfCovNew(:,2)*180/pi,360),'b-'); hold on;
ylabel('$\theta$, deg','Interpreter','latex')
subplot(4,1,4)
plot([0 t/3600],storeEkfCov(:,4)*180/pi,'b-'); grid on; hold on;
plot([revs*P/3600 t2/3600],storeEkfCovNew(:,4)*180/pi,'b-'); hold on;
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('EKF $1\sigma$ Uncertainty','Interpreter','latex');
set(gca,'gridlinestyle','--')

figure();
subplot(4,1,1)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],ekfMeans(:,1)/1e3,'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],ekfMeansNew(:,1)/1e3,'b-'); hold on;
plot([0 t3/3600],fTot(:,1)/1e3, 'r--');
ylabel('$r$, m','Interpreter','latex')
subplot(4,1,2)
plot([0 t/3600],ekfMeans(:,3)/1e3,'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],ekfMeansNew(:,3)/1e3,'b-'); hold on;
plot([0 t3/3600],fTot(:,3)/1e3, 'r--');
ylabel('$\dot{r}$, m/sec','Interpreter','latex')
subplot(4,1,3)
set(0,'DefaultAxesFontName', 'Arial'); 
plot([0 t/3600],mod(ekfMeans(:,2)*180/pi,360),'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],mod(ekfMeansNew(:,2)*180/pi,360),'b-'); hold on;
plot([0 t3/3600],mod(fTot(:,2)*180/pi,360), 'r--');
ylabel('$\theta$, deg','Interpreter','latex')
subplot(4,1,4)
plot([0 t/3600],ekfMeans(:,4)*180/pi,'g-'); grid on; hold on;
plot([revs*P/3600 t2/3600],ekfMeansNew(:,4)*180/pi,'b-'); hold on;
plot([0 t3/3600],fTot(:,4)*180/pi, 'r--');
ylabel('$\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('EKF Trajectory','Interpreter','latex');
legend('Traj Prior to Update','Traj After Update','True Traj');
set(gca,'gridlinestyle','--')

% Residuals 
PredTraj = [ekfMeans; ekfMeansNew]; % Predicted Trajectory
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
ylabel('$\Delta\dot{\theta}$, deg/sec','Interpreter','latex')
xlabel('\fontname{Times New Roman} Length of Propagation, hr');
subplot(4,1,1)
title('EKF Trajectory Error','Interpreter','latex');
set(gca,'gridlinestyle','--')