% Collin York
% Initialize Indirect Optimization
clear all; clc; close all;
colors = lines(5);

% Plotting Options
plot_opt.i = 0;             % last plot figure number
plot_opt.indirect = false;   % whether or not to plot indirect optimization results
plot_opt.actual = true;     % whether or not to plot comparison of Nav, Actual, and Target
plot_opt.nav = false;        % whether or not to plot navigation results

% Options for simulation

% Which optimization method to use
%
%   - indirect
%   - direct (TODO)
sim_opt.optim = 'indirect';

% Which estimation method to use
%
%   - ekf   Extended Kalman Filter
%   - ut    Unscented Transform + Unscented Kalman Filter
sim_opt.estim = 'ekf';

% Tolerance to check if state has reached final value
sim_opt.stateTol = 1e-6;

sim_opt.maxCount = 50;
%% Define Dimensional Initial Conditions
% Dimensional initial conditions
T = 0.005605;             % Thrust, kN; TODO - DEFINE THIS
Isp = 3000;             % Specific impulse, sec
g0 = 9.80665/1000;      % Mean Earth gravity, km/s^2
M0 = 500;               % Initial S/C Mass, kg


r0 = 42164;             % Initial orbital radius of chaser, km
theta0 = 0;             % Initial longigutde of chaser, rad
rdot0 = 0;              % Initial radial velocity of chaser, km/s

charL = r0;             % characteristic length, km; Initial Orbit Radius
charM = M0;             % characteristic mass, kg; S/C mass
muEarth = 3.986e5;      % Earth gravity parameter, km^3/s^2
charT = sqrt(charL^3/muEarth); % characteristic time, sec; sqrt(r0^3/mu)

% Compute thetadot for a circular orbit; rad/s
thetadot0 = sqrt(muEarth/r0^3);

range_err = (10/1000)/charL;            % 10 m -> nondim
rr_err = (1/1000/1000)*charT/charL;     % 1 mm/s -> nondim
            
% Chaser Nondimensional Parameters
% These will be fed in by main code
Chaser.T = T/charM/charL*charT^2;        % Non-dim thrust
Chaser.mdot = T/(Isp*g0)/charM*charT;    % Non-dim mdot
Chaser.m0 = M0/charM;                    % Non-dim mass; starts at 1; not state variable
Chaser.ts_opt = 1800/charT;               % Non-dim time-of-flight for one "leg" between observations 

% Filtered Nav Initial State for Optimizer at t0
Nav.r = 1;
Nav.theta = 0;
Nav.rdot = 0;
Nav.thetadot = 1;
Nav.P = 1e-5*eye(4);                    % Initial Covariance to test EKF
Nav.X_history = {};
Nav.t_history = {};

% Actual Initial State for Optimizer at t0
Actual.X = [1; 0; 0; 1] + sqrt(Nav.P)*randn(4,1);    %Current state [r, theta, rdot, thetadot]
Actual.X_history = {};      % Each cell holds the state history for one segment
Actual.t_history = {};      % Each cell holds the time associated with the state history
Actual.alpha_history = {};
Actual.alpha_t_history = {};

% Target Actual State at t = 0 (note: not same as updated t0)
Target.r0 = 1.01;
Target.theta0 = 0.21;
Target.rdot0 = 0;
Target.thetadot0 = sqrt(1/Target.r0^3);

% Data from Montenbruk & Gill: 
%   - At LEO (e.g. Iridium at r = 63780 + 780), most significant accel pert
%   is J_2,0 w/ magnitude 1e-5 km/s^2
%   - At GEO (r = 42000 km), most signficant accel pert are J_2,0 and Lunar 
%   gravity w/ magnitude 1e-8 km/s^2
% A = [zeros(2,4); zeros(2,2), eye(2)];
% Cov.R = A*1e-7*(charT^2/charL);   % Acceleration Process Noise (xdot = f(x,u,t) + C*w)
Cov.R = zeros(4);       % No noise in EOMs

switch(sim_opt.estim)
    case 'ekf'
        % Noise Covariance
        Cov.Z = [range_err, 0; 0, rr_err];  % Measurement noise (y = Hx + Gz)
    case 'ut'
        % Some arbitrary covariance matrix
        P0 = rand(4);
        P0 = P0*P0.' * 1e-5;        % Use small values to avoid larger errors that crash the Chaser into Earth
        Cov.P0 = P0;
        Cov.alpha = 1;
        Cov.beta = 2.0;
        Cov.dt = 0.05;    % Propagate step-size (nondimensional time)
end

t_now = 0;          % t_now is the current reoptimization time (not always 0)
tf_rel_guess = pi;  % pi is good guess when t_now = 0, will update as tf-t_now

% helpful IC, replaced by lambda_f when re-optimizing
lambda0_guess = [22; -7; 20; -5];

if(plot_opt.actual)
    plot_opt.i = plot_opt.i + 1;
    h_traj = figure(plot_opt.i);
    hold on;
    ax_traj = gca;
    title('Trajectory');
    xlabel('x, nondim');
    ylabel('y, nondim');
    set(ax_traj, 'fontsize', 12);
    
    plot_opt.i = plot_opt.i + 1;
    h_state = figure(plot_opt.i);
    subplot(4,1,1); hold on;
    ax_r = gca;
    title('State Histories');
    ylabel('r, nd');
    set(ax_r, 'fontsize', 12);
    
    subplot(4,1,2); hold on;
    ax_theta = gca;
    ylabel('\theta, rad');
    set(ax_theta, 'fontsize', 12);
    
    subplot(4,1,3); hold on;
    ax_rdot = gca;
    ylabel('$$\dot{r}$$, nd', 'interpreter', 'latex');
    set(ax_rdot, 'fontsize', 12);
    
    subplot(4,1,4); hold on;
    ax_thetadot = gca;
    ylabel('$$\dot{\theta}$$, nd', 'interpreter', 'latex');
    xlabel('Time, nd');
    set(ax_thetadot, 'fontsize', 12);
end
% Begin Mission Loop
gameover = false;
count = 0;

while(~gameover && count < sim_opt.maxCount)

    fprintf('******************\nIteration %02d\n******************\n',...
        count);
    %% Using current state est., compute optimal traj. to reach target
    %   
    %   Optimization method should spit out the FULL optimal control
    %   from the current state to the final rendezvous
    %
    %   We won't necessarily use the entire optimal solution, just a
    %   section of it that corresponds to our flight time between
    %   observations
    switch(lower(sim_opt.optim))
        case 'indirect'
            % Indirect Optimization - Collin
            % Outputs: 
            % - alpha: Control history for the entire optimal trajectory
            % - alpha_t: times associated with alpha
            % - tf: final time on the optimal trajectory, relative to
            % mission start
            % - t_seg: time at the end of "the segment", relative to
            % mission start
            % - lambda_seg: costate values at the end of "the segment", 
            %   i.e., part way through the optimal trajectory.
            if(count > 0)
                lambda0_guess = lambda_seg;         % Update
                t_now = t_seg;                      % Update starting time
                tf_rel_guess = tf - t_now ;         % Update TOF guess
                if tf_rel_guess <= 0
                    tf_rel_guess = 1.1*Chaser.ts_opt;
                end
            end
            [alpha, alpha_t, tf, t_seg, lambda_seg, plot_opt] = indirect_fcn(Chaser,...
                Target, Nav, t_now, plot_opt, lambda0_guess, tf_rel_guess);
            if t_now == 0
                tf_initial = tf;
            end
    end
    
    
    %% Propagate for one "step" until next observation
    
    % Propagate the true state with interpolated process noise
    Actual = prop_actual(Actual, Chaser, Cov, alpha, alpha_t, t_now, t_seg);
    
    switch(lower(sim_opt.estim))
        case 'ekf'
            % EKF function to propagate State covariance in P in continuous
            % time with acceleration process covariance
            Nav = prop_EKF(Nav, Chaser, Cov, alpha, alpha_t, t_now, t_seg);
        case 'ut'
            % Construct state for UT
            s0 = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
            
            % Run unscented transform
            [intMeans, intCovars, statesOut, Wm, Wc] = prop_UT(s0, Cov,...
                Chaser, alpha, alpha_t, t_now, t_seg);
            
            Nav.X_history{end+1} = intMeans;
            Nav.t_history{end+1} = [t_now:Cov.dt:t_seg, t_seg];
            
%             if(count == 28)
%                 keyboard;
%             end
            % Store propagated Sigmas
            storeCov = zeros(size(intCovars,3),4);
            for i = 1:size(intCovars,3)
                storeCov(i,:) = sqrt(diag(intCovars(1:4,1:4,i)));
            end
    end
    
    %% Make observation and update estimate
    %   
    %   * Update state estimate
    %   * Update covariance
    
    switch(lower(sim_opt.estim))
        case 'ekf'
            H = [1, 0, 0, 0; 0, 0, 1, 0];  % y = H*x + G*z; r and rdot measurements
            G = eye(2);         % y = H*x + G*z
            L_k = (Nav.P * H.') / (H*Nav.P*H.' + G*Cov.Z*G.');
            y = H*Actual.X + [sqrt(Cov.Z(1,1))*randn(1,1); sqrt(Cov.Z(2,2))*randn(1,1)];  % measurement + noise; z = [sig1*z1; sig2*z2]
            % apriori nav estimate
            nav_pre = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
            % nav estimate post observation
            nav_post = nav_pre + L_k*(y - H*nav_pre);   % update
            
            % Update Nav state after observation
            Nav.r = nav_post(1);
            Nav.theta = nav_post(2);
            Nav.rdot = nav_post(3);
            Nav.thetadot = nav_post(4);
            Nav.P = (eye(4) - L_k*H)*Nav.P;     % covariance update
        case 'ut'
            % Update propagated particles using a UKF
            truObs = [Actual.X(1); Actual.X(3)];
            
            z = [range_err; rr_err];                % Measurement Noise
            h = @(j) updatePolarMeasurement(j);                                                                              
            x_initial = statesOut(:,:,end);         % Sigma Points                   
            w = 0;                                  % Process Noise Standard Deviation                              
            obs = [truObs(1,1); truObs(2,1)];       % Single [r,rhoDot] Measurement                     
            num_iterations = 1;
            [x_update, postUpdateCov] = ukf(h, x_initial, Nav.P, w, z,...
                obs, num_iterations, Wm, Wc);
            
            % Update nav state
            x_update = mean(x_update, 2);
            Nav.r = x_update(1);
            Nav.theta = x_update(2);
            Nav.rdot = x_update(3);
            Nav.thetadot = x_update(4);
            Nav.P = postUpdateCov;
    end
    
    if(plot_opt.actual)
        % Plot Actual path (where we actually are)
        x = Actual.X_history{end}(:,1).*cos(Actual.X_history{end}(:,2));
        y = Actual.X_history{end}(:,1).*sin(Actual.X_history{end}(:,2));
        h_true = plot(ax_traj, x, y, 'linewidth', 2, 'color', colors(1,:));
        plot(ax_traj, x(end), y(end), 'ws', 'markerfacecolor', colors(1,:));
        
        % Plot Nav solution (where we think we are)
        x = Nav.X_history{end}(:,1).*cos(Nav.X_history{end}(:,2));
        y = Nav.X_history{end}(:,1).*sin(Nav.X_history{end}(:,2));
        h_nav = plot(ax_traj, x, y, 'linewidth', 2, 'color', colors(2,:));
        plot(ax_traj, x(end), y(end), 'ws', 'markerfacecolor', colors(2,:));
        
        if(count == 0)
            grid(ax_traj, 'on');
            axis(ax_traj, 'equal');
            legend([h_true, h_nav], 'Actual', 'Nav');
        end
        
        plot(ax_r, Actual.t_history{end}, Actual.X_history{end}(:,1),...
            'linewidth', 2, 'color', colors(1,:));
        plot(ax_theta, Actual.t_history{end}, Actual.X_history{end}(:,2),...
            'linewidth', 2, 'color', colors(1,:));
        plot(ax_rdot, Actual.t_history{end}, Actual.X_history{end}(:,3),...
            'linewidth', 2, 'color', colors(1,:));
        plot(ax_thetadot, Actual.t_history{end}, Actual.X_history{end}(:,4),...
            'linewidth', 2, 'color', colors(1,:));
        
        plot(ax_r, Nav.t_history{end}, Nav.X_history{end}(:,1),...
            'linewidth', 2, 'color', colors(2,:));
        plot(ax_theta, Nav.t_history{end}, Nav.X_history{end}(:,2),...
            'linewidth', 2, 'color', colors(2,:));
        plot(ax_rdot, Nav.t_history{end}, Nav.X_history{end}(:,3),...
            'linewidth', 2, 'color', colors(2,:));
        plot(ax_thetadot, Nav.t_history{end}, Nav.X_history{end}(:,4),...
            'linewidth', 2, 'color', colors(2,:));
        
        plot(ax_r, [t_now, t_seg], Target.r0*[1,1], 'k--', 'linewidth', 2);
        plot(ax_theta, [t_now, t_seg], Target.theta0*[1,1] + Target.thetadot0*[t_now, t_seg], 'k--', 'linewidth', 2);
        plot(ax_rdot, [t_now, t_seg], Target.rdot0*[1,1], 'k--', 'linewidth', 2);
        plot(ax_thetadot, [t_now, t_seg], Target.thetadot0*[1,1], 'k--', 'linewidth', 2);
    end

    %% Evaluate status
    %   Has spacecraft reached the target?
    X_targ = getTargetState(Target, t_seg);
    fprintf('Target State = [%f, %f, %f, %f]\n', X_targ);
    
    X_chaser = [Nav.r; Nav.theta; Nav.rdot; Nav.thetadot];
    fprintf('Chaser State = [%f, %f, %f, %f]\n', X_chaser);
    
    fprintf('Covariance Matrix:\n'); disp(Nav.P);
    
    if(max(abs(eig(Nav.P))) > 0)
        [bInPos, bInVel] = checkErrEllipses(Nav, Actual, Target, t_seg, plot_opt);
        gameover = bInPos && bInVel;
    else
        gameover = norm(X_targ - X_chaser) < sim_opt.stateTol ;
    end
    % TODO -  Need to add check on covariance size??
    
    
%     keyboard;
    
    if(~gameover)
        % Update Chaser state
        Chaser.m0 = Chaser.m0 - Chaser.mdot*(t_seg - t_now);
        
        if(Chaser.m0 < 0.8)
            warning('Chaser mass = %f is unrealistic\n', Chaser.m0);
        end
    end
    
    count = count + 1;
end

if(~gameover)
    warning('Process did not converge');
end

% Plot the final set of ellipses even if we didn't plot them along
% the way
if(~plot_opt.nav)
    plot_opt.nav = true;
    checkErrEllipses(Nav, Actual, Target, t_seg, plot_opt);
end

tf_delta = t_seg - tf_intial; % Actual final time - Original Opt Final Time
fprintf('Time of Flight was %f above original optimal estimate\n', tf_delta);