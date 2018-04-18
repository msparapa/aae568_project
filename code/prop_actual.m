function Actual = prop_actual(Actual, Chaser, Cov, alpha, alpha_t, t0, tf)
% PROP_ACTUAL Propagate the true dynamics
%
%   Actual = prop_actual(Actual, Chaser, Cov, alpha, alpha_t, t0, tf)
%
%   Inputs:
%
%       - Actual: Structure the stores the true Chaser state
%       - Chaser: Structure that contains information about the Chaser s/c
%       - Cov: Structure containing information about the Covariance
%       - alpha: Optimal control time-history
%       - alpha_t: Time associated with alpha
%       - t0: Epoch of the beginning of the segment, relative to mission
%       start
%       - tf: Final time on the segment, relative to mission start
%
%   Outputs:
%
%       - Actual: Updated structure with the true Chaser state after the
%       specified propagation
%
% Collin York
% AAE 508
% Final Project

% Generate a signal to perturb the true dynamics
w = randn(4,length(alpha_t));

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
odes = @(t,X)prop_actual_odes(t, X, Chaser, Cov, alpha, alpha_t, w);
[t,X] = ode113(odes, [t0, tf], Actual.X, options);

% Update the actual state
Actual.X = X(end,:)';
Actual.X_history{end+1} = X;
Actual.t_history{end+1} = t;
Actual.alpha_history{end+1} = alpha;
Actual.alpha_t_history{end+1} = alpha_t;
return