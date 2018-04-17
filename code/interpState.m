function x = interpState(x_disc, t_disc, t)
% INTERPSTATE Interpolate to find a state between two discrete points
%
%   x = interpState(x_disc, t_disc, t) interpolates between points in the
%   discretized state vector, x_disc, to find the state at time t
%
%   Inputs:
%
%       x_disc      -   Array of state vectors; each column is a state
%       t_disc      -   Times associated with each state vector (column) in
%                       x_disc
%       t           -   The time at which to find the interpolated state
%
%   Author: Andrew Cox
%   Version: April 16, 2018

    
    ix_match = find(t_disc == t, 1);
    if(~isempty(ix_match))
        x = x_disc(:,ix_match);
    else
        ix_prev = find(t_disc < t, 1, 'last');
        ix_next = find(t_disc > t, 1, 'first');
        
        if(isempty(ix_prev))
            error('No time point found before t = %f', t);
        end
        if(isempty(ix_next))
            error('No time point found after t = %f', t);
        end
        
        x = (t - t_disc(ix_prev))/(t_disc(ix_next) - t_disc(ix_prev)) * ...
            (x_disc(:,ix_next) - x_disc(:,ix_prev)) + x_disc(:,ix_prev);
    end
end