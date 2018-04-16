function [Chi,Wm,Wc] = getUTSigmaPts_SR( mean, cholCovar, options )
%
% [Chi,Wm,Wc] = getUTSigmaPts_SR( mean, cholCovar, options )
%
% ---------------------------------------------------------------------
% 
% Description:
%
%  Get the sigma points for the input Gaussian (mean and covariance) via 
%  the square-root unscented transformation 
% 
% Inputs:
%
%  mean      - nx1 vector for the mean of the input Gaussian
%  cholCovar - nxn matrix for the covariance information for the input 
%                 Gaussian
%  options - structure containing the options for the execuation of the UT.
%               options.alpha        - double, alpha parameter for the UT
%               options.beta         - double, beta parameter for the UT
% 
% Outputs:
% 
%  Chi - nxm, Matrix of sigma points where m = 2*n + 1
%  Wm  - mx1, Vector of weights for computing the mean
%  Wc  - mx1, Vector of weights for computing the covariance
%
% Assumptions/References:
%
%  Assumes that the input cholCovar is the lower-triangular form of the
%  Cholesky decomposition of the covariance matrix.
%
% Dependencies:
%
%
% Modification History:
% 
%   17jun2016     Brandon A. Jones      original version, header added
%
% ---------------------------------------------------------------------
% Copyright 2016, Brandon A. Jones, The University of Texas at Austin


nDims = length(mean);

lambda = options.alpha*options.alpha*3-nDims;
Wm0    = lambda/(nDims+lambda);
Wc0    = Wm0 + (1 - options.alpha*options.alpha + options.beta);
WmI    = ones(2*nDims,1)./(2*(nDims+lambda));
Wm     = [Wm0; WmI];
Wc     = [Wc0; WmI];

gamma = sqrt(nDims + lambda);

sqrtCov_t_Gam      = gamma.*cholCovar;
Chi                = zeros(nDims,nDims*2+1);
Chi(:,1)           = mean;
Chi(:,2:nDims+1)   = repmat(mean,1,nDims)+sqrtCov_t_Gam;
Chi(:,nDims+2:end) = repmat(mean,1,nDims)-sqrtCov_t_Gam;

end
