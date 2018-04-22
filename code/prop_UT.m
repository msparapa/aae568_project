function [intMeans,intCovars,statesOut,Wm,Wc] = prop_UT( s0, P0, options, varargin, t, dt)
%
% [outMean,outCovar] = prop_UT( mean, covar, fcn, options, varargin )
% ---------------------------------------------------------------------
%
% Description:
%
%  Perform an unscented transformation of the input PDF (parameterized by
%  the mean and covar) to an output PDF.  The function for the map of
%  single points is provided via fcn.
%
% Inputs:
%
%  mean    - nx1 vector for the mean of the input Gaussian
%  covar   - nxn matrix for the covariance information for the input
%               Gaussian
%  fcn     - Black-box function to use to map the sigma points to the
%               output domain
%  options - structure containing the options for the execuation of the UT.
%               options.utsquareroot - (optional) If 1, then assume the
%                                           input covariance is the
%                                           Cholesky decomposition and use
%                                           the square-root UT
%               options.alpha        - double, alpha parameter for the UT
%               options.beta         - double, beta parameter for the UT
%  varargin - Other arguments passed to the input fcn
%
% Outputs:
%
%  outMean  - mx1 vector for the mean of the output Gaussian
%  outCovar - mxm matrix for the output covariance information.
%
% Assumptions/References:
%
%  If options.utsquareroot is not included or is 0, then the input is
%  assumed to be a full covariance matrix.  If it is included and equals 1,
%  then the input is assumed to be the lower-triangular Cholesky
%  decomposition of the matrix.  The output covariance information will
%  match that of the input.
%
%  Rudolph Van der Merwe and Eric A. Wan. "The square-root unscented Kalman
%     filter for state and parameter-estimation." In 2001 IEEE
%     International Conference on acoustics, speech, and signal processing,
%     Volume 6, 3461--3464, Salt Lake City, Utah, May 7-11 2001.
%
% Dependencies:
%
%   getUTSigmaPts_SR()
%
% Modification History:
%
%   22jun2016     Brandon A. Jones      original version, header added
%
% ---------------------------------------------------------------------
% Copyright 2016, Brandon A. Jones, The University of Texas at Austin

%  Determine if we were provided the square root of the covariance or the
%  full matrix
if isfield( options, 'utsquareroot' )
    useSquareRoot = options.utsquareroot;
else
    useSquareRoot = 0;
end

%  If we were provided the full covariance, get its Cholesky decomposition.
if useSquareRoot == 1
    cholCov = P0;
else
    cholCov = chol( P0, 'lower' );
end

%  Get the sigma points
[Chi,Wm,Wc] = getUTSigmaPts_SR( s0, cholCov, options );

%  Map the sigma points through the input funcion
statesOut = zeros(size(Chi,1),size(Chi,2),t/dt+1);
for i = 1:size(Chi,2)
    f = twobodyPolar(Chi(:,i),[0 t],dt);
    statesOut(:,i,:) = transpose(f);
end

% Store the means and covars at each propagation time step
n = t/dt+1;
intMeans = zeros(n,4);
intCovars = zeros(4,4,n);
for h = 1:n
    
    stateOut = statesOut(:,:,h);
    nDimsOut = size( stateOut, 1 );
    
    %  Get the mean of the output PDF
    outMean = Wm(1)*stateOut(:,1) + Wm(2)*sum(stateOut(:,2:end),2);
    
    %  Compute the square-root of the covariance for the output.  This uses the
    %  algorithm in Van der Merwe and Wan cited above.
    diffs = stateOut - repmat( outMean, 1, size(stateOut,2) );
    R  = triu(qr( ( repmat(sqrt(Wc(2:end))',nDimsOut,1).*diffs(:,2:end) )' ));
    cholCovOut = cholupdate( R(1:nDimsOut,:), Wc(1).*diffs(:,1) )';
    
    %  Make sure we provide output in the same format as the input.
    if useSquareRoot == 1
        outCovar = cholCovOut;
    else
        outCovar = cholCovOut*cholCovOut';
    end
    
    intMeans(h,:) = outMean;
    intCovars(:,:,h) = outCovar;
end

end