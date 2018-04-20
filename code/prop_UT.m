function [intMeans,intCovars,statesOut,Wm,Wc] = prop_UT( s0, Cov,...
    Chaser, alpha, alpha_t, t_now, t_seg)
%
% [outMean,outCovar] = prop_UT( s0, Cov, Chaser, alpha, alpha_t, t_now, t_seg )
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
%  s0       - nx1 vector for the mean of the input Gaussian
%  Cov      - Structure containing covariance information and options for
%             the unscented transform
%               Cov.P0           - nxn covariance matrix
%               Cov.dt           - step size for integration
%               Cov.utsquareroot - (optional) If 1, then assume the
%                                           input covariance is the
%                                           Cholesky decomposition and use
%                                           the square-root UT
%               Cov.alpha        - double, alpha parameter for the UT
%               Cov.beta         - double, beta parameter for the UT
%  Chaser   - Structure containing information about the Chaser s/c
%  alpha    - Vector of control inputs for the segment
%  alpha_t  - Times associated with alpha, relative to mission start
%  t_now    - Epoch of the beginning of the propagation, relative to
%             mission start
%  t_seg    - Epoch of the end of the propagation, relative to mission
%             start
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
%   twobodyPolar()
%
% Modification History:
%
%   22jun2016     Brandon A. Jones      original version, header added
%
% ---------------------------------------------------------------------
% Copyright 2016, Brandon A. Jones, The University of Texas at Austin

%  Determine if we were provided the square root of the covariance or the
%  full matrix
if isfield( Cov, 'utsquareroot' )
    useSquareRoot = Cov.utsquareroot;
else
    useSquareRoot = 0;
end

%  If we were provided the full covariance, get its Cholesky decomposition.
if useSquareRoot == 1
    cholCov = Cov.P0;
else
    cholCov = chol( Cov.P0, 'lower' );
end

%  Get the sigma points
[Chi,Wm,Wc] = getUTSigmaPts_SR( s0, cholCov, Cov );

%  Map the sigma points through the input funcion
statesOut = zeros(size(Chi,1),size(Chi,2), length([t_now:Cov.dt:t_seg, t_seg]));
for i = 1:size(Chi,2)
    f = twobodyPolar(Chi(:,i), [t_now, t_seg], Cov.dt, Chaser,...
        Cov, alpha, alpha_t);
    statesOut(:,i,:) = f.';
end

n = size(statesOut, 3);
intMeans = zeros(n, length(s0));
intCovars = zeros(length(s0), length(s0), n);
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