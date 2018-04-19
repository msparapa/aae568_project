function [x_update,P_update] = ukf(h,x_initial,w,z,obs,num_iterations,Wm,Wc) 

[~,sigma_points]    = size(x_initial);
p1                  = size(x_initial,1);
m1                  = size(obs,1);
x_estimate          = x_initial;

Zcov = eye(m1);     % Create measurement noise covariance matrix

for j = 1:m1
  Zcov(j,j) = z(j)^2;
end

for i = 1:num_iterations  
 
   for j = 1:sigma_points
     % Create noise 
     W(:,j)          = w.*randn(p1,1);       % State Process Noise, not using here and will add random acceleration in 'accel_polar.m' instead        
     Z(:,j)          = z.*randn(m1,1);       % Measurement Noise           
     
     % Forecast Step
     forecastState   = x_estimate(:,j) + W;  % Replaced with already propagated sigma points                          
     meas            = obs + Z(:,j);         % 2-dimensional [r,rDot] vector                   
     forecastMeas    = h(forecastState)';    % Mapping from measurement space to dynamics space (which is the same in our case)        
      
     % Store estimates and measurements
     y(:,j)          = meas;                 
     y_forecast(:,j) = forecastMeas; 
     
   end
      
   % Analysis Step
   x_estimatebar     = x_estimate*Wm;   % Replaced with propagated sigma points
   ybar              = mean(y,2);       % Mean of measurement return
   y_forecastbar     = y_forecast*Wm;   % Mean of predicted measurement

   for j = 1:p1
     Ex(j,:) = [x_estimate(j,:) - x_estimatebar(j)];  % State covariance
   end

   for j = 1:m1
     Ey(j,:) = [y_forecast(j,:) - y_forecastbar(j)];  % Measurement covariance
   end
   
   Pxy            = Ex*diag(Wc)*Ey';                  % Cross-correlation
   Pyy            = Ey*diag(Wc)*Ey'+ Zcov;            % Innovation
   K              = Pxy/Pyy;                          % Kalman Gain
   x_update       = x_estimate + K*(y - y_forecast);  % State Update
   P_update       = cov(x_initial') - K*Pxy';         % Covariance Update
end
