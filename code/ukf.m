function [x_update] = ukf(h,x_initial,w,z,obs,num_iterations,Wm,Wc) 

[~,num_members]     = size(x_initial);
p1                  = size(x_initial,1);
m1                  = size(obs,1);
x_estimate          = x_initial;

Zcov = eye(m1);     % Create measurement noise covariance matrix

for j = 1:m1
  Zcov(j,j) = z(j)^2;
end

for i = 1:num_iterations  
 
   for j = 1:num_members
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
   x_estimatebar     = mean(x_estimate,2);   % Replaced with propagated sigma points
   ybar              = mean(y,2);            % Mean of measurement return
   y_forecastbar     = mean(y_forecast,2);   % Mean of predicted measurement

   for j = 1:p1
     Ex(j,:) = [x_estimate(j,:) - x_estimatebar(j)];  % State covariance
   end

   for j = 1:m1
     Ey(j,:) = [y_forecast(j,:) - y_forecastbar(j)];  % Measurement covariance
   end
   
%%
   x_estimatebar_new = zeros(p1,1); 
   y_forecastbar_new = zeros(m1,1);
   Pxy_new = zeros(p1,2);
   Pyy_new = zeros(m1,2);
 
   for m = 2:num_members
       x_estimatebar_new = x_estimatebar_new + Wm(m,1)*x_estimate(:,m);
       y_forecastbar_new = y_forecastbar_new + Wm(m,1)*y(:,m);
       Pxy_new = Pxy_new + Wc(m,1)*Ex*Ey';
       Pyy_new = Pyy_new + Wc(m,1)*Ey*Ey';
   end   
   
   K_new = Pxy_new*inv(Pyy_new);
   x_update_new = x_estimate + K_new*(y - y_forecast);
   
%%   
   Pxy            = Ex*Ey'/(num_members - 1);         % Cross-correlation
   Pyy            = Ey*Ey'/(num_members - 1) + Zcov;  % Innovation  
   K              = Pxy*inv(Pyy);                          % Kalman Gain
   x_update       = x_estimate + K*(y - y_forecast);  % State Update  
   if i == num_iterations
       x_estimatebar = mean(x_estimate,2);
   end
end
