function [x_update, x_update_vec] = ukf(h,x_initial,w,z,obs,num_iterations) 

[~,num_members]     = size(x_initial);
p1                  = size(x_initial,1);
m1                  = size(obs,1);
xvec                = [];
x_update_vec        = [];
yvec                = [];
x_estimate          = x_initial;

Zcov = eye(m1);     % create measurement noise covariance matrix

for j = 1:m1
  Zcov(j,j) = z(j)^2;
end

for i = 1:num_iterations  
 
   for j = 1:num_members
     % Create noise 
     W(:,j)          = w.*randn(p1,1);                
     Z(:,j)          = z.*randn(m1,1);                 
     
     % Forecast Step
     forecastState   = x_estimate(:,j);                             
     meas            = obs + Z(:,j);                           
     forecastMeas    = h(forecastState)';              
      
     % Store estimates and measurements
     y(:,j)          = meas;
     y_forecast(:,j) = forecastMeas; 
     
   end
   
   % Analysis Step
   x_estimatebar     = mean(x_estimate,2);    
   ybar              = mean(y,2);          
   y_forecastbar     = mean(y_forecast,2);  

   for j = 1:p1
     Ex(j,:) = [x_estimate(j,:) - x_estimatebar(j)];  
   end

   for j = 1:m1
     Ey(j,:) = [y_forecast(j,:) - y_forecastbar(j)];  
   end

   Pxy            = Ex*Ey'/(num_members - 1);
   Pyy            = Ey*Ey'/(num_members - 1) + Zcov;   
   K              = Pxy*inv(Pyy);
   x_update       = x_estimate + K*(y - y_forecast);
   x_update_vec   = [x_update_vec x_estimatebar];        
   yvec           = [yvec ybar];
   if i == num_iterations
       x_estimatebar = mean(x_estimate,2);
   end
end
