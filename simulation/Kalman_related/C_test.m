%% FPGA side (to be done in verilog)
% This will run on the FPGA. The inputs it will need are both the measured
% signal and K, which should be calculated in C on the CPU side (see
% below).

% CODE SKETCH:

%% 1: DEFINE VARIABLES
% state transition model
A =  [1               dt; ...
     -dt*trap.w_0^2   1 - dt*trap.gamma/trap.m];
% observation model
H =  [1  0];
% control input model (assuming it is Lorentz force)
B =  [0  0;...
      0  q_e/trap.m];


%% 2: INFINITE LOOP
while(1)
    % Signal is generated in particle:
    % x(:, n)
    
    % But we only detect position
    % z(n) = x(1, n) + measurement noise
    
    % Now we can do:
    
    % AHEAD OF MEASUREMENT
    % Prediction (predict state estimate)
    x_kf_predict = A*x_kf(:, n-1) + B*u(:, n-1); % 
    
    % AFTER MEASUREMENT 
    % Innovation 
    y = z(n) - H*x_kf_predict; % innovation
    % Upadate estimate
    x_kf(:, n) = x_kf_predict + K*y;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CPU side (to be done in C)
% K (aka Kalman gain) iteratively converges to constant matrix: 
% this is the only important value in the CPU. 
% Once this has converged, pass the value to the FPGA

% CODE SKETCH:

%% 1: DEFINE VARIABLES

% state transition model
A =  [1               dt; ...
     -dt*trap.w_0^2   1 - dt*trap.gamma/trap.m];
% observation model
H =  [1  0];
% process noise covariance (i.e., thermal bath noise)
R1 =  [0  0;...
      0  dt*(trap.sigma/trap.m)^2];
% observation noise covariance (i.e., measurement noise (ballanced + redpitaya))
R2 =  sigma_obs^2;
% Estimate covariance
P = [0 0;...
     0 0];

%% 2: INFINITE LOOP
while(1)
    P_predict    = A*P*A' + R1; % predict estimate covariance
    S = R2 + H*P_predict*H'; % innovation covariance
    K = P_predict*H'/S; % Optimal Kalman gain
    P = (eye(2) - K*H)*P_predict; % update estimate covariance
end
    
    
    
    
    
    
    
    