%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROUTINE DESCRIPTION:
%
% Simulate the signal from an optical tweezer, add noise to simulate
% detection noise, and reconstruct original signal with a Kalman filter
% (KF). Then do LQR feedback. Plot results in real time.
%
% The simulation is performed considering the digital processing is
% performed with a red pitaya and
%
% FOR REFERENCE: 
% 
% Check Gerard's PhD notes, notation is the same
% LQG: http://www.dt.fee.unicamp.br/~jbosco/ia856/KF%20and%20LQG_Makila.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%           n_samples    -- number of samples to be simulated
%           GLOBAL VARIABLES
%           feedback     -- boolean, activate or not feedback
%           drift        -- boolean, drift or not parameters (i.e., w_0 and
%                           gamma)
%           adaptive     -- boolean, activate or not machine learning to
%                           optimize parameters in the model
%
% OUTPUTS:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE FUNCTION CALL:
% 
% e.g. >> kalman_test(1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATE OF LAST UPDATE:
%
% GP Conangla, 8.2018
%
% COMMENTS: when the damping is very low, the discretization step may be
% too large and that can lead to errors (in a way, Kalman filter is a Euler
% method, which converges very badly in the discrete case)
%
% PENDING: 
% machine learning to get gamma and w_0. check what happens when Q2 is zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, x, x_kf] = kalman_test(n_samples, signal_handle, energy_handle, phonons_handle, phonons_real_handle, voltage_handle, parameters_handle)
global feedback;
global drift;
global adaptive;

%% Preliminary operations
% real signal -- parameters definitions
global stop_signal;
dt = 8e-9; % time step of red pitaya
trap = load_opt_tweezer;
a = 'drift_2n_sde';
b = 'diffusion_2n_sde';
T_period = 2*pi/trap.w_0;
num_T_oscilloscope = 6;
q_e = 1.60217662e-19; % charge electron
k_B = 1.38065e-23; % Boltzmann constant in J*K^-1
h_bar = 1.0545718e-34; % reduced planck constant

% observed signal parameters
v_to_el = 577; % V to V/m, voltage to electric field between electrodes
cal_factor1 = 1.77e7; % V/m, calibration factor for a particle of 235nm diameter
cal_factor2 = 58; % bit/nm
sigma_pd = 0.6089e-9; % sigma detector in m
sigma_rp =  0.5e-3/cal_factor1; % (std in volts)
sigma_obs = sqrt(sigma_pd^2 + sigma_rp^2);
n_energy = 1600;
t_energy = n_energy*dt;

% print parameters on GUI
text_GUI_1 = sprintf('f_0:  %.3d kHz', trap.w_0/(2*pi)/1000);
text_GUI_2 = sprintf('      p: %.3g mBar', trap.pressure);
text_GUI_3 = sprintf('      s_obs: %.3g', sigma_obs);
text_GUI = strcat(text_GUI_1, text_GUI_2, text_GUI_3);
set(parameters_handle, 'String', text_GUI);                                    %REMOVE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KF initialization
% Preallocate signal vectors
x = zeros(2, n_samples + 1); % original signal
z = zeros(1, n_samples + 1); % corrupted signal 
u = zeros(2, n_samples + 1); % control vector (for feedback)
x_kf = zeros(2, n_samples + 1); % reconstructed KF signal
t = linspace(0, dt*n_samples, n_samples + 1); % time vector
sk_vec = [-1 1]; % for the RK method
% energy signals
x2_energy = x(1, :);
x2_energy_T = x(1, :);
x2_energy_phonons = x(1, :);
x2_energy_phonons_real = x(1, :);

%% State space Matrices
% state transition model
A =  [1               dt; ...
     -dt*trap.w_0^2   1 - dt*trap.gamma/trap.m];
% observation model
H =  [1  0];
% process noise covariance (i.e., thermal bath noise)    <- sigma^2   NOOOOOOOOOOOOOOO===============================O?
R1 =  [0  0;...
      0  dt*(trap.sigma/trap.m)^2];
% observation noise covariance (i.e., measurement noise (ballanced + redpitaya))
R2 =  sigma_obs^2;
% feedback noise (i.e., AO redpitaya noise)
% control input model (assuming it is Lorentz force)
B =  [0  0;...
      0  q_e/trap.m];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% LQR optimization matrices
Q1 = [1 0;...
      0 0]; % we only care about x^2(t)
Q2 = 1e-20*eye(2); % -> we don't care about expending more in u(t)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
% state vector 
x(:, 1) = [6.9663e-09; 2.5248e-02]; %[normrnd(0, sqrt(k_B*trap.T/(trap.m*trap.w_0^2)));...
           %normrnd(0, sqrt(k_B*trap.T/(trap.m)))]
% observed vector
z(1) = x(1, 1) + sigma_obs*normrnd(0,1);
% estimation (assume initially known, for now)
x_kf(:, 1) = x(:, 1);
% initial feedback
u(:, 1) = -A*x(:, 1);
% Estimate covariance
P = [0 0;...
     0 0];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General loop: generate signal and reconstruct it
for n = 2:n_samples
    % Generate signal
    dW_t = normrnd(0,1,[2,1])*sqrt(dt);
    S_k = sk_vec(randi(2));
    k_1 = feval(a, t(n - 1), x(:,n - 1), trap)*dt + (dW_t - S_k*sqrt(dt)).*feval(b, t(n - 1), x(:,n - 1), trap);
    k_2 = feval(a, t(n), x(:,n - 1) + k_1, trap)*dt + (dW_t + S_k*sqrt(dt)).*feval(b, t(n - 1), x(:,n - 1), trap);
    k_mean = 1/2*(k_1 + k_2);
    x(:, n) = x(:, n - 1) + k_mean + B*u(:, n - 1); % add feedback B*u part
    
    % Generate observed signal (in this case, H = (1 0; 0 0), since we only observe x)
    z(n) = x(1, n) + sigma_obs*normrnd(0,1);
    
    %% Reconstruct signal with KF
    % Prediction (AHEAD OF TIME)
    x_kf_predict = A*x_kf(:, n-1) + B*u(:, n-1); % predict state estimate
    P_predict    = A*P*A' + R1; % predict estimate covariance
    
    % Update (not ahead of time)
    y = z(n) - H*x_kf_predict; % innovation
    S = R2 + H*P_predict*H'; % innovation covariance
    K = P_predict*H'/S; % Optimal Kalman gain
    
    correc = K*y;
    position = correc(1)/x_kf_predict(1);
    velocity = correc(2)/x_kf_predict(2);
    x_kf(:, n) = x_kf_predict + K*y; % update state estimate
    P = (eye(2) - K*H)*P_predict; % update estimate covariance
    y = z(n) - H*x_kf(:, n); % post-fit innovation
    
    %% LQR feedback (remember KF + LQR = LQQ)
    % Optimization criterion: Minimize J = \int_0^\infty x'*Q1*x + u'*Q2*u dt.
    % In this case, Q2 = 0 (we don't care about using a lot of energy in u(t)).
    
    % Solution of Ricatti equation, A'SA - S - A'SB(B'SB + Q2)^(-1)B'SA + Q1 = 0,
    % where S is the unknown, A, B from the state space equations, Q1, Q2
    % described above.
    [S , ~, ~] = dare(A, B, Q1, Q2);
    % Calculate control law
    L = (B'*S*B + Q2)\B'*S*A;
    % L
    % Control input
    if(n>1200)                                           %if(feedback)
        %L = L.*[0 0; 0 1];
        %-L*x_kf(:, n)
        u(:, n) = [0; -L(2, 1)*x_kf(1, n) - L(2, 2)*x_kf(2, n)]; % -L*x_kf(:, n); 
        %L(2, 1)*x_kf(1, n)/(L(2, 2)*x_kf(2, n))
        % add Red pitaya discretization
        u(:, n) = rp_resolution(u(:, n)/v_to_el)*v_to_el;
        %color_plot = '[0.8 0.8 0.1]';
        color_plot = '[0.2 0.8 0.1]';
    else
        u(:, n) = [0; 0];
        color_plot = '[0.3 0.8 0.1]';
    end
    % Calculate energy
    if(n > n_energy)
        x2_energy(n) = 1/t_energy*sum(x(1, n - n_energy : n).^2)*dt*trap.m*trap.w_0^2; % particle energy (J)
        x2_energy_T(n) = x2_energy(n)/k_B; % energy (K)
        x2_energy_phonons(n) = x2_energy(n)/(h_bar*trap.w_0) - 0.5; % energy (n)
        x2_energy_phonons_real(n) = 1/t_energy*sum(x_kf(1, n - n_energy : n).^2)*dt*trap.m*trap.w_0^2/(h_bar*trap.w_0) - 0.5; % energy (n)
    end
    n
    %% Plot the results in real time (like an oscilloscope)
    if n == 4000 %3 == mod(n, 100)
        % Plot signal
        h = figure(3);
        %axes(signal_handle);
        %cla reset;
        t_min = max(0, t(n) - num_T_oscilloscope*0.666*T_period);
        t_max =  max(num_T_oscilloscope*T_period, t(n) + num_T_oscilloscope*0.333*T_period);
        t_ini = find(t > t_min);
        nice_plot(t(t_ini:n), z(1, t_ini:n), 'time (s)', '$x(t)$', 'KF simulation over opt. tweezer trajectory', 'color', 'r', 'linewidth', '1');
        nice_plot(t(t_ini:n), x(1, t_ini:n), 'time (s)', '$x(t)$', 'KF simulation over opt. tweezer trajectory', 'linewidth', '3');
        nice_plot(t(t_ini:n), x_kf(1, t_ini:n), 'time (s)', '$x(t)$', 'KF simulation over opt. tweezer trajectory', 'linewidth', '3', 'color', color_plot);
        %nice_plot(t(t_ini:n), x_kf(2, t_ini:n)/(2*pi*125000), 'time (s)', '$x(t)$', 'KF simulation over opt. tweezer trajectory', 'linewidth', '3', 'color', '[0.8 0.4 0.1]');
        xlim([t_min, t(n)]);
        %xlim([t_min, t_max]);

%         % Plot voltage used for feedback
%         axes(voltage_handle);
%         cla reset;
%         % divided by v_to_el to get from electric field to Volts
%         nice_plot(t(t_ini:n), u(2, t_ini:n)/v_to_el, 'time (s)', '$V(t)$', 'Feedback voltage', 'linewidth', '3', 'color', '[0 0.6 1]');        
%         xlim([t_min, t_max]);
% 
%         % Plot also energy if enough (> n_energy) samples
%         if(n > n_energy)
%             % Plot energy in Kelvin
%             axes(energy_handle);
%             cla reset;
%             nice_plot(t(t_ini:n), x2_energy_T(t_ini:n), 'time (s)', 'Temp. (K)', 'Energy', 'color', 'r', 'linewidth', '2');
%             xlim([t_min, t_max]);
%             % Plot energy in number of estimated phonons
%             axes(phonons_handle);
%             cla reset;
%             nice_plot(t(t_ini:n), x2_energy_phonons(t_ini:n), 'time (s)', '$\langle n \rangle$', 'Phonons', 'color', 'b', 'linewidth', '2');
%             xlim([t_min, t_max]);
%             % Plot energy in number of real phonons
%             axes(phonons_real_handle);
%             cla reset;
%             nice_plot(t(t_ini:n), x2_energy_phonons_real(t_ini:n), 'time (s)', '$\langle n \rangle$', 'Est. phonons', 'color', '[0.6 0 1]', 'linewidth', '2');
%             xlim([t_min, t_max]);
%         else
%             % Plot energy in Kelvin
%             axes(energy_handle);
%             cla reset;
%             nice_plot(t(t_ini:n), zeros(1, length(t(t_ini:n))), 'time (s)', 'Temp. (K)', 'Energy', 'color', 'r', 'linewidth', '2');
%             xlim([t_min, t_max]);
%             % Plot energy in number of estimated phonons
%             axes(phonons_handle);
%             cla reset;
%             nice_plot(t(t_ini:n), zeros(1, length(t(t_ini:n))), 'time (s)', '$\langle n \rangle$', 'Phonons', 'color', 'b', 'linewidth', '2');
%             xlim([t_min, t_max]);
%             % Plot energy in number of real phonons
%             axes(phonons_real_handle);
%             cla reset;
%             nice_plot(t(t_ini:n), zeros(1, length(t(t_ini:n))), 'time (s)', '$\langle n \rangle$', 'Est. phonons', 'color', '[0.6 0 1]', 'linewidth', '2');
%             xlim([t_min, t_max]);
%         end
        pause(0.0001);
    end
    %% End executation if button pressed
    if(stop_signal)
        stop_signal = false;
        break;
    end
end

end

