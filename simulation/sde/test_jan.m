%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test expressions with simulations
%
% GP Conangla 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[] = test_jan(t_end, h, n_averages)

% Definitions: experimental parameters (like in pdf)
trap = load_paul_trap;

% Parameters as defined in Jan & Dwight's calculations
lambda_1 = trap.m*trap.eps^2/(2*trap.gamma^3);
sigma_bound = sqrt(trap.sigma^2*trap.gamma/(trap.m*trap.eps^2));
fprintf('Value of characteristic time (tau): %.3e\n', 1/lambda_1);

% Run simulation
[t, y, mean_y, var_y] = full_numerical_sde([0 t_end], h, '2nd_order', [0 0], 'runge_kutta', n_averages, 'paul_trap');

% Compare with theory
var_short_times = trap.sigma^2/trap.gamma^2*t;

%% Plot results, 'color', 'b', 'LineWidth', 3.5
% plot var(x(t)) estimation
figure;
clf;
string_title = sprintf('%d runs of Runge-Kutta, 2nd order method', n_averages);
nice_plot(t, var_y(:,1), 'time (s)', 'Estimated $\textbf{E}[x_t^2]$', string_title, 'color', 'b');
nice_plot(t, var_short_times, 'time (s)', 'Estimated $\textbf{E}[x_t^2]$', string_title, 'color', 'r', 'linestyle', '--');

% plot var(v(t)) estimation
figure;
clf;
string_title = sprintf('%d runs of Runge-Kutta, 2nd order method', n_averages);
nice_plot(t, var_y(:,2), 'time (s)', 'Estimated $\textbf{E}[v_t^2]$', string_title, 'color', 'b');

end
