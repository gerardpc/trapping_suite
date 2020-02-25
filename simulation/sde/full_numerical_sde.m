%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full_numerical_sde.m
% 
% Calculate sample path/estimate moments of SDE
% by using Runge_kutta/Euler-Maruyama + averaging 
% Richardson extrapolation is on hold: works but is irrelevant, since main 
% source of error is number of sample paths, and not dt size
% Output is plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
%
% inputs: time_interval  -- vector with [t_ini t_end]
%         h              -- time step size
%         order          -- string: '1st_order' or '2nd_order' (order of
%                           differential equation
%         initial_values -- vector with initial values, e.g. [0; 0],
%                           can also be 'rand_ini'
%         method         -- string: 'euler_maruyama' (Euler-Maruyama) order 1/0.5
%                                   'runge_kutta' (Runge kutta) order 1/1 
%         num_averages   -- Number of averages to estimate moments
%
% output: t              -- Time vector
%         y              -- vector with numerical solution (if 1st order, 
%                           only one column with position, if 2nd order, in
%                           2nd column velocity)
%         mean_y         -- E(y) estimated from num_averages runs
%         var_y          -- E(y^2) estimated from num_averages runs 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For references, see: Numerical solution of SDE through
%                      computer experiments, pg. 201
%
% GP Conangla 2015, updated 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[t, y, mean_y, var_y] = full_numerical_sde(time_interval, h, order, initial_values, method, num_averages, varargin)

%% Switch SDE ordere
switch(order)
    case '1st_order'
        a_function = 'drift_1n_sde';
        b_function = 'diffusion_1n_sde';
    case '2nd_order'
        a_function = 'drift_2n_sde';
        b_function = 'diffusion_2n_sde';
end

% Load experiment parameters
if(any(ismember(varargin, 'paul_trap')))
    trap = load_paul_trap;
elseif(any(ismember(varargin, 'optical_tweezer')))
    trap = load_opt_tweezer;
end

%% Solve numerically
[t, y, mean_y, var_y] = solve_sde(time_interval, h, a_function, b_function, initial_values, method, num_averages, trap);
 
%% Plot results
plot_results = find(ismember(varargin, 'plot'));
if(~isempty(plot_results))
    %% Compare solution with theory
    A = trap.G_norm;
    C = trap.w_0^2;
    C2 = (C - A^2/4);
    B = trap.sigma/trap.m;
    % variance E(x^2(t)) of processes starting at [0,0]
    var_short = B^2/(2*A*C)*(1 - exp(-A*t).*(C/C2 - A^2/(4*C2)*cos(2*sqrt(C2)*t) + A/(2*sqrt(C2))*sin(2*sqrt(C2)*t)));
    
    %% The plots themselves
    % plot estimated moment
    figure(1);
    clf;
    string_title = sprintf('%d runs of %s, %s SDE', num_averages, method, order);
    nice_plot(t, var_y(:,1), 'time', 'Estimated $\textbf{E}[x^2(t)]$', string_title);
    nice_plot(t, var_short, 'time', 'Estimated $\textbf{E}[x^2(t)]$', string_title, 'color', 'r');

    % plot sample paths
    figure(2);
    clf;
    string_title = sprintf('Sample path of %s, %s SDE', method, order);
    subplot(2,1,1);
    nice_plot(t, y(:,1), 'time', 'x(t)', string_title);
    subplot(2,1,2);
    nice_plot(t, y(:,2), 'time', 'v(t)', string_title, 'color', 'r'); 
end

end









