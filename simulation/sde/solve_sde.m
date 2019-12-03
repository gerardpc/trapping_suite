%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve_sde(t_interval, h, a, b, y0, method, num_averages, varargin)
%
% Numerically solve SDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
% INPUTS: 
%         t_interval      -- the interval [t0,tn]
%         h               -- the step size
%         a               -- the function f(t,y) (string name, as 'sin')
%         b               -- matrix with sigmas of noise term
%         initial_values  -- the initial values, can also be 'rand_ini'
%         method          -- 'euler_maruyama' or 'runge_kutta'
%         num_averages    -- Perform num_averages averages to estimate
%                            moments
%         varargin
%
% OUTPUTS: t, y, mean_y, var_y 
%                         -- the node, a sample process and the value of 
%                            first 2 moments of y
%
% For references, see: Numerical solution of SDE through
%                      computer experiments, pg. 201
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, updated 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y1, mean_y, var_y] = solve_sde(t_interval, h, a, b, initial_values, method, num_averages, varargin)

%% Calculate sample paths
% First of all, check if initial condition is deterministic or random. In
% case it should be random, generate y0
y0 = get_y0(initial_values);

% Now calculate first sample path (will be used to hold current moment estimation)
tic;
t_ini = toc;
[t, y1] = feval(method, a, b, t_interval, y0, h, varargin{1});

% This vector will hold the variance estimate
v_y1 = y1.^2;

% Print expected time if num_averages > 1
if num_averages > 1
    t_first_cycle = toc - t_ini;
    fprintf('Expected simulation time: %s\n\n', datestr((t_first_cycle*num_averages)/24/3600,'HH:MM:SS'));
    current_percent = 0;
end

% Calculate rest of num_averages - 1 sample paths
for i = 2:num_averages
    % print elapsed time and cycle number
    if i/num_averages*100 > current_percent
        fprintf('Cycle num. %d out of %d (%.3g %%). Elapsed time: %s\n', i, num_averages, i/num_averages*100, datestr((toc - t_ini)/24/3600,'HH:MM:SS'));
        current_percent = current_percent + 1;
    end

    % get initial conditions
    y0 = get_y0(initial_values);
    % calculate sample path
    [~, y1_sample_path] = feval(method, a, b, t_interval, y0, h, varargin{1});

    % Use the new sample path to update mean (we will divide by n later)
    y1 = y1 + y1_sample_path;

    % Same with 2nd order moment (will divide by n-1 later)
    v_y1 = v_y1 + y1_sample_path.^2;
end

%% outputs: normalize by num_averages
mean_y = y1/num_averages;
var_y = v_y1/(num_averages - 1);

end







