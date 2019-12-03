%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc_traces(time_interval, dt)
%
% Simulate traces of particle in trap, then estimate PSD from it and fit
% analytic expression to extract parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%           time_interval    -- [t(ini) t(end)], interval where trace is
%                               simulated
%           dt               -- precision of simulation and, indirectly,
%                               sampling period. This doesn't mean that it
%                               can be set to any value! For dt too large
%                               the method doesn't converge
%           num_traces       -- number of simulated traces and thus PSDs
%           varargin         -- 'figure_num', '2' -> number of figure to plot
%                            -- 'print_log', 'results.txt' -> print to log
%                               instead of on screen
%                            -- 'optical_tweezer'/'paul_trap' (to get
%                               parameters from file
%
% OUTPUTS:  
%           smooth_psd       -- Estimated PSD from averaging num_traces
%                               periodograms
%           f_period         -- freq. vector for smooth_psd
%
% e.g. >> calc_traces([0 8e-3], 1e-7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[smooth_psd, f_period] = calc_traces(time_interval, dt, num_traces, varargin)

%% Variable arguments in
figure_num = find(ismember(varargin, 'figure_num'));
if(~isempty(figure_num))
    figure_num = str2double(varargin{figure_num + 1});
else
    figure_num = 1;
end

print_log = find(ismember(varargin, 'print_log'));
if(~isempty(print_log))
    file_name = varargin{print_log + 1};
    fileID = fopen(file_name, 'w'); % fprintf writes to specified file
else
    fileID = 1; % otherwise just print at std output, i.e., command line
end

%% welcome message
fprintf(fileID, '\n');
fprintf(fileID, '#########################################################\n');
fprintf(fileID, 'Starting calculations!\n');

%% parameters
% Load experiment parameters (opt trap/paul trap/other)
if(any(ismember(varargin, 'optical_tweezer')))
    trap = load_opt_tweezer;
elseif(any(ismember(varargin, 'paul_trap')))
    trap = load_paul_trap;
end

n_fft = 2^24; % number of points in FFT
fs = 1/dt;
num_samples = round((time_interval(end) - time_interval(1))/dt) + 1;

%% preallocate matrix space
psd_vector = zeros(n_fft/2 + 1, 1); % divided by 2 because periodogram only gives one side of PSD

%% Simulate
% start_timer
t_ini = toc;
current_percent = 0;

for i = 1:num_traces
    % Calculate trace
    [~, y, ~, ~] = full_numerical_sde(time_interval, dt, '2nd_order', [0 0], 'runge_kutta', 1, varargin{:});
    
    % Calculate PSD: periodogram
    [pxx, f_period] = periodogram(y(:, 1), rectwin(num_samples), n_fft, fs); 
    psd_vector = psd_vector + pxx/2; % divided by 2 because it's double sided
    
    % Print expected time if num_traces > 1 and it's the first trace
    if num_traces > 1 && i == 1
        fprintf(fileID, '\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n');
        t_first_cycle = toc - t_ini;
        fprintf(fileID, 'Calculating sample paths.\n');
        fprintf(fileID, 'Expected simulation time: %s\n\n', datestr((t_first_cycle*num_traces)/24/3600,'dd HH:MM:SS'));
    end
    
    % print elapsed time and cycle number
    if i/num_traces*100 > current_percent
        fprintf(fileID, 'Trace num. %d out of %d (%.3g %%). Elapsed time: %s\n', i, num_traces, i/num_traces*100, datestr((toc - t_ini)/24/3600,'HH:MM:SS'));
        current_percent = current_percent + 1;
    end
end

%% Estimate PSD from simulated traces
fprintf(fileID, '\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n');
fprintf(fileID, 'Estimating PSD\n\n');
smooth_psd = psd_vector/num_traces;
var_y = std(y(:, 1))^2;
% Check if process exploded, x_t -> infty
if isnan(var_y)
    % Don't calculate PSD or fit, plot trace of process
    fprintf(fileID, '\nProcess exploded! Skipping PSD\n');
else    
    df = f_period(2) - f_period(1);
    var_y_w =  2*trapz(smooth_psd)*df;
    fprintf(fileID, 'Equipartition check: E(x^2)/expected value\n');
    fprintf(fileID, 'From time trace: %.3g\n', trap.k_B*trap.T/(trap.m*trap.w_0^2)/var_y);
    fprintf(fileID, 'From PSD:        %.3g\n', trap.k_B*trap.T/(trap.m*trap.w_0^2)/var_y_w);

    %% Fit function to PSD (assuming harmonic oscillator -> case of optical tweezer)
    fprintf(fileID, '\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n');
    fprintf(fileID, 'Fitting Lorentzian to data\n\n');
    [~, ~, ww_fit, curve_fit, curve_fit_2] = fit_lorentzian(smooth_psd(1:2^6:end), 2*pi*f_period(1:2^6:end));

    %% Plot results (PSD + fit)
    figure(figure_num);
    clf;
    nice_plot(f_period, (smooth_psd), '$f$ (Hz)', '$S_x$', 'Estimated PSD from Monte Carlo traces', 'semilogy', 'linewidth', '1');
    nice_plot(ww_fit/(2*pi), (curve_fit), '$f$ (Hz)', '$S_x$', 'Estimated PSD from Monte Carlo traces', 'color', 'r', 'semilogy');
    nice_plot(ww_fit/(2*pi), (curve_fit_2), '$f$ (Hz)', '$S_x$', 'Estimated PSD from Monte Carlo traces', 'color', 'y', 'linestyle', '--', 'semilogy');
    xlim([0 2*trap.w_0/(2*pi)]);
end

%% Finished, goodbye message
fprintf(fileID, '\n');
fprintf(fileID, '#########################################################\n');
fprintf(fileID, 'THE END\n');
if fileID ~= 1
    fclose(fileID);
end
end















