%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Richardson extrapolation to estimate SDE moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
% inputs: a          -- the function f(t,y) (string name, as 'sin')
%         b          -- matrix with sigmas of noise term
%         t_interval -- the interval [t0,tn]
%         y0         -- the initial values
%         h          -- the step size
% output: t, y, mean_y, var_y 
%                    -- the node, a sample process and the value of 
%                       first 2 moments of y
%
% For references, see: Numerical solution of SDE through
%                      computer experiments, pg. 201
%
% GP Conangla, updated 4.7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y1, mean_y, var_y] = richardson(a, b, num_averages, t_interval, y0, h)
%% start timer
tic;
%% Calculate sample paths
% first sample path (will be used to hold current moment estimation)
t_ini = toc;
[t, y1] = runge_kutta(a, b, t_interval, y0, h);
[~, y2]  = euler_maruyama(a, b, t_interval, y0, h/2);
[~, y3]  = euler_maruyama(a, b, t_interval, y0, h/4);
[~, y4]  = euler_maruyama(a, b, t_interval, y0, h/8);

v_y1 = y1.^2;
v_y2 = y2.^2;
v_y3 = y3.^2;
v_y4 = y4.^2;

if num_averages > 1
    t_first_cycle = toc - t_ini;
    fprintf('Expected simulation time: %s\n\n', datestr((t_first_cycle*num_averages)/24/3600,'HH:MM:SS'));
    current_percent = 0;
end

for i = 2:num_averages;
    % print elapsed time and cycle number
    if i/num_averages*100 > current_percent
        fprintf('Cycle num. %d out of %d (%.3g %%). Elapsed time: %s\n', i, num_averages, i/num_averages*100, datestr((toc)/24/3600,'HH:MM:SS'));
        current_percent = current_percent + 1;
    end

    % rest of the num_averages - 1 sample paths
    [~, y1_sample_path] = runge_kutta(a, b, t_interval, y0, h);
    [~, y2_sample_path]  = euler_maruyama(a, b, t_interval, y0, h/2);
    [~, y3_sample_path]  = euler_maruyama(a, b, t_interval, y0, h/4);
    [~, y4_sample_path]  = euler_maruyama(a, b, t_interval, y0, h/8);

    % Now use the new path to update mean (we will divide by n later)
    y1 = y1 + y1_sample_path;
    y2 = y2 + y2_sample_path;
    y3 = y3 + y3_sample_path;
    y4 = y4 + y4_sample_path;

    % Same with 2nd order moment (will divide by n-1 later)
    v_y1 = v_y1 + y1_sample_path.^2;
    v_y2 = v_y2 + y2_sample_path.^2;
    v_y3 = v_y3 + y3_sample_path.^2;
    v_y4 = v_y4 + y4_sample_path.^2;
end

%% Estimation of moments
% mean
m_y1 = y1/num_averages;
m_y2 = y2/num_averages;
m_y3 = y3/num_averages;
m_y4 = y4/num_averages;
% var
var_y1 = v_y1/(num_averages - 1);
var_y2 = v_y2/(num_averages - 1);
var_y3 = v_y3/(num_averages - 1);
var_y4 = v_y4/(num_averages - 1);

%% Extrapolation
% get only vector elements that coincide with y1
m_y2 = m_y2(1:2:end,:);
m_y3 = m_y3(1:4:end,:);
m_y4 = m_y4(1:8:end,:);

var_y2 = var_y2(1:2:end,:);
var_y3 = var_y3(1:4:end,:);
var_y4 = var_y4(1:8:end,:);

% Extrapolation of mean
mean_y = 1/21*(64*m_y4 - 56*m_y3 + 14*m_y2 - m_y1);

% Extrapolation of 2nd order moment
var_y = 1/21*(64*var_y4 - 56*var_y3 + 14*var_y2 - var_y1);

%% outputs
end







