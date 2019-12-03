% function[psdcell] = script_250(t_interval, dt, num_averages)

t_interval = [0 20e-3]; %20e-3
dt = 2.5e-8;
num_averages = 1000; % 1000

global xi
global f_0

%% With driving

f_0 = 10*1.60217662e-19*30/2*577;

xi = -7e12;
fprintf('\nCycle number %d\n', 2);
fprintf('\nCurrent value of xi: %.3e\n', xi);
[psdcell_50mbar{8}, f_period] = calc_traces(t_interval, dt, num_averages, 'figure_num', num2str(2), 'optical_tweezer');
pause(1);
