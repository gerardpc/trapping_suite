%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load_opt_tweezer()
%
% Read optical tweezers parameters from file 'particle_tweezer.txt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%         opt        -- struct with fields (e.g opt.freq or opt.gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[opt] = load_opt_tweezer()
% read file of optical tweezer parameters
file_ID = fopen('particle_tweezer.txt');
text_opt = textscan(file_ID, '%s', 36);
fclose(file_ID);

% put numerical values in vector
values_opt = str2double(text_opt{1});
values_opt(isnan(values_opt) == 1) = [];

% add rest of parameters
T = values_opt(1);
k_B = values_opt(2);
r = values_opt(3);
pressure = values_opt(4);
density = values_opt(5);
freq = values_opt(6);
m = 4/3*pi*r^3*density;
gamma_air = 6*pi*r*1.84e-5;
gamma = gamma_air*pressure/1010;
G_norm = gamma/m;
w_0 = 2*pi*freq;
sigma = sqrt(2*k_B*T*gamma);

% put values of parameters under struct with understandable names
opt = struct('T', T, 'k_B', k_B, 'r', r,...
    'pressure', pressure, 'density', density, 'freq', freq,...
    'm', m, 'gamma_air', gamma_air, 'gamma', gamma, 'G_norm', G_norm, 'w_0', w_0, ...
    'sigma', sigma);
end